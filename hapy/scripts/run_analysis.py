#!/usr/bin/env python
"""
Run analysis on a single galaxy cutout set (1 GNU-parallel task).


"""

import argparse
from pathlib import Path
import glob
import sys
import logging
from datetime import datetime

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table

import json
import re
import time

import hapy
from hapy.ellipse.photometry import run_ellipse_photometry
from hapy.ellipse.utils import infer_ellipse_from_r_cutout, ellipse_missing
from hapy.galfittools.rungalfit import RunGalfit
from hapy.imagetools.imutils import get_pixel_scale_from_filename
from hapy.imagetools.plotting import plot_mask_ellipse_diagnostic
from hapy.imagetools.plotting import plot_segmentation_diagnostic
from hapy.masktools.api import MaskEngine, EllipseParams
from hapy.masktools.gaia import make_gaia_mask,  get_gaia_stars, galaxy_overlaps_bright_star
from hapy.masktools.maskops import distance_to_nearest_mask, largest_mask_region, ellipse_mask_fraction
#from hapy.masktools.types import build_ell0_from_metadata
from hapy.hatools.results import write_result_row_ecsv
from hapy.geometry.adapters import pa_ccw_north_to_photutils_theta, photutils_theta_to_pa_ccw_north
from hapy.utils.paths import astromatic_dir
from hapy.hatools.utils import zp_scale_r_to_ha

#from hapy.utils.logging_utils import setup_logging

from hapy.ellipse.profile_summary import summarize_dual_profiles


def copy_image2_fields_to_row(e, row, prefix):
    fields = [
        ("SKYSTD_ADU", "sky_noise2"),
        ("SKYMED_ADU", "sky2"),
        ("SKYSTD_PHYS", "im2_skynoise"),
        ("M20", "M20_2"),
        ("ASYM", "asym2"),
        ("ASYM_ERR", "asym2_err"),
        ("SCALE_ADU_CGS", "uconversion2"),
    ]

    for outk, attr in fields:
        sv = _scalar(getattr(e, attr, None))
        if sv is not None:
            row[f"{prefix}_{outk}"] = sv


def copy_hapy_cs_fields_to_row(e, row, prefix, pixscale):
    fields = [
        ("HAPY_NPIX", "H_HAPY_NPIX"),
        ("HAPY_FILLFRAC", "H_HAPY_FILLFRAC"),
        ("HAPY_SNP_ALL", "H_HAPY_SNP_ALL"),
        ("HAPY_SNP_DET", "H_HAPY_SNP_DET"),
        ("H_HAPY_GINI_THRESHOLD", "ha_gini_threshold"),
        ("R_HAPY_GINI_THRESHOLD", "r_gini_threshold"),        
        ("HAPY_GINI", "H_HAPY_GINI"),
        ("HAPY_M20", "H_HAPY_M20"),
        ("HAPY_ASYM", "H_HAPY_ASYM"),
        ("HAPY_ASYM_ERR", "H_HAPY_ASYM_ERR"),
        ("HAPY_MTOT", "H_HAPY_MTOT"),
        ("HAPY_M20SUM", "H_HAPY_M20SUM"),
        ("HAPY_FLUX_SEG", "H_HAPY_FLUX_SEG"),
        ("HAPY_MTOT2", "H_HAPY_MTOT2"),
    ]

    for outk, attr in fields:
        sv = _scalar(getattr(e, attr, None))
        if sv is not None:
            row[f"{prefix}_{outk}"] = sv

    rmom = _scalar(getattr(e, "H_HAPY_RMOM", None))
    if rmom is not None:
        row[f"{prefix}_HAPY_RMOM_ARCSEC"] = rmom * pixscale

    row[f"{prefix}_HAPY_MORPH_OK"] = bool(getattr(e, "HAPY_MORPH_OK", False))
    row[f"{prefix}_HAPY_MORPH_FLAG"] = int(getattr(e, "HAPY_MORPH_FLAG", 0))

def ellipse_image_coverage(data, ell0_params):
    import numpy as np
    from photutils.aperture import EllipticalAperture

    ap = EllipticalAperture(
        (ell0_params.xc, ell0_params.yc),
        a=ell0_params.sma_pix,
        b=ell0_params.sma_pix * ell0_params.ba,
        theta=np.deg2rad(ell0_params.theta_deg),
    )

    aper_mask = ap.to_mask(method="center")
    in_ellipse = aper_mask.to_image(data.shape) > 0

    n_total = int(np.sum(in_ellipse))
    good = in_ellipse & np.isfinite(data)

    n_good = int(np.sum(good))
    n_missing = n_total - n_good
    frac_missing = n_missing / n_total if n_total > 0 else np.nan

    return {
        "npix_total": n_total,
        "npix_good": n_good,
        "npix_missing": n_missing,
        "missing_frac": frac_missing,
    }



def prefix_dict_keys(d, prefix):
    return {f"{prefix}_{k}": v for k, v in d.items()}
            
def init_cutout_logger(tag: str, cutdir: str | Path, level: str = "INFO",
                       log_to_console: bool = False, log_dir: str | Path | None = None):
    """
    Create a per-cutout logger writing to <cutout_dir>/<tag>.log (or log_dir).
    Safe for parallel runs because each cutout has its own log file.
    """
    cutdir = Path(cutdir)
    #cutdir = root.parent

    if log_dir is None:
        log_dir = cutdir / "logs" 
    else:
        log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir /f"{tag}.analysis.log"

    logger = logging.getLogger(f"hapy.run_analysis.{tag}")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.propagate = False  # don't double-log via root logger

    # Prevent duplicate handlers if init is called twice in same process
    if logger.handlers:
        for h in list(logger.handlers):
            logger.removeHandler(h)

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | pid=%(process)d | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(log_path, mode="w")  # overwrite each run
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    if log_to_console:
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        logger.addHandler(sh)

    logger.info(f"Log file: {log_path}")
    return logger, log_path


def _default_sex_config() -> str:
    # Preferred: resolve relative to installed/imported hapy package
    try:
        pkg_root = Path(hapy.__file__).resolve().parent
        p = pkg_root / "astromatic" / "default.sex.HDI.mask"
        if p.exists():
            return str(p)
    except Exception:
        pass

    # Secondary: try importlib.resources relative to "hapy" package root
    try:
        p = resources.files("hapy") / "astromatic" / "default.sex.HDI.mask"
        return str(p)
    except Exception:
        pass

    # Last resort: relative to this file (repo layout)
    return str(Path(__file__).resolve().parents[1] / "astromatic" / "default.sex.HDI.mask")
    
def _pick_one(pattern: str) -> str | None:
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None

def _progress_cb(stage, fraction, message=None):
    # keep stdout simple for GNU parallel logs
    if message:
        print(stage, fraction, message)
    else:
        print(stage, fraction)


def _scalar(v):
    """Best-effort conversion to JSON/ECSV-friendly scalar."""
    if v is None:
        return None
    # astropy Quantity
    try:
        v = v.value
    except Exception:
        pass
    # numpy scalar
    try:
        import numpy as np
        if isinstance(v, (np.generic,)):
            return v.item()
        if isinstance(v, (np.ndarray, list, tuple)):
            return None  # skip arrays here on purpose
    except Exception:
        pass
    # python numeric / bool / str
    if isinstance(v, (int, float, bool, str)):
        return v
    # last resort
    try:
        return float(v)
    except Exception:
        return str(v)

def _get(e, name, default=None):
    return getattr(e, name, default)

def _finite(x):
    try:
        return np.isfinite(x)
    except Exception:
        return False

def _finite_dict(d):
    return all(np.isfinite(float(v)) for v in d.values())



def _galfit_stage(rg, args, init, do_conv: bool, n_hi=4.0, logger=None):
    """
    init: dict with xobj,yobj,mag,rad,nsersic,BA,PA, first_time
    returns: (res, meta)
    """
    if do_conv:
        rg.enable_convolution()
        stage = "CV"
    else:
        if hasattr(rg, "disable_convolution"):
            rg.disable_convolution()
        stage = "NC"

    if logger:
        logger.info(f"GALFIT {stage} start init={init}")
        
    def unstable(res):
        if _scalar(res.error) != 0: return True
        if _scalar(res.comp1.numerical_error_flag) != 0: return True
        if _scalar(res.chi2nu) > 5: return True
        if _scalar(res.comp1.re) <= 0: return True
        n = _scalar(res.comp1.n)
        if n > n_hi: return True
        ba = _scalar(res.comp1.ba)
        if ba <= 0.05 or ba > 1.0: return True
        return False

    rg.set_sersic_params(
        xobj=init["xobj"], yobj=init["yobj"],
        mag=init["mag"], rad=init["rad"],
        nsersic=init["nsersic"], BA=init["BA"], PA=init["PA"],
        fitmag=1, fitcenter=1, fitrad=1, fitBA=1, fitPA=1, fitn=1,
        first_time=init.get("first_time", 0),
    )
    rg.set_sky(args.sky)
    #res = rg.run_and_parse()

    try:
        res = rg.run_and_parse()
    except Exception as e:
        if logger:
            logger.exception(f"GALFIT {stage} failed during run_and_parse: {e}")
        raise
    
    meta = {"stage": stage, "rerun_fixed_n": False, "unstable": False}

    if logger:
        try:
            logger.info(
                f"GALFIT {stage} result: chi2nu={_scalar(res.chi2nu)} "
                f"re={_scalar(res.comp1.re)} n={_scalar(res.comp1.n)} "
                f"ba={_scalar(res.comp1.ba)} pa={_scalar(res.comp1.pa)} "
                f"numerr={_scalar(res.comp1.numerical_error_flag)} err={_scalar(res.error)}"
            )
        except Exception:
            logger.info(f"GALFIT {stage} result parsed (logging fields failed)")
            
    # rerun if high n
    rerun=False
    if _scalar(res.comp1.n) > n_hi:
        meta["rerun_fixed_n"] = True
        if logger:
            logger.warning(f"GALFIT {stage}: high-n detected (n={_scalar(res.comp1.n)}); rerunning with same initial conditions except fixing n=4")
        rg.set_sersic_params(
            xobj=init["xobj"], yobj=init["yobj"],
            mag=init["mag"], rad=init["rad"],
            nsersic=4.0,
            BA=init["BA"], PA=init["PA"],
            fitmag=1, fitcenter=1, fitrad=1, fitBA=1, fitPA=1, fitn=0,
            first_time=0,
        )
        rerun = True
    if rerun:
        rg.set_sky(args.sky)
        #res = rg.run_and_parse()
        try:
            res = rg.run_and_parse()
        except Exception as e:
            if logger:
                logger.exception(f"GALFIT {stage} rerun (fixed n) failed: {e}")
            raise
        if logger:
            logger.info(
                f"GALFIT {stage} rerun result: chi2nu={_scalar(res.chi2nu)} "
                f"re={_scalar(res.comp1.re)} ba={_scalar(res.comp1.ba)} pa={_scalar(res.comp1.pa)} "
                f"numerr={_scalar(res.comp1.numerical_error_flag)} err={_scalar(res.error)}"
            )       

    meta["unstable"] = unstable(res)
    if logger:
        logger.info(f"GALFIT {stage} done: unstable={meta['unstable']} rerun_fixed_n={meta['rerun_fixed_n']}")
    return res, meta

def _store_galfit(row, res, prefix, pixscale):
    row[f"{prefix}XC"] = _scalar(res.comp1.xc)
    row[f"{prefix}XC_ERR"] = _scalar(res.comp1.xc_err)
    row[f"{prefix}YC"] = _scalar(res.comp1.yc)
    row[f"{prefix}YC_ERR"] = _scalar(res.comp1.yc_err)

    row[f"{prefix}MAG"] = _scalar(res.comp1.mag)
    row[f"{prefix}MAG_ERR"] = _scalar(res.comp1.mag_err)

    # convert to arcsec
    r = _scalar(res.comp1.re)
    if r is not None:
        r = r * pixscale
    row[f"{prefix}RE_ARCSEC"] = r

    r = _scalar(res.comp1.re_err)
    if r is not None:
        r = r * pixscale
    
    row[f"{prefix}RE_ERR_ARCSEC"] = r

    
    row[f"{prefix}N"] = _scalar(res.comp1.n)
    row[f"{prefix}N_ERR"] = _scalar(res.comp1.n_err)
    row[f"{prefix}BA"] = _scalar(res.comp1.ba)
    row[f"{prefix}BA_ERR"] = _scalar(res.comp1.ba_err)
    row[f"{prefix}PA"] = _scalar(res.comp1.pa)
    row[f"{prefix}PA_ERR"] = _scalar(res.comp1.pa_err)

    row[f"{prefix}SKY"] = _scalar(res.sky)
    row[f"{prefix}SKY_ERR"] = _scalar(res.sky_err)
    row[f"{prefix}CHISQ"] = _scalar(res.chi2nu)
    row[f"{prefix}NUMERR"] = _scalar(res.comp1.numerical_error_flag)
    row[f"{prefix}ERROR"] = _scalar(res.error)


def print_statmorph(mobj):
    for k in mobj.__dict__.keys():
        if k.startswith('_'):
            continue
        print(f"\t{k}: {mobj.__dict__[k]}")


# ---- statmorph (best-effort; only a few key fields to start) ----

def _pull_statmorph(row, prefix, mobj,pixscale):
    if mobj is None:
        return
    
    for outk, attr in [
            ("FLAG", "flag"),
            ("XCENTROID", "xc_centroid"),
            ("YCENTROID", "yc_centroid"),
            ("GINI", "gini"),
            ("M20", "m20"),
            ("C", "concentration"),
            ("A", "asymmetry"),
            ("S", "smoothness"),
            ("RPETRO_ELLIP_ARCSEC", "rpetro_ellip"),
            ("RHALF_ELLIP_ARCSEC", "rhalf_ellip"),
            ("R20_ARCSEC", "r20"),
            ("R50_ARCSEC", "r50"),
            ("R80_ARCSEC", "r80"),
            ("RMAX_CIRCLE_ARCSEC","rmax_circ"),
            ("RMAX_ELLIP_ARCSEC","rmax_ellip"),            
            ("SERSIC_AMP", "sersic_amplitude"),
            ("SERSIC_RHALF_ARCSEC","sersic_rhalf"),
            ("SERSIC_N","sersic_n"),
            ("SERSIC_XC","sersic_xc"),
            ("SERSIC_YC","sersic_yc"),
            ("SERSIC_ELLIP","sersic_ellip"),
            ("SERSIC_THETA","sersic_theta"),
            ("SERSIC_CHISQ_DOF","sersic_chi2_dof"),
            ("SERSIC_FLAG","flag_sersic"),
            ("SN_PER_PIXEL","sn_per_pixel"),
            ("SKY_MEAN","sky_mean"),
            ("SKY_MEDIAN","sky_median"),
            ("SKY_SIGMA","sky_sigma"),            
                    ]:
        #row[f"{prefix}_{outk}"] = _scalar(getattr(mobj, attr))
        try:
            r = _scalar(getattr(mobj, attr))
            if "_ARCSEC" in outk and r:
                r = r * pixscale
            row[f"{prefix}_{outk}"] = r 
        except Exception:
            pass

def cutfile(cutdir, tag, suffix):
    return cutdir / f"{tag}-{suffix}.fits"
        
def resolve_psf_path(psfdir, parent_rimage):
    """Return PSF path derived from parent R coadd filename, or None."""
    if not psfdir or not parent_rimage:
        return None
    psf_name = str(parent_rimage).replace(".fits", "-psf.fits")
    psf_path = Path(psfdir) / psf_name
    return psf_path if psf_path.exists() else None
    
def initialize_result_row():
    """Return a fully populated results-row dict with frozen schema."""

    row = {}
 

    # ---------- virgo identifiers ----------
    row["VFID"] = ""      # e.g., "VFID3084"
    row["GALNAME"] = ""   # e.g., "NGC3512" (optional but handy)
    row["OBJID"] = ""
    # ---------- coordinates ----------
    row["RA"] = np.nan
    row["DEC"] = np.nan
    row["REDSHIFT"] = np.nan
    row["VR"] = np.nan    
    # ---- meta data ---
    row["HAPY_VERSION"] = ""
    row["RUN_DATE"] = ""
    for k in [
        "TELESCOPE",
        "DATEOBS",
        "POINTING",
        "SCHEME",
        "PARENT_RIMAGE",
        "PARENT_HIMAGE",
        "HFILTER"
    ]:
        row[k] = ""

    row["X_PARENT"] = np.nan
    row["Y_PARENT"] = np.nan    
    row["PIXSCALE"] = np.nan

    row["CSZP_SOURCE"] = ""
    row["CSZP_LOCAL_SKY"] = ""
    # removing these
    not_needed = ["R_FITS", "CS_FITS"]
    # ---------- identity ----------
    for k in [
        "TAG", "CUTDIR",
         "MASK_FITS","MASK_SOURCE","PSF_FITS",
         "R_FITS", "CS_FITS","SIGMA_FITS",
         "RFILTER_FILENAME", "RFILTER_CENTER","RFILTER_WIDTH",
         "HFILTER_FILENAME", "HFILTER_CENTER","HFILTER_WIDTH",         
            ]:
        row[k] = ""
        
    #row["psf_ok"] = False

    row["PSF_SOURCE"] = ""   # "cli" | "psf_dir" | ""

    for k in ["STAGE", "STATUS"]:
        row[k] = ""

    for k in ["MASK_SEC", "PHOT_SEC","HAPY_MORPH_SEC","PROFILES_SEC","CSGR_SEC", "SM_SEC", "GAL_NC_SEC","GAL_CV_SEC", "TOTAL_SEC"]:
        row[k] = np.nan    
   
    # ---------- pipeline status ----------
    for k in ["PSF_OK", "MASK_OK", "PHOT_OK", \
                  "HAPY_MORPH_OK",\
                  "R_PROFILE_OK","H_PROFILE_OK",\
                  "R_SM_OK","H_SM_OK",\
                  "GAL_NC_OK", "GAL_CV_OK"]:#, "galfit_ok"]:
        row[k] = False

    row["GAL_CV_INIT_FROM_NC"] = False
    row["GAL_NC_RERUN_FIXEDN"] = False
    row["GAL_CV_RERUN_FIXEDN"] = False

    row["ELL_MISMATCH"] = False

    for k in [
        "R_PETRO_OK",
        "R_EXPFIT_OK",
        "R_LOGFIT_OK",
    ]:
        row[k] = False

    for k in [
        "H_PETRO_OK",
        "H_EXPFIT_OK",
        "H_LOGFIT_OK",
    ]:
        row[k] = False

    row["BRIGHT_STAR_FLAG"] = False
    row["BRIGHT_STAR_DIST_ARCSEC"] = np.nan
    row["BRIGHT_STAR_MASKRAD_ARCSEC"] = np.nan
    row["BRIGHT_STAR_MAG"] = np.nan 
    row["ELL0_MASKFRAC"] = np.nan
    row["ELL0_MASK_WARN"] = False
    row["ELL0_NMASKPIX"] = np.nan
    row["ELL0_NTOTPIX"] = np.nan
    

    # ---------- cutout properties ----------
    for k in [
        "CUTOUT_SCALE", "CUTOUT_XSIZE", "CUTOUT_YSIZE",
        "FILTER_CORRECTION", "FILTER_RATIO"
    ]:
        row[k] = np.nan

    # ---------- ELL0 ----------
    for k in [
        "ELL0_SMA_ARCSEC", "ELL0_BA", "ELL0_PA_DEG",
        "ELL0_XC", "ELL0_YC"
    ]:
        row[k] = np.nan
    row["ELL0_SOURCE"] = ""

    # # ----- ELL0 good/back pixels
    row["CUTOUT_ELL0_MISSING_FRAC_R"] = np.nan
    row["CUTOUT_ELL0_MISSING_FRAC_H"] = np.nan
    row["CUTOUT_ELL0_MISSING_FRAC_MAX"] = np.nan

    row["CUTOUT_ELL0_NPIX_TOTAL_R"] = np.nan
    row["CUTOUT_ELL0_NPIX_TOTAL_H"] = np.nan

    row["CUTOUT_ELL0_NPIX_ONIMAGE_R"] = np.nan
    row["CUTOUT_ELL0_NPIX_ONIMAGE_H"] = np.nan

    row["CUTOUT_ELL0_NPIX_GOOD_R"] = np.nan
    row["CUTOUT_ELL0_NPIX_GOOD_H"] = np.nan



    # ---------- ellipse ----------
    for k in [
        "ELLIP_XCENTROID", "ELLIP_YCENTROID",
        "ELLIP_SMA_PIX","ELLIP_SMA_ARCSEC","ELLIP_BA", "ELLIP_B_ARCSEC",
        "ELLIP_EPS", "ELLIP_THETA_RAD", "ELLIP_PA_DEG",
        "R_ELLIP_GINI","H_ELLIP_GINI", "ELLIP_SOURCE_SUM",
        "ELLIP_BA", "ELLIP_SEGMENT_FLUX",
        "ELLIP_SEGMENT_MAG",
    ]:
        row[k] = np.nan

    row["ELLIP_CENTER_METHOD"] = ""
    row["AREA_GUESS_ELLIPSE_PIX"] = np.nan
    row["AREA_GUESS_ELLIPSE_UNMASKED_PIX"] = np.nan
    row["MASKFRAC_GUESS_ELLIPSE"] = np.nan

    # ---------- R band ----------
    for k in [
        "R_FWHM_PSF","R_FWHM_SE", "R_SKYSTD_ADU", "R_SKYMED_ADU",
        "R_SKYSTD_PHYS", "R_M20", "R_ASYM", "R_ASYM_ERR"
    ]:
        row[k] = np.nan

    # ---------- H band ----------
    for k in [
        "H_FWHM_PSF", "H_FWHM_SE", "H_SKYSTD_ADU", "H_SKYMED_ADU",
        "H_SKYSTD_PHYS", "H_M20", "H_ASYM", "H_ASYM_ERR"
    ]:
        row[k] = np.nan

    row["R_SCALE_ADU_CGS"] = np.nan
    row["H_SCALE_ADU_CGS"] = np.nan

    # geometry
    row["R_HAPY_NPIX"] = np.nan
    row["H_HAPY_NPIX"] = np.nan
    row["H_HAPY_FILLFRAC"] = np.nan

    # snr
    row["R_HAPY_SNP_ALL"] = np.nan    
    row["H_HAPY_SNP_ALL"] = np.nan
    row["H_HAPY_SNP_DET"] = np.nan
    row["H_HAPY_GINI_THRESHOLD"] = np.nan
    row["R_HAPY_GINI_THRESHOLD"] = np.nan    

    # morphology
    row["R_HAPY_XC"] = np.nan # R and H have the same center
    row["R_HAPY_YC"] = np.nan
    
    row["R_HAPY_GINI"] = np.nan
    row["H_HAPY_GINI"] = np.nan
    row["R_HAPY_M20"] = np.nan
    row["H_HAPY_M20"] = np.nan
    row["R_HAPY_ASYM"] = np.nan
    row["H_HAPY_ASYM"] = np.nan
    row["R_HAPY_ASYM_ERR"] = np.nan
    row["H_HAPY_ASYM_ERR"] = np.nan
    row["R_HAPY_ASYM_XC"] = np.nan
    row["R_HAPY_ASYM_YC"] = np.nan        
    row["R_HAPY_MTOT"] = np.nan
    row["H_HAPY_MTOT"] = np.nan
    row["R_HAPY_M20SUM"] = np.nan
    row["H_HAPY_M20SUM"] = np.nan        

    row["H_HAPY_FLUX_SEG"] = np.nan        
    row["H_HAPY_MTOT2"] = np.nan        
    row["H_HAPY_RMOM_ARCSEC"] = np.nan        
    
    row["R_HAPY_FLUX_SEG"] = np.nan        
    row["R_HAPY_MTOT2"] = np.nan        
    row["R_HAPY_RMOM_ARCSEC"] = np.nan

    
    
    # HAPY_MORPH_FLAG bit meanings:
    # 1  = empty or invalid r-band morphology mask
    # 2  = no Halpha image available
    # 4  = no Halpha pixels above threshold inside r-mask
    # 8  = invalid/non-finite metric or sky noise
    # 16 = exception during HAPY morphology calculation
    row["HAPY_MORPH_FLAG"] = 0

    # ---------- mismatch ----------
    for k in ["ELL_DC_PX", "ELL_DBA", "ELL_DPA_DEG", "ELL_SMA_RATIO"]:
        row[k] = np.nan

    
    # ---------- statmorph ----------
    sm_suffixes = [
        "FLAG",
        "XCENTROID", "YCENTROID", "GINI", "M20",
        "C", "A", "S",
        "RPETRO_ELLIP_ARCSEC", "RHALF_ELLIP_ARCSEC",
        "R20_ARCSEC", "R50_ARCSEC", "R80_ARCSEC",
        "RMAX_CIRCLE_ARCSEC", "RMAX_ELLIP_ARCSEC",
        "SERSIC_AMP", 
        "SERSIC_RHALF_ARCSEC",
        "SERSIC_N",
        "SERSIC_XC",
        "SERSIC_YC",
        "SERSIC_ELLIP",
        "SERSIC_THETA",
        "SERSIC_CHISQ_DOF",
        "SERSIC_FLAG",
        "SN_PER_PIXEL",
        "SKY_MEAN",
        "SKY_MEDIAN",
        "SKY_SIGMA",        
    ]

    for band in ["R", "H"]:
        for s in sm_suffixes:
            row[f"{band}_SM_{s}"] = np.nan


    # ---------- R profile summary ----------
    for k in [
        "R_PROFILE_NGOOD",
        "R_PROFILE_MASKFRAC_MAX",
        "R25_ARCSEC", "R25_PIX",
        "R50_ARCSEC", "R50_PIX",
        "R75_ARCSEC", "R75_PIX",
        "R24_ARCSEC", "R24_ARCSEC_ERR",
        "R24_MAG", "R24_MAG_ERR",
        "R25_ISO_ARCSEC", "R25_ISO_ARCSEC_ERR",
        "R25_ISO_MAG", "R25_ISO_MAG_ERR",
        "R25P5_ARCSEC", "R25P5_ARCSEC_ERR",
        "R25P5_MAG", "R25P5_MAG_ERR",
        "R24_VEGA_ARCSEC", "R24_VEGA_ARCSEC_ERR",
        "R24_VEGA_MAG", "R24_VEGA_MAG_ERR",
        "R25_VEGA_ARCSEC", "R25_VEGA_ARCSEC_ERR",
        "R25_VEGA_MAG", "R25_VEGA_MAG_ERR",
        "R30R24_FLUX_CGS", "R30R24_FLUX_CGS_ERR",
        "R24_FLUX_CGS", "R24_FLUX_CGS_ERR",
        "R_C30", "R_C30_ERR",
        "R_PETRO_RAD_ARCSEC",
        "R_PETRO_FLUX",
        "R_PETRO_FLUX_CGS", "R_PETRO_FLUX_CGS_ERR",
        "R_PETRO_MAG",
        "R_PETRO_R50_ARCSEC",
        "R_PETRO_R90_ARCSEC",
        "R_PETRO_CON",
        "R_EXPFIT_I0", "R_EXPFIT_K", "R_EXPFIT_RE_ARCSEC",
        "R_LOGFIT_SLOPE", "R_LOGFIT_INTERCEPT", "R_LOGFIT_RE_ARCSEC",
        "R_TOT_MAG_SNR",
        "R_TOT_FLUX_CGS_SNR", "R_TOT_FLUX_CGS_SNR_ERR",
        "R_SNR_TRUNC_ARCSEC", 
        "R_PROFILE_PEAK_BIN", "R_PROFILE_PEAK_SMA"
    ]:
        row[k] = np.nan

    row["R_PROFILE_NONCENTRAL_PEAK"] = False

    # ---------- Halpha profile summary ----------
    for k in [
        "H_PROFILE_NGOOD",
        "H_PROFILE_LONGRUN",
        "H_NDET_RUNS",
        "H_PROFILE_MASKFRAC_MAX",
        "H_MAXDET_ARCSEC",
        "H_MAXDET_PIX",
        "H_TOT_FLUX_CGS",
        "H_TOT_FLUX_CGS_ERR",
        "H_SNR_TRUNC_ARCSEC",
        "H25_ARCSEC",
        "H25_PIX",
        "H50_ARCSEC",
        "H50_PIX",
        "H75_ARCSEC",
        "H75_PIX",
        "H_ISO5E17_ARCSEC",
        "H_ISO5E17_ARCSEC_ERR",
        "H_ISO5E17_FLUX_CGS",
        "H_ISO5E17_FLUX_CGS_ERR",
        "H_ISO17E18_ARCSEC",
        "H_ISO17E18_ARCSEC_ERR",
        "H_ISO17E18_FLUX_CGS",
        "H_ISO17E18_FLUX_CGS_ERR",
        "H30R24_FLUX_CGS",
        "H30R24_FLUX_CGS_ERR",
        "H_R24_FLUX_CGS",
        "H_R24_FLUX_CGS_ERR",
        "H_C30_R24",
        "H_C30_R24_ERR",
        "H_R95_R24_ARCSEC",
        "H_PETRO_RAD_ARCSEC",
        "H_PETRO_FLUX",
        "H_PETRO_FLUX_CGS",
        "H_PETRO_FLUX_CGS_ERR",
        "H_PETRO_MAG",
        "H_PETRO_R50_ARCSEC",
        "H_PETRO_R90_ARCSEC",
        "H_PETRO_CON",
        "H_EXPFIT_I0",
        "H_EXPFIT_K",
        "H_EXPFIT_RE_ARCSEC",
        "H_LOGFIT_SLOPE",
        "H_LOGFIT_INTERCEPT",
        "H_LOGFIT_RE_ARCSEC",
    ]:
        row[k] = np.nan


    # ---------- GALFIT NC ----------
    for k in [
        "GAL_XC", "GAL_XC_ERR",
        "GAL_YC", "GAL_YC_ERR",
        "GAL_MAG", "GAL_MAG_ERR",
        "GAL_RE_ARCSEC", "GAL_RE_ERR_ARCSEC",
        "GAL_N", "GAL_N_ERR",
        "GAL_BA", "GAL_BA_ERR",
        "GAL_PA", "GAL_PA_ERR",
        "GAL_SKY", "GAL_SKY_ERR",
        "GAL_CHISQ"
    ]:
        row[k] = np.nan

    for k in ["GAL_NUMERR", "GAL_ERROR"]:
        row[k] = 0



    # ---------- GALFIT CV ----------
    for k in [
        "GAL_CXC", "GAL_CXC_ERR",
        "GAL_CYC", "GAL_CYC_ERR",
        "GAL_CMAG", "GAL_CMAG_ERR",
        "GAL_CRE_ARCSEC", "GAL_CRE_ERR_ARCSEC",
        "GAL_CN", "GAL_CN_ERR",
        "GAL_CBA", "GAL_CBA_ERR",
        "GAL_CPA", "GAL_CPA_ERR",
        "GAL_CSKY", "GAL_CSKY_ERR",
        "GAL_CCHISQ"
    ]:
        row[k] = np.nan

    for k in ["GAL_CNUMERR", "GAL_CERROR"]:
        row[k] = 0
        
    #row["GAL_CV_OK"] = False

    
    return row


def init_csgr_row_defaults(row):
    """
    Initialize optional CS-gr output columns.

    These stay NaN/False/empty when no *-CS-gr.fits image exists.
    """

    defaults = {
        # file/status
        "CSGR_EXISTS": False,
        "CSGR_FITS": "",
        "CSGR_SEC": np.nan,
        "CSGR_PHOT_OK": False,
        "CSGR_HAPY_MORPH_OK": False,
        "CSGR_HAPY_MORPH_FLAG": 0,
        "CSGR_CONTSCL": np.nan,

        # image2 / ellipse photometry scalars
        "CSGR_SKYSTD_ADU": np.nan,
        "CSGR_SKYMED_ADU": np.nan,
        "CSGR_SKYSTD_PHYS": np.nan,
        "CSGR_M20": np.nan,
        "CSGR_ASYM": np.nan,
        "CSGR_ASYM_ERR": np.nan,
        "CSGR_SCALE_ADU_CGS": np.nan,

        # HAPY morphology, CS-gr image side only
        "CSGR_HAPY_NPIX": np.nan,
        "CSGR_HAPY_FILLFRAC": np.nan,
        "CSGR_HAPY_SNP_ALL": np.nan,
        "CSGR_HAPY_SNP_DET": np.nan,
        "CSGR_H_HAPY_GINI_THRESHOLD": np.nan,
        "CSGR_R_HAPY_GINI_THRESHOLD": np.nan,        
        "CSGR_HAPY_GINI": np.nan,
        "CSGR_HAPY_M20": np.nan,
        "CSGR_HAPY_ASYM": np.nan,
        "CSGR_HAPY_ASYM_ERR": np.nan,
        "CSGR_HAPY_MTOT": np.nan,
        "CSGR_HAPY_M20SUM": np.nan,
        "CSGR_HAPY_FLUX_SEG": np.nan,
        "CSGR_HAPY_MTOT2": np.nan,
        "CSGR_HAPY_RMOM_ARCSEC": np.nan,
    }

    for key, val in defaults.items():
        row.setdefault(key, val)

    csgr_profile_columns = [
        'CSGR_H25_ARCSEC', 'CSGR_H25_PIX', 'CSGR_H30R24_FLUX_CGS', 'CSGR_H30R24_FLUX_CGS_ERR',
        'CSGR_H50_ARCSEC', 'CSGR_H50_PIX', 'CSGR_H75_ARCSEC', 'CSGR_H75_PIX',
        'CSGR_H_C30_R24', 'CSGR_H_C30_R24_ERR', 'CSGR_H_EXPFIT_I0', 'CSGR_H_EXPFIT_K',
        'CSGR_H_EXPFIT_OK', 'CSGR_H_EXPFIT_RE_ARCSEC', 'CSGR_H_ISO17E18_ARCSEC',
        'CSGR_H_ISO17E18_ARCSEC_ERR', 'CSGR_H_ISO17E18_FLUX_CGS', 'CSGR_H_ISO17E18_FLUX_CGS_ERR',
        'CSGR_H_ISO5E17_ARCSEC', 'CSGR_H_ISO5E17_ARCSEC_ERR', 'CSGR_H_ISO5E17_FLUX_CGS',
        'CSGR_H_ISO5E17_FLUX_CGS_ERR', 'CSGR_H_LOGFIT_INTERCEPT', 'CSGR_H_LOGFIT_OK',
        'CSGR_H_LOGFIT_RE_ARCSEC', 'CSGR_H_LOGFIT_SLOPE', 'CSGR_H_MAXDET_ARCSEC',
        'CSGR_H_MAXDET_PIX', 'CSGR_H_NDET_RUNS', 'CSGR_H_PETRO_CON', 'CSGR_H_PETRO_FLUX',
        'CSGR_H_PETRO_FLUX_CGS', 'CSGR_H_PETRO_FLUX_CGS_ERR', 'CSGR_H_PETRO_MAG',
        'CSGR_H_PETRO_OK', 'CSGR_H_PETRO_R50_ARCSEC', 'CSGR_H_PETRO_R90_ARCSEC',
        'CSGR_H_PETRO_RAD_ARCSEC', 'CSGR_H_PROFILE_LONGRUN', 'CSGR_H_PROFILE_MASKFRAC_MAX',
        'CSGR_H_PROFILE_NGOOD', 'CSGR_H_PROFILE_OK', 'CSGR_H_R24_FLUX_CGS',
        'CSGR_H_R24_FLUX_CGS_ERR', 'CSGR_H_R95_R24_ARCSEC', 'CSGR_H_SNR_TRUNC_ARCSEC',
        'CSGR_H_TOT_FLUX_CGS', 'CSGR_H_TOT_FLUX_CGS_ERR', 'CSGR_R24_ARCSEC',
        'CSGR_R24_ARCSEC_ERR', 'CSGR_R24_FLUX_CGS', 'CSGR_R24_FLUX_CGS_ERR',
        'CSGR_R24_MAG', 'CSGR_R24_MAG_ERR', 'CSGR_R24_VEGA_ARCSEC',
        'CSGR_R24_VEGA_ARCSEC_ERR', 'CSGR_R24_VEGA_MAG', 'CSGR_R24_VEGA_MAG_ERR',
        'CSGR_R25P5_ARCSEC', 'CSGR_R25P5_ARCSEC_ERR', 'CSGR_R25P5_MAG',
        'CSGR_R25P5_MAG_ERR', 'CSGR_R25_ARCSEC', 'CSGR_R25_ISO_ARCSEC',
        'CSGR_R25_ISO_ARCSEC_ERR', 'CSGR_R25_ISO_MAG', 'CSGR_R25_ISO_MAG_ERR',
        'CSGR_R25_PIX', 'CSGR_R25_VEGA_ARCSEC', 'CSGR_R25_VEGA_ARCSEC_ERR',
        'CSGR_R25_VEGA_MAG', 'CSGR_R25_VEGA_MAG_ERR', 'CSGR_R30R24_FLUX_CGS',
        'CSGR_R30R24_FLUX_CGS_ERR', 'CSGR_R50_ARCSEC', 'CSGR_R50_PIX',
        'CSGR_R75_ARCSEC', 'CSGR_R75_PIX', 'CSGR_R_C30', 'CSGR_R_C30_ERR',
        'CSGR_R_EXPFIT_I0', 'CSGR_R_EXPFIT_K', 'CSGR_R_EXPFIT_OK',
        'CSGR_R_EXPFIT_RE_ARCSEC', 'CSGR_R_LOGFIT_INTERCEPT', 'CSGR_R_LOGFIT_OK',
        'CSGR_R_LOGFIT_RE_ARCSEC', 'CSGR_R_LOGFIT_SLOPE', 'CSGR_R_PETRO_CON',
        'CSGR_R_PETRO_FLUX', 'CSGR_R_PETRO_FLUX_CGS', 'CSGR_R_PETRO_FLUX_CGS_ERR',
        'CSGR_R_PETRO_MAG', 'CSGR_R_PETRO_OK', 'CSGR_R_PETRO_R50_ARCSEC',
        'CSGR_R_PETRO_R90_ARCSEC', 'CSGR_R_PETRO_RAD_ARCSEC',
        'CSGR_R_PROFILE_MASKFRAC_MAX', 'CSGR_R_PROFILE_NGOOD', 'CSGR_R_PROFILE_OK',
        'CSGR_R_SNR_TRUNC_ARCSEC', 'CSGR_R_TOT_FLUX_CGS_SNR',
        'CSGR_R_TOT_FLUX_CGS_SNR_ERR', 'CSGR_R_TOT_MAG_SNR',
    ]

    for col in csgr_profile_columns:
        if col.endswith("_OK") or col.endswith("_LONGRUN"):
            row.setdefault(col, False)
        elif col.endswith("_NGOOD") or col.endswith("_NDET_RUNS"):
            row.setdefault(col, -1)
        else:
            row.setdefault(col, np.nan)

    return row

def initialize_sourcecatalog_moments(row, prefix1="R", prefix2="H"):
    import numpy as np

    cols = [
        "MOMENTS_OK",
        "COV_XX", "COV_YY", "COV_XY",
        "SEMIMAJOR_SIGMA", "SEMIMINOR_SIGMA", "AREA",
        "ORIENTATION",
        "ECCENTRICITY", "ELONGATION","GINI",
        "KRON_RADIUS", "KRON_FLUX"
    ]

    # initialize for safe stacking
    for prefix in [prefix1, prefix2]:
        for col in cols:
            row[f"{prefix}_SC_{col}"] = False if col == "MOMENTS_OK" else np.nan
        row[f"{prefix}_SC_UNITS"] = ""
    return row

def add_sourcecatalog_moments(row, e, prefix1="R", prefix2="H", pixel_scale=None):
    import numpy as np

    pixscale = np.nan if pixel_scale is None else float(pixel_scale)

    def safe_float(x):
        try:
            return float(np.asarray(x).squeeze())
        except Exception:
            return np.nan

    def fill_from_sourcecat(prefix, cat, index):
        try:
            src = cat[index]
            cov = src.covariance

            sma_pix = safe_float(src.semimajor_sigma.value)
            smb_pix = safe_float(src.semiminor_sigma.value)
            kron_pix = safe_float(src.kron_radius.value)

            if np.isfinite(pixscale) and pixscale > 0:
                row[f"{prefix}_SC_COV_XX"] = safe_float(cov[0, 0]) * pixscale**2
                row[f"{prefix}_SC_COV_YY"] = safe_float(cov[1, 1]) * pixscale**2
                row[f"{prefix}_SC_COV_XY"] = safe_float(cov[0, 1]) * pixscale**2
                row[f"{prefix}_SC_SEMIMAJOR_SIGMA"] = sma_pix * pixscale
                row[f"{prefix}_SC_SEMIMINOR_SIGMA"] = smb_pix * pixscale
                row[f"{prefix}_SC_AREA"] = np.pi * sma_pix * smb_pix * pixscale**2
                row[f"{prefix}_SC_KRON_RADIUS"] = kron_pix * pixscale
                row[f"{prefix}_SC_UNITS"] = "arcsec"
            else:
                row[f"{prefix}_SC_COV_XX"] = safe_float(cov[0, 0])
                row[f"{prefix}_SC_COV_YY"] = safe_float(cov[1, 1])
                row[f"{prefix}_SC_COV_XY"] = safe_float(cov[0, 1])
                row[f"{prefix}_SC_SEMIMAJOR_SIGMA"] = sma_pix
                row[f"{prefix}_SC_SEMIMINOR_SIGMA"] = smb_pix
                row[f"{prefix}_SC_AREA"] = np.pi * sma_pix * smb_pix
                row[f"{prefix}_SC_KRON_RADIUS"] = kron_pix
                row[f"{prefix}_SC_UNITS"] = "pix"

            row[f"{prefix}_SC_KRON_FLUX"] = safe_float(src.kron_flux)
            row[f"{prefix}_SC_ORIENTATION"] = safe_float(src.orientation.value)
            row[f"{prefix}_SC_ECCENTRICITY"] = safe_float(src.eccentricity.value)
            row[f"{prefix}_SC_ELONGATION"] = safe_float(src.elongation.value)
            row[f"{prefix}_SC_GINI"] = safe_float(src.gini)
            row[f"{prefix}_SC_MOMENTS_OK"] = True

        except Exception as err:
            print(f"WARNING: could not extract SourceCatalog moments for {prefix}: {err}")

    if hasattr(e, "cat") and hasattr(e, "objectIndex"):
        fill_from_sourcecat(prefix1, e.cat, e.objectIndex)

    if hasattr(e, "cat2") and hasattr(e, "objectIndex2"):
        fill_from_sourcecat(prefix2, e.cat2, e.objectIndex2)

    return row

def _hcl_col(prefix, name):
    """
    Build a clump-analysis column name.

    Examples
    --------
    _hcl_col("HCL_", "NCLUMP")    -> "HCL_NCLUMP"
    _hcl_col("HCL_GR_", "NCLUMP") -> "HCL_GR_NCLUMP"
    """

    if prefix is None:
        prefix = ""

    if len(prefix) > 0 and not prefix.endswith("_"):
        prefix = prefix + "_"

    return prefix + name

def _pix_to_arcsec(value, pixel_scale):
    """
    Convert pixels to arcsec.

    Parameters
    ----------
    value : float
        Value in pixels.

    pixel_scale : float
        Pixel scale in arcsec/pixel.
    """

    if pixel_scale is None or not np.isfinite(pixel_scale):
        return np.nan

    if value is None:
        return np.nan

    try:
        value = float(value)
    except Exception:
        return np.nan

    if not np.isfinite(value):
        return np.nan

    return value * pixel_scale


def _pixarea_to_arcsec2(value, pixel_scale):
    """
    Convert pixel area to arcsec^2.
    """

    if pixel_scale is None or not np.isfinite(pixel_scale):
        return np.nan

    if value is None:
        return np.nan

    try:
        value = float(value)
    except Exception:
        return np.nan

    if not np.isfinite(value):
        return np.nan

    return value * pixel_scale**2

def initialize_hapy_clumps(row, prefix="HCL_"):
    """
    Initialize H-alpha clump-analysis columns in the run_analysis output row.

    Parameters
    ----------
    row : dict-like
        Output row used by run_analysis.py.

    prefix : str
        Prefix for column names. Use, e.g., "HCL_" for CS-ZP and
        "HCL_GR_" for CS-gr.
    """

    # ------------------------------------------------------------
    # Status / provenance
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "OK")] = False
    row[_hcl_col(prefix, "STATUS")] = "not_run"
    row[_hcl_col(prefix, "INPUT_IMAGE")] = ""

    # ------------------------------------------------------------
    # Detection counts
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "NCLUMP")] = -1
    row[_hcl_col(prefix, "NCLUMP_GOOD")] = -1
    row[_hcl_col(prefix, "NPEAK")] = -1
    row[_hcl_col(prefix, "NPEAK_IN_CLUMPS")] = -1
    row[_hcl_col(prefix, "NPOINTSRC")] = -1

    # ------------------------------------------------------------
    # Threshold / background
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "THRESHOLD")] = np.nan
    row[_hcl_col(prefix, "BKG_MEAN")] = np.nan
    row[_hcl_col(prefix, "BKG_MEDIAN")] = np.nan
    row[_hcl_col(prefix, "BKG_RMS")] = np.nan

    # ------------------------------------------------------------
    # Clump flux and area summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "FLUX_FOOTPRINT")] = np.nan
    row[_hcl_col(prefix, "FLUX_SUM")] = np.nan
    row[_hcl_col(prefix, "FLUX_FRAC")] = np.nan

    row[_hcl_col(prefix, "BRIGHT_FLUX")] = np.nan
    row[_hcl_col(prefix, "BRIGHT_FRAC")] = np.nan
    row[_hcl_col(prefix, "TOP2_FRAC")] = np.nan
    row[_hcl_col(prefix, "TOP3_FRAC")] = np.nan
    row[_hcl_col(prefix, "PARAM_NUCLEAR_RADIUS_FWHM")] = np.nan
    # ------------------------------------------------------------
    # Nuclear clump summaries
    # Nuclear = clump centroid within one measured H-alpha FWHM
    # of the adopted galaxy center.
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "HAS_NUCLEAR")] = False
    row[_hcl_col(prefix, "NNUCLEAR")] = -1
    row[_hcl_col(prefix, "NUCLEAR_FLUX_FRAC")] = np.nan

    # ------------------------------------------------------------
    # Clump flux dominance
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "BRIGHT_TO_SECOND_FLUX")] = np.nan

    # ------------------------------------------------------------
    # Pixel scale / angular units
    # ------------------------------------------------------------
    #row[_hcl_col(prefix, "PIXSCALE")] = np.nan  # arcsec / pixel

    # ------------------------------------------------------------
    # Clump area summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "AREA_ARCSEC2")] = np.nan
    row[_hcl_col(prefix, "FOOTPRINT_AREA_ARCSEC2")] = np.nan
    row[_hcl_col(prefix, "AREA_FRAC")] = np.nan



    # ------------------------------------------------------------
    # Peak summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "PEAK_MAX")] = np.nan
    row[_hcl_col(prefix, "PEAK_SUM")] = np.nan
    row[_hcl_col(prefix, "POINTSRC_FLUX_SUM")] = np.nan

    # ------------------------------------------------------------
    # Clump centroid summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "XCEN_FLUXWT")] = np.nan
    row[_hcl_col(prefix, "YCEN_FLUXWT")] = np.nan
    row[_hcl_col(prefix, "XCEN_BRIGHT")] = np.nan
    row[_hcl_col(prefix, "YCEN_BRIGHT")] = np.nan



    # ------------------------------------------------------------
    # Clump centroid offsets relative to galaxy center
    # Stored in arcsec, not pixels.
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "RMIN_ARCSEC")] = np.nan
    row[_hcl_col(prefix, "RMAX_ARCSEC")] = np.nan
    row[_hcl_col(prefix, "RMED_ARCSEC")] = np.nan
    row[_hcl_col(prefix, "RMEAN_ARCSEC")] = np.nan
    row[_hcl_col(prefix, "RFLUXWT_ARCSEC")] = np.nan
    row[_hcl_col(prefix, "RBRIGHT_ARCSEC")] = np.nan
    row[_hcl_col(prefix, "DX_BRIGHT_ARCSEC")] = np.nan
    row[_hcl_col(prefix, "DY_BRIGHT_ARCSEC")] = np.nan
    
    # ------------------------------------------------------------
    # Saved product paths
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "CATALOG")] = ""
    row[_hcl_col(prefix, "SEGMAP")] = ""
    row[_hcl_col(prefix, "PEAKS")] = ""
    row[_hcl_col(prefix, "POINTSRC")] = ""
    row[_hcl_col(prefix, "DIAG")] = ""

    # ------------------------------------------------------------
    # Important configuration parameters
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "PARAM_NSIGMA")] = np.nan
    row[_hcl_col(prefix, "PARAM_NPIXELS")] = -1
    row[_hcl_col(prefix, "PARAM_DEBLEND")] = False
    row[_hcl_col(prefix, "PARAM_NLEVELS")] = -1
    row[_hcl_col(prefix, "PARAM_CONTRAST")] = np.nan
    row[_hcl_col(prefix, "PARAM_MODE")] = ""

    row[_hcl_col(prefix, "PARAM_FIND_PEAKS")] = False
    row[_hcl_col(prefix, "PARAM_PEAK_BOX_SIZE")] = -1
    row[_hcl_col(prefix, "PARAM_PEAK_MIN_SEP")] = -1

    return row




def write_hapy_clumps(
    row,
    clump_result,
    prefix="HCL_",
    input_image="CS-ZP",
    config=None,
    pixel_scale=None,
    failed=False,
):
    """
    Write H-alpha clump-analysis outputs into the run_analysis output row.

    Parameters
    ----------
    row : dict-like
        Output row used by run_analysis.py.

    clump_result : ClumpAnalysisResult or None
        Result returned by analyze_halpha_clumps / measure_halpha_clumps.

    prefix : str
        Prefix for column names. Use, e.g., "HCL_" for CS-ZP and
        "HCL_GR_" for CS-gr.

    input_image : str
        Label for the H-alpha image used, e.g. "CS-ZP" or "CS-gr".

    config : ClumpDetectionConfig or None
        Configuration used for the clump analysis. If provided, key
        parameters are stored in the row.

    failed : bool
        Set True if clump analysis failed.
    """

    # Make sure columns exist even if this is called directly.
    #initialize_hapy_clumps(row, prefix=prefix)

    row[_hcl_col(prefix, "INPUT_IMAGE")] = input_image
    #row[_hcl_col(prefix, "FAILED")] = bool(failed)

    if failed or clump_result is None:
        row[_hcl_col(prefix, "OK")] = False
        row[_hcl_col(prefix, "STATUS")] = "failed"
        return row

    s = clump_result.summary

    row[_hcl_col(prefix, "OK")] = True
    row[_hcl_col(prefix, "STATUS")] = "ok"

    # ------------------------------------------------------------
    # Detection counts
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "NCLUMP")] = int(getattr(s, "n_clumps", -1))
    row[_hcl_col(prefix, "NCLUMP_GOOD")] = int(getattr(s, "n_clumps_good", -1))
    row[_hcl_col(prefix, "NPEAK")] = int(getattr(s, "n_peaks", -1))
    row[_hcl_col(prefix, "NPEAK_IN_CLUMPS")] = int(getattr(s, "n_peaks_in_clumps", -1))
    row[_hcl_col(prefix, "NPOINTSRC")] = int(getattr(s, "n_point_sources", -1))

    # ------------------------------------------------------------
    # Threshold / background
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "THRESHOLD")] = getattr(s, "threshold", np.nan)
    row[_hcl_col(prefix, "BKG_MEAN")] = getattr(s, "background_mean", np.nan)
    row[_hcl_col(prefix, "BKG_MEDIAN")] = getattr(s, "background_median", np.nan)
    row[_hcl_col(prefix, "BKG_RMS")] = getattr(s, "background_rms", np.nan)

    # ------------------------------------------------------------
    # Flux and area summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "FLUX_FOOTPRINT")] = getattr(
        s, "total_halpha_flux_footprint", np.nan
    )
    row[_hcl_col(prefix, "FLUX_SUM")] = getattr(s, "total_clump_flux", np.nan)
    row[_hcl_col(prefix, "FLUX_FRAC")] = getattr(s, "clump_flux_fraction", np.nan)

    row[_hcl_col(prefix, "BRIGHT_FLUX")] = getattr(s, "brightest_clump_flux", np.nan)
    row[_hcl_col(prefix, "BRIGHT_FRAC")] = getattr(
        s, "brightest_clump_fraction", np.nan
    )
    row[_hcl_col(prefix, "TOP2_FRAC")] = getattr(s, "top2_clump_fraction", np.nan)
    row[_hcl_col(prefix, "TOP3_FRAC")] = getattr(s, "top3_clump_fraction", np.nan)

    # ------------------------------------------------------------
    # Nuclear clump summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "HAS_NUCLEAR")] = bool(
        getattr(s, "has_nuclear", False)
        )

    row[_hcl_col(prefix, "NNUCLEAR")] = int(
        getattr(s, "n_nuclear", -1)
        )

    row[_hcl_col(prefix, "NUCLEAR_FLUX_FRAC")] = getattr(
        s, "nuclear_flux_frac", np.nan
        )
    row[_hcl_col(prefix, "PARAM_NUCLEAR_RADIUS_FWHM")] = getattr(
        config, "nuclear_radius_fwhm", np.nan
        )
    # ------------------------------------------------------------
    # Clump flux dominance
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "BRIGHT_TO_SECOND_FLUX")] = getattr(
        s, "bright_to_second_flux", np.nan
        )

    # ------------------------------------------------------------
    # Area summaries
    # Store angular areas in arcsec^2.
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "AREA_ARCSEC2")] = _pixarea_to_arcsec2(
        getattr(s, "clump_area_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "FOOTPRINT_AREA_ARCSEC2")] = _pixarea_to_arcsec2(
        getattr(s, "footprint_area_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "AREA_FRAC")] = getattr(s, "clump_area_fraction", np.nan)

    # ------------------------------------------------------------
    # Peak summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "PEAK_MAX")] = getattr(s, "peak_max", np.nan)
    row[_hcl_col(prefix, "PEAK_SUM")] = getattr(s, "peak_sum", np.nan)
    row[_hcl_col(prefix, "POINTSRC_FLUX_SUM")] = getattr(
        s, "point_source_flux_sum", np.nan
    )

    # ------------------------------------------------------------
    # Clump centroid summaries
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "XCEN_FLUXWT")] = getattr(
        s, "flux_weighted_xcentroid", np.nan
    )
    row[_hcl_col(prefix, "YCEN_FLUXWT")] = getattr(
        s, "flux_weighted_ycentroid", np.nan
    )
    row[_hcl_col(prefix, "XCEN_BRIGHT")] = getattr(s, "brightest_xcentroid", np.nan)
    row[_hcl_col(prefix, "YCEN_BRIGHT")] = getattr(s, "brightest_ycentroid", np.nan)



    # ------------------------------------------------------------
    # Clump centroid offset summaries
    # Stored in arcsec, not pixels.
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "RMIN_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "rmin_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "RMAX_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "rmax_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "RMED_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "rmed_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "RMEAN_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "rmean_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "RFLUXWT_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "rfluxwt_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "RBRIGHT_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "rbright_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "DX_BRIGHT_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "dx_bright_pix", np.nan),
        pixel_scale,
        )

    row[_hcl_col(prefix, "DY_BRIGHT_ARCSEC")] = _pix_to_arcsec(
        getattr(s, "dy_bright_pix", np.nan),
        pixel_scale,
        )


    # ------------------------------------------------------------
    # Saved product paths
    # ------------------------------------------------------------
    row[_hcl_col(prefix, "CATALOG")] = getattr(s, "catalog_path", "")
    row[_hcl_col(prefix, "SEGMAP")] = getattr(s, "segmentation_path", "")
    row[_hcl_col(prefix, "PEAKS")] = getattr(s, "peaks_path", "")
    row[_hcl_col(prefix, "POINTSRC")] = getattr(s, "point_sources_path", "")
    row[_hcl_col(prefix, "DIAG")] = getattr(s, "diagnostic_path", "")

    # ------------------------------------------------------------
    # Important configuration parameters
    # ------------------------------------------------------------
    if config is not None:
        row[_hcl_col(prefix, "PARAM_NSIGMA")] = getattr(config, "nsigma", np.nan)
        row[_hcl_col(prefix, "PARAM_NPIXELS")] = int(getattr(config, "npixels", -1))
        row[_hcl_col(prefix, "PARAM_DEBLEND")] = bool(getattr(config, "deblend", False))
        row[_hcl_col(prefix, "PARAM_NLEVELS")] = int(getattr(config, "nlevels", -1))
        row[_hcl_col(prefix, "PARAM_CONTRAST")] = getattr(config, "contrast", np.nan)
        row[_hcl_col(prefix, "PARAM_MODE")] = getattr(config, "mode", "")

        row[_hcl_col(prefix, "PARAM_FIND_PEAKS")] = bool(
            getattr(config, "find_peaks", False)
        )
        row[_hcl_col(prefix, "PARAM_PEAK_BOX_SIZE")] = int(
            getattr(config, "peak_box_size", -1)
        )

        peak_min_sep = getattr(config, "peak_min_separation", None)
        if peak_min_sep is None:
            peak_min_sep = -1
        row[_hcl_col(prefix, "PARAM_PEAK_MIN_SEP")] = int(peak_min_sep)

    return row



def _print_psf_image(p,logger):
    if logger is not None:
        logger.info(f"PSF image = {str(p)}")
            
    else:
        print(f"PSF image = {str(p)}")
    
def pick_psf_path_and_source(args, params,logger=None):
    psf_image = getattr(args, "psf_image", None)
    if psf_image:
        p = Path(psf_image)
        _print_psf_image(p,logger)
        return (str(p), "cli") if p.exists() else (None, "cli_missing")

    psfdir = getattr(args, "psf_dir", None) or getattr(args, "psfdir", None)
    parent = params.get("parent_rimage", "")
    if psfdir and parent:
        name = str(parent).replace(".fits", "-psf.fits")
        p = Path(psfdir) / name
        _print_psf_image(p,logger)
        return (str(p), "psf_dir") if p.exists() else (None, "psf_dir_missing")
    _print_psf_image(None, logger)
    return None, ""

def check_table(results_table):
    # checking table
    print()
    print("CHECKING OUTPUT TABLE")
    print()
    from astropy.table import Table
    t = Table.read(results_table, format="ascii.ecsv")
    print(t.dtype)
    print(t[0])

    t.write("tmp.ecsv", format="ascii.ecsv", overwrite=True)
    t2 = Table.read("tmp.ecsv", format="ascii.ecsv")
    assert t.colnames == t2.colnames

    #for c in ["ELLIP_MASKED_FRACTION", "SM_R_FLAG", "SM_H_FLAG"]:
    for c in ["SM_R_FLAG", "SM_H_FLAG"]:
        print(c, t[c].dtype, t[c][0])

def valid_file(path):
    p = Path(path)
    return p.is_file() and p.stat().st_size > 0


def load_gaia_table(params, args, logger):
    scheme = (params.get("scheme") or "").lower()
    if args.no_gaia or scheme == "archive":
        logger.info("Gaia masking disabled")
        return None

    parent_rimage = params.get("parent_rimage")
    logger.info(f"in load_gaia_table, parent_rimage: {parent_rimage}")
    if not parent_rimage or not args.gaia_dir:
        logger.info(f"in load_gaia_table, not loading gaia table: {parent_rimage}")
        return None

    path = Path(args.gaia_dir) / parent_rimage.replace(".fits", "-gaia.fits")
    logger.info(f"in load_gaia_table, gaia table: {path}")    
    if path.exists():
        logger.info(f"Using Gaia catalog: {path}")
        return Table.read(path)

    logger.warning(f"Gaia catalog not found: {path}")
    return None

def build_mask_for_cutout(
    *,
    cutdir,
    tag,
    r_fits,
    params,
    args,
    logger,
    row,
    results_path,
    ell0_params,
    sma_arcsec,
    pixscale,
    ra,
    dec,
    use_gaia,
):
    """
    Build a mask for the current cutout, write it to <tag>-mask.fits,
    update relevant row fields, and return:

        mask, mask_fits, row
    """
    if args.seconfig is None:
        raise ValueError("--sex-config must be set when --make-mask or --force-mask is used")

    row["STAGE"] = "mask"
    logger.info("STAGE: mask")
    t0 = time.perf_counter()

    mask_out = cutdir / f"{tag}-mask.fits"
    logger.info(f"Writing mask to {mask_out}")

    # --- Convert to pixels ---
    sma_pix = sma_arcsec / pixscale

    # convert CCW from N angle to photutils CCW from +x
    galaxy_ellipse = ell0_params

    # -- look for gaia table
    gaia_table = load_gaia_table(params, args, logger)

    engine = MaskEngine(
        image_fits=r_fits,
        sepath=args.sepath,
        config=args.seconfig,
        threshold=args.sethreshold,
        snr=args.sesnr,
        minarea=args.seminarea,
        add_gaia_stars=use_gaia,
    )

    max_fwhm = max(row["R_FWHM_PSF"], row["H_FWHM_PSF"])
    if max_fwhm > 2.5:
        #try using se fwhm
        try:
            logger.info(f" PSF FWHM is big ({max_fwhm} arcsec), checking SE value instead")
            max_fwhm = max(row["R_FWHM_SE"], row["H_FWHM_SE"])
        except:
            logger.info(f" PSF FWHM is big ({max_fwhm} arcsec), no SE value. setting max to 1.5 arcsec")
            max_fwhm = 1.5
    gaia_min_radius_arcsec = 4 * max_fwhm
    logger.info(f"Gaia min radius (arcsec) = {gaia_min_radius_arcsec}")

    gaia_min_radius_deg = gaia_min_radius_arcsec / 3600.0

    mask = engine.build_initial_mask(
        galaxy_ellipse=galaxy_ellipse,
        progress_callback=_progress_cb,
        grow_size=int(args.grow_size),
        grow_iterations=int(args.grow_iterations),
        gaia_table=gaia_table,
        gaia_min_radius=gaia_min_radius_deg,
    )

    engine.write_mask(mask_out)
    mask_fits = mask_out

    if use_gaia and gaia_table is not None:
        bright_flag, dist_arcsec, maskrad_arcsec, bright_mag = galaxy_overlaps_bright_star(
            ra,
            dec,
            gaia_table,
            mag_limit=10,
            radius_col="radius",
            min_radius_arcsec=gaia_min_radius_arcsec,
        )

        row["BRIGHT_STAR_FLAG"] = bright_flag
        row["BRIGHT_STAR_DIST_ARCSEC"] = dist_arcsec
        row["BRIGHT_STAR_MASKRAD_ARCSEC"] = maskrad_arcsec
        row["BRIGHT_STAR_MAG"] = bright_mag
    
    row["MASK_OK"] = True
    row["MASK_FITS"] = str(mask_fits)
    row["MASK_SEC"] = time.perf_counter() - t0


    return mask, mask_fits, row

def archive_existing_mask(mask_path: Path) -> Path:
    """
    Rename existing <tag>-mask.fits to the next available
    <tag>-mask-N.fits and return the archived path.
    """
    stem = mask_path.stem
    suffix = mask_path.suffix
    parent = mask_path.parent

    i = 1
    while True:
        archived = parent / f"{stem}-{i}{suffix}"
        if not archived.exists():
            mask_path.rename(archived)
            return archived
        i += 1

        
def main():

    p = argparse.ArgumentParser(
    description="Run headless analysis on one galaxy cutout directory"
    )

    # ============================================================
    # Required input
    # ============================================================
    g_req = p.add_argument_group("Required Input")

    g_req.add_argument(
        "--cutout-dir",
        required=True,
        help="Cutout directory (e.g. survey_run/cutouts/<tag>)"
        )

    # ============================================================
    # Main pipeline controls
    # ============================================================
    g_main = p.add_argument_group("Pipeline Controls")

    g_main.add_argument("--make-mask", action="store_true",
                        help="Build/write mask before photometry/galfit")
    g_main.add_argument("--csgr", action="store_true",
                        help="Compute photometry on csgr images")
    g_main.add_argument("--statmorph", action="store_true",
                        help="Compute statmorph structural parameters")
    g_main.add_argument("--galfit", action="store_true",
                        help="Run GALFIT after photometry")
    g_main.add_argument("--convflag", action="store_true", default=False,
                        help="Run GALFIT convolution stage (requires PSF)")
    g_main.add_argument("--no-diagnostic-plots", action="store_true",
                        help="Don't write diagnostic plot (R image + mask + ellipses)")


    # ============================================================
    # H-alpha clump analysis
    # ============================================================
    g_clumps = p.add_argument_group("H-alpha Clump Analysis")

    g_clumps.add_argument(
        "--clumps",
        action="store_true",
        help=(
            "Run H-alpha clump analysis inside the central R-band "
            "segmentation footprint."
        ),
    )

    g_clumps.add_argument(
        "--clump-nsigma",
        type=float,
        default=3.0,
        help=(
            "Detection threshold for H-alpha clumps in units of the local "
            "background RMS. Default is 3.0."
        ),
    )
    g_clumps.add_argument(
        "--clump-npixels",
        type=int,
        default=5,
        help=(
            "Minimum number of connected pixels required for an H-alpha "
            "clump detection. Default is 5."
        ),
    )

    g_clumps.add_argument(
        "--no-clump-deblend",
        action="store_true",
        help="Disable photutils deblending for H-alpha clump detections.",
    )
    g_clumps.add_argument(
        "--clump-nlevels",
        type=int,
        default=64,
        help="Number of deblending levels. Default is 32.",
    )
    g_clumps.add_argument(
        "--clump-deblend-mode",
        default="linear",
        choices=["exponential", "linear", "sinh"],
        help=(
            "Multi-threshold spacing mode for photutils.deblend_sources. "
            "Default is linear."
        ),
    )    
    g_clumps.add_argument(
        "--clump-contrast",
        type=float,
        default=0.001,
        help=(
            "Deblending contrast for H-alpha clumps. Smaller values split "
            "more aggressively. Default is 0.001."
        ),
    )

    g_clumps.add_argument(
        "--no-clump-peaks",
        action="store_true",
        help="Disable local peak finding within the H-alpha clump image.",
    )
    g_clumps.add_argument(
        "--clump-peak-box-size",
        type=int,
        default=5,
        help="Box size used by photutils.find_peaks. Default is 5.",
    )
    g_clumps.add_argument(
        "--clump-peak-min-separation",
        type=int,
        default=None,
        help=(
            "Minimum separation in pixels between local H-alpha peaks. "
            "Default is None."
        ),
    )

    g_clumps.add_argument(
        "--clump-background-grow-radius",
        type=float,
        default=10,
        help=(
            "Grow radius passed to calculate_background_photutils when "
            "estimating the clump detection threshold. Default is 10."
        ),
    )
    g_clumps.add_argument(
        "--clump-background-npixels",
        type=int,
        default=10,
        help=(
            "Minimum object size passed to calculate_background_photutils "
            "for background masking. Default is 10."
        ),
    )
    g_clumps.add_argument(
        "--clump-background-mask-nsigma",
        type=float,
        default=2.0,
        help=(
            "Object-mask threshold passed to calculate_background_photutils. "
            "Default is 2.0."
        ),
    )
    g_clumps.add_argument(
        "--clump-background-clip-sigma",
        type=float,
        default=3.0,
        help=(
            "Sigma-clipping threshold used by calculate_background_photutils. "
            "Default is 3.0."
        ),
    )
    g_clumps.add_argument(
        "--clump-background-clip-maxiters",
        type=int,
        default=5,
        help=(
            "Maximum sigma-clipping iterations used by "
            "calculate_background_photutils. Default is 5."
        ),
    )

    g_clumps.add_argument(
        "--clump-min-flux",
        type=float,
        default=None,
        help=(
            "Optional minimum clump flux used for summary statistics. "
            "The full clump catalog is still saved. Default is None."
        ),
    )
    g_clumps.add_argument(
        "--clump-min-area",
        type=int,
        default=None,
        help=(
            "Optional minimum clump area in pixels used for summary statistics. "
            "The full clump catalog is still saved. Default is None."
        ),
    )

    g_clumps.add_argument(
        "--no-clump-diagnostic",
        action="store_true",
        help="Do not save the H-alpha clump diagnostic image.",
    )
    g_clumps.add_argument(
        "--clump-diagnostic-format",
        default="png",
        choices=["png", "pdf"],
        help="File format for the H-alpha clump diagnostic image. Default is png.",
    )
    g_clumps.add_argument(
        "--clump-diagnostic-percent",
        type=float,
        default=99.5,
        help=(
            "Percentile stretch for the H-alpha clump diagnostic image. "
            "Default is 99.5."
        ),
    )
    g_clumps.add_argument(
        "--clump-plot-kron-apertures",
        action="store_true",
        help="Overlay SourceCatalog Kron apertures on the clump diagnostic image.",
    )

    g_clumps.add_argument(
        "--clump-point-sources",
        action="store_true",
        help=(
            "Also run optional point-source-like H-alpha knot detection. "
            "This is not used as the primary clump definition."
        ),
    )
    g_clumps.add_argument(
        "--clump-point-source-method",
        default="dao",
        choices=["dao", "iraf"],
        help="Point-source finder to use if --clump-point-sources is set.",
    )
    g_clumps.add_argument(
        "--clump-point-source-fwhm",
        type=float,
        default=None,
        help=(
            "FWHM in pixels for optional point-source detection. Required "
            "if --clump-point-sources is set."
        ),
    )
    g_clumps.add_argument(
        "--clump-point-source-threshold-nsigma",
        type=float,
        default=5.0,
        help=(
            "Point-source detection threshold in units of background RMS. "
            "Default is 5.0."
        ),
    )


    # ============================================================
    # Masking (SExtractor + Gaia)
    # ============================================================
    g_mask = p.add_argument_group("Masking Options")

    g_mask.add_argument("--sepath", default="sex",
                        help="Path to SExtractor executable")
    g_mask.add_argument("--seconfig", default=_default_sex_config(),
                        help="SExtractor config file path")
    g_mask.add_argument("--sethreshold", type=float, default=0.005,
                        help="SExtractor detection/deblend threshold. Default is 0.005.")
    g_mask.add_argument("--sesnr", type=float, default=5.0,
                        help="SExtractor SNR threshold. Default is 5.")
    g_mask.add_argument("--seminarea", type=int, default=5,
                        help="SExtractor minimum object area. Default is 5.")
    g_mask.add_argument("--grow-size", type=int, default=7,
                        help="Grow size in mask expansion. Default is 7.")
    g_mask.add_argument("--grow-iterations", type=int, default=4,
                        help="Number of mask-growth iterations. Default is 4.")
    g_mask.add_argument("--gaia-dir", default="gaia_catalogs",
                        help="Directory containing precomputed Gaia catalogs (default: gaia_catalogs)")
    g_mask.add_argument("--no-gaia", action="store_true",
                        help="Disable Gaia star masking")
    g_mask.add_argument("--force-mask", action="store_true",
                        help="Rebuild mask even if an existing mask file is present")


    # ============================================================
    # GALFIT
    # ============================================================
    g_gal = p.add_argument_group("GALFIT Options")

    g_gal.add_argument("--psf-dir", dest="psf_dir", default=None,
                       help="Directory containing PSF images (derived from parent_rimage)")
    g_gal.add_argument("--psf-image", dest="psf_image", default=None,
                       help="Explicit PSF image path (overrides --psf-dir)")
    g_gal.add_argument("--psf-oversampling", type=int, default=2,
                       help="PSF oversampling factor")
    g_gal.add_argument("--ncomp", type=int, default=1, choices=[1, 2],
                       help="Number of GALFIT components")
    g_gal.add_argument("--magzp", type=float, default=None,
                       help="Override PHOTZP passed to GALFIT")
    g_gal.add_argument("--sky", type=float, default=0.0,
                       help="Fixed sky value for GALFIT")

    # ============================================================
    # Logging
    # ============================================================
    g_log = p.add_argument_group("Logging")

    g_log.add_argument("--log-level", default="INFO",
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       help="Logging verbosity")
    g_log.add_argument("--log-to-console", action="store_true",
                       help="Also print logs to stdout")
    g_log.add_argument("--log-dir", default=None,
                       help="Optional directory for logs (default: cutout directory)")


    # ============================================================
    # Input overrides (advanced use)
    # ============================================================
    g_io = p.add_argument_group("Input Overrides (Advanced)")

    g_io.add_argument("--r", dest="r_fits", default=None,
                      help="Override R-band FITS path")
    g_io.add_argument("--cs", dest="cs_fits", default=None,
                      help="Override continuum-sub FITS path")
    g_io.add_argument("--mask", dest="mask_fits", default=None,
                      help="Override mask FITS path")
    g_io.add_argument("--sigma-image", dest="sigma_image", default=None,
                      help="Override sigma/RMS image")




    # ============================================================
    # Metadata overrides (advanced)
    # ============================================================
    g_meta = p.add_argument_group("Metadata Overrides (Advanced)")

    g_meta.add_argument("--image2-filter", dest="image2_filter", default=None,
                        help="Override image2 filter (otherwise from metadata.json)")
    g_meta.add_argument("--filter-ratio", dest="filter_ratio", type=float, default=None,
                        help="Override FLTRATIO for image2 flux calibration")
    g_meta.add_argument("--objra", type=float, default=None,
                    help="Override object RA (deg). If set, overrides metadata.json RA for WCS-derived center.")
    g_meta.add_argument("--objdec", type=float, default=None,
                    help="Override object Dec (deg). If set, overrides metadata.json Dec for WCS-derived center.")
    g_meta.add_argument("--sma-arcsec", type=float, default=None,
                        help="Override metadata ellipse SMA (arcsec)")
    g_meta.add_argument("--ba", type=float, default=None,
                        help="Override metadata ellipse b/a")
    g_meta.add_argument("--pa-deg", type=float, default=None,
                        help="Override metadata PA_DEG (CCW from North)")
    g_meta.add_argument("--fixcenter", action="store_true",
                        help="Hold ellipse center fixed during photometry")
    
    args = p.parse_args()
    
 

    cutdir = Path(args.cutout_dir)
    if not cutdir.exists():
        raise FileNotFoundError(f"Cutout directory not found: {cutdir}")

    tag = cutdir.name
    root = tag

    prefix = tag


    params_path = cutdir / "metadata.json"

    if not params_path.exists():
        raise RuntimeError(
            f"metadata.json not found for root {root}. "
            "Cutouts may be outdated or improperly generated."
        )
    
    params = json.loads(params_path.read_text())
    
    results_path = cutdir / f"{tag}-results.ecsv"

    print(f"DEBUG: cutdir={cutdir},tag={tag}, root={root}")
    # --- initialize logger
    logger, log_path = init_cutout_logger(
        tag=tag,
        cutdir=cutdir,
        level=args.loglevel if hasattr(args, "loglevel") else "INFO",
        log_to_console=args.log_to_console,
        log_dir=args.log_dir,  # usually None
        )

    
    # Auto-detect common filenames if not provided.
    # Adjust these glob patterns to match your exact suffix conventions.
    #r_fits = args.r_fits or _pick_one(str(cutdir/ f"{tag}*-R.fits")) or _pick_one(str(cutdir / f"{tag}*-r.fits"))
    r_fits = (
        args.r_fits
        or (str(cutdir / params.get("r_fits")) if params.get("r_fits") else None)
        or _pick_one(str(cutdir / f"{tag}*-R.fits"))
        or _pick_one(str(cutdir / f"{tag}*-r.fits"))
        )
    if r_fits is None:
        raise FileNotFoundError(f"Could not find R-band FITS in: {cutdir}")

    #cs_fits = args.cs_fits or _pick_one(str(cutdir / f"{tag}*-CS-ZP.fits")) or _pick_one(str(cutdir / f"{tag}*-cs.fits"))
    cs_fits = (
        args.cs_fits
        or (str(cutdir / params.get("cs_fits")) if params.get("cs_fits") else None)
        or _pick_one(str(cutdir / f"{tag}*-CS-ZP.fits"))
        or _pick_one(str(cutdir / f"{tag}*-CS.fits"))
        or _pick_one(str(cutdir / f"{tag}*-cs.fits"))
        )
    # why are we looking for a mask when we are suppose to make one?
    #mask_fits = args.mask_fits or _pick_one(str(cutdir / f"{tag}*-mask.fits"))

    
    sigma_image = args.sigma_image or _pick_one(str(cutdir / f"{tag}*-sigma.fits")) or _pick_one(str(cutdir / f"{tag}*-rms.fits"))
    psf_image = args.psf_image or _pick_one(str(cutdir / f"{tag}*-psf.fits"))


            
    row = initialize_result_row()

    row = initialize_sourcecatalog_moments(row, prefix1="R", prefix2="H")
    
    row = init_csgr_row_defaults(row)

    row = initialize_sourcecatalog_moments(row, prefix1="CSGR_R", prefix2="CSGR_H")    

    row = initialize_hapy_clumps(row, prefix="HCL_")
    row = initialize_hapy_clumps(row, prefix="CSGR_HCL_")    
    # look for CS-gr image and log it if found
    csgr_fits = (
        str(cutdir / params.get("csgr_fits")) if params.get("csgr_fits") else None
    ) or _pick_one(str(cutdir / f"{tag}*-CS-gr.fits"))


    if csgr_fits:
        logger.info(f"Found CS-gr image: {csgr_fits}")
        row["CSGR_EXISTS"] = True
        row["CSGR_FITS"] = Path(csgr_fits).name
        row["CSGR_CONTSCL"] = fits.getheader(csgr_fits).get("CONTSCL", np.nan)
    else:
        logger.info("No CS-gr image found")
    
    
    from datetime import datetime
    row["RUN_DATE"] = datetime.utcnow().strftime("%Y-%m-%d")

    from importlib.metadata import version
    row["HAPY_VERSION"] = version("hapy")
    
    row["CUTDIR"] = str(cutdir)
    row["TAG"] = Path(root).name
    row["STAGE"] = "init"
    row["STATUS"] = "running"
    row["MASK_SEC"] = 0.0
    row["PHOT_SEC"] = 0.0
    row["HAPY_MORPH_SEC"] = 0.0
    row["PROFILES_SEC"] = 0.0     
    row["CSGR_SEC"] = 0.0
    row["GAL_NC_SEC"] = 0.0
    row["GAL_CV_SEC"] = 0.0    
    row["SM_SEC"] = 0.0    
    row["TOTAL_SEC"] = 0.0

    row["R_FITS"] = r_fits
    row["CS_FITS"] = cs_fits
    row["SIGMA_FITS"] = sigma_image    

    
    # pixel scale
    #pixscale = args.pixscale
    #if pixscale is None:
    pixscale = get_pixel_scale_from_filename(r_fits)
    print(f"DEBUG: pixscale = {pixscale:.4e}")
    row["PIXSCALE"] = round(float(pixscale),4)
    # --- Load cutout image for WCS + shape ---
    r_data, hdr = fits.getdata(r_fits, header=True)
    ny, nx = r_data.shape
    wcs = WCS(hdr)
    
    magzp = args.magzp if args.magzp is not None else float(hdr.get("PHOTZP", 25.0))
    
    row["CUTOUT_XSIZE"] = nx
    row["CUTOUT_YSIZE"] = ny
    # Default center = image center
    xc = nx / 2.0
    yc = ny / 2.0


    t0_total = time.perf_counter()

    #######################################################
    # --- check for valid input ellipse
    #######################################################    
    if ellipse_missing(params):
        print("\nGetting initial ellipse estimate from photutils...\n")
        # check for gaia catalog

        # 
        # get gaia mask
        gaia_table = load_gaia_table(params, args, logger)
        if gaia_table is None:
            brightstar, star_xpix, star_ypix = get_gaia_stars(r_fits)
        else:
            brightstar = gaia_table
            from astropy.coordinates import SkyCoord
            starcoord = SkyCoord(
                brightstar["ra"],
                brightstar["dec"],
                frame="icrs",
                unit="deg",
                )

            star_xpix, star_ypix = wcs.world_to_pixel(starcoord)

        mask_array = np.zeros_like(r_data,  dtype=np.int32)
        # could add a radius_scale_factor that depends on image FWHM
        # (FWHM-1.5)
        mask_array, gaia_mask = make_gaia_mask(mask_array,star_xpix,star_ypix,pixscale/3600.,gaia_table=gaia_table,radius_scale_factor=1)

        # convert to boolean mask
        gaia_mask = gaia_mask > 0

        # get ellipse from photutils
        ell = infer_ellipse_from_r_cutout(r_data=r_data, user_mask=gaia_mask)
        if ell is not None:
            # if agc has a valid radius and BA, then keep?
            #print("DEBUG: original radius = ",params["sma_arcsec"])
            radius_scale_factor = 1.2
            #print("DEBUG: new radius = ",ell.sma_pix * pixscale)
            
            # we tested this scale factor to see what looks reasonable
            # settled on radius_scale_factor = 1.2
            
            params["sma_arcsec"] = float(ell.sma_pix * pixscale * radius_scale_factor)
            params["ba"] = float(ell.ba)
            params["pa_deg"] = float(photutils_theta_to_pa_ccw_north(ell.theta_deg))
            params["ell0_source"] = "quick_photutils"
            params["ell0_ok"] = True
        else:
            params["ell0_source"] = "quick_photutils_failed"
            params["ell0_ok"] = False

        # write back (development mode)
        params_path.write_text(json.dumps(params, indent=2))

    # --- START OF PASTE
    # --- update row with other info from metadata.json
    scheme = (params.get("scheme") or "").lower()

    use_gaia = (not args.no_gaia) and (scheme != "archive")    

    row["TELESCOPE"] = params.get("telescope", "")
    row["DATEOBS"]   = params.get("dateobs", "")
    row["POINTING"]  = params.get("pointing", "")
    row["SCHEME"]    = params.get("scheme", "")
    row["PARENT_RIMAGE"]  = params.get("parent_rimage", "")
    row["PARENT_HIMAGE"] = params.get("parent_haimage", "")
    row["X_PARENT"]  = params.get("x_parent", np.nan)
    row["Y_PARENT"] = params.get("y_parent", np.nan)
    row["CSZP_SOURCE"] = params.get("cszp_source", "")
    row["CSZP_LOCAL_SKY"] = bool(params.get("cszp_local_sky", False))

    # These are more survey/workflow-specific; leave archive rows at initialized defaults
    if scheme != "archive":
        row["HFILTER"] = params.get("hafilter")
        row["CUTOUT_SCALE"] = params.get("cutout_scale")
        row["FILTER_CORRECTION"] = params.get("filter_correction")

        # TODO get this information from ZP ratio

        # Prefer FILTER_RATIO calculated from PHOTZP values.
        # Fall back to metadata only if PHOTZP values are unavailable/non-finite.
        rheader = fits.getheader(r_fits)
        csheader = fits.getheader(cs_fits)
        zp_r = rheader.get("PHOTZP", np.nan)
        zp_h = csheader.get("PHOTZP", np.nan)

        filter_ratio = zp_scale_r_to_ha(zp_h, zp_r, logger=logger)

        if np.isfinite(filter_ratio):
            logger.info(
                "Calculated FILTER_RATIO from PHOTZP values: zp_h=%.4f zp_r=%.4f filter_ratio=%.6g",
                zp_h, zp_r, filter_ratio,
            )
        else:
            filter_ratio = params.get("filter_ratio", np.nan)

            try:
                filter_ratio = float(filter_ratio)
            except Exception:
                filter_ratio = np.nan

            if np.isfinite(filter_ratio):
                logger.warning(
                    "Using FILTER_RATIO from metadata because PHOTZP ratio could not be calculated: %.6g",
                    filter_ratio,
                )
            else:
                logger.warning(
                    "FILTER_RATIO could not be calculated from PHOTZP and is missing from metadata; "
                    "physical flux calibration will be NaN."
                )

        # row["FILTER_RATIO"] = filter_ratio


        # filter_ratio = params.get("filter_ratio", None)
        
        # if filter_ratio is None:
        #     filter_ratio = np.nan
        #     logger.warning("FLTRATIO missing from metadata; physical flux calibration will be NaN.")
        # row["FILTER_RATIO"] = filter_ratio

        # --- Construct the name of the psf image
        psf_path, psf_source = pick_psf_path_and_source(args, params,logger=logger)
        row["PSF_FITS"] = str(psf_path) if psf_path else ""
        row["PSF_OK"] = bool(psf_path)
        row["PSF_SOURCE"] = psf_source
        psf_ok = row["PSF_OK"]
    else:
        psf_path = None
        psf_ok = False
        
    # --- Get ellipse parameters ---
    sma_arcsec = float(params["sma_arcsec"])
    ba = float(params["ba"])
    pa_deg = float(params["pa_deg"])  # CCW from N, from input catalog

    # Try WCS-based centering using stored RA/DEC
    objid = params.get("objid", Path(root).name)
    ra = args.objra if args.objra is not None else params.get("ra")
    dec = args.objdec if args.objdec is not None else params.get("dec")
    row["RA"] = ra
    row["DEC"] = dec
    row["OBJID"] = objid
    row["REDSHIFT"] = params.get("redshift")
    row["VR"] = params.get("vr")

    # Optional floats: only set if present
    val = params.get("rimage_fwhm_se_arcsec")
    if val is not None:
        row["R_FWHM_SE"] = float(val)

    val = params.get("rimage_fwhm_psf_arcsec")
    if val is not None:
        row["R_FWHM_PSF"] = float(val)

    val = params.get("himage_fwhm_se_arcsec")
    if val is not None:
        row["H_FWHM_SE"] = float(val)

    val = params.get("himage_fwhm_psf_arcsec")
    if val is not None:
        row["H_FWHM_PSF"] = float(val)

    # Optional strings
    val = params.get("rfilter_name")
    if val is not None:
        row["RFILTER_FILENAME"] = val

    val = params.get("hafilter_name")
    if val is not None:
        row["HFILTER_FILENAME"] = val

    # Optional filter metadata
    val = params.get("rfilter_center_A")
    if val is not None:
        row["RFILTER_CENTER"] = float(val)

    val = params.get("rfilter_width_A")
    if val is not None:
        row["RFILTER_WIDTH"] = float(val)

    val = params.get("hafilter_center_A")
    if val is not None:
        row["HFILTER_CENTER"] = float(val)

    val = params.get("hafilter_width_A")
    if val is not None:
        row["HFILTER_WIDTH"] = float(val)

    # --- END OF METADATA TRANSFER  

    if ra is not None and dec is not None:
        try:
            # note this is only as good as the header wcs, aka not that good for archive sample!
            xw, yw = wcs.world_to_pixel_values(float(ra), float(dec))
            if np.isfinite(xw) and np.isfinite(yw):
                xc, yc = float(xw), float(yw)
        except Exception:
            pass
    print(f"DEBUG: xc={xc:.1f},yc={yc:.1f},\nra={ra:.6f},dec={dec:.6f}")
    # -- let user input ellipse geometry that is different from what is in metadata.json
    if args.sma_arcsec is not None and args.ba is not None and args.pa_deg is not None:
        sma_arcsec = float(args.sma_arcsec)
        ba = float(args.ba)
        pa_deg = float(args.pa_deg)


    objid = row.get("OBJID", "") or ""

    # VFID is the survey ID prefix for Virgo objects: "VFID####"
    m = re.match(r"^(VFID\d+)", objid)
    if m:
        row["VFID"] = m.group(1)

    # Optional: split out the NED name part from "VFID####-NEDname"
    if "-" in objid:
        row["GALNAME"] = objid.split("-", 1)[1]

    # --- Store the initial ellipse used as input to masking
    row["ELL0_SMA_ARCSEC"] = sma_arcsec
    row["ELL0_BA"] = ba
    row["ELL0_PA_DEG"] = pa_deg # CCW from N, where N is +y axis, and W is +x axis
    row["ELL0_XC"] = xc
    row["ELL0_YC"] = yc
    row["ELL0_SOURCE"] = "metadata.json" if params_path.exists() else "args"
    print(f"DEBUG: xc={xc:.1f},yc={yc:.1f},\nra={ra:.6f},dec={dec:.6f}")
    ell0_params = EllipseParams(
        xc = xc,
        yc = yc,
        ba = ba,
        sma_pix = sma_arcsec/pixscale,
        theta_deg = photutils_theta_to_pa_ccw_north(pa_deg)
        )

    # --- store coverage information about initial ellipse
    r_cov = ellipse_image_coverage(r_data, ell0_params)

    cs_data, hhdr = fits.getdata(cs_fits, header=True)    
    h_cov = ellipse_image_coverage(cs_data, ell0_params)

    row["CUTOUT_ELL0_NPIX_TOTAL_R"] = r_cov["npix_total"]
    row["CUTOUT_ELL0_NPIX_TOTAL_H"] = h_cov["npix_total"]

    row["CUTOUT_ELL0_NPIX_GOOD_R"] = r_cov["npix_good"]
    row["CUTOUT_ELL0_NPIX_GOOD_H"] = h_cov["npix_good"]

    row["CUTOUT_ELL0_MISSING_FRAC_R"] = r_cov["missing_frac"]
    row["CUTOUT_ELL0_MISSING_FRAC_H"] = h_cov["missing_frac"]

    row["CUTOUT_ELL0_MISSING_FRAC_MAX"] = np.nanmax([
        r_cov["missing_frac"],
        h_cov["missing_frac"],
        ])

    # cutout_map = {
    #     "cutout_ell0_missing_frac_r": "CUTOUT_ELL0_MISSING_FRAC_R",
    #     "cutout_ell0_missing_frac_h": "CUTOUT_ELL0_MISSING_FRAC_H",
    #     "cutout_ell0_missing_frac_max": "CUTOUT_ELL0_MISSING_FRAC_MAX",
    #     "cutout_ell0_npix_total_r": "CUTOUT_ELL0_NPIX_TOTAL_R",
    #     "cutout_ell0_npix_total_h": "CUTOUT_ELL0_NPIX_TOTAL_H",
    #     "cutout_ell0_npix_onimage_r": "CUTOUT_ELL0_NPIX_ONIMAGE_R",
    #     "cutout_ell0_npix_onimage_h": "CUTOUT_ELL0_NPIX_ONIMAGE_H",
    #     "cutout_ell0_npix_good_r": "CUTOUT_ELL0_NPIX_GOOD_R",
    #     "cutout_ell0_npix_good_h": "CUTOUT_ELL0_NPIX_GOOD_H",
    # }

    # for pkey, rkey in cutout_map.items():
    #     val = params.get(pkey)
    #     if val is not None:
    #         row[rkey] = float(val)
        

    ################################################################
    # Mask block
    ################################################################

    manual_mask = cutdir / f"{tag}-mask-manual.fits"

    mask_fits = (
        args.mask_fits
        or (cutdir / params["mask_fits"] if params.get("mask_fits") else None)
        or (manual_mask if manual_mask.exists() else None)
        or _pick_one(str(cutdir / f"{tag}*-mask.fits"))
    )

    mask_fits = Path(mask_fits) if mask_fits is not None else None

    mask_out = cutdir / f"{tag}-mask.fits"

    if args.force_mask:
        logger.info("Force-mask enabled: rebuilding mask")

        if mask_out.exists():
            archived = archive_existing_mask(mask_out)
            logger.info(f"Archived existing mask to {archived}")

        mask, mask_fits, row = build_mask_for_cutout(
            cutdir=cutdir,
            tag=tag,
            r_fits=r_fits,
            params=params,
            args=args,
            logger=logger,
            row=row,
            results_path=results_path,
            ell0_params=ell0_params,
            sma_arcsec=sma_arcsec,
            pixscale=pixscale,
            ra=ra,
            dec=dec,
            use_gaia=use_gaia,
        )
        row["MASK_SOURCE"] = "rebuilt"

    elif mask_fits is not None and mask_fits.exists():
        logger.info(f"Using existing mask: {mask_fits}")
        mask = fits.getdata(mask_fits)
        row["MASK_OK"] = True
        row["MASK_FITS"] = str(mask_fits)

        if args.mask_fits:
            row["MASK_SOURCE"] = "cli"
        elif params.get("mask_fits"):
            row["MASK_SOURCE"] = "params"
        elif mask_fits == manual_mask:
            row["MASK_SOURCE"] = "manual"
        else:
            row["MASK_SOURCE"] = "auto"

    elif args.make_mask:
        logger.info("No existing mask found; building mask")

        mask, mask_fits, row = build_mask_for_cutout(
            cutdir=cutdir,
            tag=tag,
            r_fits=r_fits,
            params=params,
            args=args,
            logger=logger,
            row=row,
            results_path=results_path,
            ell0_params=ell0_params,
            sma_arcsec=sma_arcsec,
            pixscale=pixscale,
            ra=ra,
            dec=dec,
            use_gaia=use_gaia,
        )
        row["MASK_SOURCE"] = "built"

    else:
        logger.info("No mask provided and mask building not requested")
        mask = None
        row["MASK_OK"] = False
        row["MASK_SOURCE"] = "none"

    if mask is not None:
        res = ellipse_mask_fraction(mask, ell0_params)
        print("DEBUG ellipse_mask_fraction type:", type(res))
        print("DEBUG ellipse_mask_fraction value:", res)
        #import inspect
        #print("DEBUG ellipse_mask_fraction module:", ellipse_mask_fraction.__module__)
        #print("DEBUG ellipse_mask_fraction file:", inspect.getsourcefile(ellipse_mask_fraction))
        #print("DEBUG ellipse_mask_fraction source:")
        #print(inspect.getsource(ellipse_mask_fraction))
        row["ELL0_MASKFRAC"] = res.frac_masked
        row["ELL0_MASK_WARN"] = res.frac_masked > 0.3
        row["ELL0_NMASKPIX"] = res.n_masked
        row["ELL0_NTOTPIX"] = res.n_total

    write_result_row_ecsv(results_path, row)

    
 


    ################################################################
    # phot block
    ################################################################
    row["STAGE"] = "phot"
    logger.info("STAGE: phot")
    
    t0 = time.perf_counter()

    hafilter = row["HFILTER"]
    if args.image2_filter is not None:
        hafilter = args.image2_filter
        row["HFILTER"] = hafilter

    # filter ratio should come from photometric zps
    
    filter_ratio = row["FILTER_RATIO"]

    print(f"DEBUG: xc={xc:.1f},yc={yc:.1f},\nra={ra:.6f},dec={dec:.6f}")
    e = run_ellipse_photometry(
        r_fits=r_fits,
        cs_fits=cs_fits,
        mask_fits=mask_fits,
        image2_filter=hafilter,
        filter_ratio=filter_ratio,
        objra=ra,
        objdec=dec,
        fixcenter=args.fixcenter,
        logger=logger,
        #run_statmorph=args.statmorph,
        #write_prefix=prefix,
    )

    # ---- photometry summary (scalar-only; arrays stay in the photometry table files) ----
    row["PHOT_SEC"] = _scalar(time.perf_counter() - t0)
    row["PHOT_OK"] = True

    # ---- core ellipse / detection-derived quantities ----
    FIELDS = [
        ("ELLIP_XCENTROID", "xcenter_fit"),
        ("ELLIP_YCENTROID", "ycenter_fit"),
        ("ELLIP_SMA_ARCSEC", "sma_fit"),
        ("ELLIP_SMA_PIX", "sma_fit"),
        ("ELLIP_B_ARCSEC", "b"),
        ("ELLIP_EPS", "eps_fit"),
        ("ELLIP_THETA_RAD", "pa_fit"),
        ("R_ELLIP_GINI", "gini"),
        ("H_ELLIP_GINI", "gini2"),        
        ("ELLIP_SOURCE_SUM", "source_sum"),
        ("ELLIP_SEGMENT_FLUX", "photutils_segment_flux"),
        ("ELLIP_SEGMENT_MAG","photutils_segment_mag"),
        ("R_SKYSTD_ADU", "sky_noise"),
        ("R_SKYMED_ADU", "sky"),
        ("R_SKYSTD_PHYS", "im1_skynoise"),
        ("R_M20", "M20_1"),
        ("R_ASYM", "asym"),
        ("R_ASYM_ERR", "asym_err"),
        ("R_SCALE_ADU_CGS", "uconversion1"),
        ("AREA_GUESS_ELLIPSE_PIX", "area_guess_ellipse_pix"),
        ("AREA_GUESS_ELLIPSE_UNMASKED_PIX", "area_guess_ellipse_unmasked_pix"),
        ("MASKFRAC_GUESS_ELLIPSE", "maskfrac_guess_ellipse"),
        ("R_PROFILE_NONCENTRAL_PEAK","r_profile_noncentral_peak"),
        ("R_PROFILE_PEAK_BIN","r_profile_peak_bin"),
        ("R_PROFILE_PEAK_SMA","r_profile_peak_sma"),
        ]


    for outk, attr in FIELDS:
        v = getattr(e, attr, None)

        sv = _scalar(v)
        
        if sv is not None:
            if "_ARCSEC" in outk:
                #print(f"Converting to arcsec for ",outk, attr)
                sv = sv * pixscale
            row[outk] = sv  # leave as np.nan if missing/array/etc.
        
    # convert theta to PA and store
    if row["ELLIP_THETA_RAD"] is not None:
        # convert to PA in deg and store result
        # ELLIP_THETA_RAD measured from +x axis
        phot_theta_deg = (np.degrees(float(row["ELLIP_THETA_RAD"])) % 180.0)
        phot_pa_deg = photutils_theta_to_pa_ccw_north(phot_theta_deg)  # inverse of your adapter
        row["ELLIP_PA_DEG"] = phot_pa_deg
        
    # add photutils B/A
    row["ELLIP_BA"] = 1. - float(row["ELLIP_EPS"])
    
    row["ELLIP_CENTER_METHOD"] = e.center_method_best
    # JSON field (stable schema)
    # mf = getattr(e, "masked_fraction", None)
    # if mf is not None:
    #     try:
    #         row["ELLIP_MASKED_FRACTION"] = json.dumps(mf)
    #     except TypeError:
    #         # if mf contains numpy types etc.
    #         row["ELLIP_MASKED_FRACTION"] = json.dumps(mf, default=_scalar)
        

    # If image2 exists (e.g., continuum-sub / HA), capture analogous scalars
    FIELDS2 = [
                ("H_SKYSTD_ADU", "sky_noise2"),
                ("H_SKYMED_ADU", "sky2"),
                ("H_SKYSTD_PHYS", "im2_skynoise"),
                ("H_M20", "M20_2"),
                ("H_ASYM", "asym2"),
                ("H_ASYM_ERR", "asym2_err"),
                ("H_SCALE_ADU_CGS","uconversion2"),
                
            ]

    if getattr(e, "image2", None) is not None:
        for outk, attr in FIELDS2:
            v = getattr(e, attr, None)
            sv = _scalar(v)
            if sv is not None:
                row[outk] = sv



    add_sourcecatalog_moments(row, e, prefix1="R", prefix2="H",pixel_scale = float(pixscale))
    
    try:
        phot_xc = float(row["ELLIP_XCENTROID"])
        phot_yc = float(row["ELLIP_YCENTROID"])
        phot_sma_pix = float(row["ELLIP_SMA_ARCSEC"])/pixscale
        phot_ba = 1.0 - float(row["ELLIP_EPS"])
        #phot_pa_deg = (np.degrees(float(row["ELLIP_THETA_RAD"])) % 180.0)
        #phot_pa_deg = photutils_theta_to_pa_ccw_north(theta_phot_deg)  # inverse of your adapter


        # save this for later
        
        dx = phot_xc - float(row["ELL0_XC"])
        dy = phot_yc - float(row["ELL0_YC"])
        dc = float(np.hypot(dx, dy))

        dba = abs(phot_ba - float(row["ELL0_BA"]))

        dpa = phot_pa_deg - float(row["ELL0_PA_DEG"])

        row["ELL_DC_PX"] = dc
        row["ELL_DBA"] = float(dba)
        row["ELL_DPA_DEG"] = float(dpa)

        # size ratio: prefer arcsec if pixscale known, else pixels
        if "pixscale" in locals() and pixscale:
            phot_sma_arcsec = phot_sma_pix * float(pixscale)
            row["ELL_SMA_RATIO"] = phot_sma_arcsec / float(row["ELL0_SMA_ARCSEC"])
        else:
            # fallback: compare in pixels if you have ELL0_SMA_ARCSEC only -> skip ratio
            row["ELL_SMA_RATIO"] = np.nan

        sma_ratio = row["ELL_SMA_RATIO"]
        row["ELL_MISMATCH"] = bool(
            (dc > 10.0) or (dba > 0.2) or (np.abs(dpa) > 10.0)
            #or (np.isfinite(sma_ratio) and ((sma_ratio < 0.5) or (sma_ratio > 2.0)))
            )
    except Exception:
        # if any missing keys, just don't set mismatch fields
        pass
    phot_xc = float(row["ELLIP_XCENTROID"])
    phot_yc = float(row["ELLIP_YCENTROID"])
    phot_sma_pix = float(row["ELLIP_SMA_ARCSEC"])/pixscale
    phot_ba = 1.0 - float(row["ELLIP_EPS"])
    # ELLIP_THETA_RAD measured from +x axis
    phot_theta_deg = (np.degrees(float(row["ELLIP_THETA_RAD"])) % 180.0)
    ellphot_params = EllipseParams(
        xc = phot_xc,
        yc = phot_yc,
        sma_pix = phot_sma_pix,
        ba = phot_ba,
        theta_deg = phot_theta_deg
    )

    
    outfile = cutdir / f"{tag}-diagnostic.png"
    outfile2 = cutdir / f"{tag}-seg-diagnostic.png"
    se_seg_fits = cutdir / f"{tag}-R-segmentation.fits"
    phot_seg_fits = cutdir / f"{tag}-R-phot-segmentation.fits"    

    plot_segmentation_diagnostic(
        r_fits=str(r_fits),
        se_seg_fits=se_seg_fits,
        mask_fits=str(mask_fits),
        phot_seg_fits=phot_seg_fits,
        e0=ell0_params,
        eph=ellphot_params,
        outfile=str(outfile2),
        row=row,
        )
    plot_mask_ellipse_diagnostic(
        r_fits=str(r_fits),
        mask_fits=str(mask_fits),
        e0=ell0_params,
        eph=ellphot_params,
        outfile=str(outfile),
        row=row,
        )

    # calculate hapy gini
    t0 = time.perf_counter()

    # could pass in conversions from ADU to cgs
    e.run_hapy_morphology(nsigma=3, RSKY_SIGMA_FLOOR_SB=4e-16, HSKY_SIGMA_FLOOR_SB=5e-17)


    print(f"DEBUG: ASYM {e.R_HAPY_ASYM}, center={e.R_HAPY_ASYM_CENTER}")
    row["R_HAPY_NPIX"] = e.R_HAPY_NPIX     
    row["H_HAPY_NPIX"] = e.H_HAPY_NPIX
    row["H_HAPY_FILLFRAC"] = e.H_HAPY_FILLFRAC
    row["R_HAPY_SNP_ALL"] = e.R_HAPY_SNP_ALL
    row["H_HAPY_SNP_ALL"] = e.H_HAPY_SNP_ALL
    row["H_HAPY_SNP_DET"] = e.H_HAPY_SNP_DET    
    row["H_HAPY_GINI_THRESHOLD"] = e.ha_gini_threshold
    row["R_HAPY_GINI_THRESHOLD"] = e.r_gini_threshold        
    row["R_HAPY_XC"] = e.R_HAPY_XC
    row["R_HAPY_YC"] = e.R_HAPY_YC
    row["R_HAPY_GINI"] = e.R_HAPY_GINI
    row["H_HAPY_GINI"] = e.H_HAPY_GINI
    row["R_HAPY_M20"] = e.R_HAPY_M20
    row["H_HAPY_M20"] = e.H_HAPY_M20
    row["R_HAPY_ASYM"] = e.R_HAPY_ASYM
    row["H_HAPY_ASYM"] = e.H_HAPY_ASYM
    row["R_HAPY_ASYM_ERR"] = e.R_HAPY_ASYM_ERR
    row["H_HAPY_ASYM_ERR"] = e.H_HAPY_ASYM_ERR
    row["R_HAPY_ASYM_XC"] = e.R_HAPY_ASYM_CENTER[1]
    row["R_HAPY_ASYM_YC"] = e.R_HAPY_ASYM_CENTER[0]
    row["R_HAPY_MTOT"] = e.R_HAPY_MTOT
    row["H_HAPY_MTOT"] = e.H_HAPY_MTOT
    row["R_HAPY_M20SUM"] = e.R_HAPY_M20SUM
    row["H_HAPY_M20SUM"] = e.H_HAPY_M20SUM
    
    row["H_HAPY_FLUX_SEG"] = e.H_HAPY_FLUX_SEG
    row["H_HAPY_MTOT2"] = e.H_HAPY_MTOT2
    row["H_HAPY_RMOM_ARCSEC"] = e.H_HAPY_RMOM * pixscale
    
    row["R_HAPY_FLUX_SEG"] = e.R_HAPY_FLUX_SEG
    row["R_HAPY_MTOT2"] = e.R_HAPY_MTOT2
    row["R_HAPY_RMOM_ARCSEC"] = e.R_HAPY_RMOM * pixscale
    # set flags
    row["HAPY_MORPH_OK"] = e.HAPY_MORPH_OK       
    row["HAPY_MORPH_FLAG"] = e.HAPY_MORPH_FLAG   
    # TODONE add HAPY SNP for H and R

    row["HAPY_MORPH_SEC"] = _scalar(time.perf_counter() - t0)        
    

    # ---- RUN HAPY CLUMP ANALYSIS  ----------- #

    if args.clumps:
        print("DEBUG: trying to run clump analysis")
        from hapy.ellipse.clumps import ClumpDetectionConfig

        clump_config = ClumpDetectionConfig(
            nsigma=args.clump_nsigma,
            npixels=args.clump_npixels,

            deblend=not args.no_clump_deblend,
            nlevels=args.clump_nlevels,
            contrast=args.clump_contrast,
            mode=args.clump_deblend_mode,
            
            find_peaks=not args.no_clump_peaks,
            peak_box_size=args.clump_peak_box_size,
            peak_min_separation=args.clump_peak_min_separation,

            background_grow_radius=args.clump_background_grow_radius,
            background_npixels=args.clump_background_npixels,
            background_mask_nsigma=args.clump_background_mask_nsigma,
            background_clip_sigma=args.clump_background_clip_sigma,
            background_clip_maxiters=args.clump_background_clip_maxiters,

            min_flux=args.clump_min_flux,
            min_area=args.clump_min_area,

            save_diagnostic=not args.no_clump_diagnostic,
            diagnostic_format=args.clump_diagnostic_format,
            diagnostic_percent=args.clump_diagnostic_percent,
            plot_kron_apertures=args.clump_plot_kron_apertures,

            find_point_sources=args.clump_point_sources,
            point_source_method=args.clump_point_source_method,
            point_source_fwhm=args.clump_point_source_fwhm,
            point_source_threshold_nsigma=args.clump_point_source_threshold_nsigma,
        )

        # temporarily moving this out of try/except for debugging
        clump_result = e.measure_halpha_clumps(
            config=clump_config,
            output_dir=cutdir,
            basename=tag,
            overwrite=True,
            update_results=False,   # keep this False while testing
        )
        

        try:
            clump_result = e.measure_halpha_clumps(
                config=clump_config,
                output_dir=cutdir,
                basename=tag,
                overwrite=True,
                update_results=False,   # keep this False while testing
            )

            write_hapy_clumps(
                row,
                clump_result,
                prefix="HCL_",
                input_image="CS-ZP",
                config=clump_config,
                pixel_scale=pixscale,
                failed=False,
                )
            
            logger.info(
                "H-alpha clump analysis complete: NCLUMP=%d, NPEAK=%d",
                clump_result.summary.n_clumps,
                getattr(clump_result.summary, "n_peaks", 0),
            )

        except Exception as err:
            logger.warning("H-alpha clump analysis failed: %s", err)
            print("H-alpha clump analysis failed: %s", err)

            if hasattr(e, "results") and isinstance(e.results, dict):
                e.results["HCL_OK"] = False
                e.results["HCL_STATUS"] = "ok"
            write_hapy_clumps(
                row,
                None,
                prefix="HCL_",
                input_image="CS-ZP",
                config=clump_config,
                failed=True,
                )
                


    # ---- FIT PROFILES!  ----------- #    

    if valid_file(e.photfile) and valid_file(e.photfile2):
        row["PHOT_OK"] = True
        t0 = time.perf_counter()    
        rtab = Table.read(e.photfile)
        hatab = Table.read(e.photfile2)

        profile_results = summarize_dual_profiles(
            rtab=rtab,
            hatab=hatab,
            r_magzp=magzp,
        )

        row.update(profile_results)
        row["PROFILES_SEC"] = _scalar(time.perf_counter() - t0)        

    # Write/update per-galaxy results row
    write_result_row_ecsv(results_path, row)
    
    if not args.no_diagnostic_plots:
        print("making diagnostic plots...")
        e.plot_fancy_profiles()
        e.draw_phot_results_mpl()

    ################################################################
    # block for cs-gr if the image exists
    ################################################################

    if csgr_fits and args.csgr:
        t0 = time.perf_counter()
        logger.info("Running optional CS-gr ellipse photometry")

        row["CSGR_EXISTS"] = True
        row["CSGR_FITS"] = Path(csgr_fits).name

        e_gr = run_ellipse_photometry(
            r_fits=r_fits,
            cs_fits=csgr_fits,
            mask_fits=mask_fits,
            image2_filter=hafilter,
            filter_ratio=filter_ratio,
            objra=ra,
            objdec=dec,
            fixcenter=args.fixcenter,
            logger=logger,
            fileid="csgr"
        )

        copy_image2_fields_to_row(e_gr, row, "CSGR")

        e_gr.run_hapy_morphology()
        copy_hapy_cs_fields_to_row(e_gr, row, "CSGR", pixscale=pixscale)

        if valid_file(e_gr.photfile) and valid_file(e_gr.photfile2):
            row["CSGR_PHOT_OK"] = True

            rtab = Table.read(e_gr.photfile)
            hatab = Table.read(e_gr.photfile2)

            profile_results_gr = summarize_dual_profiles(
                rtab=rtab,
                hatab=hatab,
                r_magzp=magzp,
            )

            row.update(prefix_dict_keys(profile_results_gr, "CSGR"))
            
        add_sourcecatalog_moments(row, e_gr, prefix1="CSGR_R", prefix2="CSGR_H",pixel_scale = float(pixscale))
        row["CSGR_SEC"] = _scalar(time.perf_counter() - t0)
    

    # ---- RUN HAPY CLUMP ANALYSIS  ----------- #

    if args.clumps and csgr_fits and args.csgr:

        print("DEBUG: trying to run clump analysis")
        from hapy.ellipse.clumps import ClumpDetectionConfig

        clump_config = ClumpDetectionConfig(
            nsigma=args.clump_nsigma,
            npixels=args.clump_npixels,

            deblend=not args.no_clump_deblend,
            nlevels=args.clump_nlevels,
            contrast=args.clump_contrast,
            mode=args.clump_deblend_mode,
            
            find_peaks=not args.no_clump_peaks,
            peak_box_size=args.clump_peak_box_size,
            peak_min_separation=args.clump_peak_min_separation,

            background_grow_radius=args.clump_background_grow_radius,
            background_npixels=args.clump_background_npixels,
            background_mask_nsigma=args.clump_background_mask_nsigma,
            background_clip_sigma=args.clump_background_clip_sigma,
            background_clip_maxiters=args.clump_background_clip_maxiters,

            min_flux=args.clump_min_flux,
            min_area=args.clump_min_area,

            save_diagnostic=not args.no_clump_diagnostic,
            diagnostic_format=args.clump_diagnostic_format,
            diagnostic_percent=args.clump_diagnostic_percent,
            plot_kron_apertures=args.clump_plot_kron_apertures,

            find_point_sources=args.clump_point_sources,
            point_source_method=args.clump_point_source_method,
            point_source_fwhm=args.clump_point_source_fwhm,
            point_source_threshold_nsigma=args.clump_point_source_threshold_nsigma,
        )

        # temporarily moving this out of try/except for debugging
        # csgr_clump_result = e_gr.measure_halpha_clumps(
        #     config=clump_config,
        #     output_dir=cutdir,
        #     basename=tag,
        #     overwrite=True,
        #     update_results=False,   # keep this False while testing
        # )
        

        try:
            csgr_clump_result = e_gr.measure_halpha_clumps(
                config=clump_config,
                output_dir=cutdir,
                basename=tag,
                overwrite=True,
                update_results=False,   # keep this False while testing
            )

            write_hapy_clumps(
                row,
                csgr_clump_result,
                prefix="CSGR_HCL_",
                input_image="CS-GR",
                config=clump_config,
                pixel_scale=pixscale,
                failed=False,
                )
            
            logger.info(
                "H-alpha clump analysis complete: NCLUMP=%d, NPEAK=%d",
                csgr_clump_result.summary.n_clumps,
                getattr(csgr_clump_result.summary, "n_peaks", 0),
            )

        except Exception as err:
            logger.warning("H-alpha clump analysis failed: %s", err)
            print("H-alpha clump analysis failed: %s", err)

            if hasattr(e, "results") and isinstance(e.results, dict):
                e_gr.results["CSGR_HCL_OK"] = False
                e_gr.results["CSGR_HCL_STATUS"] = "failed"
            write_hapy_clumps(
                row,
                None,
                prefix="CSGR_HCL_",
                input_image="CS-GR",
                config=clump_config,
                failed=True,
                )
                
        
    ################################################################
    # statmorph block
    ################################################################

    if args.statmorph:
        t0 = time.perf_counter()
        logger.info("STAGE: statmorph")
        e.run_statmorph_supervisor()
        if e.statmorph_flag:
            #_pull_statmorph(row,"R_SM", getattr(e, "morph", None), pixscale)
            try:
                _pull_statmorph(row,"R_SM", getattr(e, "morph", None), pixscale)
                # statmorph sets flag == 1 for a problem, so need to negate it
                row["R_SM_OK"] = True
            except Exception:
                pass

            try:
                _pull_statmorph(row,"H_SM", getattr(e, "morph2", None), pixscale)
                row["H_SM_OK"] = True
            except Exception:
                pass
        # write table after statmorph
        row["SM_SEC"] = _scalar(time.perf_counter() - t0)        
        write_result_row_ecsv(results_path, row)
 

        
    ################################################################
    # galfit block
    ################################################################

    if args.galfit:
        print("starting galfit ...")
        #print("DEBUG: cutdir = ",root)
        #print("DEBUG: tag = ",tag)        
        row["STAGE"] = "galfit_nc"
        logger.info("STAGE: galfit NC")
        t0 = time.perf_counter()
        #galname = Path(root).name  # no .fits; matches your test
        galname = tag  # no .fits; matches your test
        pscale = get_pixel_scale_from_filename(r_fits)

        r_data, hdr = fits.getdata(r_fits, header=True)
        ny, nx = r_data.shape
        xminfit, xmaxfit = 1, nx
        yminfit, ymaxfit = 1, ny


        convflag = bool(args.convflag)

        # set convolution box size
        # GALFIT manual recommends a box size of >=20 seeing diameters
        nconvolution_scale = 20

        h_fwhm = params.get("himage_fwhm_psf_arcsec", None)
        r_fwhm = params.get("rimage_fwhm_psf_arcsec", None)

        # fallbacks for metadata-light / archive / VESTIGE-style metadata
        generic_fwhm = params.get("fwhm_psf_arcsec", None)
        ha_finaliq = params.get("ha_finaliq", None)
        r_finaliq = params.get("r_finaliq", None)

        if h_fwhm is not None and np.isfinite(float(h_fwhm)):
            seeing_arcsec = float(h_fwhm)
        elif r_fwhm is not None and np.isfinite(float(r_fwhm)):
            seeing_arcsec = float(r_fwhm)
        elif generic_fwhm is not None and np.isfinite(float(generic_fwhm)):
            seeing_arcsec = float(generic_fwhm)
        elif ha_finaliq is not None and np.isfinite(float(ha_finaliq)):
            seeing_arcsec = float(ha_finaliq)
        elif r_finaliq is not None and np.isfinite(float(r_finaliq)):
            seeing_arcsec = float(r_finaliq)
        else:
            seeing_arcsec = 2.0
            logger.info("no FWHM in metadata.json - assuming seeing = 2 arcsec")

        convolution_size = nconvolution_scale * seeing_arcsec / pixscale

        # GALFIT wants an integer convolution box size
        convolution_size = int(np.ceil(convolution_size))

        # keep it within image bounds
        convolution_size = min(convolution_size, nx)

        logger.info(
            "GALFIT convolution box: seeing=%.3f arcsec, pixscale=%.4f arcsec/pix, size=%d pix",
            seeing_arcsec,
            pixscale,
            convolution_size,
        )

        # going back to original convolutionsize
        #convolution_size = min(nx, ny)
        if psf_path is not None:
            psf_image = str(psf_path)
        else:
            psf_image = None
        rg = RunGalfit(
            galname=galname,
            image=r_fits,
            sigma_image=sigma_image,
            psf_image=psf_image,
            psf_oversampling=args.psf_oversampling,
            mask_image=mask_fits,
            xminfit=xminfit,
            yminfit=yminfit,
            xmaxfit=xmaxfit,
            ymaxfit=ymaxfit,
            convolution_size=convolution_size, # this is the full image
            magzp=magzp,
            pscale=pscale,
            convflag=convflag,
            fitallflag=False,
            ncomp=args.ncomp,
            asym=False,
            workdir=cutdir,
            localize_files=True,
            psf_local_name="psf.fits",
        )

        t0 = time.perf_counter()

        #xc = nx / 2
        #yc = ny / 2

        # try to get a more sensible initial radius for galfit
        sma_pix = sma_arcsec / pixscale
        rad_init = max(sma_pix, 30)
        
        init0 = dict(xobj=xc, yobj=yc, mag=10.0, rad=rad_init, nsersic=2.0, BA=0.7, PA=0.0, first_time=1)
        
        # --- No convolution ---
        try:
            res_nc, meta_nc = _galfit_stage(rg, args, init0, do_conv=False, logger=logger)
            _store_galfit(row, res_nc, "GAL_", pixscale)
            row["GAL_NC_RERUN_FIXEDN"] = meta_nc["rerun_fixed_n"]
            row["GAL_NC_OK"] = not meta_nc["unstable"]
            row["GAL_NC_SEC"] = time.perf_counter() - t0
        except Exception as e:
            logger.exception(f"GALFIT NC failed: {e}")
            row["GAL_NC_OK"] = False
        write_result_row_ecsv(results_path, row)

        # print this by setting --log-to-console at command line
        # print(
        #     f"GALFIT NC: chi2nu={_scalar(res_nc.chi2nu):.3f} "
        #     f"re={_scalar(res_nc.comp1.re):.2f} "
        #     f"n={_scalar(res_nc.comp1.n):.2f} "
        #     f"ba={_scalar(res_nc.comp1.ba):.2f} "
        #     f"pa={_scalar(res_nc.comp1.pa):.1f}"
        #     )
        
        if args.convflag and not psf_ok:
            logger.warning("convflag requested but PSF not available; skipping convolution.")
        
        if psf_ok and args.convflag and row["GAL_NC_OK"]:
            # --- Convolution (init from NC) ---
            t0 = time.perf_counter()
            row["STAGE"] = "galfit_cv"
            logger.info("STAGE: galfit CV")

            init_from_nc = dict(
                xobj=_scalar(res_nc.comp1.xc), yobj=_scalar(res_nc.comp1.yc),
                mag=_scalar(res_nc.comp1.mag), rad=_scalar(res_nc.comp1.re),
                nsersic=_scalar(res_nc.comp1.n), BA=_scalar(res_nc.comp1.ba), PA=_scalar(res_nc.comp1.pa),
                first_time=0,
                )


            use_nc_init = (
                (not meta_nc["unstable"])
                and (_scalar(res_nc.error) == 0)
                and (_scalar(res_nc.comp1.numerical_error_flag) == 0)
                and _finite_dict(init_from_nc)
                and (_scalar(res_nc.comp1.re) > 3.)
                )

            if use_nc_init:
                init_cv = init_from_nc
                row["GAL_CV_INIT_FROM_NC"] = True

            else:
                init_cv = dict(init0)
                init_cv["first_time"] = 0
                row["GAL_CV_INIT_FROM_NC"] = False
                logger.warning("GALFIT CV init: NC results flagged; using fallback init0")
                
            try:
                res_cv, meta_cv = _galfit_stage(rg, args, init_cv, do_conv=True, logger=logger)
                _store_galfit(row, res_cv, "GAL_C", pixscale)
                row["GAL_CV_RERUN_FIXEDN"] = meta_cv["rerun_fixed_n"]
                row["GAL_CV_OK"] = not meta_cv["unstable"]

                row["GAL_CV_SEC"] = time.perf_counter() - t0
                #row["galfit_ok"] = row["GAL_CV_OK"]  # or (NC_OK and CV_OK) if you prefer
            except Exception as e:
                logger.exception(f"GALFIT CV failed: {e}")
                row["GAL_CV_OK"] = False
                # keep NC results, continue

            # print(
            #     f"GALFIT CV: chi2nu={_scalar(res_cv.chi2nu):.3f} "
            #     f"re={_scalar(res_cv.comp1.re):.2f} "
            #     f"n={_scalar(res_cv.comp1.n):.2f} "
            #     f"ba={_scalar(res_cv.comp1.ba):.2f} "
            #     f"pa={_scalar(res_cv.comp1.pa):.1f}"
            #     )
        
        write_result_row_ecsv(results_path, row)
    row["TOTAL_SEC"] = time.perf_counter() - t0_total
    row["STAGE"] = "done"
    row["STATUS"] = "ok"

    write_result_row_ecsv(results_path, row)
    print(f"Wrote results: {results_path}")
    #return results_path
    return 0
    #return e

        
if __name__ == "__main__":
    #e = main()

    raise SystemExit(main())
    
    # checking table - comment after check
    #check_table(results_table)
