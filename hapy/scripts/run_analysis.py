#!/usr/bin/env python
"""
Run analysis on a single galaxy cutout set (1 GNU-parallel task).

Typical inputs are the *root prefix* created by build_cutout_name(), e.g.
cutouts/VF1234-HDI-20200226-p012/VF1234-HDI-20200226-p012

Examples:
  python scripts/run_analysis.py --root cutouts/.../VF1234-HDI-20200226-p012
  ls cutouts/*/* | grep -v '\.fits$' | parallel python scripts/run_analysis.py --root {}


Running in parallel:

In two steps.  Step 1/2
```
for d in cutouts/*/; do
  b=$(basename "$d")
  echo "${d}${b}"
done > roots.txt
```

Step 2/2
```
parallel -j 8 python scripts/run_analysis.py --root {} --make-mask --galfit :::: roots.txt
```

One Step, Step 1/1:
```
for d in cutouts/*/; 
  b=$(basename "$d")
  echo "${d}${b}"
done | parallel -j 8 python scripts/run_analysis.py --root {} --make-mask --galfit
```


TESTING:

python ~/github/hapy/scripts/run_analysis.py --root cutouts/<tag>/<tag> --make-mask --galfit --convflag 1

cat cutouts/<tag>/<tag>-results.ecsv


from coadds directory where you just made cutouts:

python ~/github/hapy/scripts/run_analysis.py --root cutouts/VFID3084-NGC3512-HDI-20200226-p012/VFID3084-NGC3512-HDI-20200226-p012 --make-mask --convflag 0 --psf-image VF-165.869+28.044-HDI-20200226-p012-r-psf.fits --statmorph --image2-filter 4 --galfit

To make use of --cutout-dir instead of --root:

python ~/github/hapy/scripts/run_analysis.py --cutout-dir cutouts/VFID3084-NGC3512-HDI-20200226-p012 --make-mask --convflag 0 --psf-image VF-165.869+28.044-HDI-20200226-p012-r-psf.fits --statmorph --image2-filter 4 --galfit


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
from hapy.masktools.api import MaskEngine, EllipseParams
from hapy.masktools.gaia import make_gaia_mask,  get_gaia_stars, galaxy_overlaps_bright_star
from hapy.masktools.maskops import distance_to_nearest_mask, largest_mask_region, ellipse_mask_fraction
#from hapy.masktools.types import build_ell0_from_metadata
from hapy.hatools.results import write_result_row_ecsv
from hapy.geometry.adapters import pa_ccw_north_to_photutils_theta, photutils_theta_to_pa_ccw_north
from hapy.utils.paths import astromatic_dir 
#from hapy.utils.logging_utils import setup_logging

from hapy.ellipse.profile_summary import summarize_dual_profiles



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

def _store_galfit(row, res, prefix):
    row[f"{prefix}XC"] = _scalar(res.comp1.xc)
    row[f"{prefix}XC_ERR"] = _scalar(res.comp1.xc_err)
    row[f"{prefix}YC"] = _scalar(res.comp1.yc)
    row[f"{prefix}YC_ERR"] = _scalar(res.comp1.yc_err)

    row[f"{prefix}MAG"] = _scalar(res.comp1.mag)
    row[f"{prefix}MAG_ERR"] = _scalar(res.comp1.mag_err)
    row[f"{prefix}RE"] = _scalar(res.comp1.re)
    row[f"{prefix}RE_ERR"] = _scalar(res.comp1.re_err)
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

def _pull_statmorph(row, prefix, mobj):
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
            ("RPETRO_ELLIP", "rpetro_ellip"),
            ("RHALF_ELLIP", "rhalf_ellip"),
            ("R20", "r20"),
            ("R50", "r50"),
            ("R80", "r80"),
            ("RMAX_CIRCLE","rmax_circ"),
            ("RMAX_ELLIP","rmax_ellip"),            
            ("SERSIC_AMP", "sersic_amplitude"),
            ("SERSIC_RHALF","sersic_rhalf"),
            ("SERSIC_N","sersic_n"),
            ("SERSIC_XC","sersic_xc"),
            ("SERSIC_YC","sersic_yc"),
            ("SERSIC_ELLIP","sersic_ellip"),
            ("SERSIC_THETA","sersic_theta"),
            ("SERSIC_CHISQ_DOF","sersic_chi2_dof"),
            ("SERSIC_FLAG","flag_sersic"),
            ("SN_PER_PIXEL","sn_per_pixel"),
                    ]:
        row[f"{prefix}_{outk}"] = _scalar(getattr(mobj, attr))
        try:
            row[f"{prefix}_{outk}"] = _scalar(getattr(mobj, attr))
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


    # removing these
    not_needed = ["R_FITS", "CS_FITS"]
    # ---------- identity ----------
    for k in [
        "TAG", "CUTDIR",
         "MASK_FITS","PSF_FITS",
         "R_FITS", "CS_FITS","SIGMA_FITS",
         "RFILTER_FILENAME", "RFILTER_CENTER","RFILTER_WIDTH",
         "HFILTER_FILENAME", "HFILTER_CENTER","HFILTER_WIDTH",         
            ]:
        row[k] = ""
        
    #row["psf_ok"] = False

    row["PSF_SOURCE"] = ""   # "cli" | "psf_dir" | ""

    for k in ["STAGE", "STATUS"]:
        row[k] = ""

    for k in ["MASK_SEC", "PHOT_SEC", "GALFIT_SEC", "TOTAL_SEC"]:
        row[k] = np.nan    
   
    # ---------- pipeline status ----------
    for k in ["PSF_OK", "MASK_OK", "PHOT_OK", \
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
    row["BRIGHT_STAR_DIST"] = np.nan
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

    # ---------- ellipse ----------
    for k in [
        "ELLIP_XCENTROID", "ELLIP_YCENTROID",
        "ELLIP_SMA_PIX", "ELLIP_B_PIX",
        "ELLIP_EPS", "ELLIP_THETA_RAD",
        "ELLIP_GINI_DET", "ELLIP_SOURCE_SUM",
        "ELLIP_BA", "ELLIP_SEGMENT_FLUX",
        "ELLIP_SEGMENT_MAG"
    ]:
        row[k] = np.nan


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

    # ---------- mismatch ----------
    for k in ["ELL_DC_PX", "ELL_DBA", "ELL_DPA_DEG", "ELL_SMA_RATIO"]:
        row[k] = np.nan

    # ---------- statmorph ----------
    sm_suffixes = [
        "FLAG",
        "XCENTROID", "YCENTROID", "GINI", "M20",
        "C", "A", "S",
        "RPETRO_ELLIP", "RHALF_ELLIP",
        "R20", "R50", "R80",
        "RMAX_CIRCLE", "RMAX_ELLIP",
        "SERSIC_AMP", 
        "SERSIC_RHALF",
        "SERSIC_N",
        "SERSIC_XC",
        "SERSIC_YC",
        "SERSIC_ELLIP",
        "SERSIC_THETA",
        "SERSIC_CHISQ_DOF",
        "SERSIC_FLAG",
        "SN_PER_PIXEL",
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
    ]:
        row[k] = np.nan


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
        "GAL_RE", "GAL_RE_ERR",
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
        "GAL_CRE", "GAL_CRE_ERR",
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
    g_main.add_argument("--statmorph", action="store_true",
                        help="Compute statmorph structural parameters")
    g_main.add_argument("--galfit", action="store_true",
                        help="Run GALFIT after photometry")
    g_main.add_argument("--convflag", action="store_true", default=False,
                        help="Run GALFIT convolution stage (requires PSF)")
    g_main.add_argument("--no-diagnostic-plots", action="store_true",
                        help="Don't write diagnostic plot (R image + mask + ellipses)")


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
                        help="SExtractor SNR threshold.  Default is 5.")
    g_mask.add_argument("--seminarea", type=int, default=5,
                        help="SExtractor minimum object area. Default is 7.")
    g_mask.add_argument("--grow-size", default=7,
                        help="Grow size in mask expansion.  Default is 5.")
    g_mask.add_argument("--grow-iterations", default=4,
                        help="Grow size in mask expansion.  Default is 4.")
    #g_mask.add_argument("--gaiapath", default=None,
    #                    help="Path to Gaia catalog file")
    g_mask.add_argument("--gaia-dir", default="gaia_catalogs",
                        help="Directory containing precomputed Gaia catalogs (default: gaia_catalogs)")
    g_mask.add_argument("--no-gaia", action="store_true",
                        help="Disable Gaia star masking")
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
    r_fits = args.r_fits or _pick_one(str(cutdir/ f"{tag}*-R.fits")) or _pick_one(str(cutdir / f"{tag}*-r.fits"))
    if r_fits is None:
        raise FileNotFoundError(f"Could not find R-band FITS in: {cutdir}")

    cs_fits = args.cs_fits or _pick_one(str(cutdir / f"{tag}*-CS-ZP.fits")) or _pick_one(str(cutdir / f"{tag}*-cs.fits"))

    # why are we looking for a mask when we are suppose to make one?
    mask_fits = args.mask_fits or _pick_one(str(cutdir / f"{tag}*-mask.fits"))

    
    sigma_image = args.sigma_image or _pick_one(str(cutdir / f"{tag}*-sigma.fits")) or _pick_one(str(cutdir / f"{tag}*-rms.fits"))
    psf_image = args.psf_image or _pick_one(str(cutdir / f"{tag}*-psf.fits"))


            
    row = initialize_result_row()

    
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
    row["GALFIT_SEC"] = 0.0
    row["TOTAL_SEC"] = 0.0

    row["R_FITS"] = r_fits
    row["CS_FITS"] = cs_fits
    row["SIGMA_FITS"] = sigma_image    

    
    # pixel scale
    #pixscale = args.pixscale
    #if pixscale is None:
    pixscale = get_pixel_scale_from_filename(r_fits)


    # --- Load cutout image for WCS + shape ---
    data, hdr = fits.getdata(r_fits, header=True)
    ny, nx = data.shape
    wcs = WCS(hdr)
    
    magzp = args.magzp if args.magzp is not None else float(hdr.get("PHOTZP", 25.0))
    
    row["CUTOUT_XSIZE"] = nx
    row["CUTOUT_YSIZE"] = ny
    # Default center = image center
    xc = nx / 2.0
    yc = ny / 2.0


    t0_total = time.perf_counter()

    params_path = cutdir / "metadata.json"

    if not params_path.exists():
        raise RuntimeError(
            f"metadata.json not found for root {root}. "
            "Cutouts may be outdated or improperly generated."
        )
    
    params = json.loads(params_path.read_text())

    # --- check for valid input ellipse
    if ellipse_missing(params):
        print("\nGetting initial ellipse estimate from photutils...\n")
        # get gaia mask
        brightstar, star_xpix, star_ypix = get_gaia_stars(r_fits)
        mask_array = np.zeros_like(data,  dtype=np.int32)
        # could add a radius_scale_factor that depends on image FWHM
        # (FWHM-1.5)
        mask_array, gaia_mask = make_gaia_mask(mask_array,star_xpix,star_ypix,pixscale/3600.,gaia_table=brightstar,radius_scale_factor=1)

        # convert to boolean mask
        gaia_mask = gaia_mask > 0

        # get ellipse from photutils
        ell = infer_ellipse_from_r_cutout(r_data=data, user_mask=gaia_mask)
        if ell is not None:
            # if agc has a valid radius and BA, then keep?
            #print("DEBUG: original radius = ",params["sma_arcsec"])
            radius_scale_factor = 1.2
            #print("DEBUG: new radius = ",ell.sma_pix * pixscale)
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


    # --- update row with other info from metadata.json
    row["TELESCOPE"] = params.get("telescope", "")
    row["DATEOBS"]   = params.get("dateobs", "")
    row["POINTING"]  = params.get("pointing", "")
    row["SCHEME"]    = params.get("scheme", "")
    row["PARENT_RIMAGE"]  = params.get("parent_rimage", "")
    row["PARENT_HIMAGE"] = params.get("parent_haimage", "")

    row["HFILTER"] = params.get("hafilter")
    row["CUTOUT_SCALE"] = params.get("cutout_scale")
    row["FILTER_CORRECTION"] = params.get("filter_correction")


    # TODO get this information from ZP ratio
    filter_ratio = params.get("filter_ratio", None)
    if filter_ratio is None:
        filter_ratio = np.nan
        logger.warning("FLTRATIO missing from metadata; physical flux calibration will be NaN.")

    row["FILTER_RATIO"] = filter_ratio           

    # --- Get ellipse parameters ---    
    sma_arcsec = float(params["sma_arcsec"])
    ba = float(params["ba"])
    pa_deg = float(params["pa_deg"]) # CCW from N, from input catalog
    #xc = float(params["xc"])
    #yc = float(params["yc"])    

    # Try WCS-based centering using stored RA/DEC
    #ra = params.get("ra", None)
    #dec = params.get("dec", None)
    objid = params.get("objid", Path(root).name)
    ra = args.objra if args.objra is not None else params.get("ra")
    dec = args.objdec if args.objdec is not None else params.get("dec")
    row["RA"] = ra
    row["DEC"] = dec
    row["OBJID"] = objid

    row["R_FWHM_SE"] = float(params.get("rimage_fwhm_se_arcsec"))
    row["R_FWHM_PSF"] = float(params.get("rimage_fwhm_psf_arcsec"))    
    row["H_FWHM_SE"] = float(params.get("himage_fwhm_se_arcsec"))
    row["H_FWHM_PSF"] = float(params.get("himage_fwhm_psf_arcsec"))    
    row["RFILTER_FILENAME"] = params.get("rfilter_name")
    row["RFILTER_CENTER"] = float(params.get("rfilter_center_A"))
    row["RFILTER_WIDTH"] = float(params.get("rfilter_width_A"))    
    row["HFILTER_FILENAME"] = params.get("hafilter_name")
    row["HFILTER_CENTER"] = float(params.get("hafilter_center_A"))
    row["HFILTER_WIDTH"] = float(params.get("hafilter_width_A"))    


    if ra is not None and dec is not None:
        try:
            xw, yw = wcs.world_to_pixel_values(float(ra), float(dec))
            if np.isfinite(xw) and np.isfinite(yw):
                xc, yc = float(xw), float(yw)
        except Exception:
            pass

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

    ell0_params = EllipseParams(
        xc = xc,
        yc = yc,
        ba = ba,
        sma_pix = sma_arcsec/pixscale,
        theta_deg = photutils_theta_to_pa_ccw_north(pa_deg)
        )

    # --- Construct the name of the psf image
    psf_path, psf_source = pick_psf_path_and_source(args, params,logger=logger)
    row["PSF_FITS"] = str(psf_path) if psf_path else ""
    row["PSF_OK"] = bool(psf_path)
    row["PSF_SOURCE"] = psf_source
    psf_ok = row["PSF_OK"]

    #print("TESTING: psf_path = ",psf_path)
    #sys.exit()


    
    if args.make_mask:

        if args.seconfig is None:
            raise ValueError("--sex-config must be set when --make-mask is used")
        
        row["STAGE"] = "mask"
        logger.info("STAGE: mask")
        
        t0 = time.perf_counter()

        # choose output mask name if not provided/found
        mask_out = mask_fits or (cutdir / f"{tag}-mask.fits")

        print(f"DEBUG: mask_out={mask_out}")


        


        #row["sma_arcsec"] = sma_arcsec
        #row["ba"] = ba
        #row["pa_deg"] = pa_deg        

        # --- Convert to pixels ---
        sma_pix = sma_arcsec / pixscale

        # convert CCW from N angle to photutils CCW from +x
        galaxy_ellipse = ell0_params
        #theta_deg = pa_ccw_north_to_photutils_theta(pa_deg)
        #galaxy_ellipse = EllipseParams(
        #    xc=xc,
        #    yc=yc,
        #    sma_pix=sma_pix,
        #    ba=ba,
        #    theta_deg=theta_deg,
        #)

        # -- look for gaia table from parent image
        parent_rimage = params.get("parent_rimage", None)        
        gaia_table_path = None
        if parent_rimage:
            gaia_dir = Path(args.gaia_dir)
            gaia_table_path = gaia_dir / parent_rimage.replace(".fits", "-gaia.fits")
    

        if not args.no_gaia and parent_rimage:
            gaia_table_path = Path(args.gaia_dir) / parent_rimage.replace(".fits", "-gaia.fits")
            if gaia_table_path.exists():
                gaia_table = Table.read(gaia_table_path)
                logger.info(f"Using local Gaia catalog: {gaia_table_path}")
            else:
                logger.warning(f"Local Gaia catalog not found: {gaia_table_path}")
        
        engine = MaskEngine(
            image_fits=r_fits,
            sepath=args.sepath,
            #gaiapath=args.gaiapath,
            config=args.seconfig,
            threshold=args.sethreshold,
            snr=args.sesnr,
            minarea=args.seminarea,
            add_gaia_stars=(not args.no_gaia),
        )
        # calculate the min radius to use for gaia stars
        max_fwhm = max(row["R_FWHM_PSF"], row["H_FWHM_PSF"])
        gaia_min_radius_arcsec = 4 * max_fwhm
        logger.info(f"Gaia min radius (arcsec) = {gaia_min_radius_arcsec}")

        # convert gaia min radius to get
        gaia_min_radius_deg = gaia_min_radius_arcsec/3600.
        mask = engine.build_initial_mask(
            galaxy_ellipse=galaxy_ellipse,
            progress_callback=_progress_cb,
            grow_size=int(args.grow_size),
            grow_iterations=int(args.grow_iterations),
            gaia_table = gaia_table,
            gaia_min_radius = gaia_min_radius_deg,
        )

        #mask_out = mask_fits or (root + "-mask.fits")
        engine.write_mask(mask_out)
        mask_fits = mask_out



        row["MASK_OK"] = True
        row["MASK_FITS"] = str(mask_fits)

        row["MASK_SEC"] = time.perf_counter() - t0
        #row["mask_ok"] = True
        write_result_row_ecsv(results_path, row)

        bright_flag, dist_arcsec, maskrad_arcsec, bright_mag = galaxy_overlaps_bright_star(
            ra,
            dec,
            gaia_table,
            mag_limit=10,
            radius_col="radius",
            min_radius_arcsec = gaia_min_radius_arcsec
            )
        
        row["BRIGHT_STAR_FLAG"] = bright_flag
        row["BRIGHT_STAR_DIST"] = dist_arcsec
        row["BRIGHT_STAR_MASKRAD_ARCSEC"] = maskrad_arcsec
        row["BRIGHT_STAR_MAG"] = bright_mag        


        res = ellipse_mask_fraction(mask, ell0_params)
        row["ELL0_MASKFRAC"] = res.frac_masked
        row["ELL0_MASK_WARN"] = res.frac_masked > 0.5        
        row["ELL0_NMASKPIX"] = res.n_masked
        row["ELL0_NTOTPIX"] = res.n_total

        # under development
        #ellipse_pixels = aper_image > 0

        #largest_blob = largest_mask_region(mask_image, ellipse_pixels)

        #row["ELL_LARGEST_MASK"] = largest_blob
        #largest_blob = largest_mask_region(mask_image, ellipse_pixels)

        #row["ELL0_LARGEST_MASK"] = largest_blob
        
        # dmask = distance_to_nearest_mask(mask_image, ell0.xc, ell0.yc)

        # row["NEAR_MASK_DIST_PIX"] = dmask
        # row["NEAR_MASK_WARN"] = dmask < (2 * ell0.sma_pix)
        


        
    row["STAGE"] = "phot"
    logger.info("STAGE: phot")
    
    t0 = time.perf_counter()

    hafilter = row["HFILTER"]
    if args.image2_filter is not None:
        hafilter = args.image2_filter
        row["HFILTER"] = hafilter
        
    e = run_ellipse_photometry(
        r_fits=r_fits,
        cs_fits=cs_fits,
        mask_fits=mask_fits,
        image2_filter=hafilter,
        filter_ratio=filter_ratio,
        objra=ra,
        objdec=dec,
        fixcenter=args.fixcenter,
        #run_statmorph=args.statmorph,
        #write_prefix=prefix,
    )

    # ---- photometry summary (scalar-only; arrays stay in the photometry table files) ----
    row["PHOT_SEC"] = _scalar(time.perf_counter() - t0)
    row["PHOT_OK"] = True

    # ---- core ellipse / detection-derived quantities ----
    FIELDS = [
        ("ELLIP_XCENTROID", "xcenter"),
        ("ELLIP_YCENTROID", "ycenter"),
        ("ELLIP_SMA_PIX", "sma"),
        ("ELLIP_B_PIX", "b"),
        ("ELLIP_EPS", "eps"),
        ("ELLIP_THETA_RAD", "theta"),
        ("ELLIP_GINI_DET", "gini"),
        ("ELLIP_SOURCE_SUM", "source_sum"),
        ("ELLIP_SEGMENT_FLUX", "photutils_segment_flux"),
        ("ELLIP_SEGMENT_MAG","photutils_segment_mag"),
        ("R_SKYSTD_ADU", "sky_noise"),
        ("R_SKYMED_ADU", "sky"),
        ("R_SKYSTD_PHYS", "im1_skynoise"),
        ("R_M20", "M20_1"),
        ("R_ASYM", "asym"),
        ("R_ASYM_ERR", "asym_err"),
        ("AREA_GUESS_ELLIPSE_PIX", "area_guess_ellipse_pix"),
        ("AREA_GUESS_ELLIPSE_UNMASKED_PIX", "area_guess_ellipse_unmasked_pix"),
        ("MASKFRAC_GUESS_ELLIPSE", "maskfrac_guess_ellipse"),
        ]


    for outk, attr in FIELDS:
        v = getattr(e, attr, None)
        sv = _scalar(v)
        if sv is not None:
            row[outk] = sv  # leave as np.nan if missing/array/etc.

    # add photutils B/A
    row["ELLIP_BA"] = 1. - float(row["ELLIP_EPS"])

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
            ]

    if getattr(e, "image2", None) is not None:
        for outk, attr in FIELDS2:
            v = getattr(e, attr, None)
            sv = _scalar(v)
            if sv is not None:
                row[outk] = sv

        
    try:
        phot_xc = float(row["ELLIP_XCENTROID"])
        phot_yc = float(row["ELLIP_YCENTROID"])
        phot_sma_pix = float(row["ELLIP_SMA_PIX"])
        phot_ba = 1.0 - float(row["ELLIP_EPS"])
        #phot_pa_deg = (np.degrees(float(row["ELLIP_THETA_RAD"])) % 180.0)
        #phot_pa_deg = photutils_theta_to_pa_ccw_north(theta_phot_deg)  # inverse of your adapter

        # ELLIP_THETA_RAD measured from +x axis
        phot_theta_deg = (np.degrees(float(row["ELLIP_THETA_RAD"])) % 180.0)
        phot_pa_deg = photutils_theta_to_pa_ccw_north(phot_theta_deg)  # inverse of your adapter

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

    # ---- FIT PROFILES!  ----------- #

    if valid_file(e.photfile) and valid_file(e.photfile2):
        row["PHOT_OK"] = True

        rtab = Table.read(e.photfile)
        hatab = Table.read(e.photfile2)

        profile_results = summarize_dual_profiles(
            rtab=rtab,
            hatab=hatab,
            r_magzp=magzp,
        )

        row.update(profile_results)

    # Write/update per-galaxy results row
    write_result_row_ecsv(results_path, row)


    if args.statmorph:
        logger.info("STAGE: statmorph")
        e.run_statmorph_supervisor()
        if e.statmorph_flag:
            #_pull_statmorph(row,"R_SM", getattr(e, "morph", None))
            try:
                _pull_statmorph(row,"R_SM", getattr(e, "morph", None))
                # statmorph sets flag == 1 for a problem, so need to negate it
                row["R_SM_OK"] = True
            except Exception:
                pass

            try:
                _pull_statmorph(row,"H_SM", getattr(e, "morph2", None))
                row["H_SM_OK"] = True
            except Exception:
                pass
        # write table after statmorph
        write_result_row_ecsv(results_path, row)
 
    
    if not args.no_diagnostic_plots:
        print("making diagnostic plots...")
        e.plot_fancy_profiles()
        e.draw_phot_results_mpl()
        
        phot_xc = float(row["ELLIP_XCENTROID"])
        phot_yc = float(row["ELLIP_YCENTROID"])
        phot_sma_pix = float(row["ELLIP_SMA_PIX"])
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
        plot_mask_ellipse_diagnostic(
            r_fits=str(r_fits),
            mask_fits=str(mask_fits),
            e0=ell0_params,
            eph=ellphot_params,
            outfile=str(outfile),
            row=row,
            )
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

        data, hdr = fits.getdata(r_fits, header=True)
        ny, nx = data.shape
        xminfit, xmaxfit = 1, nx
        yminfit, ymaxfit = 1, ny

        

        convflag = bool(args.convflag)

        # set convolution box size
        nconvolution_scale = 20 # galfit manual says use box size of 20 or more seeing diameters
        if params['himage_fwhm_psf_arcsec'] is not None:
            convolution_size = nconvolution_scale * float(params['himage_fwhm_psf_arcsec'])/pixscale
        elif params['rimage_fwhm_psf_arcsec'] is not None:
            convolution_size = nconvolution_scale * float(params['rimage_fwhm_psf_arcsec'])/pixscale
        else:
            # set to the number of pixels with
            # assume seeing = 2 arcsec, and 0.4"/pixels
            logger.info("no FWHM for in metadata.json - assuming 2 arcsec")
            convolution_size = nconvolution_scale * 2/0.4

        
        convolution_size = min(convolution_size, nx)
        # going back to original convolutionsize
        #convolution_size = min(nx, ny)
        rg = RunGalfit(
            galname=galname,
            image=r_fits,
            sigma_image=sigma_image,
            psf_image=str(psf_path),
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
            _store_galfit(row, res_nc, "GAL_")
            row["GAL_NC_RERUN_FIXEDN"] = meta_nc["rerun_fixed_n"]
            row["GAL_NC_OK"] = not meta_nc["unstable"]

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
                _store_galfit(row, res_cv, "GAL_C")
                row["GAL_CV_RERUN_FIXEDN"] = meta_cv["rerun_fixed_n"]
                row["GAL_CV_OK"] = not meta_cv["unstable"]

                row["GALFIT_SEC"] = time.perf_counter() - t0
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

        
if __name__ == "__main__":
    #results_table = main()
    raise SystemExit(main())
    
    # checking table - comment after check
    #check_table(results_table)
