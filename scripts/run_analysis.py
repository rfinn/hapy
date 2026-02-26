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
for d in cutouts/*/; do
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

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import json
import re
import time

from hapy.ellipse.photometry import run_ellipse_photometry
from hapy.galfittools.rungalfit import RunGalfit
from hapy.imagetools.imutils import get_pixel_scale_from_filename
from hapy.masktools.api import MaskEngine, EllipseParams
from hapy.hatools.results import write_result_row_ecsv


def _default_sex_config() -> str:
    # hapy/astromatic/default.sex.HDI.mask (adjust if package name differs)
    try:
        with resources.as_file(resources.files("hapy.astromatic") / "default.sex.HDI.mask") as p:
            return str(p)
    except Exception:
        # fallback relative to repo
        return str(Path(__file__).resolve().parents[1] / "hapy" / "astromatic" / "default.sex.HDI.mask")

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

def resolve_psf_path(psfdir, parent_rimage):
    """Return PSF path derived from parent R coadd filename, or None."""
    if not psfdir or not parent_rimage:
        return None
    psf_name = str(parent_rimage).replace(".fits", "-psf.fits")
    psf_path = Path(psfdir) / psf_name
    return psf_path if psf_path.exists() else None
    
def _galfit_stage(rg, args, init, do_conv: bool, n_hi=8.0):
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
    res = rg.run_and_parse()

    meta = {"stage": stage, "rerun_fixed_n": False, "unstable": False}

    # rerun if high n
    if _scalar(res.comp1.n) > n_hi:
        meta["rerun_fixed_n"] = True
        rg.set_sersic_params(
            xobj=_scalar(res.comp1.xc), yobj=_scalar(res.comp1.yc),
            mag=_scalar(res.comp1.mag), rad=_scalar(res.comp1.re),
            nsersic=4.0, BA=_scalar(res.comp1.ba), PA=_scalar(res.comp1.pa),
            fitmag=1, fitcenter=1, fitrad=1, fitBA=1, fitPA=1, fitn=0,
            first_time=0,
        )
        rg.set_sky(args.sky)
        res = rg.run_and_parse()

    meta["unstable"] = unstable(res)
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

def initialize_result_row():
    """Return a fully populated results-row dict with frozen schema."""

    row = {}

    for k in [
        "TELESCOPE",
        "DATEOBS",
        "POINTING",
        "SCHEME",
        "PARENT_RIMAGE",
        "PARENT_HAIMAGE",
    ]:
        row[k] = ""

    # ---------- virgo identifiers ----------
    row["VFID"] = ""      # e.g., "VFID3084"
    row["GALNAME"] = ""   # e.g., "NGC3512" (optional but handy)

    
    # ---------- identity ----------
    for k in [
        "objid", "tag", "root",
        "r_fits", "cs_fits", "mask_fits",
        "psf_fits", "sigma_fits", "scheme"
    ]:
        row[k] = ""
        
    row["psf_ok"] = False

    row["PSF_SOURCE"] = ""   # "cli" | "psf_dir" | ""
    
    # ---------- pipeline status ----------
    for k in ["mask_ok", "phot_ok", "galfit_ok"]:
        row[k] = False

    for k in ["STAGE", "STATUS"]:
        row[k] = ""

    for k in ["MASK_SEC", "PHOT_SEC", "GALFIT_SEC", "TOTAL_SEC"]:
        row[k] = np.nan

    # ---------- coordinates ----------
    row["ra"] = np.nan
    row["dec"] = np.nan

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
        "ELLIP_GINI_DET", "ELLIP_SOURCE_SUM"
    ]:
        row[k] = np.nan

    row["ELLIP_MASKED_FRACTION"] = ""

    # ---------- R band ----------
    for k in [
        "R_FWHM", "R_SKYSTD_ADU", "R_SKYMED_ADU",
        "R_SKYSTD_PHYS", "R_M20", "R_ASYM", "R_ASYM_ERR"
    ]:
        row[k] = np.nan

    # ---------- H band ----------
    for k in [
        "H_SKYSTD_ADU", "H_SKYMED_ADU",
        "H_SKYSTD_PHYS", "H_M20", "H_ASYM", "H_ASYM_ERR"
    ]:
        row[k] = np.nan

    # ---------- mismatch ----------
    for k in ["ELL_DC_PX", "ELL_DBA", "ELL_DPA_DEG", "ELL_SMA_RATIO"]:
        row[k] = np.nan
    row["ELL_MISMATCH"] = False

    # ---------- statmorph ----------
    sm_suffixes = [
        "XCENTROID", "YCENTROID", "GINI", "M20",
        "C", "A", "S",
        "RPETRO_ELLIP", "RHALF_ELLIP",
        "R20", "R50", "R80"
    ]

    for band in ["R", "H"]:
        for s in sm_suffixes:
            row[f"{band}_SM_{s}"] = np.nan
        row[f"{band}_SM_FLAG"] = False

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

    row["GAL_NC_RERUN_FIXEDN"] = False
    row["GAL_NC_OK"] = False

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

    row["GAL_CV_RERUN_FIXEDN"] = False
    row["GAL_CV_OK"] = False

    return row

def pick_psf_path_and_source(args, params):
    psf_image = getattr(args, "psf_image", None)
    if psf_image:
        p = Path(psf_image)
        return (str(p), "cli") if p.exists() else (None, "cli_missing")

    psfdir = getattr(args, "psf_dir", None) or getattr(args, "psfdir", None)
    parent = params.get("parent_rimage", "")
    if psfdir and parent:
        name = str(parent).replace(".fits", "-psf.fits")
        p = Path(psfdir) / name
        return (str(p), "psf_dir") if p.exists() else (None, "psf_dir_missing")

    return None, ""

def main():
    p = argparse.ArgumentParser(description="Run headless analysis on one galaxy cutout set")
    p.add_argument("--root", help="Cutout root prefix (no extension) with relative path. e.g. --root cutouts/VFID3084-NGC3512-HDI-20200226-p012/VFID3084-NGC3512-HDI-20200226-p012.  See --cutout-dir for easier interface for Virgo or UAT surveys.")
    p.add_argument("--cutout-dir", default=None, help="Cutout directory containing the root basename (optional).  Example --cutout-dir VFID3084-NGC3512-HDI-20200226-p012.  This will construct the root filename as cutouts/VFID3084-NGC3512-HDI-20200226-p012/VFID3084-NGC3512-HDI-20200226-p012.")    
    p.add_argument("--r", dest="r_fits", default=None, help="Override R-band FITS path")
    p.add_argument("--cs", dest="cs_fits", default=None, help="Override CS FITS path")
    p.add_argument("--mask", dest="mask_fits", default=None, help="Override mask FITS path")

    # photometry knobs (passed through)
    p.add_argument("--image2-filter", dest="image2_filter", default=None)
    p.add_argument("--filter-ratio", dest="filter_ratio", type=float, default=None)
    p.add_argument("--objra", type=float, default=None)
    p.add_argument("--objdec", type=float, default=None)
    p.add_argument("--fixcenter", action="store_true")
    p.add_argument("--statmorph", action="store_true")

    p.add_argument("--prefix", default=None, help="Output prefix tag (default: root basename)")
    p.add_argument("--no-plots", action="store_true", help="Skip profile/diagnostic plots")

    # MASK
    p.add_argument("--make-mask", action="store_true", help="Build/write mask before photometry/galfit")
    p.add_argument("--sepath", default="sex")
    p.add_argument("--gaiapath", default=None)
    p.add_argument("--sex-config", dest="sex_config", default=_default_sex_config(), help="SExtractor config file path (default: hapy/astromatic/default.sex.HDI.mask)")
    
    p.add_argument("--threshold", type=float, default=0.005)
    p.add_argument("--snr", type=float, default=10.0)
    p.add_argument("--minarea", type=int, default=5)
    p.add_argument("--no-gaia", action="store_true", help="Disable Gaia star masking")
    p.add_argument("--pixscale", type=float, default=None, help="Override pixel scale (arcsec/pix)")
    p.add_argument("--sma-arcsec", type=float, default=None, help="Ellipse semi-major axis in arcsec (optional)")
    p.add_argument("--ba", type=float, default=None, help="Ellipse b/a (optional)")
    p.add_argument("--pa-deg", type=float, default=None, help="Ellipse PA deg (optional)")
    
    # GALFIT
    p.add_argument("--galfit", action="store_true", help="Run GALFIT after photometry")
    p.add_argument("--sigma-image", dest="sigma_image", default=None, help="Override sigma/RMS image (optional)")
    p.add_argument("--psf-dir", dest="psf_dir", default=None, help="Directory containing PSF image.  Assumes the PSF image in the r-band coadd image (from metadata.json) rfilename.replace('.fits','-psf.fits')")    
    p.add_argument("--psf-image", dest="psf_image", default=None, help="Override PSF image in metadata.json (optional)")
    p.add_argument("--psf-oversampling", type=int, default=2)
    p.add_argument("--convflag", type=int, default=1, help="1 convolve with PSF, 0 otherwise")
    p.add_argument("--ncomp", type=int, default=1, choices=[1, 2])
    p.add_argument("--magzp", type=float, default=None)
    p.add_argument("--sky", type=float, default=0.0)
    
    args = p.parse_args()

    if args.cutout_dir is not None:
        d = Path(args.cutout_dir)
        # allow passing either "cutouts/<tag>" or "<tag>"
        if not d.exists():
            d = Path("cutouts") / d
        if not d.exists():
            raise FileNotFoundError(f"cutout dir not found: {args.cutout_dir}")
        args.root = str(d / d.name)   # matches your current convention
    else:
        # if user passed cutouts/<tag>/<tag>, collapse to cutouts/<tag>/<tag>
        rp = Path(args.root)
        if rp.parent.name == rp.name:
            args.root = str(rp)  # already canonical
        elif (rp.parent / rp.parent.name) == rp:
            args.root = str(rp)  # already canonical
        # otherwise leave as-is
    
    root = args.root
    root_base = Path(root).name
    prefix = args.prefix or root_base


    cutdir = Path(root).parent
    tag = Path(root).name
    results_path = cutdir / f"{tag}-results.ecsv"
    
    # Auto-detect common filenames if not provided.
    # Adjust these glob patterns to match your exact suffix conventions.
    r_fits = args.r_fits or _pick_one(root + "*-R.fits") or _pick_one(root + "*-r.fits")
    if r_fits is None:
        raise FileNotFoundError(f"Could not find R-band FITS for root: {root}")

    cs_fits = args.cs_fits or _pick_one(root + "*-CS-ZP.fits") or _pick_one(root + "*-cs.fits") or _pick_one(root + "*-cs.fits")
    mask_fits = args.mask_fits or _pick_one(root + "*-mask.fits")
    sigma_image = args.sigma_image or _pick_one(root + "*-sigma.fits") or _pick_one(root + "*-rms.fits")
    psf_image = args.psf_image or _pick_one(root + "*-psf.fits")



            
    row = initialize_result_row()
    row["root"] = str(root)
    row["tag"] = Path(root).name
    row["STAGE"] = "init"
    row["STATUS"] = "running"
    row["MASK_SEC"] = 0.0
    row["PHOT_SEC"] = 0.0
    row["GALFIT_SEC"] = 0.0
    row["TOTAL_SEC"] = 0.0
    
    # pixel scale
    pixscale = args.pixscale
    if pixscale is None:
        pixscale = get_pixel_scale_from_filename(r_fits)


    # --- Load cutout image for WCS + shape ---
    data, hdr = fits.getdata(r_fits, header=True)
    ny, nx = data.shape
    wcs = WCS(hdr)

    # Default center = image center
    xc = nx / 2.0
    yc = ny / 2.0








    t0_total = time.perf_counter()

    params_path = Path(root).parent / "metadata.json"

    if not params_path.exists():
        raise RuntimeError(
            f"metadata.json not found for root {root}. "
            "Cutouts may be outdated or improperly generated."
        )
    
    params = json.loads(params_path.read_text())
    row["TELESCOPE"] = params.get("telescope", "")
    row["DATEOBS"]   = params.get("dateobs", "")
    row["POINTING"]  = params.get("pointing", "")
    row["SCHEME"]    = params.get("scheme", "")
    row["PARENT_RIMAGE"]  = params.get("parent_rimage", "")
    row["PARENT_HAIMAGE"] = params.get("parent_haimage", "")


    # --- Get ellipse parameters ---    
    sma_arcsec = float(params["sma_arcsec"])
    ba = float(params["ba"])
    pa_deg = float(params["pa_deg"])
    #xc = float(params["xc"])
    #yc = float(params["yc"])    

    # Try WCS-based centering using stored RA/DEC
    ra = params.get("ra", None)
    dec = params.get("dec", None)
    objid = params.get("objid", Path(root).name)

    row["ra"] = ra
    row["dec"] = dec
    row["objid"] = objid

    if ra is not None and dec is not None:
        try:
            xw, yw = wcs.world_to_pixel_values(float(ra), float(dec))
            if np.isfinite(xw) and np.isfinite(yw):
                xc, yc = float(xw), float(yw)
        except Exception:
            pass

    # -- let use input an ellipse geometry that is different from what is in metadata.json
    if args.sma_arcsec is not None and args.ba is not None and args.pa_deg is not None:
        sma_arcsec = float(args.sma_arcsec)
        ba = float(args.ba)
        pa_deg = float(args.pa_deg)


    objid = row.get("objid", "") or ""

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
    row["ELL0_PA_DEG"] = pa_deg
    row["ELL0_XC"] = xc
    row["ELL0_YC"] = yc
    row["ELL0_SOURCE"] = "metadata.json" if params_path.exists() else "args"


    # --- Construct the name of the psf image
    psf_path, psf_source = pick_psf_path_and_source(args, params)
    row["psf_fits"] = str(psf_path) if psf_path else ""
    row["psf_ok"] = bool(psf_path)
    row["PSF_SOURCE"] = psf_source
    psf_ok = row["psf_ok"]

    print("TESTING: psf_path = ",psf_path)
    #sys.exit()
    if args.make_mask:

        if args.sex_config is None:
            raise ValueError("--sex-config must be set when --make-mask is used")
        
        row["STAGE"] = "mask"
        t0 = time.perf_counter()

        # choose output mask name if not provided/found
        mask_out = mask_fits or (root + "-mask.fits")






        #row["sma_arcsec"] = sma_arcsec
        #row["ba"] = ba
        #row["pa_deg"] = pa_deg        

        # --- Convert to pixels ---
        sma_pix = sma_arcsec / pixscale

        galaxy_ellipse = EllipseParams(
            xc=xc,
            yc=yc,
            sma_pix=sma_pix,
            ba=ba,
            pa_deg=pa_deg,
        )
        engine = MaskEngine(
            image_fits=r_fits,
            sepath=args.sepath,
            gaiapath=args.gaiapath,
            config=args.sex_config,
            threshold=args.threshold,
            snr=args.snr,
            minarea=args.minarea,
            add_gaia_stars=(not args.no_gaia),
        )

        mask = engine.build_initial_mask(
            galaxy_ellipse=galaxy_ellipse,
            progress_callback=_progress_cb,
        )

        mask_out = mask_fits or (root + "-mask.fits")
        engine.write_mask(mask_out)
        mask_fits = mask_out



        row["mask_ok"] = True
        row["mask_fits"] = str(mask_fits)

        row["MASK_SEC"] = time.perf_counter() - t0
        row["mask_ok"] = True
        write_result_row_ecsv(results_path, row)


    row["STAGE"] = "phot"
    t0 = time.perf_counter()
    
    e = run_ellipse_photometry(
        r_fits=r_fits,
        cs_fits=cs_fits,
        mask_fits=mask_fits,
        image2_filter=args.image2_filter,
        filter_ratio=args.filter_ratio,
        objra=args.objra,
        objdec=args.objdec,
        fixcenter=args.fixcenter,
        run_statmorph=args.statmorph,
        write_prefix=prefix,
    )

    # ---- photometry summary (scalar-only; arrays stay in the photometry table files) ----
    row["PHOT_SEC"] = _scalar(time.perf_counter() - t0)
    row["phot_ok"] = True

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
        ("R_FWHM", "fwhm"),
        ("R_SKYSTD_ADU", "sky_noise"),
        ("R_SKYMED_ADU", "sky"),
        ("R_SKYSTD_PHYS", "im1_skynoise"),
        ("R_M20", "M20_1"),
        ("R_ASYM", "asym"),
        ("R_ASYM_ERR", "asym_err"),
        ]


    for outk, attr in FIELDS:
        v = getattr(e, attr, None)
        sv = _scalar(v)
        if sv is not None:
            row[outk] = sv  # leave as np.nan if missing/array/etc.
                
    # JSON field (stable schema)
    mf = getattr(e, "masked_fraction", None)
    if mf is not None:
        try:
            row["ELLIP_MASKED_FRACTION"] = json.dumps(mf)
        except TypeError:
            # if mf contains numpy types etc.
            row["ELLIP_MASKED_FRACTION"] = json.dumps(mf, default=_scalar)
        

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

    # compute difference between input ellipse and photutils ellipse
    def _angle_diff_deg(a, b):
        d = (a - b + 90.0) % 180.0 - 90.0
        return abs(d)

    try:
        phot_xc = float(row["ELLIP_XCENTROID"])
        phot_yc = float(row["ELLIP_YCENTROID"])
        phot_sma_pix = float(row["ELLIP_SMA_PIX"])
        phot_ba = 1.0 - float(row["ELLIP_EPS"])
        phot_pa_deg = (np.degrees(float(row["ELLIP_THETA_RAD"])) % 180.0)

        dx = phot_xc - float(row["ELL0_XC"])
        dy = phot_yc - float(row["ELL0_YC"])
        dc = float(np.hypot(dx, dy))

        dba = abs(phot_ba - float(row["ELL0_BA"]))
        dpa = _angle_diff_deg(phot_pa_deg, float(row["ELL0_PA_DEG"]))

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
            (dc > 5.0) or (dba > 0.2) or (dpa > 25.0) or
            (np.isfinite(sma_ratio) and ((sma_ratio < 0.5) or (sma_ratio > 2.0)))
            )
    except Exception:
        # if any missing keys, just don't set mismatch fields
        pass

    # ---- statmorph (best-effort; only a few key fields to start) ----
    def _pull_statmorph(prefix, mobj):
        if mobj is None:
            return
        for outk, attr in [
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
            ("FLAG", "flag"),
        ]:
            try:
                row[f"{prefix}_{outk}"] = _scalar(getattr(mobj, attr))
            except Exception:
                pass

    try:
        _pull_statmorph("R_SM", getattr(e, "morph", None))
        row["R_SM_FLAG"] = _scalar(getattr(e, "statmorph_flag", None))
    except Exception:
        pass

    try:
        _pull_statmorph("H_SM", getattr(e, "morph2", None))
        row["H_SM_FLAG"] = _scalar(getattr(e, "statmorph_flag2", None))
    except Exception:
        pass

    # Write/update per-galaxy results row
    write_result_row_ecsv(results_path, row)


    if not args.no_plots:
        e.plot_fancy_profiles()
        e.draw_phot_results_mpl()

    if args.galfit:
        row["STAGE"] = "galfit"
        t0 = time.perf_counter()
        galname = root  # no .fits; matches your test
        pscale = get_pixel_scale_from_filename(r_fits)

        data, hdr = fits.getdata(r_fits, header=True)
        ny, nx = data.shape
        xminfit, xmaxfit = 1, nx
        yminfit, ymaxfit = 1, ny

        
        magzp = args.magzp if args.magzp is not None else float(hdr.get("PHOTZP", 25.0))
        convflag = bool(args.convflag)

        
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
            convolution_size=min(nx, ny),
            magzp=magzp,
            pscale=pscale,
            convflag=convflag,
            fitallflag=False,
            ncomp=args.ncomp,
            asym=False,
        )

        t0 = time.perf_counter()

        #xc = nx / 2
        #yc = ny / 2
        init0 = dict(xobj=xc, yobj=yc, mag=15.0, rad=10.0, nsersic=2.0, BA=0.7, PA=0.0, first_time=1)

        # --- No convolution ---
        res_nc, meta_nc = _galfit_stage(rg, args, init0, do_conv=False)
        _store_galfit(row, res_nc, "GAL_")
        row["GAL_NC_RERUN_FIXEDN"] = meta_nc["rerun_fixed_n"]
        row["GAL_NC_OK"] = not meta_nc["unstable"]

        if psf_ok:
            # --- Convolution (init from NC) ---
            init_cv = dict(
                xobj=_scalar(res_nc.comp1.xc), yobj=_scalar(res_nc.comp1.yc),
                mag=_scalar(res_nc.comp1.mag), rad=_scalar(res_nc.comp1.re),
                nsersic=_scalar(res_nc.comp1.n), BA=_scalar(res_nc.comp1.ba), PA=_scalar(res_nc.comp1.pa),
                first_time=0,
                )
            res_cv, meta_cv = _galfit_stage(rg, args, init_cv, do_conv=True)
            _store_galfit(row, res_cv, "GAL_C")
            row["GAL_CV_RERUN_FIXEDN"] = meta_cv["rerun_fixed_n"]
            row["GAL_CV_OK"] = not meta_cv["unstable"]

            row["GALFIT_SEC"] = time.perf_counter() - t0
            row["galfit_ok"] = row["GAL_CV_OK"]  # or (NC_OK and CV_OK) if you prefer

        
        write_result_row_ecsv(results_path, row)
    row["TOTAL_SEC"] = time.perf_counter() - t0_total
    row["STAGE"] = "done"
    row["STATUS"] = "ok"

    write_result_row_ecsv(results_path, row)
    print(f"Wrote results: {results_path}")
    return results_path

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

    for c in ["ELLIP_MASKED_FRACTION", "SM_R_FLAG", "SM_H_FLAG"]:
        print(c, t[c].dtype, t[c][0])
        
if __name__ == "__main__":
    results_table = main()

    # checking table - comment after check
    #check_table(results_table)
