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
"""

import argparse
from pathlib import Path
import glob

from astropy.io import fits
from astropy.wcs import WCS

import json


from hapy.ellipse.photometry import run_ellipse_photometry
from hapy.galfittools.rungalfit import RunGalfit
from hapy.imagetools.imutils import get_pixel_scale_from_filename
from hapy.masktools.api import MaskEngine, EllipseParams
from hapy.hatools.results import write_result_row_ecsv


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
    
def main():
    p = argparse.ArgumentParser(description="Run headless analysis on one galaxy cutout set")
    p.add_argument("--root", required=True, help="Cutout root prefix (no extension)")
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
    p.add_argument("--sex-config", dest="sex_config", default=None, help="e.g. default.sex.HDI.mask")
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
    p.add_argument("--psf-image", dest="psf_image", default=None, help="Override PSF image (optional)")
    p.add_argument("--psf-oversampling", type=int, default=2)
    p.add_argument("--convflag", type=int, default=1, help="1 convolve with PSF, 0 otherwise")
    p.add_argument("--ncomp", type=int, default=1, choices=[1, 2])
    p.add_argument("--magzp", type=float, default=None)
    p.add_argument("--sky", type=float, default=0.0)
    
    args = p.parse_args()

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

    cs_fits = args.cs_fits or _pick_one(root + "*-CS.fits") or _pick_one(root + "*-cs.fits")
    mask_fits = args.mask_fits or _pick_one(root + "*-mask.fits")
    sigma_image = args.sigma_image or _pick_one(root + "*-sigma.fits") or _pick_one(root + "*-rms.fits")
    psf_image = args.psf_image or _pick_one(root + "*-psf.fits")


    row = dict(
        objid=tag,                # you can replace with params["objid"] if you prefer
        tag=tag,
        root=str(root),
        r_fits=str(r_fits),
        cs_fits=str(cs_fits) if cs_fits else "",
        mask_fits=str(mask_fits) if mask_fits else "",
        psf_fits=str(psf_image) if psf_image else "",
        sigma_fits=str(sigma_image) if sigma_image else "",
        scheme="",                # optional: pass in from CLI if you want
        mask_ok=False,
        phot_ok=False,
        galfit_ok=False,
    )
    
    if args.make_mask:
        # choose output mask name if not provided/found
        mask_out = mask_fits or (root + "-mask.fits")

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

        # --- Get ellipse parameters ---
        params_path = Path(root).parent / "mask_params.json"

        if params_path.exists():
            params = json.loads(params_path.read_text())

            sma_arcsec = float(params["sma_arcsec"])
            ba = float(params["ba"])
            pa_deg = float(params["pa_deg"])

            # Try WCS-based centering using stored RA/DEC
            ra = params.get("ra", None)
            dec = params.get("dec", None)
            objid = params.get("objid", Path(root).name)
            
            if ra is not None and dec is not None:
                try:
                    xw, yw = wcs.world_to_pixel_values(float(ra), float(dec))
                    if np.isfinite(xw) and np.isfinite(yw):
                        xc, yc = float(xw), float(yw)
                except Exception:
                    pass

        elif args.sma_arcsec is not None and args.ba is not None and args.pa_deg is not None:
            sma_arcsec = float(args.sma_arcsec)
            ba = float(args.ba)
            pa_deg = float(args.pa_deg)

        else:
            raise ValueError(
                "Masking requires ellipse params. Provide mask_params.json or "
                "--sma-arcsec/--ba/--pa-deg."
            )

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
        write_result_row_ecsv(results_path, row)


        
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
    row["phot_ok"] = True

    # Core ellipse / detection-derived quantities
    for outk, attr in [
        ("ELLIP_XCENTROID", "xcenter"),
        ("ELLIP_YCENTROID", "ycenter"),
        ("ELLIP_SMA_GUESS_PIX", "sma"),
        ("ELLIP_B_GUESS_PIX", "b"),
        ("ELLIP_EPS_GUESS", "eps"),
        ("ELLIP_THETA_GUESS_RAD", "theta"),
        ("ELLIP_GINI_DET", "gini"),
        ("ELLIP_SOURCE_SUM", "source_sum"),
        ("ELLIP_MASKED_FRACTION", "masked_fraction"),
        ("R_FWHM", "fwhm"),
        ("R_SKYNOISE", "sky_noise"),
        ("R_SKY", "sky"),
        ("M20_R", "M20_1"),
        ("ASYM_R", "asym"),
        ("ASYM_R_ERR", "asym_err"),
    ]:
        try:
            row[outk] = _scalar(getattr(e, attr))
        except Exception:
            pass

    # If image2 exists (e.g., continuum-sub / HA), capture analogous scalars
    try:
        if getattr(e, "image2_flag", False):
            for outk, attr in [
                ("H_SKYNOISE", "sky_noise2"),
                ("H_SKY", "sky2"),
                ("M20_H", "M20_2"),
                ("ASYM_H", "asym2"),
                ("ASYM_H_ERR", "asym2_err"),
            ]:
                try:
                    row[outk] = _scalar(getattr(e, attr))
                except Exception:
                    pass
    except Exception:
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
        _pull_statmorph("SM_R", getattr(e, "morph", None))
        row["SM_R_FLAG"] = _scalar(getattr(e, "statmorph_flag", None))
    except Exception:
        pass

    try:
        _pull_statmorph("SM_H", getattr(e, "morph2", None))
        row["SM_H_FLAG"] = _scalar(getattr(e, "statmorph_flag2", None))
    except Exception:
        pass

    # Write/update per-galaxy results row
    write_result_row_ecsv(results_path, row)


    if not args.no_plots:
        e.plot_fancy_profiles()
        e.draw_phot_results_mpl()

    if args.galfit:
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
            psf_image=psf_image,
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

        xc = nx / 2
        yc = ny / 2
        rg.set_sersic_params(
            xobj=xc, yobj=yc,
            mag=15.0,
            rad=10.0,
            nsersic=2.0,
            BA=0.7,
            PA=0.0,
            fitmag=1, fitcenter=1, fitrad=1, fitBA=1, fitPA=1, fitn=1,
            first_time=1,
        )
        rg.set_sky(args.sky)
        rg.run_and_parse()
        row["galfit_ok"] = True
        write_result_row_ecsv(results_path, row)
if __name__ == "__main__":
    main()
