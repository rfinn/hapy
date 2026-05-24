#!/usr/bin/env python

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

from hapy.ellipse.photometry import run_ellipse_photometry
from hapy.ellipse.profile_summary import summarize_dual_profiles
from hapy.hatools.utils import zp_scale_r_to_ha
#from hapy.utils.results_table import _scalar


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
def valid_file(path):
    return path is not None and Path(path).is_file()


def find_one(cutdir, patterns):
    for pat in patterns:
        matches = sorted(cutdir.glob(pat))
        if matches:
            return matches[0]
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Measure ellipse profiles on Halpha-with-continuum images using r-band geometry."
    )
    parser.add_argument("--cutout-dir", required=True, help="HAPY cutout directory")
    parser.add_argument("--image2-filter", default=None, help="Override Halpha filter name")
    parser.add_argument("--fixcenter", action="store_true", help="Fix center during ellipse fitting")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output table")
    args = parser.parse_args()

    cutdir = Path(args.cutout_dir)
    tag = cutdir.name

    results_path = cutdir / f"{tag}-results-ha-with-continuum.ecsv"
    if results_path.exists() and not args.overwrite:
        print(f"Skipping existing output: {results_path}")
        return

    r_fits = cutdir / f"{tag}-R.fits"
    h_fits = cutdir / f"{tag}-Ha.fits"

    if not r_fits.exists():
        raise FileNotFoundError(f"Missing r-band cutout: {r_fits}")
    if not h_fits.exists():
        raise FileNotFoundError(f"Missing Halpha cutout: {h_fits}")

    manual_mask = cutdir / f"{tag}-mask-manual.fits"
    auto_mask = cutdir / f"{tag}-mask.fits"
    mask_fits = manual_mask if manual_mask.exists() else auto_mask
    if not mask_fits.exists():
        mask_fits = None

    params_path = cutdir / "metadata.json"
    if not params_path.exists():
        raise FileNotFoundError(f"Missing metadata.json: {params_path}")

    params = json.loads(params_path.read_text())

    ra = params.get("ra")
    dec = params.get("dec")
    hafilter = args.image2_filter or params.get("hafilter")

    if ra is None or dec is None:
        raise ValueError(f"Missing ra/dec in metadata.json for {tag}")

    rheader = fits.getheader(r_fits)
    hheader = fits.getheader(h_fits)

    pixscale = abs(float(rheader.get("PIXSCALE", params.get("pixscale", np.nan))))
    if not np.isfinite(pixscale):
        # fallback from CD matrix, deg/pix -> arcsec/pix
        pixscale = abs(float(rheader.get("CD1_1", np.nan))) * 3600.0

    magzp = float(rheader.get("PHOTZP", params.get("r_photzp", np.nan)))

    zp_r = rheader.get("PHOTZP", np.nan)
    zp_h = hheader.get("PHOTZP", np.nan)
    filter_ratio = zp_scale_r_to_ha(zp_h, zp_r)

    row = {
        "TAG": tag,
        "R_FITS": r_fits.name,
        "H_FITS": h_fits.name,
        "MASK_FITS": "" if mask_fits is None else Path(mask_fits).name,
        "HFILTER": "" if hafilter is None else hafilter,
        "FILTER_RATIO": filter_ratio,
        "RA": float(ra),
        "DEC": float(dec),
        "PHOT_OK": False,
        "PHOT_SEC": np.nan,
    }

    t0 = time.perf_counter()

    e = run_ellipse_photometry(
        r_fits=str(r_fits),
        cs_fits=str(h_fits),
        mask_fits=None if mask_fits is None else str(mask_fits),
        image2_filter=hafilter,
        filter_ratio=filter_ratio,
        objra=ra,
        objdec=dec,
        fixcenter=args.fixcenter,
        logger=None,
        fileid="HaCont",
    )

    row["PHOT_SEC"] = _scalar(time.perf_counter() - t0)

    fields = [
        ("ELLIP_XCENTROID", "xcenter_fit"),
        ("ELLIP_YCENTROID", "ycenter_fit"),
        ("ELLIP_SMA_ARCSEC", "sma_fit"),
        ("ELLIP_SMA_PIX", "sma_fit"),
        ("ELLIP_B_ARCSEC", "b"),
        ("ELLIP_EPS", "eps_fit"),
        ("ELLIP_THETA_RAD", "pa_fit"),
        ("R_SKYSTD_ADU", "sky_noise"),
        ("R_SKYMED_ADU", "sky"),
        ("R_SKYSTD_PHYS", "im1_skynoise"),
        ("H_SKYSTD_ADU", "sky_noise2"),
        ("H_SKYMED_ADU", "sky2"),
        ("H_SKYSTD_PHYS", "im2_skynoise"),
        ("H_SCALE_ADU_CGS", "uconversion2"),
    ]

    for outk, attr in fields:
        sv = _scalar(getattr(e, attr, None))
        if sv is not None:
            if outk.endswith("_ARCSEC"):
                sv = sv * pixscale
            row[outk] = sv

    if row.get("ELLIP_THETA_RAD", None) is not None:
        row["ELLIP_THETA_DEG"] = np.degrees(float(row["ELLIP_THETA_RAD"])) % 180.0

    if row.get("ELLIP_EPS", None) is not None:
        row["ELLIP_BA"] = 1.0 - float(row["ELLIP_EPS"])

    if valid_file(e.photfile) and valid_file(e.photfile2):
        rtab = Table.read(e.photfile)
        htab = Table.read(e.photfile2)

        profile_results = summarize_dual_profiles(
            rtab=rtab,
            hatab=htab,
            r_magzp=magzp,
        )

        row.update(profile_results)
        row["PHOT_OK"] = True

    outtab = Table(rows=[row])
    outtab.write(results_path, format="ascii.ecsv", overwrite=True)
    print(f"Wrote {results_path}")


if __name__ == "__main__":
    main()
