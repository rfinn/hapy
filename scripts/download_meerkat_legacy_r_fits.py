#!/usr/bin/env python

"""
Download Legacy Survey r-band FITS cutouts for the Meerkat interesting sample.

The default catalog is:
    ~/research/JWST/Meerkat_interesting_v2.csv

Required input columns:
    VFID_1, RA_1, DEC_1
"""

import argparse
import csv
import os
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlretrieve

import numpy as np
from astropy.io import fits


DEFAULT_CSV = "~/research/JWST/Meerkat_interesting_v2.csv"
REQUIRED_COLS = ["VFID_1", "RA_1", "DEC_1"]
RETRIABLE_HTTP_CODES = {429, 500, 502, 503, 504}


def expand_path(path):
    """
    Expand ~ and environment variables in a user-supplied path.
    """
    return Path(os.path.expandvars(os.path.expanduser(str(path))))


def clean_galid(galid):
    """
    Make a galaxy id safe for use in output filenames.
    """
    galid = str(galid).strip()
    galid = re.sub(r"[^\w.-]+", "_", galid)
    return galid or "unknown"


def read_catalog(csvfile):
    """
    Yield catalog rows with validated VFID, RA, and DEC values.
    """
    with open(csvfile, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [col for col in REQUIRED_COLS if col not in reader.fieldnames]

        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        for lineno, row in enumerate(reader, start=2):
            vfid = str(row["VFID_1"]).strip()

            if not vfid:
                raise ValueError(f"Missing VFID_1 on line {lineno}")

            try:
                ra = float(row["RA_1"])
                dec = float(row["DEC_1"])
            except ValueError as err:
                raise ValueError(f"Bad RA_1/DEC_1 on line {lineno}: {err}")

            yield {
                "VFID_1": vfid,
                "RA_1": ra,
                "DEC_1": dec,
                "lineno": lineno,
            }


def legacy_r_fits_url(ra, dec, imsize, pixscale=0.262, layer="ls-dr9", band="r"):
    """
    Build the Legacy Survey FITS cutout URL.
    """
    return (
        "https://www.legacysurvey.org/viewer/cutout.fits?"
        + urlencode(
            {
                "ra": ra,
                "dec": dec,
                "layer": layer,
                "size": int(imsize),
                "pixscale": pixscale,
                "bands": band,
            }
        )
    )


def urlretrieve_with_retries(url, filename, retries=8, sleep0=15, verbose=False):
    """
    Retrieve a URL with simple backoff for transient server errors.
    """
    last_err = None

    for attempt in range(1, retries + 1):
        try:
            if verbose:
                print(f"download attempt {attempt}/{retries}: {filename}")
            return urlretrieve(url, filename)

        except HTTPError as err:
            last_err = err
            if err.code not in RETRIABLE_HTTP_CODES:
                raise

        except URLError as err:
            last_err = err

        wait = sleep0 * attempt
        if verbose:
            print(f"WARNING: download failed: {last_err}; sleeping {wait}s")
        time.sleep(wait)

    raise last_err


def get_first_image_data(hdul):
    """
    Return image data from either extension 1 or the primary HDU.
    """
    if len(hdul) > 1 and hdul[1].data is not None:
        return hdul[1].data
    return hdul[0].data


def validate_fits_image(fits_name):
    """
    Confirm that the downloaded FITS file opens and contains nonzero image data.
    """
    with fits.open(fits_name) as hdul:
        data = get_first_image_data(hdul)

        if data is None:
            raise ValueError("no image data found")

        data = np.asarray(data)
        if data.size == 0:
            raise ValueError("empty image data array")

        if np.all(data == 0):
            raise ValueError("image data are all zero")


def download_one(
    row,
    outdir,
    imsize,
    pixscale=0.262,
    layer="ls-dr9",
    band="r",
    overwrite=False,
    dry_run=False,
    verbose=False,
):
    """
    Download one Legacy Survey FITS cutout.
    """
    galid = clean_galid(row["VFID_1"])
    fits_name = outdir / f"{galid}-legacy-{int(imsize)}-{band}.fits"
    url = legacy_r_fits_url(
        row["RA_1"],
        row["DEC_1"],
        imsize=imsize,
        pixscale=pixscale,
        layer=layer,
        band=band,
    )

    result = {
        "VFID_1": row["VFID_1"],
        "RA_1": row["RA_1"],
        "DEC_1": row["DEC_1"],
        "fits_name": str(fits_name),
        "status": "",
        "message": "",
        "url": url,
    }

    if dry_run:
        result["status"] = "dry-run"
        result["message"] = "not downloaded"
        print(url)
        return result

    if fits_name.exists() and not overwrite:
        try:
            validate_fits_image(fits_name)
            result["status"] = "exists"
            result["message"] = "valid existing file"
            return result
        except Exception as err:
            if verbose:
                print(f"WARNING: existing file failed validation; re-downloading: {err}")

    tmp_name = fits_name.with_suffix(fits_name.suffix + ".part")

    if verbose:
        print(f"downloading {row['VFID_1']} -> {fits_name}")
        print(url)

    try:
        urlretrieve_with_retries(url, tmp_name, verbose=verbose)
        tmp_name.replace(fits_name)
        validate_fits_image(fits_name)
        result["status"] = "downloaded"
        result["message"] = "ok"

    except Exception as err:
        if tmp_name.exists():
            tmp_name.unlink()
        result["status"] = "failed"
        result["message"] = str(err)

    return result


def write_manifest(results, manifest):
    """
    Write a CSV summary of download results.
    """
    fieldnames = ["VFID_1", "RA_1", "DEC_1", "fits_name", "status", "message", "url"]

    with open(manifest, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


def download_meerkat_legacy_r_fits(
    csvfile=DEFAULT_CSV,
    outdir=None,
    cutout_arcmin=10.0,
    pixscale=0.262,
    layer="ls-dr9",
    overwrite=False,
    dry_run=False,
    verbose=False,
):
    """
    Download 10 arcmin by 10 arcmin Legacy Survey r-band FITS images.
    """
    csvfile = expand_path(csvfile)

    if outdir is None:
        outdir = csvfile.parent / "legacy-r-fits"
    else:
        outdir = expand_path(outdir)

    imsize = int(round(cutout_arcmin * 60.0 / pixscale))

    if verbose:
        print(f"catalog = {csvfile}")
        print(f"outdir = {outdir}")
        print(f"cutout_arcmin = {cutout_arcmin}")
        print(f"pixscale = {pixscale} arcsec/pixel")
        print(f"Legacy cutout size = {imsize} pixels")

    rows = list(read_catalog(csvfile))

    if not dry_run:
        outdir.mkdir(parents=True, exist_ok=True)

    results = []
    for row in rows:
        result = download_one(
            row,
            outdir=outdir,
            imsize=imsize,
            pixscale=pixscale,
            layer=layer,
            band="r",
            overwrite=overwrite,
            dry_run=dry_run,
            verbose=verbose,
        )
        results.append(result)

        if result["status"] == "failed":
            print(f"WARNING: {row['VFID_1']} failed: {result['message']}")
        elif verbose:
            print(f"{row['VFID_1']}: {result['status']}")

    manifest = outdir / "legacy_r_fits_manifest.csv"
    if not dry_run:
        write_manifest(results, manifest)

    n_downloaded = sum(r["status"] == "downloaded" for r in results)
    n_exists = sum(r["status"] == "exists" for r in results)
    n_failed = sum(r["status"] == "failed" for r in results)

    print(f"galaxies: {len(results)}")
    print(f"downloaded: {n_downloaded}")
    print(f"already present: {n_exists}")
    print(f"failed: {n_failed}")

    if not dry_run:
        print(f"manifest: {manifest}")
        print(f"FITS files: {outdir}")

    return results


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Download Legacy Survey r-band FITS cutouts for the Meerkat "
            "interesting sample."
        )
    )

    parser.add_argument(
        "csvfile",
        nargs="?",
        default=DEFAULT_CSV,
        help=f"Input CSV with columns {', '.join(REQUIRED_COLS)}. Default: {DEFAULT_CSV}",
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: legacy-r-fits next to the input CSV.",
    )

    parser.add_argument(
        "--cutout-arcmin",
        type=float,
        default=10.0,
        help="Cutout width/height in arcmin. Default is 10.",
    )

    parser.add_argument(
        "--pixscale",
        type=float,
        default=0.262,
        help="Legacy cutout pixel scale in arcsec/pixel. Default is 0.262.",
    )

    parser.add_argument(
        "--layer",
        default="ls-dr9",
        help="Legacy Survey layer. Default is ls-dr9.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download FITS files even if they already exist.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print URLs without downloading files.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print download details.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    download_meerkat_legacy_r_fits(
        csvfile=args.csvfile,
        outdir=args.outdir,
        cutout_arcmin=args.cutout_arcmin,
        pixscale=args.pixscale,
        layer=args.layer,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
