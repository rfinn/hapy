#!/usr/bin/env python

"""
fetch_legacy_cutouts.py

Download Legacy Survey grz FITS cutouts and jpg images for HAPY cutout directories.

Outputs are written to:
    cutouts/<tag>/legacy/
"""

import argparse
import json
from pathlib import Path

from astropy.io import fits

from hapy.imagetools.downloads import get_legacy_images


def read_metadata(cutout_dir):
    meta_file = Path(cutout_dir) / "metadata.json"
    if not meta_file.exists():
        raise FileNotFoundError(f"Missing metadata.json in {cutout_dir}")

    with open(meta_file) as fh:
        meta = json.load(fh)

    # adjust these keys if your metadata uses different names
    ra = meta.get("ra", None)
    dec = meta.get("dec", None)

    if ra is None or dec is None:
        raise ValueError(f"Could not find ra/dec in {meta_file}")

    return float(ra), float(dec)


def get_cutout_imsize(cutout_dir, tag, fallback=512):
    """
    Determine Legacy cutout size in pixels based on the HAPY R-band cutout.
    Uses the larger image dimension.
    """
    rfiles = sorted(Path(cutout_dir).glob(f"{tag}*-R.fits"))
    if len(rfiles) == 0:
        rfiles = sorted(Path(cutout_dir).glob(f"{tag}*-r.fits"))

    if len(rfiles) == 0:
        return fallback

    data = fits.getdata(rfiles[0])
    ny, nx = data.shape
    return int(max(nx, ny))


def fetch_one(cutout_dir, pixscale=1.0, layer="ls-dr10", verbose=False):
    cutout_dir = Path(cutout_dir).resolve()
    tag = cutout_dir.name
    legacy_dir = cutout_dir / "legacy"
    legacy_dir.mkdir(parents=True, exist_ok=True)

    ra, dec = read_metadata(cutout_dir)
    imsize = get_cutout_imsize(cutout_dir, tag)

    if verbose:
        print(f"\n{tag}")
        print(f"  ra, dec = {ra:.6f}, {dec:.6f}")
        print(f"  imsize  = {imsize}")
        print(f"  outdir  = {legacy_dir}")

    for band in ["g", "r", "z"]:
        get_legacy_images(
            ra=ra,
            dec=dec,
            galid=tag,
            pixscale=pixscale,
            imsize=imsize,
            band=band,
            makeplots=False,
            subfolder=str(legacy_dir),
            verbose=verbose,
            layer=layer,
        )


def iter_cutout_dirs(cutout_root):
    cutout_root = Path(cutout_root)
    for subdir in sorted(cutout_root.iterdir()):
        if not subdir.is_dir():
            continue
        if (subdir / "metadata.json").exists():
            yield subdir


def main():
    parser = argparse.ArgumentParser(
        description="Download Legacy Survey images for HAPY cutout directories."
    )
    parser.add_argument("--cutout-dir", default=None,
                        help="Single cutout directory")
    parser.add_argument("--cutout-root", default=None,
                        help="Root directory containing cutout subdirectories")
    parser.add_argument("--pixscale", type=float, default=0.262,
                        help="Legacy cutout pixel scale in arcsec/pixel.  Default is 0.262")
    parser.add_argument("--layer", default="ls-dr10",
                        help="Legacy Survey viewer layer")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.cutout_dir is not None:
        fetch_one(args.cutout_dir, pixscale=args.pixscale,
                  layer=args.layer, verbose=args.verbose)

    elif args.cutout_root is not None:
        for cutout_dir in iter_cutout_dirs(args.cutout_root):
            try:
                fetch_one(cutout_dir, pixscale=args.pixscale,
                          layer=args.layer, verbose=args.verbose)
            except Exception as e:
                print(f"WARNING: failed for {cutout_dir}: {e}")

    else:
        raise SystemExit("Provide --cutout-dir or --cutout-root")


if __name__ == "__main__":
    main()
