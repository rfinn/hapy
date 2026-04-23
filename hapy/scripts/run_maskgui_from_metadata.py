#!/usr/bin/env python

import argparse
import json
import subprocess
from pathlib import Path


def pick_one(*paths):
    for p in paths:
        if p is not None and Path(p).exists():
            return str(p)
    return None


def resolve_gaia_catalog(meta, args, verbose=False):
    """
    Resolve a direct Gaia catalog path using the same logic as run_analysis.

    Returns
    -------
    str or None
        Path to Gaia catalog FITS if found, else None.
    """
    scheme = (meta.get("scheme") or "").lower()
    if args.no_gaia or scheme == "archive":
        if verbose:
            print("Gaia masking disabled")
        return None

    parent_rimage = meta.get("parent_rimage")
    if not parent_rimage or not args.gaia_dir:
        if verbose:
            print(f"Not loading Gaia catalog: parent_rimage={parent_rimage}, gaia_dir={args.gaia_dir}")
        return None

    gaia_path = Path(args.gaia_dir) / parent_rimage.replace(".fits", "-gaia.fits")
    if gaia_path.exists():
        if verbose:
            print(f"Using Gaia catalog: {gaia_path}")
        return str(gaia_path)

    if verbose:
        print(f"WARNING: Gaia catalog not found: {gaia_path}")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Launch run_maskgui using metadata.json from a cutout directory."
    )
    parser.add_argument(
        "--cutout-dir",
        required=True,
        help="Cutout directory containing metadata.json and FITS files",
    )
    parser.add_argument(
        "--gaia-dir",
        default=None,
        help="Directory containing precomputed Gaia catalogs",
    )
    parser.add_argument(
        "--no-gaia",
        action="store_true",
        help="Disable Gaia masking",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print command before running",
    )

    args, extra = parser.parse_known_args()

    cutdir = Path(args.cutout_dir).resolve()
    meta_path = cutdir / "metadata.json"

    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata.json: {meta_path}")

    with open(meta_path) as fh:
        meta = json.load(fh)

    tag = cutdir.name

    rimage = pick_one(
        cutdir / f"{tag}-R.fits",
        cutdir / f"{tag}-r.fits",
    )
    if rimage is None:
        raise FileNotFoundError(f"Could not find r-band image in {cutdir}")

    haimage = pick_one(
        cutdir / f"{tag}-CS-ZP.fits",
        cutdir / f"{tag}-CS.fits",
        cutdir / f"{tag}-Ha.fits",
    )

    weightim = pick_one(
        cutdir / f"{tag}-R.weight.fits",
        cutdir / f"{tag}-r.weight.fits",
        cutdir / f"{tag}-weight.fits",
    )

    gaia_catalog = resolve_gaia_catalog(meta, args, verbose=args.verbose)

    cmd = [
        "run_maskgui",
        "--image", rimage,
        "--title", tag,
        "--objra", str(meta["ra"]),
        "--objdec", str(meta["dec"]),
        "--objsma", str(meta["sma_arcsec"]),
        "--objba", str(meta["ba"]),
        "--objpa", str(meta["pa_deg"]),
    ]

    if haimage is not None:
        cmd += ["--haimage", haimage]

    if weightim is not None:
        cmd += ["--weightim", weightim]

    if gaia_catalog is not None:
        cmd += ["--gaia-catalog", gaia_catalog]
        max_fwhm = max(meta["rimage_fwhm_psf_arcsec"], meta["himage_fwhm_psf_arcsec"])
        gaia_min_radius_arcsec = 4 * max_fwhm
        cmd += ["--gaia-min-radius", gaia_min_radius_arcsec]
    elif args.no_gaia:
        cmd += ["--no-gaia"]

    # pass through any additional args directly to run_maskgui
    cmd += extra

    if args.verbose:
        print("Running:")
        print(" ".join(cmd))

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

    
