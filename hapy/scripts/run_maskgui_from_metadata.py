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
        help="Optional Gaia catalog directory to pass through",
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

    if args.gaia_dir is not None:
        cmd += ["--gaia-dir", args.gaia_dir]

    if args.no_gaia:
        cmd += ["--no-gaia"]

    # pass through any extra args directly to run_maskgui
    cmd += extra

    if args.verbose:
        print("Running:")
        print(" ".join(cmd))

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
