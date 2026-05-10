#!/usr/bin/env python

"""
Reproject Legacy Survey g/r/z images onto the Halpha cutout footprint.

Examples
--------
Single cutout:

    python ~/github/hapy/hapy/scripts/make_legacy_reprojections.py cutouts/<tag>

From the directory containing cutouts/:

    find cutouts -mindepth 1 -maxdepth 1 -type d -name 'VFID*' | sort > cutout_list.txt

    mkdir -p logs_legacy_reproject

    parallel --bar -j 12 'python ~/github/hapy/hapy/scripts/make_legacy_reprojections.py {} > logs_legacy_reproject/{/}.log 2>&1' :::: cutout_list.txt


"""

import argparse
import sys

from hapy.hatools.reproject_legacy import make_legacy_reprojections


def main():
    parser = argparse.ArgumentParser(
        description="Reproject Legacy g/r/z images onto the Halpha cutout footprint."
    )

    parser.add_argument(
        "cutdir",
        help="HAPY cutout directory, e.g. cutouts/VFID0377-IC1210-BOK-20210414-VFID0422",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing *-ha.fits reprojection products.",
    )

    args = parser.parse_args()

    try:
        outputs = make_legacy_reprojections(
            args.cutdir,
            overwrite=args.overwrite,
        )
    except Exception as err:
        print(f"FAILED {args.cutdir}: {err}", file=sys.stderr)
        sys.exit(1)

    for outfile in outputs:
        print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
