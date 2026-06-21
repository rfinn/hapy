#!/usr/bin/env python

import argparse
import os
from pathlib import Path

import numpy as np
from astropy.table import Table

from hapy.hatools.sfr import make_sfr_map


def get_velocity_from_scheme(cutout_dir, scheme="generic", tabledir=None):
    cutout_dir = Path(cutout_dir)
    vfid = cutout_dir.name.split("-")[0]

    if scheme == "virgo":
        if tabledir is None:
            tabledir = Path(os.getenv("HOME")) / "research/Virgo/tables-north/v2"
        else:
            tabledir = Path(tabledir)

        envfile = tabledir / "vf_v2_environment.fits"
        env = Table.read(envfile)

        match = np.where(env["VFID"] == vfid)[0]
        if len(match) == 0:
            raise ValueError(f"Could not find {vfid} in {envfile}")

        return float(env["Vcosmic"][match[0]])

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Make linear SFR map from HAPY continuum-subtracted Halpha image."
    )

    parser.add_argument("cutout_dir", help="HAPY cutout directory")

    parser.add_argument(
        "--scheme",
        choices=["generic", "virgo", "agc"],
        default="generic",
        help="Dataset scheme; virgo uses vf_v2_environment.fits Vcosmic.",
    )

    parser.add_argument(
        "--tabledir",
        default=None,
        help="Directory containing Virgo tables; used for --scheme virgo.",
    )

    parser.add_argument(
        "--velocity",
        type=float,
        default=None,
        help="Velocity for distance, km/s. Overrides metadata unless --scheme virgo is used.",
    )

    parser.add_argument(
        "--distance-mpc",
        type=float,
        default=None,
        help="Distance in Mpc. Overrides velocity-based distance.",
    )

    parser.add_argument(
        "--prefer",
        choices=["csgr", "cszp"],
        default="csgr",
        help="Preferred continuum-subtracted image.",
    )

    parser.add_argument(
        "--output-suffix",
        default="",
        help="Optional output suffix, e.g. csgr or cszp.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output.",
    )

    args = parser.parse_args()

    velocity = args.velocity

    scheme_velocity = get_velocity_from_scheme(
        args.cutout_dir,
        scheme=args.scheme,
        tabledir=args.tabledir,
    )

    if scheme_velocity is not None:
        velocity = scheme_velocity

    out = make_sfr_map(
        args.cutout_dir,
        velocity_kms=velocity,
        distance_mpc=args.distance_mpc,
        prefer=args.prefer,
        output_suffix=args.output_suffix,
        overwrite=args.overwrite,
    )

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
