#!/usr/bin/env python

import argparse

from hapy.hatools.mstar import make_mstar_map
from pathlib import Path
import os
import numpy as np
from astropy.table import Table


def get_velocity_from_scheme(cutout_dir, scheme="generic", tabledir=None):
    cutout_dir = Path(cutout_dir)
    vfid = cutout_dir.name.split("-")[0]

    if scheme == "virgo":
        if tabledir is None:
            tabledir = Path(os.getenv("HOME")) / "research/Virgo/tables-north/v2"
        else:
            tabledir = Path(tabledir)

        env = Table.read(tabledir / "vf_v2_environment.fits")
        match = np.where(env["VFID"] == vfid)[0]

        if len(match) == 0:
            raise ValueError(f"Could not find {vfid} in vf_v2_environment.fits")

        return float(env["Vcosmic"][match[0]])

    return None

def main():
    parser = argparse.ArgumentParser(
        description="Make stellar-mass map from Legacy r and g-r images."
    )

    parser.add_argument("cutout_dir", help="HAPY cutout directory")
    parser.add_argument("--velocity", type=float, default=None, help="Velocity for distance, km/s")
    parser.add_argument("--distance-mpc", type=float, default=None, help="Distance in Mpc")
    parser.add_argument("--smoothing", type=int, default=15, help="Boxcar smoothing size")
    parser.add_argument("--output-suffix", default="vcosmic", help="Output filename suffix")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output")
    parser.add_argument(
        "--scheme",
        choices=["generic", "virgo", "agc"],
        default="generic",
        help="Filename parsing scheme for coadd images.",
    )

    parser.add_argument(
        "--tabledir",
        default=None,
        help="Directory containing Virgo tables; used for --scheme virgo.",
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

    

    out = make_mstar_map(
        args.cutout_dir,
        velocity_kms=args.velocity,
        distance_mpc=args.distance_mpc,
        smoothing=args.smoothing,
        output_suffix=args.output_suffix,
        overwrite=args.overwrite,
    )

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
