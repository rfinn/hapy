#!/usr/bin/env python

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


def make_simple_mask(infile, outfile=None, overwrite=False):
    infile = Path(infile)

    if outfile is None:
        outfile = infile.with_name(infile.name.replace(".fits", ".weight.fits"))
    else:
        outfile = Path(outfile)

    data, header = fits.getdata(infile, header=True)

    mask = np.zeros(data.shape, dtype=np.uint8)
    mask[np.isfinite(data) & (data != 0)] = 1

    header["MASKTYPE"] = ("SIMPLE", "1=finite nonzero science pixel")
    header["MASKSRC"] = (infile.name, "Source coadd image")
    header["MASKVAL"] = ("1 good, 0 bad", "Mask convention")

    fits.PrimaryHDU(data=mask, header=header).writeto(outfile, overwrite=overwrite)

    print(f"Wrote {outfile}")
    print(f"Good fraction: {mask.mean():.4f}")

    return outfile


def main():
    parser = argparse.ArgumentParser(
        description="Make a simple valid-data mask from a coadd image."
    )
    parser.add_argument("infile", help="Input coadd FITS image")
    parser.add_argument(
        "-o", "--outfile",
        default=None,
        help="Output mask FITS file. Default: <input>.simple-mask.fits",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file",
    )

    args = parser.parse_args()

    make_simple_mask(
        infile=args.infile,
        outfile=args.outfile,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
