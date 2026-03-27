#!/usr/bin/env python

from pathlib import Path
import numpy as np
from astropy.io import fits


def read_mask(path):
    """Read FITS mask and return boolean array."""
    with fits.open(path) as hdul:
        data = hdul[0].data

    if data is None:
        raise ValueError(f"No data found in {path}")

    # Treat any nonzero value as masked
    mask = np.asarray(data) != 0
    return mask


def combine_masks(mask1_path, mask2_path, output_path):
    """Combine two masks with logical OR and write to FITS."""

    mask1 = read_mask(mask1_path)
    mask2 = read_mask(mask2_path)

    if mask1.shape != mask2.shape:
        raise ValueError(
            f"Mask shapes do not match: {mask1.shape} vs {mask2.shape}"
        )

    combined = mask1 | mask2  # logical OR

    # Convert to int (0/1) for FITS output
    combined_int = combined.astype(np.uint8)

    # Copy header from first mask (optional but useful for WCS)
    with fits.open(mask1_path) as hdul:
        header = hdul[0].header

    fits.writeto(output_path, combined_int, header=header, overwrite=True)

    print(f"Wrote combined mask: {output_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Combine two FITS masks using logical OR"
    )
    parser.add_argument("mask1", help="First mask FITS file")
    parser.add_argument("mask2", help="Second mask FITS file")
    parser.add_argument(
        "-o", "--output",
        help="Output mask filename (default: derived name)",
        default=None
    )

    args = parser.parse_args()

    mask1 = Path(args.mask1)
    mask2 = Path(args.mask2)

    if args.output is None:
        # Example: n4567.mask.fits
        root = mask1.name.replace("ha.mask.fits", "").replace("r.mask.fits", "")
        output = mask1.with_name(f"{root}.mask.fits")
    else:
        output = Path(args.output)

    combine_masks(mask1, mask2, output)
