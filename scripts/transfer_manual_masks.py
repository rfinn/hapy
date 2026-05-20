#!/usr/bin/env python

"""
Hydrid INT coadds - I created a new hapy-output-20260519

- need to copy or reproject the existing manual masks from hapy-output-20260417

"""

from pathlib import Path
import argparse
import shutil

import numpy as np
from astropy.io import fits
from reproject import reproject_interp


INSTRUMENTS_COPY = ("BOK", "HDI", "MOS")


def find_one(pattern, directory):
    matches = sorted(directory.glob(pattern))
    if len(matches) == 0:
        return None
    if len(matches) > 1:
        print(f"WARNING: multiple matches for {pattern} in {directory}; using {matches[0].name}")
    return matches[0]


def reproject_mask_to_ha(old_mask, new_ha, outmask, overwrite=False, dry_run=False):
    print(f"REPROJECT INT: {old_mask} -> {outmask}")
    print(f"  target grid: {new_ha}")

    if dry_run:
        return "dryrun"

    if outmask.exists() and not overwrite:
        print(f"  exists, skipping: {outmask}")
        return "exists"

    mask_data, mask_hdr = fits.getdata(old_mask, header=True)
    target_hdr = fits.getheader(new_ha)

    reproj, footprint = reproject_interp(
        (mask_data, mask_hdr),
        target_hdr,
        order="nearest-neighbor",
    )

    #newmask = np.zeros(reproj.shape, dtype=np.uint8)
    #newmask[np.isfinite(reproj) & (reproj > 0)] = 1

    newmask = np.zeros(reproj.shape, dtype=np.int16)

    good = np.isfinite(reproj) & (reproj > 0)
    newmask[good] = np.rint(reproj[good]).astype(np.int16)

    out_hdr = target_hdr.copy()
    out_hdr["MASKMAN"] = (True, "Manual mask")
    out_hdr["MASKREP"] = (True, "Manual mask reprojected")
    out_hdr["MASKSRC"] = (old_mask.name, "Source manual mask")
    out_hdr["MASKTARG"] = (new_ha.name, "Target image grid")
    out_hdr["MASKTYPE"] = ("manual", "Mask type")
    out_hdr["MASKVAL"] = ("0=unmasked, >0=masked IDs", "Mask convention")
    out_hdr["MASKLAB"] = (True, "Mask preserves object/region labels")

    fits.PrimaryHDU(data=newmask, header=out_hdr).writeto(outmask, overwrite=True)

    return "reprojected"


def copy_mask(old_mask, outmask, overwrite=False, dry_run=False):
    print(f"COPY: {old_mask} -> {outmask}")

    if dry_run:
        return "dryrun"

    if outmask.exists() and not overwrite:
        print(f"  exists, skipping: {outmask}")
        return "exists"

    outmask.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(old_mask, outmask)

    return "copied"


def process_tag(src_tag_dir, dst_root, overwrite=False, dry_run=False):
    tag = src_tag_dir.name
    dst_tag_dir = dst_root / tag

    old_mask = find_one("*-mask-manual.fits", src_tag_dir)
    if old_mask is None:
        return tag, "no_manual_mask"

    if not dst_tag_dir.exists():
        return tag, "missing_dst_tag_dir"

    outmask = dst_tag_dir / old_mask.name

    if "INT" in tag:
        new_ha = (
            find_one("*-Ha.fits", dst_tag_dir)
            or find_one("*-Halpha.fits", dst_tag_dir)
            or find_one("*-Ha6657.fits", dst_tag_dir)
        )

        if new_ha is None:
            return tag, "missing_new_ha"

        status = reproject_mask_to_ha(
            old_mask=old_mask,
            new_ha=new_ha,
            outmask=outmask,
            overwrite=overwrite,
            dry_run=dry_run,
        )
        return tag, status

    if any(inst in tag for inst in INSTRUMENTS_COPY):
        status = copy_mask(
            old_mask=old_mask,
            outmask=outmask,
            overwrite=overwrite,
            dry_run=dry_run,
        )
        return tag, status

    return tag, "unknown_instrument"


def main():
    parser = argparse.ArgumentParser(
        description="Transfer manual masks from old HAPY cutouts to new cutouts."
    )
    parser.add_argument(
        "--src-root",
        default="/data-pool/Halpha/hapy-output-20260417/cutouts",
        help="Source cutouts root",
    )
    parser.add_argument(
        "--dst-root",
        default="/data-pool/Halpha/hapy-output-20260519/cutouts",
        help="Destination cutouts root",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing destination manual masks",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned actions without writing files",
    )

    args = parser.parse_args()

    src_root = Path(args.src_root)
    dst_root = Path(args.dst_root)

    src_dirs = sorted([p for p in src_root.iterdir() if p.is_dir()])

    counts = {}

    for src_tag_dir in src_dirs:
        tag, status = process_tag(
            src_tag_dir,
            dst_root=dst_root,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )
        counts[status] = counts.get(status, 0) + 1
        print(f"{tag}: {status}")

    print("\nSummary")
    for key in sorted(counts):
        print(f"  {key}: {counts[key]}")


if __name__ == "__main__":
    main()
