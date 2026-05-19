#!/usr/bin/env python3
"""
Copy PSF-related header cards from psf images to rebuilt INT coadds.

Use case:
- new rebuilt INT images live in /data-pool/Halpha/coadds-v20260518/
- psf images live in /data-pool/Halpha/psf-images-v20260518/
- filenames are assumed to match exactly

Cards copied:
    FWHM
    SEFWHM
    STD
    PSFSTD
    OVERSAMP
"""

from pathlib import Path
from astropy.io import fits
import argparse
import shutil
import sys

PSF_KEYS = ["FWHM", "SEFWHM", "STD", "PSFSTD", "OVERSAMP"]


def copy_header_cards(
    new_file: Path,
    old_file: Path,
    keys=PSF_KEYS,
    dry_run: bool = False,
    backup: bool = False,
    verbose: bool = True,
):
    """Copy selected header cards from old_file to new_file."""
    if not old_file.exists():
        print(f"WARNING: missing source file: {old_file}")
        return False

    with fits.open(old_file) as hdul_old:
        old_header = hdul_old[0].header

    found = {}
    missing = []
    for key in keys:
        if key in old_header:
            found[key] = (old_header[key], old_header.comments[key])
        else:
            missing.append(key)

    if verbose:
        print(f"\n{new_file.name}")
        print(f"  source: {old_file}")
        if found:
            print(f"  copying: {', '.join(found.keys())}")
        if missing:
            print(f"  missing in source: {', '.join(missing)}")

    if dry_run:
        return True

    if backup:
        backup_path = new_file.with_suffix(new_file.suffix + ".bak")
        if not backup_path.exists():
            shutil.copy2(new_file, backup_path)
            if verbose:
                print(f"  backup: {backup_path}")

    with fits.open(new_file, mode="update") as hdul_new:
        hdr = hdul_new[0].header
        for key, (value, comment) in found.items():
            hdr[key] = (value, comment)
        hdul_new.flush()

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Copy PSF header cards from old INT coadds to rebuilt INT coadds."
    )
    parser.add_argument(
        "--newdir",
        required=True,
        help="Directory containing coadds (e.g. /data-pool/Halpha/coadds-v20260518/)",
    )
    parser.add_argument(
        "--psfdir",
        required=True,
        help="Directory containing matching psf images (e.g. /data-pool/Halpha/psf-images-v20260518/)",
    )
    parser.add_argument(
        "--pattern",
        default="*INT*.fits",
        help='Glob pattern for files to update (default: "*INT*.fits")',
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be updated without modifying files",
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Write a .bak copy of each updated file before editing",
    )
    args = parser.parse_args()

    newdir = Path(args.newdir)
    psfdir = Path(args.psfdir)

    if not newdir.is_dir():
        print(f"ERROR: newdir does not exist: {newdir}")
        sys.exit(1)
    if not psfdir.is_dir():
        print(f"ERROR: olddir does not exist: {psfdir}")
        sys.exit(1)

    files = sorted(newdir.glob(args.pattern))
    if not files:
        print(f"No files matched {args.pattern} in {newdir}")
        sys.exit(0)

    n_ok = 0
    n_missing = 0

    for new_file in files:
        psf_file = psfdir / new_file.name.replace(".fits","-psf.fits")
        ok = copy_header_cards(
            new_file,
            psf_file,
            dry_run=args.dry_run,
            backup=args.backup,
        )
        if psf_file.exists():
            n_ok += int(ok)
        else:
            n_missing += 1

    print("\nSummary")
    print(f"  matched files updated: {n_ok}")
    print(f"  missing source files:  {n_missing}")
    print(f"  dry run:               {args.dry_run}")


if __name__ == "__main__":
    main()
