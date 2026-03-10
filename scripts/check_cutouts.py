#!/usr/bin/env python

"""
check_cutouts.py

Verify that get_cutouts produced expected cutout directories
and that required files exist.

Matching logic
--------------
Coadds and cutout directories are linked using the last three fields:

    <INSTRUMENT>-<DATE>-<TARGET>

Examples
--------
Coadd:
    VF-126.291+27.988-HDI-20180313-p004-R.fits
    -> HDI-20180313-p004

Cutout dir:
    VFID2995-NGC4793-BOK-20210416-VFID2997
    -> BOK-20210416-VFID2997

This is robust even when the NED name contains extra '-'.

Usage
-----
check_cutouts ready_coadds.txt cutouts

or

python check_cutouts.py ready_coadds.txt cutouts
"""

from __future__ import annotations

import argparse
from pathlib import Path


def coadd_key_from_filename(path: str | Path) -> str | None:
    """
    Build a matching key from a coadd image filename.

    Example
    -------
    VF-126.291+27.988-BOK-20210416-VFID2997-R.fits
    -> BOK-20210416-VFID2997
    """
    stem = Path(path).stem
    parts = stem.split("-")

    # remove trailing band token if present
    if len(parts) > 0 and parts[-1] in {"r", "R"}:
        parts = parts[:-1]

    if len(parts) < 3:
        return None

    return "-".join(parts[-3:])


def cutout_key_from_dirname(dirname: str | Path) -> str | None:
    """
    Build a matching key from a cutout directory name.

    Example
    -------
    VFID2995-NGC4793-BOK-20210416-VFID2997
    -> BOK-20210416-VFID2997

    Works even if NED name contains extra '-'.
    """
    parts = Path(dirname).name.split("-")

    if len(parts) < 4:
        return None

    return "-".join(parts[-3:])


def check_cutouts(coadd_list: str | Path, cutout_dir: str | Path) -> None:
    coadd_list = Path(coadd_list)
    cutout_dir = Path(cutout_dir)

    if not coadd_list.exists():
        raise FileNotFoundError(f"Could not find coadd list: {coadd_list}")

    if not cutout_dir.exists():
        raise FileNotFoundError(f"Could not find cutout directory: {cutout_dir}")

    # --------------------------------------------------
    # Read coadd list
    # --------------------------------------------------
    with open(coadd_list) as f:
        coadds = [line.strip() for line in f if line.strip()]

    # --------------------------------------------------
    # Collect cutout directories
    # --------------------------------------------------
    cutdirs = sorted([d for d in cutout_dir.iterdir() if d.is_dir()])

    # map cutout key -> list of matching directories
    cutdir_map: dict[str, list[Path]] = {}
    bad_cutdir_names: list[Path] = []

    for d in cutdirs:
        key = cutout_key_from_dirname(d.name)
        if key is None:
            bad_cutdir_names.append(d)
            continue
        cutdir_map.setdefault(key, []).append(d)

    # --------------------------------------------------
    # Check which coadds produced no cutout directories
    # --------------------------------------------------
    no_cutouts: list[str] = []
    bad_coadd_names: list[str] = []

    for f in coadds:
        key = coadd_key_from_filename(f)
        if key is None:
            bad_coadd_names.append(f)
            continue

        if key not in cutdir_map:
            no_cutouts.append(f)

    # --------------------------------------------------
    # Check cutout dirs for key files
    # --------------------------------------------------
    missing_r: list[Path] = []
    missing_cs: list[Path] = []

    for d in cutdirs:
        tag = d.name
        rfile = d / f"{tag}-R.fits"
        csfile = d / f"{tag}-CS-ZP.fits"

        if not rfile.exists():
            missing_r.append(d)

        if not csfile.exists():
            missing_cs.append(d)

    # --------------------------------------------------
    # Summary report
    # --------------------------------------------------
    print("\nCUTOUT SUMMARY")
    print("--------------")
    print(f"Input coadds:              {len(coadds)}")
    print(f"Cutout directories:        {len(cutdirs)}")
    print(f"Coadds with no cutouts:    {len(no_cutouts)}")
    print(f"Cutout dirs missing R:     {len(missing_r)}")
    print(f"Cutout dirs missing CS:    {len(missing_cs)}")
    print(f"Bad coadd names:           {len(bad_coadd_names)}")
    print(f"Bad cutout dir names:      {len(bad_cutdir_names)}")

    if no_cutouts:
        print("\nFirst few coadds with no cutouts:")
        for f in no_cutouts[:10]:
            print(f"  {f}")

    if missing_r:
        print("\nFirst few cutout dirs missing R:")
        for d in missing_r[:10]:
            print(f"  {d}")

    if missing_cs:
        print("\nFirst few cutout dirs missing CS:")
        for d in missing_cs[:10]:
            print(f"  {d}")

    if bad_coadd_names:
        print("\nFirst few coadd filenames that could not be parsed:")
        for f in bad_coadd_names[:10]:
            print(f"  {f}")

    if bad_cutdir_names:
        print("\nFirst few cutout directory names that could not be parsed:")
        for d in bad_cutdir_names[:10]:
            print(f"  {d}")

    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check cutout production from get_cutouts"
    )

    parser.add_argument(
        "coadd_list",
        help="Text file listing coadd images used for get_cutouts",
    )

    parser.add_argument(
        "cutout_dir",
        help="Directory containing cutout subdirectories",
    )

    args = parser.parse_args()

    check_cutouts(args.coadd_list, args.cutout_dir)


if __name__ == "__main__":
    main()
