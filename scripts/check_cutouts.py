#!/usr/bin/env python

"""
check_cutouts.py

Verify that get_cutouts produced expected cutout directories
and that required files exist.

Example
-------
check_cutouts ready_coadds.txt cutouts
"""

import argparse
from pathlib import Path

from pathlib import Path

def coadd_key_from_filename(path):
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
    if parts[-1] in {"r", "R"}:
        parts = parts[:-1]

    if len(parts) < 3:
        return None

    return "-".join(parts[-3:])


def cutout_key_from_dirname(dirname):
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

def check_cutouts(coadd_list, cutout_dir):

    coadd_list = Path(coadd_list)
    cutout_dir = Path(cutout_dir)

    if not coadd_list.exists():
        raise FileNotFoundError(coadd_list)

    if not cutout_dir.exists():
        raise FileNotFoundError(cutout_dir)

    # --------------------------------------------------
    # Read coadd list
    # --------------------------------------------------
    with open(coadd_list) as f:
        coadds = [line.strip() for line in f if line.strip()]

    # --------------------------------------------------
    # Collect cutout directories
    # --------------------------------------------------
    cutdirs = [d for d in cutout_dir.iterdir() if d.is_dir()]

    # --------------------------------------------------
    # Check coadds with no cutouts
    # --------------------------------------------------
    no_cutouts = []

    for f in coadds:
        base = Path(f).stem
        #matches = [d for d in cutdirs if d.name.startswith(base)]
        cutdir_map = {}
        for d in cutdirs:
            key = cutout_key_from_dirname(d.name)
            if key is not None:
                cutdir_map.setdefault(key, []).append(d)

        no_cutouts = []

        for f in coadds:
            key = coadd_key_from_filename(f)
            if key is None or key not in cutdir_map:
                no_cutouts.append(f)

        
        if len(matches) == 0:
            no_cutouts.append(f)

    # --------------------------------------------------
    # Check required files
    # --------------------------------------------------
    missing_r = []
    missing_cs = []

    for d in cutdirs:

        tag = d.name
        rfile = d / f"{tag}-R.fits"
        csfile = d / f"{tag}-CS.fits"

        if not rfile.exists():
            missing_r.append(d)

        if not csfile.exists():
            missing_cs.append(d)

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------
    print("\nCUTOUT SUMMARY")
    print("--------------")
    print(f"Input coadds:           {len(coadds)}")
    print(f"Cutout directories:     {len(cutdirs)}")
    print(f"Coadds with no cutouts: {len(no_cutouts)}")
    print(f"Dirs missing R:         {len(missing_r)}")
    print(f"Dirs missing CS:        {len(missing_cs)}")

    if no_cutouts:
        print("\nExample missing coadds:")
        for f in no_cutouts[:5]:
            print("  ", f)

    if missing_r:
        print("\nExample dirs missing R:")
        for d in missing_r[:5]:
            print("  ", d)

    if missing_cs:
        print("\nExample dirs missing CS:")
        for d in missing_cs[:5]:
            print("  ", d)


def main():

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
