#!/usr/bin/env python

"""

RUN COMMAND:

python ~/github/hapy/scripts/make_vfs_hapy_rowmatched.py merged_results_virgo_20260514_with_best_duplicate.fits

OUTPUT:
merged_results_virgo_20260514_with_best_duplicate_vfs_rowmatched.fits

"""

import argparse
import numpy as np
from pathlib import Path
from astropy.table import Table
from hapy.utils.results_table import prepare_analysis_table

def safe_str(x):
    try:
        return str(x).strip()
    except Exception:
        return ""


def infer_vfid_from_tag(tag):
    """
    Extract VFID0000-style ID from a HAPY TAG.
    """
    parts = safe_str(tag).split("-")
    for part in parts:
        if part.startswith("VFID"):
            return part
    return ""


def get_vfid_col(tab):
    """
    Prefer DUP_GALID if present, otherwise VFID, otherwise parse TAG.
    """
    if "DUP_GALID" in tab.colnames:
        return np.array([safe_str(x) for x in tab["DUP_GALID"]])

    if "VFID" in tab.colnames:
        return np.array([safe_str(x) for x in tab["VFID"]])

    if "TAG" in tab.colnames:
        return np.array([infer_vfid_from_tag(x) for x in tab["TAG"]])

    raise ValueError("Could not find DUP_GALID, VFID, or TAG column.")


def make_vfs_hapy_rowmatched(tab, n_vfs=6780):
    """
    Make 6780-row HAPY table row-matched to VFS VFID order.

    Keeps only rows with USE_FOR_SCIENCE=True if that column exists.
    HAPY column names are preserved.
    Adds:
      HAPY_HAS_OBS
      HAPY_NOBS
    """
    all_vfids = get_vfid_col(tab)

    if "USE_FOR_SCIENCE" in tab.colnames:
        use_tab = tab[np.array(tab["USE_FOR_SCIENCE"], dtype=bool)]
    else:
        print("WARNING: USE_FOR_SCIENCE not found; using all rows.")
        use_tab = tab

    use_vfids = get_vfid_col(use_tab)

    rows = []

    # columns to copy from HAPY table
    copy_cols = list(tab.colnames)

    for i in range(n_vfs):
        vfid = f"VFID{i:04d}"

        all_matches = np.where(all_vfids == vfid)[0]
        use_matches = np.where(use_vfids == vfid)[0]

        out = {
            "VFID": vfid,
            "HAPY_HAS_OBS": len(all_matches) > 0,
            "HAPY_NOBS": len(all_matches),
        }

        if len(use_matches) == 1:
            row = use_tab[use_matches[0]]
            for col in copy_cols:
                out[col] = row[col]

        elif len(use_matches) > 1:
            print(f"WARNING: {vfid} has {len(use_matches)} USE_FOR_SCIENCE rows; using first.")
            row = use_tab[use_matches[0]]
            for col in copy_cols:
                out[col] = row[col]

        rows.append(out)

    return Table(rows=rows)


def main():
    parser = argparse.ArgumentParser(
        description="Create a VFS-row-matched HAPY table ordered by VFID."
    )

    parser.add_argument(
        "hapy_table",
        help="HAPY merged_results table, preferably *_with_best_duplicate.fits",
    )

    parser.add_argument(
        "--outfile",
        default=None,
        help="Output FITS file. Default: <input_stem>_vfs_rowmatched.fits",
    )

    parser.add_argument(
        "--n-vfs",
        type=int,
        default=6780,
        help="Number of VFS rows/VFIDs. Default: 6780.",
    )

    args = parser.parse_args()

    infile = Path(args.hapy_table)
    tab = Table.read(infile)

    # add additional columns
    tab = prepare_analysis_table(tab)

    # get best duplicate    
    if 'BEST_DUPLICATE' in tab.colnames:
        tab = tab[tab['BEST_DUPLICATE']]
    else:
        print("WARNING: BEST_DUPLICATE column was not found in", args.hapy_table)
        print("PLEASE SELECT BEST DUPLICATE")
    


    
    outtab = make_vfs_hapy_rowmatched(tab, n_vfs=args.n_vfs)

    from datetime import datetime
    today = datetime.now().strftime("%Y%m%d")

    name_parts = str(infile.stem).split('_')
    for t in name_parts:
        if t.startswith('2026'):
            merge_date = t
            break
    #merge_date = t[-1].replace(".fits","")
    if args.outfile is None:
        outfile = infile.parent / f"vf_v2_hapy_v{merge_date}_{today}.fits"
    else:
        outfile = Path(args.outfile)

    outtab.write(outfile, overwrite=True)

    print(f"Wrote row-matched HAPY table:")
    print(f"  {outfile}")
    print(f"Rows: {len(outtab)}")
    print(f"HAPY detections: {np.sum(outtab['HAPY_HAS_OBS'])}")


if __name__ == "__main__":
    main()
