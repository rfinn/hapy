#!/usr/bin/env python

"""
summarize_run.py

Print summary statistics for a HAPY pipeline run
using the merged results table.

Usage
-----
python summarize_run.py merged_results.fits
"""
import argparse
import sys
import numpy as np
from astropy.table import Table
from collections import Counter

def find_status_columns(tab):
    """Return boolean-style pipeline status columns."""
    cols = []
    for name in tab.colnames:
        if name.endswith("_OK") or name.endswith("_FLAG"):
            cols.append(name)
    return cols

def summarize(tablefile, scheme):

    tab = Table.read(tablefile)

    n = len(tab)

    print()
    print("HAPY RUN SUMMARY")
    print("----------------")
    print(f"Total galaxies: {n}")
    if scheme == "virgo":
        if "VFID" in tab.colnames:
            print(f"Unique galaxies: {len(set(tab['VFID']))}")
        print()
    elif scheme == "agc":
        if "OBJID" in tab.colnames:
            print(f"Unique galaxies: {len(set(tab['OBJID']))}")
        print()

    # ---------- pipeline flags ----------
    print("Pipeline completion")
    print("-------------------")    
    flags = find_status_columns(tab)


    for f in flags:
        try:
            col = np.array(tab[f], dtype=bool)
        except Exception:
            print(f"{f:18s}: could not interpret as boolean")
            continue

        if ('SM_FLAG' in f) or ('SM_SERSIC_FLAG' in f):
            ntrue = np.sum(col)
            nfalse = np.sum(~col)
            pct = 100 * nfalse / n
            
        else:
            ntrue = np.sum(col)
            nfalse = np.sum(~col)
            pct = 100 * ntrue / n

        print(f"{f:18s}: {ntrue:4d} OK  | {nfalse:4d} FAIL  ({pct:5.1f}%)")

    print()



    print()
    if "R_PROFILE_OK" in tab.colnames and "H_PROFILE_OK" in tab.colnames:
        both = np.logical_and(tab["R_PROFILE_OK"], tab["H_PROFILE_OK"])
        pct = 100 * np.sum(both) / n
        print(f"{'PROFILES_BOTH':18s}: {np.sum(both):4d} ({pct:5.1f}%)")
        print()

    if "R_SM_FLAG" in tab.colnames and "H_SM_FLAG" in tab.colnames:
        both = np.logical_and(tab["R_SM_FLAG"], tab["H_SM_FLAG"])
        print(f"{'STATMORPH_BOTH':18s}: {np.sum(both):4d}")

    print()
    # ---------- STATUS ----------
    if "STATUS" in tab.colnames:

        print("STATUS counts")
        print("-------------")

        counts = Counter(tab["STATUS"])

        for k, v in counts.items():
            print(f"{str(k):12s}: {v}")

        print()

    # ---------- STAGE ----------
    if "STAGE" in tab.colnames:

        print("STAGE counts")
        print("------------")

        counts = Counter(tab["STAGE"])

        for k, v in counts.items():
            print(f"{str(k):12s}: {v}")

        print()

    # ---------- runtime ----------
    runtime_cols = [
        "MASK_SEC",
        "PHOT_SEC",
        "GALFIT_SEC",
        "TOTAL_SEC",
    ]

    print("Runtime medians (sec)")
    print("---------------------")

    for col in runtime_cols:
        if col in tab.colnames:
            vals = tab[col]
            vals = vals[np.isfinite(vals)]

            if len(vals) > 0:
                print(f"{col:12s}: {np.median(vals):.2f}")

    print()
    bad = tab[~tab["PHOT_OK"]]
    print(f"Number with bad phot = {len(bad)}")

def main():
    parser = argparse.ArgumentParser(
        description="Merge per-galaxy ECSV results into a single FITS table."
    )
    parser.add_argument(
        "--infile",
        required=True,
        help="file to pass in, e.g. merged_results.fits."
    )

    parser.add_argument(
        "--scheme",
        choices=["virgo", "agc"],
        required=True,
        help="Pipeline stage whose results should be merged."
    )

    args = parser.parse_args()
    summarize(args.infile, args.scheme)


if __name__ == "__main__":
    main()
