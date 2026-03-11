#!/usr/bin/env python

"""
summarize_run.py

Print summary statistics for a HAPY pipeline run
using the merged results table.

Usage
-----
python summarize_run.py merged_results.fits
"""

import sys
import numpy as np
from astropy.table import Table
from collections import Counter


def summarize(tablefile):

    tab = Table.read(tablefile)

    n = len(tab)

    print()
    print("HAPY RUN SUMMARY")
    print("----------------")
    print(f"Total galaxies: {n}")
    if "VFID" in tab.colnames:
        print(f"Unique galaxies: {len(set(tab['VFID']))}")
    print()

    # ---------- pipeline flags ----------
    flags = [
        "MASK_OK",
        "PHOT_OK",
        "PSF_OK",
        "GAL_NC_OK",
        "GAL_CV_OK",
        "R_PROFILE_OK",
        "HA_PROFILE_OK",
    ]

    print("Pipeline completion")
    print("-------------------")

    for f in flags:
        if f in tab.colnames:
            col = np.array(tab[f], dtype=bool)

            ntrue = np.sum(col)
            nfalse = np.sum(~col)

            pct = 100 * ntrue / n

            print(f"{f:18s}: {ntrue:4d} OK  | {nfalse:4d} FAIL  ({pct:5.1f}%)")
        


    print()
    if "R_PROFILE_OK" in tab.colnames and "HA_PROFILE_OK" in tab.colnames:
        both = np.logical_and(tab["R_PROFILE_OK"], tab["HA_PROFILE_OK"])
        print(f"{'PROFILES_BOTH':18s}: {np.sum(both):4d}")

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

    if len(sys.argv) < 2:
        print("Usage: summarize_run.py merged_results.fits")
        sys.exit()

    summarize(sys.argv[1])


if __name__ == "__main__":
    main()
