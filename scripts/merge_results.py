#!/usr/bin/env python
"""
merge_results.py

Merge per-galaxy ECSV result files into a single survey-level FITS table.

Recursively searches subdirectories under --indir for *-results.ecsv files.

Designed for Virgo / AGC survey structure:

    cutouts/
        VFIDxxxx-.../
            VFIDxxxx-...-results.ecsv

Output:
    merged_results.fits (binary FITS table)

Author: Rose Finn
"""

import argparse
from pathlib import Path
from astropy.table import Table, vstack
import sys

def find_result_files(indir, pattern="*-results.ecsv"):
    """Recursively locate result files under indir."""
    files = sorted(Path(indir).rglob(pattern))
    if not files:
        raise RuntimeError(f"No files matching '{pattern}' found in {indir}")
    return files


def validate_schema(tables):
    """Ensure all tables share identical column names."""
    reference = tables[0].colnames
    for i, t in enumerate(tables[1:], start=2):
        if t.colnames != reference:
            raise RuntimeError(
                f"Schema mismatch detected in table #{i}.\n"
                f"Expected columns:\n{reference}\n\n"
                f"Found columns:\n{t.colnames}"
            )


def check_duplicate_objids(tables):
    """Detect duplicate objid entries across tables."""
    objids = []
    for t in tables:
        if "objid" not in t.colnames:
            raise RuntimeError("Column 'objid' not found in results table.")
        objids.append(t["objid"][0])

    duplicates = set([x for x in objids if objids.count(x) > 1])
    if duplicates:
        raise RuntimeError(f"Duplicate objid detected: {duplicates}")


def merge_tables(files, output):
    """Read, validate, merge, and write output FITS table."""
    print(f"Found {len(files)} result files.")
    print("Reading tables...")

    tables = [Table.read(f, format="ascii.ecsv") for f in files]

    print("Validating schema...")
    validate_schema(tables)



    print("Stacking tables...")
    merged = vstack(tables, metadata_conflicts="silent")




    # Add explicit observation ID column
    obs_ids = [Path(r).name for r in merged["root"]]
    merged["obs_id"] = obs_ids
    
    
    print(f"Writing merged table → {output}")
    merged.write(output, format="fits", overwrite=True)

    print("Done.")
    print(f"Final table rows: {len(merged)}")
    print(f"Final table columns: {len(merged.colnames)}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-galaxy ECSV results into a single FITS table."
    )
    parser.add_argument(
        "--indir",
        required=True,
        help="Root directory containing galaxy subdirectories."
    )
    parser.add_argument(
        "--pattern",
        default="*-results.ecsv",
        help="Filename pattern to search for (default: *-results.ecsv)"
    )
    parser.add_argument(
        "--out",
        default="merged_results.fits",
        help="Output FITS filename (default: merged_results.fits)"
    )

    args = parser.parse_args()

    files = find_result_files(args.indir, args.pattern)


    
    merge_tables(files, args.out)


if __name__ == "__main__":
    main()
