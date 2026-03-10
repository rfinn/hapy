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
import numpy as np
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


def validate_schema(tables, filenames):
    """Ensure all tables share identical column names."""
    keepflag = np.ones(len(tables), dtype=bool)
    reference = tables[0].colnames

    for i, t in enumerate(tables[1:], start=1):
        if t.colnames != reference:
            print(f"WAIT!!! Problem with table {filenames[i]}!!!")
            print("Schema mismatch detected.\n")
            print(f"Expected columns:\n{reference}\n")
            print(f"Found columns:\n{t.colnames}\n")
            keepflag[i] = False

    return keepflag

def _coerce_bool_col(tab, name, default=False):
    if name not in tab.colnames:
        tab[name] = np.full(len(tab), default, dtype=bool)
        return

    col = tab[name]

    # masked -> fill default
    try:
        if hasattr(col, "filled"):
            col = col.filled(default)
    except Exception:
        pass

    # coerce object/mixed to bool safely
    if getattr(col, "dtype", None) == object:
        def asbool(v):
            if v is None:
                return default
            if isinstance(v, (bool, np.bool_)):
                return bool(v)
            s = str(v).strip().lower()
            if s in ("true", "t", "1", "yes", "y"):
                return True
            if s in ("false", "f", "0", "no", "n", "", "none", "nan"):
                return False
            return default

        tab[name] = np.array([asbool(v) for v in col], dtype=bool)
    else:
        tab[name] = np.array(col, dtype=bool)
        
def merge_tables(files, output, mode):
    """Read, validate, merge, and write output FITS table."""
    print(f"Found {len(files)} result files.")
    print("Reading tables...")

    tables = [Table.read(f, format="ascii.ecsv") for f in files]

    if mode == "run_analysis":
        for t in tables:
            _coerce_bool_col(t, "R_SM_FLAG", default=False)
            _coerce_bool_col(t, "H_SM_FLAG", default=False)
        
    print("Validating schema...")
    keepflag = validate_schema(tables,files)

    print(f"\tvalidated {np.sum(keepflag)}/{len(keepflag)} tables")
    #tables = tables[keepflag]

    goodtables = []
    for i in range(len(tables)):
        if keepflag[i]:
            goodtables.append(tables[i])

    tables = goodtables

    if not tables:
        raise RuntimeError("No valid tables remain after schema validation.")
    
    print("Stacking tables...")
    merged = vstack(tables, metadata_conflicts="silent")

    if mode == "run_analysis":
        if "OBJID" in merged.colnames:
            merged["obs_id"] = [Path(str(r)).name for r in merged["OBJID"]]

    if "SIGMA_FITS" in merged.colnames:
        merged.remove_column("SIGMA_FITS")
    
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
        "--outdir",
        type=str,
        default=None,
        help="Directory where merged table will be written (default: current directory)."
    )
    parser.add_argument(
    "--pattern",
    default=None,
    help="Optional filename pattern to override the mode-specific default."
    )
    

    parser.add_argument(
        "--out",
        default="merged_results.fits",
        help="Output FITS filename (default: merged_results.fits)"
    )
    parser.add_argument(
        "--mode",
        choices=["run_analysis", "get_cutouts"],
        required=True,
        help="Pipeline stage whose results should be merged."
    )

    args = parser.parse_args()

    if args.pattern is not None:
        pattern = args.pattern
    elif args.mode == "get_cutouts":
        pattern = "cutouts_summary*.ecsv"
    else:
        pattern = "*results.ecsv"
    
    files = find_result_files(args.indir, pattern)

    if args.outdir:
        outpath = Path(args.outdir).resolve() / args.out
    else:
        outpath = Path(args.out).resolve()
    
    merge_tables(files, outpath, args.mode)


if __name__ == "__main__":
    main()
