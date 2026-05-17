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
from astropy.table import Table, vstack, Column
import sys

from collections import defaultdict

from hapy.utils.results_table import get_excluded_tags, add_review_columns



STRING_COLS = [
    "TAG",
    "OBJID",
    "VFID",
    "NEDNAME",
    "R_FITS",
    "HA_FITS",
    "CS_FITS",
    "MASK_FITS",
    "PHOTFILE",
    "PHOTFILE2",
]


def normalize_string_columns(tables, string_cols=STRING_COLS):
    for tab in tables:
        for col in string_cols:
            if col not in tab.colnames:
                continue

            values = []
            for val in tab[col]:
                if val is None:
                    values.append("")
                elif isinstance(val, bytes):
                    values.append(val.decode("utf-8", errors="ignore"))
                else:
                    values.append(str(val))

            tab[col] = Column(values, name=col, dtype="U512")

    return tables


def coerce_bool_columns(tab, columns=None):
    """
    Force selected columns to boolean to avoid merge dtype conflicts.
    """
    import numpy as np

    if columns is None:
        columns = [
            "HAPY_MORPH_OK",
        ]

    for col in columns:
        if col not in tab.colnames:
            continue

        vals = np.asarray(tab[col])

        # Handles bool, int 0/1, strings like True/False, masked values
        out = np.zeros(len(tab), dtype=bool)

        for i, v in enumerate(vals):
            if np.ma.is_masked(v):
                out[i] = False
            elif isinstance(v, str):
                out[i] = v.strip().lower() in ("true", "t", "1", "yes")
            else:
                out[i] = bool(v)

        tab[col] = out

    return tab

def keep_latest_cutout_summaries(files):
    latest = {}

    for f in files:
        f = Path(f)
        name = f.stem

        if not name.startswith("cutouts_summary-"):
            print(f"WARNING: unexpected filename, skipping: {f}")
            continue

        remainder = name[len("cutouts_summary-"):]
        parts = remainder.rsplit("-", 2)

        if len(parts) != 3:
            print(f"WARNING: could not parse summary filename, skipping: {f}")
            continue

        tag, user, ts = parts

        if tag in latest:
            old_ts, old_file = latest[tag]
            if ts > old_ts:
                latest[tag] = (ts, f)
            elif ts == old_ts:
                print(f"WARNING: duplicate timestamp for tag {tag}:")
                print(f"    {old_file}")
                print(f"    {f}")
        else:
            latest[tag] = (ts, f)

    return sorted(v[1] for v in latest.values())



def find_result_files(indir, pattern="*-results.ecsv"):
    """Recursively locate result files under indir."""
    print("Searching for files ",pattern)
    files = sorted(Path(indir).rglob(pattern))
    if not files:
        raise RuntimeError(f"No files matching '{pattern}' found in {indir}")
    return files




def validate_schema(tables, filenames, reorder=True):
    """
    Ensure all tables share the same column names, ignoring column order.

    Parameters
    ----------
    tables : list
        List of astropy Tables.
    filenames : list
        Corresponding filenames.
    reorder : bool
        If True, reorder columns in-place to match the first table.

    Returns
    -------
    keepflag : np.ndarray
        Boolean array marking tables with compatible schemas.
    """

    keepflag = np.ones(len(tables), dtype=bool)
    reference = list(tables[0].colnames)
    reference_set = set(reference)

    for i, t in enumerate(tables[1:], start=1):
        this_cols = list(t.colnames)
        this_set = set(this_cols)

        missing = sorted(reference_set - this_set)
        extra = sorted(this_set - reference_set)

        if missing or extra:
            print(f"WAIT!!! Problem with table {filenames[i]}!!!")
            print("Schema mismatch detected.\n")
            if missing:
                print(f"Missing columns:\n{missing}\n")
            if extra:
                print(f"Extra columns:\n{extra}\n")
            keepflag[i] = False
            continue

        if this_cols != reference:
            print(f"Column order differs for {filenames[i]}; reordering to match reference table.")
            if reorder:
                tables[i] = t[reference]

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


def coerce_object_columns_for_fits(tab):
    """
    Convert object dtype columns into FITS-safe numeric or string columns.
    """

    numeric_like_cols = {
        "REDSHIFT",
        "vr",
        "velocity",
        "distance",
    }

    for col in tab.colnames:
        if tab[col].dtype.kind != "O":
            continue

        if col in numeric_like_cols:
            vals = []
            for x in tab[col]:
                if x is None or str(x).strip() in ["", "None", "nan", "--"]:
                    vals.append(np.nan)
                else:
                    vals.append(float(x))
            tab[col] = np.array(vals, dtype=float)
        else:
            tab[col] = np.array(
                ["" if x is None else str(x) for x in tab[col]],
                dtype=f"U{max(1, max(len(str(x)) for x in tab[col]))}",
            )

    return tab

def merge_tables(files, output, mode, review_csv=None):
    """Read, validate, merge, and write output FITS table."""
    print(f"Merging {len(files)} result files.")
    print("Reading tables...")

    tables = []
    for f in files:
        t = Table.read(f, format="ascii.ecsv")
        ok_cols = [c for c in t.colnames if c.endswith("_OK")]
        t = coerce_bool_columns(t, columns=ok_cols)
        tables.append(t)
    
    # tables = [Table.read(f, format="ascii.ecsv") for f in files]

    # # adding protection for HAPY_MORPH_OK and other _OK columns
    # for t in tables:
    #     ok_cols = [c for c in t.colnames if c.endswith("_OK")]
    #     t = coerce_bool_columns(t, columns=ok_cols)
    
    #if mode == "run_analysis":
        #for t in tables:
        #    _coerce_bool_col(t, "R_SM_FLAG", default=False)
        #    _coerce_bool_col(t, "H_SM_FLAG", default=False)
        
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

    print("Normalizing string columns...")
    tables = normalize_string_columns(tables)

    print("Stacking tables...")
    merged = vstack(tables, metadata_conflicts="silent")

    #if mode == "run_analysis":
    #    if "OBJID" in merged.colnames:
    #        merged["obs_id"] = [Path(str(r)).name for r in merged["OBJID"]]

    if "SIGMA_FITS" in merged.colnames:
        merged.remove_column("SIGMA_FITS")

    if review_csv is not None:
        merged = add_review_columns(merged, review_csv)



    # FITS cannot write object dtype columns.
    # Convert known numeric metadata columns that may contain None/string/NaN.


    merged = coerce_object_columns_for_fits(merged)

    print(f"Writing merged table → {output}")
    merged.write(output, format="fits", overwrite=True)

    print("Done.")
    print(f"Final table rows: {len(merged)}")
    print(f"Final table columns: {len(merged.colnames)}")



def infer_scheme_from_result_files(files, max_files=None):
    """
    Infer scheme label from input result tables.

    Returns
    -------
    scheme : str or None
        Single scheme if all readable files agree, "mixed" if multiple
        schemes are found, or None if no scheme could be determined.
    """
    schemes = set()
    use_files = files if max_files is None else files[:max_files]

    for f in use_files:
        try:
            t = Table.read(f, format="ascii.ecsv")
        except Exception:
            continue

        if "SCHEME" not in t.colnames or len(t) == 0:
            continue

        val = t["SCHEME"][0]
        if hasattr(val, "item"):
            try:
                val = val.item()
            except Exception:
                pass

        if isinstance(val, str):
            val = val.strip()

        if isinstance(val, str) and val not in ("", "None", "nan"):
            schemes.add(val)

    if len(schemes) == 1:
        return next(iter(schemes))
    elif len(schemes) > 1:
        return "mixed"
    else:
        return None

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
        default=None,
        help="Output FITS filename"
    )

    parser.add_argument(
        "--mode",
        choices=["run_analysis", "get_cutouts"],
        required=True,
        help="Pipeline stage whose results should be merged."
    )
    parser.add_argument(
        "--latest-only",
        action="store_true",
        help="For get_cutouts mode, keep only the newest summary file per tag."
    )

    parser.add_argument(
    "--review-csv",
    default=None,
    help="Optional review CSV; rows with CATALOG_USE == EXCLUDE are skipped"
    )
    
    args = parser.parse_args()

    if args.pattern is not None:
        pattern = args.pattern
    elif args.mode == "get_cutouts":
        pattern = "cutouts_summary*.ecsv"
    else:
        pattern = "*results.ecsv"
    
    files = find_result_files(args.indir, pattern)
    print(f"Found {len(files)} result files.")
    
    if args.review_csv is not None:
        excluded_tags = get_excluded_tags(args.review_csv)

        kept = []
        skipped = []

        for f in files:
            tag = Path(f).parent.name
            if tag in excluded_tags:
                skipped.append(f)
            else:
                kept.append(f)

        print(f"Skipping {len(skipped)} files with CATALOG_USE == EXCLUDE")
        files = kept
    

    

    if args.mode == "get_cutouts" and args.latest_only:
        files = keep_latest_cutout_summaries(files)
        print(f"Keeping latest summary per tag: {len(files)} files")

    from datetime import datetime
    today = datetime.now().strftime("%Y%m%d")

    if args.out is None:
        scheme = infer_scheme_from_result_files(files, max_files=1)
        if scheme is not None:
            args.out = f"merged_results_{scheme}_{today}.fits"
        else:
            args.out = f"merged_results_{today}.fits"

    if args.outdir:
        outpath = Path(args.outdir).resolve() / args.out
    else:
        outpath = Path(args.out).resolve()

    merge_tables(files, outpath, args.mode, review_csv=args.review_csv)


 


if __name__ == "__main__":
    main()
