from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from hapy.utils.plotting import raincloud_by_group
from hapy.utils.results_table import safe_bool_array, safe_float_array, safe_str_array
from hapy.utils.results_table import first_existing_col, first_populated_col#, build_row_qc_flags, 
from hapy.utils.results_table import ensure_dir, median_and_mad
from hapy.utils.results_table import prepare_analysis_table


def main():


    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Survey-level QC for merged HAPY results.")
    parser.add_argument("table", help="Merged HAPY results table, e.g. merged_results.fits")
    parser.add_argument("--outdir", default="qc", help="Output directory")
    parser.add_argument(
        "--scheme",
        choices=["virgo", "agc"],
        required=True,
        help="Pipeline stage whose results should be merged."
    )
    parser.add_argument(
        "--max-ha-filter-correction",
        type=float,
        default=1.2,
        help="Threshold above which FILTER_CORRECTION is flagged as a warning for Halpha science",
    )
    args = parser.parse_args()




    tab = Table.read(args.table)
    print(f"Read {len(tab)} rows from {args.table}")

    # -- add columns for qc and science
    tab = prepare_analysis_table(tab)
