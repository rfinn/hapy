#!/usr/bin/env python
"""
USAGE:

python scripts/merge_results.py --indir cutouts --out all-results.ecsv

"""
from __future__ import annotations

import argparse
from pathlib import Path
import glob

from hapy.hatools.results import merge_result_rows_ecsv


def main():
    p = argparse.ArgumentParser(description="Merge per-galaxy ECSV result files into one table")
    p.add_argument("--indir", default="cutouts", help="Base cutouts directory (default: cutouts)")
    p.add_argument("--pattern", default="*-results.ecsv", help="Glob pattern (default: *-results.ecsv)")
    p.add_argument("--out", default="merged-results.ecsv", help="Output ECSV filename")
    args = p.parse_args()

    paths = sorted(glob.glob(str(Path(args.indir) / "*" / args.pattern)))
    outpath = Path(args.out)

    merge_result_rows_ecsv(paths, outpath, overwrite=True)
    print(f"Merged {len(paths)} files -> {outpath}")


if __name__ == "__main__":
    main()
