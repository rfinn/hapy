#!/usr/bin/env python

"""
GOAL:
- make a list of cutouts directories to input into run_analysis, but I want to skip galaxies that have CATALOG_USE == EXCLUDE in review_sample_20260514.csv.

USAGE:

python ~/github/hapy/scripts/make_run_analysis_list.py --cutout-dir cutouts --review review_sample_20260514.csv --outfile cutout_run_analysis_list.txt

"""
import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cutout-dir", default="cutouts")
    parser.add_argument("--review", default="review_sample_20260514.csv")
    parser.add_argument("--outfile", default="cutout_run_analysis_list.txt")
    parser.add_argument("--exclude-value", default="EXCLUDE")
    args = parser.parse_args()

    review = pd.read_csv(args.review, comment="#")

    if "TAG" not in review.columns:
        raise ValueError("Expected TAG column in review CSV.")
    if "CATALOG_USE" not in review.columns:
        raise ValueError("Expected CATALOG_USE column in review CSV.")

    exclude_tags = set(
        review.loc[
            review["CATALOG_USE"].astype(str).str.upper().str.strip() == args.exclude_value,
            "TAG",
        ].astype(str)
    )

    cutout_dirs = sorted(
        p for p in Path(args.cutout_dir).iterdir()
        if p.is_dir() and p.name != "cutouts_summary"
    )

    keep = [
        str(p)
        for p in cutout_dirs
        if p.name not in exclude_tags
    ]

    with open(args.outfile, "w") as f:
        for item in keep:
            f.write(item + "\n")

    print(f"Found {len(cutout_dirs)} cutout directories")
    print(f"Excluded {len(cutout_dirs) - len(keep)} directories")
    print(f"Wrote {len(keep)} directories to {args.outfile}")


if __name__ == "__main__":
    main()
