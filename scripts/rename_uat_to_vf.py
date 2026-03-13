#!/usr/bin/env python3

import os
import re
import argparse
from pathlib import Path


def convert_name(fname):
    """
    Convert UAT filename to VF naming scheme.
    """

    name = fname.name

    # Replace prefix
    name = name.replace("UAT-", "VF-", 1)

    # Replace last -NN- pattern with _NN-
    name = re.sub(r'-(\d+)-', r'_\1-', name, count=1)

    return name


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Directory containing FITS files")
    parser.add_argument("--apply", action="store_true",
                        help="Actually rename files (otherwise dry run)")
    args = parser.parse_args()

    d = Path(args.directory)

    files = sorted(d.glob("UAT-*.fits"))

    for f in files:

        newname = convert_name(f)

        newpath = f.with_name(newname)

        print(f"{f.name}  →  {newname}")

        if args.apply:
            f.rename(newpath)


if __name__ == "__main__":
    main()
