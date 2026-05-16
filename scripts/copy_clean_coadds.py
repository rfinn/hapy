#!/usr/bin/env python

from pathlib import Path
import shutil
import argparse


def clean_name(path):
    return path.name.replace("-shifted.fits", ".fits").replace("-shifted.weight.fits", ".weight.fits")


def copy_clean(src_dir, out_dir):
    src_dir = Path(src_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for infile in sorted(src_dir.glob("*.fits")):

        outfile = out_dir / clean_name(infile)

        if outfile.exists():
            print(f"SKIP exists: {outfile}")
            continue

        print(f"{infile.name} -> {outfile.name}")
        shutil.copy2(infile, outfile)


def main():
    parser = argparse.ArgumentParser(
        description="Copy coadds to a clean directory, removing '-shifted' from FITS filenames."
    )
    parser.add_argument("src_dir", help="Input coadd directory")
    parser.add_argument("out_dir", help="Output clean coadd directory")
    args = parser.parse_args()

    copy_clean(args.src_dir, args.out_dir)


if __name__ == "__main__":
    main()
