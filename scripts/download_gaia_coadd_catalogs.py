#!/usr/bin/env python
import argparse
from pathlib import Path
from hapy.masktools.gaia import get_gaia_stars

def main(args):
    # get list of coadds
    patterns = []
    if args.pre2025:
        patterns.append("*BOK*r-shifted.fits")
        patterns.append("*INT*r-shifted.fits")
        patterns.append("*HDI*R.fits")
        patterns.append("*HDI*r.fits")
        # skipping mosaic b/c nothing new, and hardly any duplicates...
    else:
        if args.prefix is not None:
            patterns.append(f"{args.prefix}*R.fits")
            patterns.append(f"{args.prefix}*r.fits")
        else:
            patterns.append("*R.fits")
            patterns.append("*r.fits")
        

    flist = []
    for pat in patterns:
        flist.extend(Path(".").glob(pat))
    flist = sorted(set(flist))
    
    
    # make output directory gaia_catalogs
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # for each coadd,
    import time
    import requests

    failed = []

    for image_name in flist:
        image_path = Path(image_name)
        outfile = outdir / f"{image_path.stem}-gaia.fits"

        if outfile.exists() and not args.overwrite:
            print(f"Skipping existing {outfile}")
            continue

        success = False

        for attempt in range(5):
            try:
                print(f"Querying Gaia for {image_path} (attempt {attempt+1}/5)")
                gaiatab, xpixel, ypixel = get_gaia_stars(str(image_path), use_cache=False)
                gaiatab.write(outfile, format="fits", overwrite=True)
                success = True
                break

            except requests.exceptions.HTTPError as e:
                print(f"HTTPError for {image_path}: {e}")
                if attempt < 4:
                    time.sleep(10 * (attempt + 1))

            except Exception as e:
                print(f"ERROR for {image_path}: {e}")
                break

        if not success:
            failed.append(str(image_path))
            print(f"FAILED: {image_path}")

        if args.testing:
            break

    if failed:
        failfile = outdir / "failed_gaia_queries.txt"
        with open(failfile, "w") as fh:
            for name in failed:
                fh.write(name + "\n")
        print(f"Wrote failed query list to {failfile}")

    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Download gaia catalogs for a set of FITS images."
    )

    g_io = p.add_argument_group("Input / Output")
    g_io.add_argument(
        "--prefix",
        default=None,
        help="Prefix for images to glob. Default: all *R.fits and *r.fits files in the current directory."
    )
    g_io.add_argument(
        "--outdir",
        default="gaia_catalogs",
        help="Output directory for Gaia catalogs. Default: gaia_catalogs"
    )
    g_io.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing Gaia catalogs.  Default is False."
    )

    p.add_argument(
        "--testing",
        action="store_true",
        help="Run on one image only."
    )
    p.add_argument(
        "--pre2025",
        action="store_true",
        help="Run on pre 2025 coadds (different naming convention)."
    )
    


    args = p.parse_args()
    raise SystemExit(main(args))
