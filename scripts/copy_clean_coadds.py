#!/usr/bin/env python
"""
python ~/github/hapy/scripts/copy_clean_coadds.py /data-pool/Halpha/coadds-pre2025/all-virgo-coadds /data-pool/Halpha/coadds-pre2025-hapy


python ~/github/hapy/scripts/copy_clean_coadds.py /data-pool/Halpha/coadds-pre2025/all-virgo-coadds /data-pool/Halpha/coadds-pre2025-hapy --clobber
"""
from pathlib import Path
import shutil
import argparse
from astropy.io import fits


HEADER_IMAGE_KEYS = [
    "HAIMAGE",
    "RIMAGE",
    "WEIGHT",
    "WTIMAGE",
    "WEIGHTIMG",
    "RWEIGHT",
    "HAWEIGHT",
    "BKGIMAGE",
]


def clean_name(name):
    return (
        str(name)
        .replace("-shifted.weight.fits", ".weight.fits")
        .replace("-shifted.fits", ".fits")
        .replace("-shifted", "")
    )


def copy_clean(src_dir, out_dir, clobber=False):
    src_dir = Path(src_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    copied = []

    for infile in sorted(src_dir.glob("*.fits")):
        outfile = out_dir / clean_name(infile.name)

        if outfile.exists() and not clobber:
            print(f"SKIP exists: {outfile.name}")
        else:
            print(f"COPY {infile.name} -> {outfile.name}")
            shutil.copy2(infile, outfile)

        copied.append(outfile)

    return copied


def rewrite_shifted_header_values(fitsfile):
    changed = False

    with fits.open(fitsfile, mode="update") as hdul:
        hdr = hdul[0].header

        for key in HEADER_IMAGE_KEYS:
            if key in hdr and isinstance(hdr[key], str) and "-shifted" in hdr[key]:
                old = hdr[key]
                hdr[key] = clean_name(old)
                print(f"  {fitsfile.name}: {key}: {old} -> {hdr[key]}")
                changed = True

        if changed:
            hdul.flush()

    return changed


def get_header_value(fitsfile, key):
    with fits.open(fitsfile) as hdul:
        return hdul[0].header.get(key)


def set_header_value(fitsfile, key, value, comment=None):
    with fits.open(fitsfile, mode="update") as hdul:
        if comment is None:
            hdul[0].header[key] = value
        else:
            hdul[0].header[key] = (value, comment)
        hdul.flush()


def is_rband_image(fitsfile):
    with fits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        return "HAIMAGE" in hdr


def backfill_rimage_haimage(out_dir):
    """
    For each r-band image with HAIMAGE:
      - clean HAIMAGE name
      - add HAIMAGE to r image if needed/cleaned
      - add RIMAGE to matching Halpha image
      - add HAIMAGE to Halpha image too, for symmetry/useful bookkeeping
    """
    out_dir = Path(out_dir)

    for rfile in sorted(out_dir.glob("*.fits")):
        if not is_rband_image(rfile):
            continue

        ha_name = clean_name(get_header_value(rfile, "HAIMAGE"))
        ha_file = out_dir / Path(ha_name).name

        if not ha_file.exists():
            print(f"WARNING: {rfile.name} has HAIMAGE={ha_name}, but file not found")
            continue

        # ensure r-band points to clean Halpha filename
        set_header_value(
            rfile,
            "HAIMAGE",
            ha_file.name,
            "Matching Halpha image",
        )

        # ensure Halpha points back to clean r-band filename
        set_header_value(
            ha_file,
            "RIMAGE",
            rfile.name,
            "Matching r-band image",
        )

        # optionally also keep HAIMAGE in Halpha header
        set_header_value(
            ha_file,
            "HAIMAGE",
            ha_file.name,
            "This Halpha image filename",
        )

        print(f"LINK {rfile.name} <-> {ha_file.name}")


def main():
    parser = argparse.ArgumentParser(
        description="Copy pre-2025 coadds to clean HAPY-ready names and repair image header links."
    )
    parser.add_argument("src_dir", help="Input coadd directory")
    parser.add_argument("out_dir", help="Output clean coadd directory")
    parser.add_argument(
        "--clobber",
        action="store_true",
        help="Overwrite existing output files",
    )

    args = parser.parse_args()

    copied = copy_clean(args.src_dir, args.out_dir, clobber=args.clobber)

    print("\nRewriting header values containing '-shifted'...")
    for f in copied:
        if f.exists():
            rewrite_shifted_header_values(f)

    print("\nBackfilling RIMAGE/HAIMAGE links...")
    backfill_rimage_haimage(args.out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
