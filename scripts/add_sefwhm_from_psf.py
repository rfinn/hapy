#!/usr/bin/env python

import argparse
from pathlib import Path
from astropy.io import fits


def find_psf(image_path, psf_dir):
    """
    Find matching PSF image for a coadd image.

    Tries common conventions:
      image.fits -> image-psf.fits
      image.fits -> image_psf.fits
      image.fits -> image*psf*.fits
    """
    image_path = Path(image_path)
    psf_dir = Path(psf_dir)

    stem = image_path.stem
    stem = stem.replace("-shifted","")

    candidates = [
        psf_dir / f"{stem}-psf.fits",
        psf_dir / f"{stem}_psf.fits",
        psf_dir / f"{stem}.psf.fits",
    ]

    for c in candidates:
        if c.exists():
            return c

    glob_matches = sorted(psf_dir.glob(f"{stem}*psf*.fits"))

    if len(glob_matches) == 1:
        return glob_matches[0]

    if len(glob_matches) > 1:
        raise RuntimeError(
            f"Multiple PSF matches for {image_path.name}:\n"
            + "\n".join(str(g) for g in glob_matches)
        )

    raise FileNotFoundError(f"No PSF found for {image_path.name} in {psf_dir}")


def get_sefwhm_from_psf(psf_path):
    with fits.open(psf_path) as hdul:
        hdr = hdul[0].header
        if "SEFWHM" not in hdr:
            raise KeyError(f"SEFWHM not found in PSF header: {psf_path}")
        return hdr["SEFWHM"]


def add_sefwhm_to_image(image_path, psf_dir, overwrite=True):

    image_path = Path(image_path)
    psf_path = find_psf(image_path, psf_dir)
    sefwhm = get_sefwhm_from_psf(psf_path)

    print(f"{image_path.name}: PSF={psf_path.name}, SEFWHM={sefwhm}")

    with fits.open(image_path, mode="update") as hdul:
        hdul[0].header["SEFWHM"] = (
            sefwhm,
            "SExtractor FWHM copied from matching PSF image",
        )
        hdul[0].header["PSFSEFWH"] = (
            psf_path.name,
            "PSF image used for SEFWHM",
        )
        hdul.flush()

    return sefwhm, psf_path


def main():
    parser = argparse.ArgumentParser(
        description="Copy SEFWHM from matching PSF headers into r-band and Halpha coadds."
    )
    parser.add_argument("rimage", help="Input r-band FITS image")
    parser.add_argument(
        "--psf-dir",
        default="/data-pool/Halpha/psf-images-pre2025/",
        help="Directory containing PSF images",
    )
    parser.add_argument(
        "--ha-key",
        default="HAIMAGE",
        help="Header keyword in r-band image giving matching Halpha image",
    )

    args = parser.parse_args()

    
    rimage = Path(rimage)

    # update r-band image
    add_sefwhm_to_image(rimage, args.psf_dir)

    # find matching Halpha image from r-band header
    with fits.open(rimage) as hdul:
        rhdr = hdul[0].header
        if args.ha_key not in rhdr:
            raise KeyError(f"{args.ha_key} not found in r-band header: {rimage}")

        ha_name = rhdr[args.ha_key]

    haimage = Path(ha_name)
    if not haimage.is_absolute():
        haimage = rimage.parent / haimage

    if not haimage.exists():
        raise FileNotFoundError(f"Halpha image from {args.ha_key} does not exist: {haimage}")

    # update Halpha image
    add_sefwhm_to_image(haimage, args.psf_dir)


if __name__ == "__main__":
    main()
