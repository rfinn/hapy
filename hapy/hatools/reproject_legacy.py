"""
Utilities to reproject Legacy Survey images onto the Halpha cutout footprint.
"""

from pathlib import Path

from astropy.io import fits
from reproject import reproject_adaptive


def reproject_image_to_reference(infile, reffile, outname, overwrite=False):
    """
    Reproject infile to the WCS/pixel grid of reffile.
    """

    infile = Path(infile)
    reffile = Path(reffile)
    outname = Path(outname)

    if outname.exists() and not overwrite:
        print(f"reprojected image exists - not redoing it: {outname}")
        return outname

    with fits.open(infile) as hinfile, fits.open(reffile) as href:
        outim, footprint = reproject_adaptive(
            hinfile,
            href[0].header,
            conserve_flux=True,
        )

        fits.writeto(
            outname,
            outim,
            href[0].header,
            overwrite=True,
        )

    return outname


def find_halpha_reference(cutdir):
    """
    Find the Halpha-frame reference image for a HAPY cutout.
    Prefer CS-ZP if present, otherwise Ha.
    """

    cutdir = Path(cutdir)
    tag = cutdir.name

    candidates = [
        cutdir / f"{tag}-CS-ZP.fits",
        cutdir / f"{tag}-CS.fits",
        cutdir / f"{tag}-Ha.fits",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        f"Could not find Halpha reference image in {cutdir}. "
        f"Tried: {[str(c) for c in candidates]}"
    )


def find_legacy_images(cutdir):
    """
    Find original Legacy g/r/z images, excluding already-reprojected products.
    """

    cutdir = Path(cutdir)
    legacy_dir = cutdir / "legacy"

    if not legacy_dir.exists():
        raise FileNotFoundError(f"Missing legacy directory: {legacy_dir}")

    images = []

    for band in ["g", "r", "z"]:
        images.extend(legacy_dir.glob(f"*-{band}.fits"))

    # Avoid reprocessing outputs such as *r-ha.fits
    images = [
        image for image in images
        if not image.name.endswith("-ha.fits")
    ]

    return sorted(images)


def make_legacy_reprojections(cutdir, overwrite=False):
    """
    Reproject Legacy g/r/z images onto the Halpha image grid.
    """

    cutdir = Path(cutdir)
    reffile = find_halpha_reference(cutdir)
    legacy_images = find_legacy_images(cutdir)

    if len(legacy_images) == 0:
        print(f"No Legacy images found in {cutdir / 'legacy'}")
        return []

    outputs = []

    for infile in legacy_images:
        outname = infile.with_name(infile.stem + "-ha.fits")
        outfile = reproject_image_to_reference(
            infile,
            reffile,
            outname,
            overwrite=overwrite,
        )
        outputs.append(outfile)

    return outputs
