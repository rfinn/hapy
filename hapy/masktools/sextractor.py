import os
import subprocess
from pathlib import Path


SEXTRACTOR_FILES = [
    "default.sex.HDI.mask",
    "default.param",
    "default.conv",
    "default.nnw",
]


def link_sextractor_files(sepath, workdir="."):
    """Create symbolic links to required SExtractor config files."""
    for fname in SEXTRACTOR_FILES:
        src = Path(sepath) / fname
        dst = Path(workdir) / fname

        if dst.exists() or dst.is_symlink():
            dst.unlink()

        dst.symlink_to(src)


def clean_sextractor_links(workdir="."):
    """Remove SExtractor symbolic links."""
    for fname in SEXTRACTOR_FILES:
        path = Path(workdir) / fname
        if path.exists() or path.is_symlink():
            path.unlink()


import numpy as np
from astropy.io import fits


def run_sextractor(
    image_name,
    config,
    threshold,
    snr,
    snr_analysis,
    minarea,
    weight_image=None,
    weight_threshold=1,
):
    """
    Run SExtractor and return segmentation array + catalog filename.
    """

    catname = image_name.replace(".fits", ".cat")
    segmentation = image_name.replace(".fits", "-segmentation.fits")

    cmd = [
        "sex",
        image_name,
        "-c", config,
        "-CATALOG_NAME", catname,
        "-CATALOG_TYPE", "FITS_1.0",
        "-DEBLEND_MINCONT", str(threshold),
        "-DETECT_THRESH", str(snr),
        "-ANALYSIS_THRESH", str(snr_analysis),
        "-CHECKIMAGE_NAME", segmentation,
        "-DETECT_MINAREA", str(minarea),
    ]

    if weight_image is not None:
        cmd += [
            "-WEIGHT_TYPE", "MAP_WEIGHT",
            "-WEIGHT_IMAGE", weight_image,
            "-WEIGHT_THRESH", str(weight_threshold),
        ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    segdata = fits.getdata(segmentation)

    return segdata, catname, segmentation


def read_se_catalog(catname, xc=None, yc=None):
    """
    Read SExtractor catalog and optionally identify central object.
    
    Returns:
        catalog_table, central_object_number (or None)
    """

    sexout = fits.getdata(catname)

    if len(sexout) == 0:
        return sexout, None

    xsex = sexout["XWIN_IMAGE"]
    ysex = sexout["YWIN_IMAGE"]

    objnumb = None

    if xc is not None and yc is not None:
        dist = np.sqrt((yc - ysex) ** 2 + (xc - xsex) ** 2)
        obj_index = np.argmin(dist)
        objnumb = sexout["NUMBER"][obj_index]

    return sexout, objnumb


