#import importlib.resources as pkg_resources
import os
import numpy as np

from astropy.wcs import WCS

# utils.py
from pathlib import Path
#import pkg_resources  # legacy, still works for pip-installed packages
import importlib.resources as resources  # modern alternative
from hapy import hatools


def get_filter_file(name):
    """
    Return the full path to a filter file in hatools/filter_traces.

    Works both for development installs (editable) and normal installs.
    """
    # Option A: Try modern importlib.resources (preferred)
    try:
        # hatools.filter_traces must have __init__.py for this to work
        with resources.path("hatools.filter_traces", name) as p:
            return str(p)
    except (ModuleNotFoundError, FileNotFoundError, AttributeError):
        # fallback: use path relative to hatools package (editable install)
        base_dir = Path(hatools.__file__).parent
        fpath = base_dir / "filter_traces/" / name
        print('fpath = ',fpath)
        if fpath.exists():
            return str(fpath)
        else:
            raise FileNotFoundError(f"Filter file {name} not found in hatools/filter_traces/")


def get_pixel_scale_chatgpt(wcs):
    """
    Compute pixel scale from WCS header.
    Assumes square pixels.
    """
    # cdelt in degrees, convert to arcsec
    if wcs.wcs.has_cd():
        cd = wcs.wcs.cd
        scale = np.sqrt(np.abs(cd[0,0]*cd[1,1] - cd[0,1]*cd[1,0])) * 3600.0
    else:
        scale = np.abs(self.wcs.wcs.cdelt[0]) * 3600.0
    return scale


def get_pixel_scale(header):
    hwcs = WCS(header)
    pixelscale = None
    try:
        pixelscale = get_pixel_scale_chatgpt(hwcs)
    except:
        try:
            pixelscale = np.abs(float(header['CD1_1']))*3600. # convert deg/pix to arcsec/pixel
        except KeyError:
            pixelscale = np.abs(float(header['PC1_1']))*3600. # Siena pipeline from astronometry.net
    return pixelscale


