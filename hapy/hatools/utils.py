#import importlib.resources as pkg_resources
import os
import numpy as np

from astropy.wcs import WCS

# utils.py
from pathlib import Path
#import pkg_resources  # legacy, still works for pip-installed packages
#import importlib.resources as resources  # modern alternative
from hapy import hatools


from importlib import resources


def zp_scale_r_to_ha(zp_ha, zp_r, logger=None):
    """Scale factor alpha so that CS = Ha - alpha * R."""
    if zp_ha is None or zp_r is None:
        return np.nan

    try:
        zp_ha = float(zp_ha)
        zp_r = float(zp_r)
    except Exception:
        if logger:
            logger.warning("Could not calculate ZP scale: non-numeric ZPs")
        return np.nan

    if not (np.isfinite(zp_ha) and np.isfinite(zp_r)):
        if logger:
            logger.warning("Could not calculate ZP scale: non-finite ZPs")
        return np.nan

    return float(10 ** (-0.4 * (zp_r - zp_ha)))

def get_filter_file(name: str) -> str:
    """
    Return the full path to a filter file in hatools/filter_traces.

    Works both for development installs (editable) and normal installs.
    """
    try:
        p = resources.files("hapy.hatools.filter_traces").joinpath(name)
        if p.is_file():
            return str(p)
    except Exception:
        pass

    # editable/dev fallback
    base_dir = Path(hatools.__file__).parent
    fpath = base_dir / "filter_traces" / name
    if fpath.exists():
        return str(fpath)

    raise FileNotFoundError(f"Filter file {name} not found in hatools/filter_traces/")



def parse_agc_name(image_name: str) -> dict:
    """
    Parse AGC/UAT Groups Survey-style coadd name.

    Expected patterns (by splitting basename on '-'):
      - Positive dec: len(parts)==7  -> UAT-RA+DEC-telescope-dateobs-pointingA-pointingB-filter
      - Negative dec: len(parts)==8  -> UAT-RA-DEC-telescope-dateobs-pointingA-pointingB-filter
    """

    name = Path(image_name).name
    parts = name.split("-")

    if len(parts) == 7:
        ra, dec = parts[1].split("+")
        telescope = parts[2]
        dateobs = parts[3]
        pointing = parts[4] + "-" + parts[5]

    elif len(parts) == 8:
        ra = parts[1]
        dec = f"-{parts[2]}"
        telescope = parts[3]
        dateobs = parts[4]
        pointing = parts[5] + "-" + parts[6]

    else:
        raise ValueError(f"Cannot parse agc coadd name: {image_name}")

    return {
        "telescope": telescope,
        "dateobs": dateobs,
        "pointing": pointing,
        "ra": float(ra),
        "dec": float(dec),
    }





def parse_virgo_name(image_name: str) -> dict:
    """
    Parse Virgo-style coadd name.

    Returns:
        dict with keys:
            telescope, dateobs, pointing, ra, dec
    """

    name = Path(image_name).name
    parts = name.split("-")

    if len(parts) == 6:
        # UAT-RA+DEC-telescope-date-pointing-filter
        ra, dec = parts[1].split("+")
        telescope = parts[2]
        dateobs = parts[3]
        pointing = parts[4]

    elif len(parts) == 7:
        # negative declination
        ra = parts[1]
        dec = f"-{parts[2]}"
        telescope = parts[3]
        dateobs = parts[4]
        pointing = parts[5]

    else:
        raise ValueError(f"Cannot parse Virgo coadd name: {image_name}")

    return {
        "telescope": telescope,
        "dateobs": dateobs,
        "pointing": pointing,
        "ra": float(ra),
        "dec": float(dec),
    }

def parse_coadd_name(path: str, scheme: str = "generic") -> dict:

    if scheme == "virgo":
        return parse_virgo_name(path)

    elif scheme == "agc":
        return parse_agc_name(path)

    else:
        return {"basename": Path(path).stem}

def build_cutout_name(tokens, galname, outdir):
    """
    Returns a root path prefix for cutouts under:
      <outdir>/<galname>/<galname>-<telescope>-<dateobs>-<pointing>

    This mirrors the old layout: one folder per galaxy.
    """
    outdir = Path(outdir)

    if "telescope" not in tokens:
        root = galdir / str(galname)
        return str(root)

    t = tokens["telescope"]
    d = tokens["dateobs"]
    p = tokens["pointing"]


    tag = f"{galname}-{t}-{d}-{p}"

    galdir = outdir / tag
    galdir.mkdir(parents=True, exist_ok=True)

    root = galdir / tag
    return str(root)


def get_survey_vectors(gcat, scheme: str):
    """
    Return (redshift_full, objid_full) arrays aligned with gcat.cat.

    objid is assumed to already exist (ensure_objid invariant).
    """

    objid = gcat.cat["objid"]

    if scheme == "virgo":
        redshift = gcat.cat["vr"] / 3.0e5
        return redshift, objid

    if scheme == "agc":
        redshift = gcat.cat["vopt"] / 3.0e5
        flag0 = gcat.cat["vopt"] == 0
        redshift[flag0] = gcat.cat["v21"][flag0] / 3.0e5
        return redshift, objid

    # generic
    return None, objid
