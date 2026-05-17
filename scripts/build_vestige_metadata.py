#!/usr/bin/env python

"""
build_vestige_metadata.py

Create HAPY cutout directories + metadata.json files for VESTIGE/CFHT images.

Expected input naming:
    VCC226_r.fits
    VCC226_ha.fits

The *_ha.fits image is assumed to be continuum-subtracted Halpha.
"""

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
import sys

def get_wcs_center(fitsfile):
    with fits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        ny, nx = hdul[0].data.shape

    w = WCS(hdr)
    ra, dec = w.pixel_to_world_values(nx / 2, ny / 2)
    return float(ra), float(dec), hdr


def match_vfs(ra, dec, v, max_sep_arcsec=60):
    imgcoord = SkyCoord(ra, dec, unit="deg")
    catcoord = SkyCoord(v.main["RA"], v.main["DEC"], unit="deg")

    idx, sep2d, _ = imgcoord.match_to_catalog_sky(catcoord)

    if sep2d.arcsec > max_sep_arcsec:
        return None, None, sep2d.arcsec

    vfid = v.main["VFID"][idx]

    # match ephot row by VFID
    eidx = np.where(v.ephot["VFID"] == vfid)[0]
    erow = v.ephot[eidx[0]] if len(eidx) > 0 else None


    idx = int(np.asarray(idx).ravel()[0])
    sep_arcsec = float(np.asarray(sep2d.arcsec).ravel()[0])

    if sep_arcsec > max_sep_arcsec:
        return None, None, sep_arcsec

    vfid = v.main["VFID"][idx]

    eidx = np.where(v.ephot["VFID"] == vfid)[0]
    erow = v.ephot[int(eidx[0])] if len(eidx) > 0 else None

    return v.main[idx], erow, sep_arcsec




def clean_name(x):
    return str(x).strip().replace(" ", "").replace("/", "-")


def get_dateobs(header):
    """
    Use DATE-OBS first, falling back to DATE1.
    Return YYYYMMDD.
    """
    date = header.get("DATE-OBS", header.get("DATE1", "UNKNOWN"))
    return str(date).split("T")[0].replace("-", "")


def make_tag(mainrow, header):
    vfid = clean_name(mainrow["VFID"])
    nedname = clean_name(mainrow["NEDname"])
    telescope = "CFHT"
    dateobs = get_dateobs(header)
    target = clean_name(header.get("OBJECT", "VESTIGE"))

    return f"{vfid}-{nedname}-{telescope}-{dateobs}-{target}"



def make_metadata(tag, rfile, hafile, mainrow, erow, rhdr, hahdr, sep_arcsec):

    pixscale = abs(float(rhdr.get("CD1_1", np.nan))) * 3600.0

    r_finaliq = float(rhdr.get("FINALIQ", np.nan))
    ha_finaliq = float(hahdr.get("FINALIQ", np.nan))

    r_fwhm_pix = (
        r_finaliq / pixscale
        if np.isfinite(r_finaliq) and np.isfinite(pixscale) and pixscale > 0
        else np.nan
    )

    ha_fwhm_pix = (
        ha_finaliq / pixscale
        if np.isfinite(ha_finaliq) and np.isfinite(pixscale) and pixscale > 0
        else np.nan
    )

    metadata = {
        "objid": str(mainrow["VFID"]),
        "VFID": str(mainrow["VFID"]),
        "NEDname": str(mainrow["NEDname"]),
        "tag": tag,

        "ra": float(mainrow["RA"]),
        "dec": float(mainrow["DEC"]),

        "telescope": "CFHT",
        "instrument": "MegaCam",
        "dateobs": get_dateobs(rhdr),
        "pointing": str(rhdr.get("OBJECT", "VESTIGE")),
        "scheme": "vestige",

        "parent_rimage": Path(rfile).name,
        "parent_haimage": Path(hafile).name,

        # ha image is already continuum-subtracted
        "parent_csimage": Path(hafile).name,

        "hafilter": "VESTIGE",
        "filter_correction": 1.0,
        "filter_ratio": 1.0,

        "match_sep_arcsec": float(sep_arcsec),

        "pixscale": pixscale,
        "pixelscale": pixscale,

        "r_exptime": float(rhdr.get("EXPTIME", np.nan)),
        "ha_exptime": float(hahdr.get("EXPTIME", np.nan)),
        "r_photzp": float(rhdr.get("PHOTZP", np.nan)),
        "ha_photzp": float(hahdr.get("PHOTZP", np.nan)),

        "r_finaliq": r_finaliq,
        "ha_finaliq": ha_finaliq,

        "rimage_fwhm_psf_arcsec": r_finaliq,
        "himage_fwhm_psf_arcsec": ha_finaliq,
        "rimage_fwhm_psf_pix": r_fwhm_pix,
        "himage_fwhm_psf_pix": ha_fwhm_pix,

        "fwhm_psf_arcsec": ha_finaliq,
        "fwhm_psf_pix": ha_fwhm_pix,

        "r_maglim": float(rhdr.get("MAGLIM", np.nan)),
        "ha_maglim": float(hahdr.get("MAGLIM", np.nan)),
    }
    return metadata

        



def process_one(rfile, args, v):
    rfile = Path(rfile)
    hafile = rfile.with_name(rfile.name.replace("_r.fits", "_ha.fits"))

    if not hafile.exists():
        print(f"SKIP: no matching Halpha image for {rfile}")
        return

    ra, dec, rhdr = get_wcs_center(rfile)
    _, _, hahdr = get_wcs_center(hafile)

    mainrow, erow, sep = match_vfs(ra, dec, v, args.max_sep)

    if mainrow is None:
        print(f"NO MATCH: {rfile} sep={sep:.1f} arcsec")
        return

    tag = make_tag(mainrow, rhdr)
    cutdir = Path(args.outdir) / tag
    cutdir.mkdir(parents=True, exist_ok=True)

    metadata_file = cutdir / "metadata.json"

    if metadata_file.exists() and not args.overwrite:
        print(f"SKIP existing: {metadata_file}")
        return

    metadata = make_metadata(
        tag=tag,
        rfile=rfile,
        hafile=hafile,
        mainrow=mainrow,
        erow=erow,
        rhdr=rhdr,
        hahdr=hahdr,
        sep_arcsec=sep,
    )

    r_out = cutdir / f"{tag}-R.fits"
    cs_out = cutdir / f"{tag}-CS-ZP.fits"

    shutil.copy2(rfile, r_out)
    shutil.copy2(hafile, cs_out)

    metadata["parent_rimage"] = r_out.name
    metadata["parent_haimage"] = cs_out.name
    metadata["parent_csimage"] = cs_out.name


    #shutil.copy2(rfile, cutdir / rfile.name)
    #shutil.copy2(hafile, cutdir / hafile.name)

    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=4)

    print(f"OK: {rfile.name} -> {cutdir} sep={sep:.1f} arcsec")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rglob", default="*_r.fits")
    parser.add_argument("--outdir", default="cutouts")
    parser.add_argument("--max-sep", type=float, default=60.0)
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    # assumes Virgo catalog object is importable the usual way
    homedir = os.getenv("HOME")
    sys.path.append(os.path.join(homedir,'github/Virgo/programs/'))
    from readtablesv2 import vtables
    
    tabledir= os.path.join(homedir,'research/Virgo/tables-north/v2/')
    tableprefix = 'vf_v2_'
                        
    v = vtables(tabledir, tableprefix)
    v.read_all()
    
    rfiles = sorted(Path(".").glob(args.rglob))

    for rfile in rfiles:
        process_one(rfile, args, v)


if __name__ == "__main__":
    main()
