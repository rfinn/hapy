"""
Purpose:

Build minimal metadata.json files for archival Virgo Hα galaxy cutouts so they can run through run_analysis.py.


Required for Pipeline operation:

metadata.json
{
  "objid": "n4354",
  "ra": 185.734,
  "dec": 13.121,
  "sma_arcsec": 55.2,
  "ba": 0.62,
  "pa_deg": 127.3
}

Easy additions:
{
  "scheme": "archive",
  "has_bad_wcs": true,

  "r_fits": "r.fits",
  "cs_fits": "ha_cs.fits",
  "mask_fits": "mask.fits",

  "rimage_fwhm_arcsec": 1.7,
  "himage_fwhm_arcsec": 1.7
}

"""

from pathlib import Path
import json
from astropy.io import fits
from astropy.table import Table

def build_metadata_archive(root_dir, virgo_catalog):

    catalog = Table.read(virgo_catalog)

    for gal_dir in Path(root_dir).iterdir():

        if not gal_dir.is_dir():
            continue

        objid = gal_dir.name

        r_image = find_r_image(gal_dir)
        ha_image = find_ha_image(gal_dir)

        mask = find_mask(gal_dir)

        row = lookup_catalog(objid, catalog)

        meta = dict(
            objid=objid,
            ra=float(row["RA"]),
            dec=float(row["DEC"]),
            sma_arcsec=float(row["SMA"]),
            ba=float(row["BA"]),
            pa_deg=float(row["PA"]),

            scheme="archive",
            has_bad_wcs=True,

            r_fits=r_image.name,
            cs_fits=ha_image.name,
            mask_fits=mask.name if mask else None,

            rimage_fwhm_arcsec=None,
            himage_fwhm_arcsec=None
        )

        with open(gal_dir / "metadata.json", "w") as f:
            json.dump(meta, f, indent=2)
