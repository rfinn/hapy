#!/usr/bin/env python

import argparse
from pathlib import Path
import json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from hapy.imagetools.plotting import display_image
from hapy.geometry.adapters import pa_ccw_north_to_photutils_theta
#from hapy.geometry.adapters import pa_ccw_north_to_photutils

"""
USAGE:
python ~/github/hapy/scripts/plot_mask_diagnostic.py --root cutouts/VFID3084-NGC3512-HDI-20200226-p012/VFID3084-NGC3512-HDI-20200226-p012

"""

def ellipse_patch(xc, yc, sma_pix, ba, pa_deg, **kwargs):
    # matplotlib Ellipse needs angle relative to +x axis
    return Ellipse(
        (xc, yc),
        width=2 * sma_pix,
        height=2 * sma_pix * ba,
        angle=pa_deg,
        fill=False,
        **kwargs,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--params", default=None)
    args = p.parse_args()

    root = args.root
    r_fits = root + "-R.fits"
    mask_fits = root + "-mask.fits"

    r_data, r_hdr = fits.getdata(r_fits, header=True)
    m_data = fits.getdata(mask_fits)

    params_path = args.params or Path(root).parent / "metadata.json"
    params = json.loads(Path(params_path).read_text())

    # world → pixel
    w = WCS(r_hdr)
    sc = SkyCoord(params["ra"] * u.deg, params["dec"] * u.deg)
    xc, yc = w.world_to_pixel(sc)

    pixscale = float(r_hdr.get("PIXSCALE", 0.426))  # fallback if needed
    sma_pix = params["sma_arcsec"] / pixscale

    theta_photutils = pa_ccw_north_to_photutils_theta(float(params["pa_deg"]))
    print(f"checking PA angles: pa_deg from metadata = {float(params['pa_deg']):.1f}, theta_photutils = {theta_photutils:.1f}")
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    plt.sca(ax[0])
    display_image(r_data)
    #ax[0].imshow(r_data, origin="lower")
    ax[0].add_patch(ellipse_patch(xc, yc, sma_pix, params["ba"], theta_photutils, edgecolor="cyan", linewidth=2))
    ax[0].set_title("R Image")

    ax[1].imshow(m_data, origin="lower")
    ax[1].add_patch(ellipse_patch(xc, yc, sma_pix, params["ba"], theta_photutils, edgecolor="cyan", linewidth=2))
    ax[1].set_title("Mask")

    plt.tight_layout()
    plt.savefig(root + "-mask_diag.png", dpi=150)
    print("Wrote:", root + "-mask_diag.png")


if __name__ == "__main__":
    main()
