#!/usr/bin/env python

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS


def display_limits(data):
    finite = np.isfinite(data)
    if not np.any(finite):
        return 0.0, 1.0
    _, med, std = sigma_clipped_stats(data[finite], sigma=3.0, maxiters=5)
    vmin = med - 1.0 * std
    vmax = med + 5.0 * std
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin, vmax = np.nanmin(data[finite]), np.nanmax(data[finite])
    return vmin, vmax


def get_gaia_xy(tab, wcs):
    """
    Try common Gaia column names and convert RA/Dec to pixel coordinates.
    """
    ra_col = None
    dec_col = None

    for cand in ["ra", "RA", "ra_epoch2000", "RA_ICRS"]:
        if cand in tab.colnames:
            ra_col = cand
            break

    for cand in ["dec", "DEC", "dec_epoch2000", "DE_ICRS", "DEC_ICRS"]:
        if cand in tab.colnames:
            dec_col = cand
            break

    if ra_col is None or dec_col is None:
        raise ValueError(
            f"Could not find RA/Dec columns in Gaia table. Found: {tab.colnames}"
        )

    x, y = wcs.world_to_pixel_values(tab[ra_col], tab[dec_col])
    return np.asarray(x), np.asarray(y)


def plot_gaia_overlay(image_fits, gaia_fits, outfile):
    data, hdr = fits.getdata(image_fits, header=True)
    wcs = WCS(hdr)

    tab = Table.read(gaia_fits)
    xg, yg = get_gaia_xy(tab, wcs)

    ny, nx = data.shape
    inside = (xg >= 0) & (xg < nx) & (yg >= 0) & (yg < ny)

    vmin, vmax = display_limits(data)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax, interpolation="nearest")

    # all Gaia sources
    ax.plot(xg, yg, "co", ms=7, mfc="none", mew=1.2, alpha=0.8, label="Gaia catalog")

    # highlight those inside image footprint
    ax.plot(
        xg[inside], yg[inside],
        "r.", ms=4, alpha=0.9, label="Inside image"
    )

    ax.set_title(Path(image_fits).name)
    ax.set_xlabel("x [pix]")
    ax.set_ylabel("y [pix]")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Wrote: {outfile}")


import sys
from pathlib import Path

def main():
    gaia_dir = Path("gaia_catalogs")
    outdir = Path("gaia_diagnostic")
    outdir.mkdir(exist_ok=True)

    if len(sys.argv) > 1:
        image_fits = Path(sys.argv[1])

        if image_fits.suffix != ".fits":
            image_fits = Path(f"{image_fits}.fits")

        root = image_fits.name.replace(".fits", "")
        gaia_fits = gaia_dir / f"{root}-gaia.fits"

        if not image_fits.exists():
            print(f"Missing image: {image_fits}")
            return
        if not gaia_fits.exists():
            print(f"Missing Gaia catalog: {gaia_fits}")
            return

        outfile = outdir / f"{root}-gaia-check.png"
        plot_gaia_overlay(image_fits, gaia_fits, outfile)
        return

    gaia_files = sorted(gaia_dir.glob("*-gaia.fits"))
    if not gaia_files:
        print("No Gaia catalogs found in gaia_catalogs/")
        return

    for gaia_fits in gaia_files:
        root = gaia_fits.name.replace("-gaia.fits", "")
        image_fits = Path(f"{root}.fits")

        if not image_fits.exists():
            print(f"Missing image for {gaia_fits.name}: expected {image_fits}")
            continue

        outfile = outdir / f"{root}-gaia-check.png"
        try:
            plot_gaia_overlay(image_fits, gaia_fits, outfile)
        except Exception as e:
            print(f"ERROR processing {root}: {e}")




if __name__ == "__main__":
    main()
