#!/usr/bin/env python

import os
import math
import argparse
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlretrieve

import pandas as pd
import numpy as np
from PIL import Image

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



from hapy.imagetools.downloads import get_legacy_jpeg



def clean_name(x):
    """
    Make a safe string for plotting labels.
    """
    if pd.isna(x):
        return ""
    return str(x)


def make_panel_label(row, radius=False):
    """
    Label each postage stamp with VFID and NEDname.
    """
    vfid = clean_name(row["VFID"])
    nedname = clean_name(row["NEDname"])

    if len(nedname) > 32:
        nedname = nedname[:29] + "..."
    s = f"{vfid}\n{nedname}"
    return s


def make_legacy_pdf(
    csvfile,
    outfile="legacy_jpegs.pdf",
    jpg_dir="legacy-jpegs",
    cutout_arcmin=10.0,
    pixscale=0.262,
    layer="ls-dr9",
    n_per_page=20,
    ncols=5,
    overwrite=False,
    radius=False,
    radius_scale=2.,
    verbose=False,
):
    """
    Read a CSV/table of galaxies, download Legacy JPEG cutouts,
    and create a multi-page PDF.

    Required input columns:
        VFID, RA_1, DEC_1, NEDname
    """
    csvfile = Path(csvfile)
    jpg_dir = Path(jpg_dir)
    jpg_dir.mkdir(exist_ok=True, parents=True)

    # sep=None allows pandas to infer comma vs tab vs whitespace-delimited files.
    tab = pd.read_csv(csvfile, sep=None, engine="python")

    required_cols = ["VFID", "RA_1", "DEC_1", "NEDname"]
    missing = [c for c in required_cols if c not in tab.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # 10 arcmin = 600 arcsec.
    # Legacy image size is in pixels.
    imsize = int(round(cutout_arcmin * 60.0 / pixscale))

    if verbose:
        print(f"cutout_arcmin = {cutout_arcmin}")
        print(f"pixscale = {pixscale} arcsec/pixel")
        print(f"Legacy cutout size = {imsize} pixels")
        print(f"number of galaxies = {len(tab)}")

    # Download images and store filenames in the table.
    jpeg_names = []

    for i, row in tab.iterrows():
        vfid = str(row["VFID"])
        ra = float(row["RA_1"])
        dec = float(row["DEC_1"])
        

        try:
            if radius:
                rad_arcsec = float(row['radius'])
                imsize = int(float(radius_scale) * float(rad_arcsec) / float(pixscale))
                print("imsize = ",imsize)
            print("DEBUG: downloading image for ",vfid)
            jpeg_name = get_legacy_jpeg(
                ra=ra,
                dec=dec,
                galid=vfid,
                pixscale=pixscale,
                imsize=imsize,
                subfolder=str(jpg_dir),
                verbose=verbose,
                layer=layer,
                overwrite=overwrite,
            )
        except Exception as err:
            print(f"WARNING: could not download {vfid}: {err}")
            jpeg_name = None

        jpeg_names.append(jpeg_name)

    tab["jpeg_name"] = jpeg_names

    # Make PDF contact sheets.
    nrows = int(math.ceil(n_per_page / ncols))
    npages = int(math.ceil(len(tab) / n_per_page))

    with PdfPages(outfile) as pdf:
        for page in range(npages):
            start = page * n_per_page
            end = min((page + 1) * n_per_page, len(tab))
            page_tab = tab.iloc[start:end]

            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(15, 12),
                constrained_layout=True,
            )

            axes = np.asarray(axes).reshape(-1)

            for ax in axes:
                ax.axis("off")

            for ax, (_, row) in zip(axes, page_tab.iterrows()):
                jpeg_name = row["jpeg_name"]
                label = make_panel_label(row)
                if radius:
                    #print(radius_scale,row['radius'])
                    try:
                        #print("test", float(radius_scale)*float(row['radius'])/60.)
                        size_arcmin = float(radius_scale)*float(row['radius'])/60.
                        #print(label)
                        label = label +  f" size={size_arcmin:.2f} arcmin"
                    except Exception as err:
                        print("warning: could not append size to label: ",err)


                ax.axis("off")

                if jpeg_name is None or not os.path.exists(jpeg_name):
                    ax.text(
                        0.5,
                        0.5,
                        f"{label}\n\nmissing image",
                        ha="center",
                        va="center",
                        fontsize=9,
                    )
                    continue

                try:
                    img = Image.open(jpeg_name)

                    # Downsample for the PDF page display.
                    # The original downloaded JPG remains unchanged on disk.
                    img.thumbnail((700, 700))

                    ax.imshow(img, origin="upper")
                    ax.set_title(label, fontsize=8)

                except Exception as err:
                    ax.text(
                        0.5,
                        0.5,
                        f"{label}\n\nproblem reading image\n{err}",
                        ha="center",
                        va="center",
                        fontsize=8,
                    )

            fig.suptitle(
                f"Legacy Survey JPEG cutouts: page {page + 1} of {npages}",
                fontsize=14,
            )

            pdf.savefig(fig)
            plt.close(fig)

    print(f"wrote {outfile}")
    print(f"downloaded jpgs are in {jpg_dir}")

    return tab


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download Legacy JPEG cutouts and make a multi-page PDF."
    )

    parser.add_argument(
        "csvfile",
        help="Input CSV/table with columns VFID, RA_1, DEC_1, NEDname.",
    )

    parser.add_argument(
        "--outfile",
        default="legacy_jpegs.pdf",
        help="Output multi-page PDF name.",
    )

    parser.add_argument(
        "--jpg-dir",
        default="legacy-jpegs",
        help="Directory where downloaded JPEGs will be saved.",
    )

    parser.add_argument(
        "--cutout-arcmin",
        type=float,
        default=5.0,
        help="Cutout size in arcmin. Default is 10 arcmin.",
    )

    parser.add_argument(
        "--pixscale",
        type=float,
        default=0.262,
        help="Legacy cutout pixel scale in arcsec/pixel.",
    )

    parser.add_argument(
        "--layer",
        default="ls-dr9",
        help="Legacy Survey layer. Default is ls-dr9.",
    )

    parser.add_argument(
        "--n-per-page",
        type=int,
        default=20,
        help="Number of galaxies per PDF page.",
    )

    parser.add_argument(
        "--ncols",
        type=int,
        default=5,
        help="Number of columns per PDF page.",
    )

    parser.add_argument(
        "--radius",
        default=False,
        action='store_true',
        help="Set this if catalog has radius column to use when making cutouts.  Should be in arcsec.",
    )

    parser.add_argument(
        "--radius-scale",
        default=4.,
        help="Multiplicative factor for radius.",
    )
    
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download JPEGs even if they already exist.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print download messages.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    make_legacy_pdf(
        csvfile=args.csvfile,
        outfile=args.outfile,
        jpg_dir=args.jpg_dir,
        cutout_arcmin=args.cutout_arcmin,
        pixscale=args.pixscale,
        layer=args.layer,
        n_per_page=args.n_per_page,
        ncols=args.ncols,
        overwrite=args.overwrite,
        radius=args.radius,
        radius_scale = args.radius_scale,
        verbose=args.verbose,
    )
