#!/usr/bin/env python

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy.io import fits
from astropy.visualization import ImageNormalize, AsinhStretch, ZScaleInterval


TEL_PRIORITY = {
    "BOK": 0,
    "INT": 0,
    "HDI": 1,
    "MOS": 1,
    "MOSAIC": 1,
}


def short_tag(tag):
    parts = str(tag).split("-")
    if len(parts) > 2:
        return "-".join(parts[2:])
    return str(tag)


def full_galname_from_tag(tag):
    parts = str(tag).split("-")
    if len(parts) >= 2:
        return "-".join(parts[:2])
    return str(tag)

def get_col(row, names, default=np.nan):
    for name in names:
        if name in row.colnames:
            return row[name]
    return default

def get_display_limits_from_row(row, shape, buffer_pix=100):
    """
    Fallback display crop using ellipse size if segmentation map is unavailable.
    """
    ny, nx = shape

    xc = safe_float(get_col(row, ["ELLIP_XCENTROID", "GAL_XC", "XC", "xcenter"]))
    yc = safe_float(get_col(row, ["ELLIP_YCENTROID", "GAL_YC", "YC", "ycenter"]))

    sma = safe_float(get_col(row, ["ELLIP_SMA_PIX", "SMA_PIX", "sma_pix"]))

    if not np.isfinite(xc):
        xc = nx / 2
    if not np.isfinite(yc):
        yc = ny / 2
    if not np.isfinite(sma):
        return None

    halfsize = sma + buffer_pix

    xmin = max(0, int(xc - halfsize))
    xmax = min(nx - 1, int(xc + halfsize))
    ymin = max(0, int(yc - halfsize))
    ymax = min(ny - 1, int(yc + halfsize))

    return (xmin, xmax), (ymin, ymax)

def infer_galid(row):
    for col in ["VFID", "OBJID", "objid", "GALID", "galid"]:
        if col in row.colnames:
            val = str(row[col])
            if val and val != "nan":
                return val.split("-")[0]
    tag = str(row["TAG"])
    return tag.split("-")[0]


def safe_float(x):
    try:
        if np.ma.is_masked(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def telescope_rank(row):
    tel = str(get_col(row, ["TELESCOPE", "telescope", "INSTRUMENT", "instrument"], ""))
    tag = str(get_col(row, ["TAG"], ""))

    text = f"{tel} {tag}".upper()
    for key, rank in TEL_PRIORITY.items():
        if key in text:
            return rank
    return 2


def score_duplicate(row):
    r_fwhm = safe_float(get_col(row, ["R_FWHM_PSF", "R_FWHM_PSF_ARCSEC"]))
    h_fwhm = safe_float(get_col(row, ["H_FWHM_PSF", "H_FWHM_PSF_ARCSEC"]))
    r_sky = safe_float(get_col(row, ["R_SKYSTD_PHYS"]))
    h_sky = safe_float(get_col(row, ["H_SKYSTD_PHYS"]))
    fcorr = safe_float(get_col(row, ["FILTER_CORRECTION"]))

    score = 0.0

    # FWHM should dominate
    for val, weight in [
        (r_fwhm, 10.0),
        (h_fwhm, 10.0),
        (r_sky, 1.0),
        (h_sky, 1.0),
    ]:
        if np.isfinite(val):
            score += weight * val
        else:
            score += 999.0

    if np.isfinite(fcorr):
        if fcorr >= 1.2:
            score += 100.0 + 50.0 * (fcorr - 1.2)
    else:
        score += 50.0

    score += 0.5 * telescope_rank(row)

    return score


def find_image(cutout_dir, tag, suffixes):
    cdir = Path(cutout_dir) / tag
    for suffix in suffixes:
        path = cdir / f"{tag}{suffix}"
        if path.exists():
            return path

    # fallback glob
    for suffix in suffixes:
        matches = sorted(cdir.glob(f"*{suffix}"))
        if matches:
            return matches[0]

    return None


def read_image(path):
    if path is None:
        return None
    with fits.open(path) as hdul:
        return hdul[0].data.astype(float)


def image_norm(data, mode):
    if data is None:
        return None

    if mode == "asinh":
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(data[np.isfinite(data)])
        return ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch())

    if mode == "zscale":
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(data[np.isfinite(data)])
        return ImageNormalize(vmin=vmin, vmax=vmax)

    return None


def add_panel_text(ax, text):
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="white",
        bbox=dict(facecolor="black", alpha=0.55, edgecolor="none"),
    )


def mark_best_panel(ax):
    for spine in ax.spines.values():
        spine.set_edgecolor("lime")
        spine.set_linewidth(4)


def plot_duplicate_group(rows, best_idx, cutout_dir, outdir, galid):
    n = len(rows)
    fig, axes = plt.subplots(
        2,
        n,
        figsize=(4.2 * n, 8),
        squeeze=False,
        constrained_layout=True,
    )

    for j, row in enumerate(rows):
        tag = str(row["TAG"])

        r_path = find_image(cutout_dir, tag, ["-R.fits", "_R.fits", "-r.fits"])
        cs_path = find_image(
            cutout_dir,
            tag,
            ["-CS-ZP.fits", "_CS-ZP.fits", "-CS.fits", "_CS.fits"],
        )

        r_img = read_image(r_path)
        cs_img = read_image(cs_path)

        # Top row: r-band, asinh
        ax = axes[0, j]
        if r_img is not None:
            ax.imshow(r_img, origin="lower", cmap="gray", norm=image_norm(r_img, "asinh"))
        else:
            ax.text(0.5, 0.5, "missing R image", ha="center", va="center")
        #ax.set_title(tag, fontsize=10)
        ax.set_title(short_tag(tag), fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])

        r_fwhm = safe_float(get_col(row, ["R_FWHM_PSF", "R_FWHM_PSF_ARCSEC"]))
        r_sky = safe_float(get_col(row, ["R_SKYSTD_PHYS"]))
        add_panel_text(ax, f"R FWHM={r_fwhm:.2f}\nR sky={r_sky:.3g}")

        limits = get_display_limits_from_row(row, r_img.shape, buffer_pix=125)

        if limits is not None:
            xlim, ylim = limits
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
        
        # Bottom row: CS-ZP, zscale
        ax = axes[1, j]
        if cs_img is not None:
            ax.imshow(cs_img, origin="lower", cmap="gray", norm=image_norm(cs_img, "zscale"))
        else:
            ax.text(0.5, 0.5, "missing CS-ZP image", ha="center", va="center")
        ax.set_xticks([])
        ax.set_yticks([])

        h_fwhm = safe_float(get_col(row, ["H_FWHM_PSF", "H_FWHM_PSF_ARCSEC"]))
        h_sky = safe_float(get_col(row, ["H_SKYSTD_PHYS"]))
        fcorr = safe_float(get_col(row, ["FILTER_CORRECTION"]))
        add_panel_text(
            ax,
            f"H FWHM={h_fwhm:.2f}\nH sky={h_sky:.3g}\nfilter corr={fcorr:.2f}",
        )

        limits = get_display_limits_from_row(row, r_img.shape, buffer_pix=125)

        if limits is not None:
            xlim, ylim = limits
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            
        if j == best_idx:
            mark_best_panel(axes[0, j])
            mark_best_panel(axes[1, j])
            axes[0, j].text(
                0.5,
                1.08,
                "BEST",
                transform=axes[0, j].transAxes,
                ha="center",
                va="bottom",
                fontsize=13,
                color="lime",
                fontweight="bold",
            )

    #fig.suptitle(f"Duplicate comparison: {galid}", fontsize=16)
    full_galname = full_galname_from_tag(str(rows[0]["TAG"]))
    fig.suptitle(f"Duplicate comparison: {full_galname}", fontsize=16)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)
    outfile = outdir / f"{galid}_duplicate_comparison.png"
    fig.savefig(outfile, dpi=150)
    plt.close(fig)

    return outfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("merged_results", help="merged_results_*.fits file")
    parser.add_argument(
        "--cutout-dir",
        default="cutouts",
        help="Directory containing cutouts/<TAG>/",
    )
    parser.add_argument(
        "--outdir",
        default="duplicate_comparison",
        help="Output directory for PNGs and table",
    )
    parser.add_argument(
        "--min-dups",
        type=int,
        default=2,
        help="Minimum number of observations required",
    )
    args = parser.parse_args()

    tab = Table.read(args.merged_results)

    if "TAG" not in tab.colnames:
        raise ValueError("Expected a TAG column in merged_results table.")

    galids = np.array([infer_galid(row) for row in tab])
    tab["DUP_GALID"] = galids

    best_rows = []
    pngs = []

    for galid in sorted(set(galids)):
        idx = np.where(galids == galid)[0]
        if len(idx) < args.min_dups:
            continue

        rows = tab[idx]
        scores = np.array([score_duplicate(row) for row in rows])
        best_local = int(np.nanargmin(scores))
        best_global = idx[best_local]

        for k, global_i in enumerate(idx):
            tab["BEST_DUPLICATE"] = False if "BEST_DUPLICATE" not in tab.colnames else tab["BEST_DUPLICATE"]

        best_rows.append(
            {
                "DUP_GALID": galid,
                "N_DUP": len(idx),
                "BEST_TAG": str(tab[best_global]["TAG"]),
                "BEST_SCORE": scores[best_local],
                "ALL_TAGS": ",".join(str(tab[i]["TAG"]) for i in idx),
            }
        )

        png = plot_duplicate_group(
            rows=rows,
            best_idx=best_local,
            cutout_dir=args.cutout_dir,
            outdir=args.outdir,
            galid=galid,
        )
        pngs.append(png)
        print(f"{galid}: best = {tab[best_global]['TAG']} -> {png}")

    best_tab = Table(rows=best_rows)
    out_table = Path(args.outdir) / "best_duplicates.fits"
    Path(args.outdir).mkdir(exist_ok=True, parents=True)
    best_tab.write(out_table, overwrite=True)

    print(f"\nWrote {len(best_tab)} duplicate selections to {out_table}")
    print(f"Wrote {len(pngs)} PNGs to {args.outdir}/")


if __name__ == "__main__":
    main()
