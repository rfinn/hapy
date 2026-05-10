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

MANUAL_BEST_TAG = {
    # "DUP_GALID": "TAG to use"
    "VFID1234": "VFID1234-NGC4567-INT-20210414-p001",
    "VFID5678": "VFID5678-NGC9999-BOK-20220312-VFID0000",
}

def finite_median(tab, col):
    if col not in tab.colnames:
        return np.nan
    vals = np.array([safe_float(x) for x in tab[col]])
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) == 0:
        return np.nan
    return np.nanmedian(vals)




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


def get_display_limits_from_row(row, shape, buffer_pix=100, halfsize_pix=None):
    ny, nx = shape

    xc = safe_float(get_col(row, ["ELLIP_XCENTROID", "GAL_XC", "XC", "xcenter"]))
    yc = safe_float(get_col(row, ["ELLIP_YCENTROID", "GAL_YC", "YC", "ycenter"]))

    sma = safe_float(get_col(row, ["ELLIP_SMA_PIX", "SMA_PIX", "sma_pix"]))

    if not np.isfinite(xc):
        xc = nx / 2
    if not np.isfinite(yc):
        yc = ny / 2

    if halfsize_pix is not None:
        halfsize = halfsize_pix
    elif np.isfinite(sma):
        halfsize = sma + buffer_pix
    else:
        return None

    xmin = max(0, int(xc - halfsize))
    xmax = min(nx - 1, int(xc + halfsize))
    ymin = max(0, int(yc - halfsize))
    ymax = min(ny - 1, int(yc + halfsize))

    return (xmin, xmax), (ymin, ymax)

def get_display_limits_from_row_v0(row, shape, buffer_pix=100):
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


def get_group_display_halfsize_arcsec(rows, min_buffer_arcsec=60, scale=1.2):
    """
    Pick one angular half-size for all duplicate panels.

    Uses the largest available SMA-like radius, then adds buffer.
    """
    sizes = []

    for row in rows:
        sma_arcsec = safe_float(
            get_col(
                row,
                [
                    "SMA_ARCSEC",
                    "sma_arcsec",
                    "ELLIP_SMA_ARCSEC",
                    "R25_ARCSEC",
                    "R24_ARCSEC",
                ],
            )
        )

        if np.isfinite(sma_arcsec) and sma_arcsec > 0:
            sizes.append(scale * sma_arcsec + min_buffer_arcsec)

    if len(sizes) == 0:
        return None

    return np.nanmax(sizes)

def infer_galid(row):
    for col in ["VFID", "OBJID", "objid", "GALID", "galid"]:
        if col in row.colnames:
            val = str(row[col])
            if val and val != "nan":
                return val.split("-")[0]
    tag = str(row["TAG"])
    return tag.split("-")[0]


def normalized_value(val, key, norms):
    med = norms.get(key, np.nan)
    if np.isfinite(val) and np.isfinite(med) and med > 0:
        return val / med
    return np.nan

def safe_float(x):
    try:
        if np.ma.is_masked(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan
    
def short_tag(tag):
    parts = str(tag).split("-")
    if len(parts) > 2:
        return "-".join(parts[2:])
    return str(tag)

def telescope_rank(row):
    tel = str(get_col(row, ["TELESCOPE", "telescope", "INSTRUMENT", "instrument"], ""))
    tag = str(get_col(row, ["TAG"], ""))

    text = f"{tel} {tag}".upper()
    for key, rank in TEL_PRIORITY.items():
        if key in text:
            return rank
    return 2


def score_duplicate(row, norms=None):
    """
    Lower score is better.

    Approximate inverse S/N cost:

        cost ~ FILTER_CORRECTION * FWHM^2 * SKYSTD

    Halpha term is weighted more strongly because duplicate choice is mainly
    driven by Halpha morphology/extent quality.
    """
    if norms is None:
        norms = {}

    def ratio_norm(val, key):
        med = norms.get(key, np.nan)
        if np.isfinite(val) and np.isfinite(med) and med > 0 and val > 0:
            return val / med
        return np.nan

    r_fwhm = safe_float(get_col(row, ["R_FWHM_PSF", "R_FWHM_PSF_ARCSEC"]))
    h_fwhm = safe_float(get_col(row, ["H_FWHM_PSF", "H_FWHM_PSF_ARCSEC"]))
    r_sky  = safe_float(get_col(row, ["R_SKYSTD_PHYS"]))
    h_sky  = safe_float(get_col(row, ["H_SKYSTD_PHYS"]))
    fcorr  = safe_float(get_col(row, ["FILTER_CORRECTION"]))

    r_fwhm_n = ratio_norm(r_fwhm, "R_FWHM_PSF")
    h_fwhm_n = ratio_norm(h_fwhm, "H_FWHM_PSF")
    r_sky_n  = ratio_norm(r_sky,  "R_SKYSTD_PHYS")
    h_sky_n  = ratio_norm(h_sky,  "H_SKYSTD_PHYS")

    # MOS sky values are artificially low after convolution.
    # Do not let MOS win because of unrealistically low SKYSTD.
    is_mos = "MOS" in str(row["TAG"]).upper() or "MOSAIC" in str(row["TAG"]).upper()
    if is_mos:
        r_sky_n = 1.0
        h_sky_n = 1.0

    # Missing values should lose.
    if not np.isfinite(r_fwhm_n):
        r_fwhm_n = 99.0
    if not np.isfinite(h_fwhm_n):
        h_fwhm_n = 99.0
    if not np.isfinite(r_sky_n):
        r_sky_n = 99.0
    if not np.isfinite(h_sky_n):
        h_sky_n = 99.0

    if not np.isfinite(fcorr) or fcorr <= 0:
        fcorr = 99.0

    r_cost = (r_fwhm_n ** 2) * r_sky_n
    h_cost = fcorr * (h_fwhm_n ** 2) * h_sky_n

    score = 0.3 * r_cost + 0.7 * h_cost

    # Strong extra penalty for large filter correction.
    # This keeps fcorr in the score continuously, but still flags risky cases.
    if fcorr >= 1.2:
        score += 10.0 * (fcorr - 1.2)

    # Mild telescope tie-breaker.
    score += 0.05 * telescope_rank(row)

    return score

def score_duplicate_v1(row, norms=None):
    """
    Lower score is better.

    Uses log(value / median) so:
      0      = typical
      < 0    = better than median
      > 0    = worse than median

    FWHM is weighted more strongly than sky.
    """
    if norms is None:
        norms = {}

    def log_norm(val, key):
        med = norms.get(key, np.nan)
        if np.isfinite(val) and np.isfinite(med) and med > 0 and val > 0:
            return np.log(val / med)
        return np.nan

    r_fwhm = safe_float(get_col(row, ["R_FWHM_PSF", "R_FWHM_PSF_ARCSEC"]))
    h_fwhm = safe_float(get_col(row, ["H_FWHM_PSF", "H_FWHM_PSF_ARCSEC"]))
    r_sky  = safe_float(get_col(row, ["R_SKYSTD_PHYS"]))
    h_sky  = safe_float(get_col(row, ["H_SKYSTD_PHYS"]))
    fcorr  = safe_float(get_col(row, ["FILTER_CORRECTION"]))

    r_fwhm_n = log_norm(r_fwhm, "R_FWHM_PSF")
    h_fwhm_n = log_norm(h_fwhm, "H_FWHM_PSF")
    r_sky_n  = log_norm(r_sky,  "R_SKYSTD_PHYS")
    h_sky_n  = log_norm(h_sky,  "H_SKYSTD_PHYS")

    score = 0.0

    # FWHM dominates; sky matters but less

    for val, weight in [
        (2.0 * r_fwhm_n, 10.0),
        (2.0 * h_fwhm_n, 10.0),
        (r_sky_n,       1.0),
        (h_sky_n,       1.0),
    ]:    

        if np.isfinite(val):
            score += weight * val
        else:
            score += 999.0

    # Strong penalty if filter correction is too large
    if np.isfinite(fcorr):
        if fcorr >= 1.2:
            score += 100.0 + 50.0 * (fcorr - 1.2)
    else:
        score += 50.0

    # Mild tie-breaker: BOK/INT preferred over HDI/MOS
    score += 0.3 * telescope_rank(row)

    return score

def score_duplicate_v0(row, norms=None):
    if norms is None:
        norms = {}

    r_fwhm = safe_float(get_col(row, ["R_FWHM_PSF", "R_FWHM_PSF_ARCSEC"]))
    h_fwhm = safe_float(get_col(row, ["H_FWHM_PSF", "H_FWHM_PSF_ARCSEC"]))
    r_sky  = safe_float(get_col(row, ["R_SKYSTD_PHYS"]))
    h_sky  = safe_float(get_col(row, ["H_SKYSTD_PHYS"]))
    fcorr  = safe_float(get_col(row, ["FILTER_CORRECTION"]))

    # --- normalize ---
    def norm(val, key):
        med = norms.get(key, np.nan)
        if np.isfinite(val) and np.isfinite(med) and med > 0:
            return val / med
        return np.nan

    r_fwhm_n = norm(r_fwhm, "R_FWHM_PSF")
    h_fwhm_n = norm(h_fwhm, "H_FWHM_PSF")
    r_sky_n  = norm(r_sky,  "R_SKYSTD_PHYS")
    h_sky_n  = norm(h_sky,  "H_SKYSTD_PHYS")

    score = 0.0

    # --- weights (now dimensionless) ---
    for val, weight in [
        (r_fwhm_n, 5.0),
        (h_fwhm_n, 5.0),
        (r_sky_n,  1.5),
        (h_sky_n,  1.5),
    ]:
        if np.isfinite(val):
            score += weight * val
        else:
            score += 999.0

    # --- filter correction penalty ---
    if np.isfinite(fcorr):
        if fcorr >= 1.2:
            score += 100.0 + 50.0 * (fcorr - 1.2)
    else:
        score += 50.0

    # --- telescope tie-break ---
    score += 0.3 * telescope_rank(row)

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
        fontsize=11,
        color="white",
        bbox=dict(facecolor="black", alpha=0.55, edgecolor="none"),
    )


def mark_best_panel(ax):
    for spine in ax.spines.values():
        spine.set_edgecolor("lime")
        spine.set_linewidth(4)


def plot_duplicate_group(rows, best_idx, cutout_dir, outdir, galid, norms=None):
    """
    Make duplicate comparison plot.

    Top row: R-band image, asinh stretch
    Bottom row: Halpha CS-ZP image, zscale stretch

    Best duplicate is outlined in green.
    """
    if norms is None:
        norms = {}

    n = len(rows)

    fig, axes = plt.subplots(
        2,
        n,
        figsize=(4.2 * n, 8),
        squeeze=False,
        constrained_layout=True,
    )

    full_galname = full_galname_from_tag(str(rows[0]["TAG"]))

    group_halfsize_arcsec = get_group_display_halfsize_arcsec(
        rows,
        min_buffer_arcsec=60,
        scale=1.2,
        )
    for j, row in enumerate(rows):
        tag = str(row["TAG"])

        r_path = find_image(
            cutout_dir,
            tag,
            ["-R.fits", "_R.fits", "-r.fits"],
        )

        cs_path = find_image(
            cutout_dir,
            tag,
            ["-CS-ZP.fits", "_CS-ZP.fits", "-CS.fits", "_CS.fits"],
        )

        r_img = read_image(r_path)
        cs_img = read_image(cs_path)

        # ------------------------------------------------------------
        # Display limits: prefer ellipse/row-based zoom over full cutout
        # ------------------------------------------------------------

        
        #limits = None
        #if r_img is not None:
        #    limits = get_display_limits_from_row(row, r_img.shape, buffer_pix=125)


        limits = None

        if r_img is not None and group_halfsize_arcsec is not None:
            pixscale = safe_float(get_col(row, ["PIXSCALE", "pixscale"]))

            if np.isfinite(pixscale) and pixscale > 0:
                buffer_pix = group_halfsize_arcsec / pixscale
                limits = get_display_limits_from_row(
                    row,
                    r_img.shape,
                    buffer_pix=buffer_pix,
                )

        if limits is None and r_img is not None:
            limits = get_display_limits_from_row(row, r_img.shape, buffer_pix=125)

        # ============================================================
        # Top row: R band
        # ============================================================
        ax = axes[0, j]

        if r_img is not None:
            ax.imshow(
                r_img,
                origin="lower",
                cmap="gray",
                norm=image_norm(r_img, "asinh"),
            )
        else:
            ax.text(
                0.5,
                0.5,
                "missing R image",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )

        ax.set_title(short_tag(tag), fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])

        if limits is not None:
            xlim, ylim = limits
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)

        r_fwhm = safe_float(get_col(row, ["R_FWHM_PSF", "R_FWHM_PSF_ARCSEC"]))
        r_sky = safe_float(get_col(row, ["R_SKYSTD_PHYS"]))

        r_fwhm_n = normalized_value(r_fwhm, "R_FWHM_PSF", norms)
        r_sky_n = normalized_value(r_sky, "R_SKYSTD_PHYS", norms)

        add_panel_text(
            ax,
            f"R FWHM={r_fwhm:.2f} ({r_fwhm_n:.2f}x)\n"
            f"R sky={r_sky:.2g} ({r_sky_n:.2f}x)",
        )

        # ============================================================
        # Bottom row: continuum-subtracted Halpha
        # ============================================================
        ax = axes[1, j]

        if cs_img is not None:
            ax.imshow(
                cs_img,
                origin="lower",
                cmap="gray",
                norm=image_norm(cs_img, "zscale"),
            )
        else:
            ax.text(
                0.5,
                0.5,
                "missing CS-ZP image",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )

        ax.set_xticks([])
        ax.set_yticks([])

        if limits is not None:
            xlim, ylim = limits
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)

        h_fwhm = safe_float(get_col(row, ["H_FWHM_PSF", "H_FWHM_PSF_ARCSEC"]))
        h_sky = safe_float(get_col(row, ["H_SKYSTD_PHYS"]))
        fcorr = safe_float(get_col(row, ["FILTER_CORRECTION"]))

        h_fwhm_n = normalized_value(h_fwhm, "H_FWHM_PSF", norms)
        h_sky_n = normalized_value(h_sky, "H_SKYSTD_PHYS", norms)

        add_panel_text(
            ax,
            f"H FWHM={h_fwhm:.2f} ({h_fwhm_n:.2f}x)\n"
            f"H sky={h_sky:.3g} ({h_sky_n:.2f}x)\n"
            f"filter corr={fcorr:.2f}",
        )

        # ------------------------------------------------------------
        # Mark best duplicate
        # ------------------------------------------------------------
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

    fig.suptitle(f"Duplicate comparison: {full_galname}", fontsize=16)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    outfile = outdir / f"{full_galname}_duplicate_comparison.png"
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

    norms = {
        "R_FWHM_PSF": finite_median(tab, "R_FWHM_PSF"),
        "H_FWHM_PSF": finite_median(tab, "H_FWHM_PSF"),
        "R_SKYSTD_PHYS": finite_median(tab, "R_SKYSTD_PHYS"),
        "H_SKYSTD_PHYS": finite_median(tab, "H_SKYSTD_PHYS"),
        }

        
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
        scores = np.array([score_duplicate(row, norms=norms) for row in rows])

        best_local = int(np.nanargmin(scores))
        best_global = idx[best_local]

        manual_tag = MANUAL_BEST_TAG.get(galid, None)

        if manual_tag is not None:
            matches = np.where(np.array([str(tab[i]["TAG"]) for i in idx]) == manual_tag)[0]

            if len(matches) == 1:
                best_local = int(matches[0])
                best_global = idx[best_local]
                manual_override = True
                override_note = "hard-coded manual override"
            else:
                manual_override = False
                override_note = f"manual override tag not found: {manual_tag}"
        else:
            manual_override = False
            override_note = ""

        
        for k, global_i in enumerate(idx):
            tab["BEST_DUPLICATE"] = False if "BEST_DUPLICATE" not in tab.colnames else tab["BEST_DUPLICATE"]

    
            
        best_rows.append(
            {
                "DUP_GALID": galid,
                "N_DUP": len(idx),
                "BEST_TAG": str(tab[best_global]["TAG"]),
                "BEST_SCORE": scores[best_local],
                "ALL_TAGS": ",".join(str(tab[i]["TAG"]) for i in idx),
                "ALL_SCORES": ",".join(f"{s:.4f}" for s in scores),
                "USE_TAG": str(tab[best_global]["TAG"]),
                "MANUAL_OVERRIDE": manual_override,
                "NOTES": override_note,
            }
        )


        png = plot_duplicate_group(
            rows=rows,
            best_idx=best_local,
            cutout_dir=args.cutout_dir,
            outdir=args.outdir,
            galid=galid,
            norms=norms,
        )


        pngs.append(png)
        print(f"{galid}: best = {tab[best_global]['TAG']} -> {png}")

    best_tab = Table(rows=best_rows)

    # Add manual override columns
    best_tab["USE_TAG"] = best_tab["BEST_TAG"]
    best_tab["MANUAL_OVERRIDE"] = False
    best_tab["NOTES"] = ""


    

    out_table = Path(args.outdir) / "best_duplicates.ecsv"
    Path(args.outdir).mkdir(exist_ok=True, parents=True)


    best_tab = Table(rows=best_rows)

    # Add override columns
    best_tab["USE_TAG"] = best_tab["BEST_TAG"]
    best_tab["MANUAL_OVERRIDE"] = False
    best_tab["NOTES"] = ""

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    # --- ECSV (pipeline-safe) ---
    ecsv_file = outdir / "best_duplicates.ecsv"
    best_tab.write(ecsv_file, format="ascii.ecsv", overwrite=True)

    # --- CSV (human-editable) ---
    csv_file = outdir / "best_duplicates.csv"
    best_tab.write(csv_file, format="ascii.csv", overwrite=True)

    print(f"\nWrote {len(best_tab)} rows:")
    print(f"  {ecsv_file}")
    print(f"  {csv_file}")




    print(f"\nWrote {len(best_tab)} duplicate selections to {out_table}")
    print(f"Wrote {len(pngs)} PNGs to {args.outdir}/")


if __name__ == "__main__":
    main()
