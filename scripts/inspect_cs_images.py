#!/usr/bin/env python
"""
Inspect HAPY continuum-subtracted images and select best duplicate observations.

Creates per-galaxy comparison PNGs for single and duplicate observations.
For duplicates, computes a quality score and writes best_duplicates.csv/ecsv.


To just write the duplicates table:
python ~/github/hapy/scripts/inspect_cs_images.py make-table merged_results_virgo_20260514.fits --outdir cs_image_inspection --min-dups 1



To create an input list for running in parallel:
python  ~/github/hapy/scripts/inspect_cs_images.py list-groups cs_image_inspection/cs_image_inspection_groups.ecsv > cs_group_list.txt


To build the plots in parallel:
parallel --bar -j 16 --joblog cs_image_plot.joblog --results cs_image_plot_logs python  ~/github/hapy/scripts/inspect_cs_images.py plot-one cs_image_inspection/cs_image_inspection_groups.ecsv {} --cutout-dir cutouts --outdir cs_image_inspection :::: cs_group_list.txt



To make VFS row-matched table
python ~/github/hapy/scripts/inspect_cs_images.py make-vfs-rowmatched merged_results_virgo_20260514_with_best_duplicate.fits --out hapy_vfs_rowmatched.fits

"""


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
    
    score = (0.3 * r_cost + 0.7 * h_cost)/1e-17

    # Strong extra penalty for large filter correction.
    # This keeps fcorr in the score continuously, but still flags risky cases.
    #if fcorr >= 1.2:
    #    score += 10.0 * (fcorr - 1.2)

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

def write_group_table(results, outfile):
    # read merged table
    # group by GALID / VFID / full_galname
    # determine best_idx / BEST flag
    # write one row per observation
    pass


def plot_group_from_table(group_table, galid, cutout_dir, outdir):
    tab = Table.read(group_table)

    rows = tab[tab["GALID"] == galid]

    best = np.where(rows["BEST"])[0]
    best_idx = int(best[0]) if len(best) > 0 else 0

    return plot_observation_group(
        rows=rows,
        best_idx=best_idx,
        cutout_dir=cutout_dir,
        outdir=outdir,
        galid=galid,
        norms=None,
        mark_best=True,
    )


def list_groups(group_table):
    tab = Table.read(group_table)
    for galid in np.unique(tab["GALID"]):
        print(galid)
        
def plot_observation_group(rows, best_idx, cutout_dir, outdir, galid, norms=None, mark_best=True):
    """
    Duplicate / continuum-subtraction comparison plot.

    Left diagnostic column:
      1. Legacy JPG image
      2. smoothed g-r image
      3. delta_mag = Halpha_mag - R_mag

    Observation columns:
      1. R-band image, asinh stretch
      2. Halpha CS-ZP image, shared zscale with CS-gr
      3. Halpha CS-gr image, shared zscale with CS-ZP, only if any CS-gr image exists

    Best duplicate is outlined in green.
    """

    if norms is None:
        norms = {}

    nobs = len(rows)
    ncols = nobs + 1   # extra left diagnostic column
    full_galname = full_galname_from_tag(str(rows[0]["TAG"]))

    # ------------------------------------------------------------
    # Check whether any duplicate has a CS-gr image
    # ------------------------------------------------------------
    csgr_paths = []
    for row in rows:
        tag = str(row["TAG"])
        csgr_path = find_image(
            cutout_dir,
            tag,
            [
                "-CS-gr.fits",
                "_CS-gr.fits",
                "-CS-GR.fits",
                "_CS-GR.fits",
                "-CS-g-r.fits",
                "_CS-g-r.fits",
            ],
        )
        csgr_paths.append(csgr_path)

    has_csgr = any(p is not None for p in csgr_paths)
    nrows = 3 if has_csgr else 2

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.0 * nrows),
        squeeze=False,
        constrained_layout=True,
    )

    # ------------------------------------------------------------
    # Use first row/tag for left diagnostic column
    # ------------------------------------------------------------
    first_tag = str(rows[0]["TAG"])
    first_cutdir = Path(cutout_dir) / first_tag

    legacy_jpg = None
    legacy_dir = first_cutdir / "legacy"
    if legacy_dir.exists():
        jpgs = sorted(list(legacy_dir.glob("*.jpg")) + list(legacy_dir.glob("*.jpeg")) + list(legacy_dir.glob("*.png")))
        if len(jpgs) > 0:
            legacy_jpg = jpgs[0]

    gr_path = find_image(
        legacy_dir,
        first_tag,
        [
            "-gr-ha-smooth.fits",
            "_gr-ha-smooth.fits",
            "-gr-smooth.fits",
            "_gr-smooth.fits",
        ],
    )

    # fallback: search directly in legacy directory
    gr_path = None
    if legacy_dir.exists():
        matches = sorted(legacy_dir.glob("*gr-smooth.fits"))
        if len(matches) > 0:
            gr_path = matches[0]

    delta_path = find_image(
        cutout_dir,
        first_tag,
        [
            "-CS-gr-delta-mag.fits",
            "_CS-gr-delta-mag.fits",
        ],
    )

    # fallback: exact current naming convention
    if delta_path is None:
        candidate = first_cutdir / f"{first_tag}-CS-gr-delta-mag.fits"
        if candidate.exists():
            delta_path = candidate

    gr_img = read_image(gr_path)
    delta_img = read_image(delta_path)

    # ------------------------------------------------------------
    # Left diagnostic column
    # ------------------------------------------------------------

    # Row 1: Legacy JPG
    ax = axes[0, 0]
    if legacy_jpg is not None:
        try:
            legacy_img = plt.imread(legacy_jpg)
            ax.imshow(legacy_img, origin="upper")
        except Exception:
            ax.text(
                0.5,
                0.5,
                "could not read\nLegacy JPG",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
    else:
        ax.text(
            0.5,
            0.5,
            "missing\nLegacy JPG",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )

    ax.set_title("Diagnostics", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    add_panel_text(ax, "Legacy JPG")

    # Row 2: smoothed g-r
    ax = axes[1, 0]
    if gr_img is not None:
        im = ax.imshow(
            gr_img,
            origin="lower",
            cmap="viridis",
            norm=image_norm(gr_img, "zscale"),
            #vmin=-0.5,
            #vmax=1.2,
        )
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    else:
        ax.text(
            0.5,
            0.5,
            "missing\ng-r image",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )

    ax.set_xticks([])
    ax.set_yticks([])
    add_panel_text(ax, "smoothed g-r")

    # Row 3: delta_mag
    if has_csgr:
        ax = axes[2, 0]
        if delta_img is not None:
            im = ax.imshow(
                delta_img,
                origin="lower",
                cmap="viridis",
                norm=image_norm(delta_img, "zscale"),
                #vmin=-0.25,
                #vmax=0.15,
            )
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        else:
            ax.text(
                0.5,
                0.5,
                "missing\nCS-gr delta-mag",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )

        ax.set_xticks([])
        ax.set_yticks([])
        add_panel_text(ax, r"$\Delta m = H\alpha - R$")

    # ------------------------------------------------------------
    # Read CS images first so CS-ZP and CS-gr can share zscale
    # ------------------------------------------------------------
    cszp_imgs = []
    csgr_imgs = []

    for j, row in enumerate(rows):
        tag = str(row["TAG"])

        cszp_path = find_image(
            cutout_dir,
            tag,
            ["-CS-ZP.fits", "_CS-ZP.fits", "-CS.fits", "_CS.fits"],
        )

        cszp_imgs.append(read_image(cszp_path))
        csgr_imgs.append(read_image(csgr_paths[j]))

    # Shared CS display norm across all CS-ZP and CS-gr images in this group
    cs_values = []
    for img in cszp_imgs + csgr_imgs:
        if img is None:
            continue
        good = np.isfinite(img)
        if np.any(good):
            cs_values.append(img[good])

    # if len(cs_values) > 0:
    #     cs_sample = np.concatenate(cs_values)
    #     cs_norm = image_norm(cs_sample, "zscale")
    # else:
    #     cs_norm = None

    # ------------------------------------------------------------
    # Observation columns start at column 1
    # ------------------------------------------------------------
    for j, row in enumerate(rows):

    
        col = j + 1
        tag = str(row["TAG"])

        r_path = find_image(
            cutout_dir,
            tag,
            ["-R.fits", "_R.fits", "-r.fits"],
        )

        cszp_path = find_image(
            cutout_dir,
            tag,
            ["-CS-ZP.fits", "_CS-ZP.fits", "-CS.fits", "_CS.fits"],
        )

        csgr_path = csgr_paths[j]

        r_img = read_image(r_path)
        cszp_img = cszp_imgs[j]
        csgr_img = csgr_imgs[j]


        cs_norm = None
        cs_pair_values = []

        for img in [cszp_img, csgr_img]:
            if img is None:
                continue
            good = np.isfinite(img)
            if np.any(good):
                cs_pair_values.append(img[good])

        if len(cs_pair_values) > 0:
            cs_pair_sample = np.concatenate(cs_pair_values)
            cs_norm = image_norm(cs_pair_sample, "zscale")
        

        # ------------------------------------------------------------
        # Shared display limits for this observation
        # ------------------------------------------------------------
        limits = None
        if r_img is not None:
            limits = get_display_limits_from_row(row, r_img.shape, buffer_pix=125)

        # ============================================================
        # Row 1: R band
        # ============================================================
        ax = axes[0, col]

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

        r_fwhm_n = safe_float(get_col(row, ["R_FWHM_PSF_NORM"]))
        r_sky_n = safe_float(get_col(row, ["R_SKYSTD_PHYS_NORM"]))

        #r_fwhm_n = normalized_value(r_fwhm, "R_FWHM_PSF", norms)
        #r_sky_n = normalized_value(r_sky, "R_SKYSTD_PHYS", norms)

        add_panel_text(
            ax,
            f"R\n"
            f"FWHM={r_fwhm:.2f} ({r_fwhm_n:.2f}x)\n"
            f"sky={r_sky_n:.2f}x",
        )

        # ============================================================
        # Row 2: CS-ZP
        # ============================================================
        ax = axes[1, col]

        if cszp_img is not None:
            ax.imshow(
                cszp_img,
                origin="lower",
                cmap="gray",
                norm=cs_norm if cs_norm is not None else image_norm(cszp_img, "zscale"),
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

        #h_fwhm_n = normalized_value(h_fwhm, "H_FWHM_PSF", norms)
        #h_sky_n = normalized_value(h_sky, "H_SKYSTD_PHYS", norms)

        h_fwhm_n = safe_float(get_col(row, ["H_FWHM_PSF_NORM"]))
        h_sky_n = safe_float(get_col(row, ["H_SKYSTD_PHYS_NORM"]))

        add_panel_text(
            ax,
            f"CS-ZP\n"
            f"H FWHM={h_fwhm:.2f} ({h_fwhm_n:.2f}x)\n"
            f"H sky={h_sky_n:.2f}x\n"
            f"fcorr={fcorr:.2f}",
        )

        # ============================================================
        # Row 3: CS-gr, optional
        # ============================================================
        if has_csgr:
            ax = axes[2, col]

            if csgr_img is not None:
                ax.imshow(
                    csgr_img,
                    origin="lower",
                    cmap="gray",
                    norm=cs_norm if cs_norm is not None else image_norm(csgr_img, "zscale"),
                )
            else:
                ax.text(
                    0.5,
                    0.5,
                    "missing CS-gr image",
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

            add_panel_text(
                ax,
                f"CS-gr\n"
                f"H FWHM={h_fwhm:.2f} ({h_fwhm_n:.2f}x)\n"
                f"H sky={h_sky_n:.2f}x\n"
                f"fcorr={fcorr:.2f}",
            )

        # ------------------------------------------------------------
        # Mark best duplicate
        # ------------------------------------------------------------
        if mark_best and j == best_idx:
            for rownum in range(nrows):
                mark_best_panel(axes[rownum, col])

            axes[0, col].text(
                0.5,
                1.08,
                "BEST",
                transform=axes[0, col].transAxes,
                ha="center",
                va="bottom",
                fontsize=13,
                color="lime",
                fontweight="bold",
            )

    if len(rows) == 1:
        fig.suptitle(f"Continuum subtraction inspection: {full_galname}", fontsize=16)
    else:
        fig.suptitle(f"Observation comparison: {full_galname}", fontsize=16)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    outfile = outdir / f"nobs{len(rows)}_{full_galname}_observation_comparison.png"

    fig.savefig(outfile, dpi=150)
    plt.close(fig)

    return outfile
        
#def plot_duplicate_group(rows, best_idx, cutout_dir, outdir, galid, norms=None):
def plot_observation_group_v1(rows, best_idx, cutout_dir, outdir, galid, norms=None, mark_best=True):

    """
    Duplicate comparison plot.

    Rows:
      1. R-band image, asinh stretch
      2. Halpha CS-ZP image, zscale stretch
      3. Halpha CS-gr image, zscale stretch, only if any CS-gr image exists

    Best duplicate is outlined in green.
    """
    if norms is None:
        norms = {}

    n = len(rows)
    full_galname = full_galname_from_tag(str(rows[0]["TAG"]))

    # Check whether any duplicate has a CS-gr image
    csgr_paths = []
    for row in rows:
        tag = str(row["TAG"])
        csgr_path = find_image(
            cutout_dir,
            tag,
            [
                "-CS-gr.fits",
                "_CS-gr.fits",
                "-CS-GR.fits",
                "_CS-GR.fits",
                "-CS-g-r.fits",
                "_CS-g-r.fits",
            ],
        )
        csgr_paths.append(csgr_path)

    has_csgr = any(p is not None for p in csgr_paths)
    nrows = 3 if has_csgr else 2

    fig, axes = plt.subplots(
        nrows,
        n,
        figsize=(4.2 * n, 4.0 * nrows),
        squeeze=False,
        constrained_layout=True,
    )

    for j, row in enumerate(rows):
        tag = str(row["TAG"])

        r_path = find_image(
            cutout_dir,
            tag,
            ["-R.fits", "_R.fits", "-r.fits"],
        )

        cszp_path = find_image(
            cutout_dir,
            tag,
            ["-CS-ZP.fits", "_CS-ZP.fits", "-CS.fits", "_CS.fits"],
        )

        csgr_path = csgr_paths[j]

        r_img = read_image(r_path)
        cszp_img = read_image(cszp_path)
        csgr_img = read_image(csgr_path)

        # ------------------------------------------------------------
        # Shared display limits for this observation
        # ------------------------------------------------------------
        limits = None
        if r_img is not None:
            limits = get_display_limits_from_row(row, r_img.shape, buffer_pix=125)

        # ============================================================
        # Row 1: R band
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
            f"R\n"
            f"FWHM={r_fwhm:.2f} ({r_fwhm_n:.2f}x)\n"
            f"sky={r_sky_n:.2f}x",
        )

        # ============================================================
        # Row 2: CS-ZP
        # ============================================================
        ax = axes[1, j]

        if cszp_img is not None:
            ax.imshow(
                cszp_img,
                origin="lower",
                cmap="gray",
                norm=image_norm(cszp_img, "zscale"),
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
            f"CS-ZP\n"
            f"H FWHM={h_fwhm:.2f} ({h_fwhm_n:.2f}x)\n"
            f"H sky={h_sky_n:.2f}x\n"
            f"fcorr={fcorr:.2f}",
        )

        # ============================================================
        # Row 3: CS-gr, optional
        # ============================================================
        if has_csgr:
            ax = axes[2, j]

            if csgr_img is not None:
                ax.imshow(
                    csgr_img,
                    origin="lower",
                    cmap="gray",
                    norm=image_norm(csgr_img, "zscale"),
                )
            else:
                ax.text(
                    0.5,
                    0.5,
                    "missing CS-gr image",
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

            add_panel_text(
                ax,
                f"CS-gr\n"
                f"H FWHM={h_fwhm:.2f} ({h_fwhm_n:.2f}x)\n"
                f"H sky={h_sky_n:.2f}x\n"
                f"fcorr={fcorr:.2f}",
            )

        # ------------------------------------------------------------
        # Mark best duplicate
        # ------------------------------------------------------------
        if mark_best and j == best_idx:
            for rownum in range(nrows):
                mark_best_panel(axes[rownum, j])

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

    #fig.suptitle(f"Duplicate comparison: {full_galname}", fontsize=16)
    if len(rows) == 1:
        fig.suptitle(f"Continuum subtraction inspection: {full_galname}", fontsize=16)
    else:
        fig.suptitle(f"Observation comparison: {full_galname}", fontsize=16)
    
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    #outfile = outdir / f"{full_galname}_observation_comparison.png"
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)
    outfile = outdir / f"nobs{len(rows)}_{full_galname}_observation_comparison.png"
    #outfile = outdir / f"{full_galname}_nobs{len(rows)}_comparison.png"

    fig.savefig(outfile, dpi=150)
    plt.close(fig)

    return outfile


def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest="mode", required=True)

    # ------------------------------------------------------------
    # make-table mode
    # ------------------------------------------------------------
    p_table = subparsers.add_parser(
        "make-table",
        help="Write duplicate-selection/group tables from merged_results.",
    )

    p_table.add_argument("merged_results", help="merged_results_*.fits file")
    p_table.add_argument("--outdir", default="cs_image_inspection")
    p_table.add_argument("--min-dups", type=int, default=1)
    p_table.add_argument("--testing", action="store_true")

    # ------------------------------------------------------------
    # plot-all mode
    # ------------------------------------------------------------
    p_plot_all = subparsers.add_parser(
        "plot-all",
        help="Plot all groups from saved group table.",
    )

    p_plot_all.add_argument("group_table", help="Group table written by make-table")
    p_plot_all.add_argument("--cutout-dir", default="cutouts")
    p_plot_all.add_argument("--outdir", default="cs_image_inspection")
    p_plot_all.add_argument("--testing", action="store_true")

    # ------------------------------------------------------------
    # plot-one mode
    # ------------------------------------------------------------
    p_plot_one = subparsers.add_parser(
        "plot-one",
        help="Plot one group from saved group table.",
    )

    p_plot_one.add_argument("group_table", help="Group table written by make-table")
    p_plot_one.add_argument("galid", help="DUP_GALID to plot")
    p_plot_one.add_argument("--cutout-dir", default="cutouts")
    p_plot_one.add_argument("--outdir", default="cs_image_inspection")

    # ------------------------------------------------------------
    # list-groups mode
    # ------------------------------------------------------------
    p_list = subparsers.add_parser(
        "list-groups",
        help="Print unique DUP_GALID values from saved group table.",
    )

    p_list.add_argument("group_table", help="Group table written by make-table")

    args = parser.parse_args()


    # ============================================================
    # make-table
    # ============================================================
    if args.mode == "make-table":
        tab = Table.read(args.merged_results)

        if "TAG" not in tab.colnames:
            raise ValueError("Expected a TAG column in merged_results table.")

        norms = {
            "R_FWHM_PSF": finite_median(tab, "R_FWHM_PSF"),
            "H_FWHM_PSF": finite_median(tab, "H_FWHM_PSF"),
            "R_SKYSTD_PHYS": finite_median(tab, "R_SKYSTD_PHYS"),
            "H_SKYSTD_PHYS": finite_median(tab, "H_SKYSTD_PHYS"),
        }

        galids = np.array([infer_galid(row) for row in tab])
        tab["DUP_GALID"] = galids

        # ------------------------------------------------------------
        # Initialize merged-results duplicate columns
        # ------------------------------------------------------------
        nrows_total = len(tab)

        tab["BEST_DUPLICATE"] = np.zeros(nrows_total, dtype=bool)
        tab["USE_FOR_SCIENCE"] = np.zeros(nrows_total, dtype=bool)
        tab["DUP_SCORE"] = np.full(nrows_total, np.nan)
        tab["N_DUP"] = np.zeros(nrows_total, dtype=int)

        tab["BEST_TAG"] = np.array([""] * nrows_total, dtype="U120")
        tab["USE_TAG"] = np.array([""] * nrows_total, dtype="U120")
        tab["DUP_NOTES"] = np.array([""] * nrows_total, dtype="U200")
        tab["MANUAL_OVERRIDE"] = np.zeros(nrows_total, dtype=bool)

        best_rows = []
        group_rows = []

        for galid in sorted(set(galids)):
            idx = np.where(galids == galid)[0]

            if len(idx) < args.min_dups:
                continue

            rows = tab[idx]

            if len(idx) > 1:
                scores = np.array([score_duplicate(row, norms=norms) for row in rows])
                best_local = int(np.nanargmin(scores))
                mark_best = True
            else:
                scores = np.array([np.nan])
                best_local = 0
                mark_best = False

            best_global = idx[best_local]

            manual_tag = MANUAL_BEST_TAG.get(galid, None)

            if manual_tag is not None:
                matches = np.where(
                    np.array([str(tab[i]["TAG"]) for i in idx]) == manual_tag
                )[0]

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

            best_tag = str(tab[best_global]["TAG"])

            # ------------------------------------------------------------
            # Update merged-results table
            # ------------------------------------------------------------
            for k, global_i in enumerate(idx):
                tab["BEST_DUPLICATE"][global_i] = (k == best_local)
                tab["USE_FOR_SCIENCE"][global_i] = (k == best_local)

                tab["DUP_SCORE"][global_i] = scores[k]
                tab["N_DUP"][global_i] = len(idx)

                tab["BEST_TAG"][global_i] = best_tag
                tab["USE_TAG"][global_i] = best_tag

                tab["MANUAL_OVERRIDE"][global_i] = manual_override
                tab["DUP_NOTES"][global_i] = override_note

            best_rows.append(
                {
                    "DUP_GALID": galid,
                    "N_DUP": len(idx),
                    "BEST_TAG": best_tag,
                    "BEST_SCORE": scores[best_local],
                    "ALL_TAGS": ",".join(str(tab[i]["TAG"]) for i in idx),
                    "ALL_SCORES": ",".join(
                        f"{s:.4f}" if np.isfinite(s) else "nan"
                        for s in scores
                    ),
                    "USE_TAG": best_tag,
                    "MANUAL_OVERRIDE": manual_override,
                    "NOTES": override_note,
                }
            )

            for k, global_i in enumerate(idx):
                group_row = {
                    "DUP_GALID": galid,
                    "N_DUP": len(idx),
                    "TAG": str(tab[global_i]["TAG"]),
                    "SCORE": scores[k],
                    "BEST_DUPLICATE": k == best_local,
                    "BEST_TAG": best_tag,
                    "USE_TAG": best_tag,
                    "MANUAL_OVERRIDE": manual_override,
                    "NOTES": override_note,
                }

                # Keep columns needed by plot_observation_group()
                plot_cols = [
                    "R_FWHM_PSF",
                    "R_FWHM_PSF_ARCSEC",
                    "H_FWHM_PSF",
                    "H_FWHM_PSF_ARCSEC",
                    "R_SKYSTD_PHYS",
                    "H_SKYSTD_PHYS",
                    "FILTER_CORRECTION",
                    "ELLIP_XCENTROID",
                    "ELLIP_YCENTROID",
                    "ELL0_XCENTROID",
                    "ELL0_YCENTROID",
                    "GAL_XC",
                    "GAL_YC",
                    "GAL_CXC",
                    "GAL_CYC",
                    "R24_PIX",
                    "R25_PIX",
                    "SMA_PIX",
                    "BA",
                    "PA",
                    "PIXSCALE",
                ]

                for col in plot_cols:
                    if col in tab.colnames:
                        group_row[col] = tab[global_i][col]

                # ------------------------------------------------------------
                # Store normalized values for plotting
                # ------------------------------------------------------------
                norm_cols = [
                    ("R_FWHM_PSF", "R_FWHM_PSF_NORM"),
                    ("H_FWHM_PSF", "H_FWHM_PSF_NORM"),
                    ("R_SKYSTD_PHYS", "R_SKYSTD_PHYS_NORM"),
                    ("H_SKYSTD_PHYS", "H_SKYSTD_PHYS_NORM"),
                ]

                for raw_col, norm_col in norm_cols:
                    if raw_col in tab.colnames:
                        val = safe_float(tab[global_i][raw_col])
                        med = safe_float(norms.get(raw_col, np.nan))

                        group_row[norm_col] = (
                            val / med
                            if np.isfinite(val)
                            and np.isfinite(med)
                            and med != 0
                            else np.nan
                        )

                group_rows.append(group_row)

            print(f"{galid}: best = {best_tag}")

            if args.testing:
                break

        best_tab = Table(rows=best_rows)
        group_tab = Table(rows=group_rows)

        outdir = Path(args.outdir)
        outdir.mkdir(exist_ok=True, parents=True)

        best_ecsv = outdir / "best_duplicates.ecsv"
        best_csv = outdir / "best_duplicates.csv"

        group_ecsv = outdir / "cs_image_inspection_groups.ecsv"
        group_csv = outdir / "cs_image_inspection_groups.csv"

        infile = Path(args.merged_results)
        merged_outfile = infile.parent / f"{infile.stem}_with_best_duplicate{infile.suffix}"
        #merged_outfile = outdir / "merged_results_with_best_duplicate.fits"

        # ------------------------------------------------------------
        # Write outputs
        # ------------------------------------------------------------
        best_tab.write(best_ecsv, format="ascii.ecsv", overwrite=True)
        best_tab.write(best_csv, format="ascii.csv", overwrite=True)

        group_tab.write(group_ecsv, format="ascii.ecsv", overwrite=True)
        group_tab.write(group_csv, format="ascii.csv", overwrite=True)

        tab.write(merged_outfile, overwrite=True)

        print(f"\nWrote {len(best_tab)} best-duplicate rows:")
        print(f"  {best_ecsv}")
        print(f"  {best_csv}")

        print(f"\nWrote {len(group_tab)} group rows:")
        print(f"  {group_ecsv}")
        print(f"  {group_csv}")

        print(f"\nWrote merged table with duplicate selection:")
        print(f"  {merged_outfile}")

        return

    
 

    # ============================================================
    # list-groups
    # ============================================================
    if args.mode == "list-groups":
        group_tab = Table.read(args.group_table)

        if "DUP_GALID" not in group_tab.colnames:
            raise ValueError("Expected DUP_GALID column in group table.")

        for galid in sorted(set(np.array(group_tab["DUP_GALID"]).astype(str))):
            print(galid)

        return

    # ============================================================
    # plot-one / plot-all
    # ============================================================
    group_tab = Table.read(args.group_table)

    if "DUP_GALID" not in group_tab.colnames:
        raise ValueError("Expected DUP_GALID column in group table.")

    if args.mode == "plot-one":
        galids_to_plot = [args.galid]
    else:
        galids_to_plot = sorted(set(np.array(group_tab["DUP_GALID"]).astype(str)))

    pngs = []

    for galid in galids_to_plot:
        rows = group_tab[np.array(group_tab["DUP_GALID"]).astype(str) == str(galid)]

        if len(rows) == 0:
            print(f"WARNING: no rows found for {galid}")
            continue

        if "BEST_DUPLICATE" in rows.colnames:
            best = np.where(np.array(rows["BEST_DUPLICATE"], dtype=bool))[0]
            best_idx = int(best[0]) if len(best) > 0 else 0
            mark_best = len(best) > 0
        else:
            best_idx = 0
            mark_best = False

        png = plot_observation_group(
            rows=rows,
            best_idx=best_idx,
            cutout_dir=args.cutout_dir,
            outdir=args.outdir,
            galid=galid,
            norms=None,
            mark_best=mark_best,
        )

        pngs.append(png)
        print(f"{galid}: wrote {png}")

        if getattr(args, "testing", False):
            break

    print(f"\nWrote {len(pngs)} plots")
    
 

 


if __name__ == "__main__":
    main()
