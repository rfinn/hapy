#!/usr/bin/env python

"""
qc_results.py

Survey-level quality control for merged HAPY results.

This script:
- reads merged_results.fits
- summarizes pipeline flags
- defines useful QC/sample masks
- writes clean/problem subsets
- generates a compact set of QC plots

Example
-------
python qc_results.py merged_results.fits --outdir qc
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from hapy.utils.plotting import raincloud_by_group
from qc_helpers import safe_bool_array, safe_float_array, first_existing_col, first_populated_col
from qc_helpers import build_row_qc_flags, ensure_dir, median_and_mad

from hapy.utils.plotting import QC_TIER_ORDER, QC_TIER_PALETTE

# ----------------------------------------------------------------------
# some science
# ----------------------------------------------------------------------
def get_ratio(a,b):
    return a/b

def add_science_columns(tab):
    """
    Add core science-ready derived columns to merged HAPY results table.

    Columns added:
      - H50_R50_RATIO
      - H75_R75_RATIO
      - H_MAXDET_R25_RATIO
      - H_PETRO_R50_RATIO
      - DELTA_GINI
      - DELTA_M20

    Safe against missing columns and divide-by-zero.
    """

    import numpy as np

    def safe_float(colname):
        if colname not in tab.colnames:
            return np.full(len(tab), np.nan)
        out = np.full(len(tab), np.nan)
        for i, v in enumerate(tab[colname]):
            try:
                out[i] = float(v)
            except Exception:
                pass
        return out

    def safe_ratio(num, den):
        out = np.full(len(num), np.nan)
        good = np.isfinite(num) & np.isfinite(den) & (den > 0)
        out[good] = num[good] / den[good]
        return out

    # --- load columns ---
    R50 = safe_float("R50_ARCSEC")
    H50 = safe_float("H50_ARCSEC")

    R75 = safe_float("R75_ARCSEC")
    H75 = safe_float("H75_ARCSEC")

    R25 = safe_float("R25_ARCSEC")
    HMAX = safe_float("H_MAXDET_ARCSEC")

    R_P50 = safe_float("R_PETRO_R50_ARCSEC")
    H_P50 = safe_float("H_PETRO_R50_ARCSEC")

    R_GINI = safe_float("R_HAPY_GINI")
    H_GINI = safe_float("H_HAPY_GINI")

    R_M20 = safe_float("R_HAPY_M20")
    H_M20 = safe_float("H_HAPY_M20")

    # --- compute ratios ---
    tab["H50_R50_RATIO"] = safe_ratio(H50, R50)
    tab["H75_R75_RATIO"] = safe_ratio(H75, R75)
    tab["H_MAXDET_R25_RATIO"] = safe_ratio(HMAX, R25)
    tab["H_PETRO_R50_RATIO"] = safe_ratio(H_P50, R_P50)

    # --- compute differences ---
    tab["DELTA_GINI"] = H_GINI - R_GINI
    tab["DELTA_M20"] = H_M20 - R_M20

    return tab

def add_qc_flags(tab):
    """
    Add QC warning and usability flags to HAPY merged results table.

    Adds:
      - WARN_MASK
      - WARN_BRIGHT_STAR
      - WARN_FILTER
      - WARN_ELLIPSE
      - WARN_WEAK_HA

      - USE_R_STRUCTURE   (0/1/2)
      - USE_HA_EXTENT     (0/1/2)
      - USE_HA_MORPH      (0/1/2)
      - USE_GALFIT        (0/1/2)

      - QC_TIER           (string: F/D/C/B/A)
    """

    import numpy as np

    n = len(tab)

    def safe_bool(colname):
        if colname not in tab.colnames:
            return np.zeros(n, dtype=bool)
        out = np.zeros(n, dtype=bool)
        for i, v in enumerate(tab[colname]):
            if isinstance(v, (bool, np.bool_)):
                out[i] = v
            else:
                s = str(v).strip().lower()
                out[i] = s in ("true", "t", "1", "yes", "y")
        return out

    def safe_float(colname):
        if colname not in tab.colnames:
            return np.full(n, np.nan)
        out = np.full(n, np.nan)
        for i, v in enumerate(tab[colname]):
            try:
                out[i] = float(v)
            except Exception:
                pass
        return out

    # -----------------------------
    # load columns
    # -----------------------------
    mask_warn = safe_bool("ELL0_MASK_WARN")
    bright_star = safe_bool("BRIGHT_STAR_FLAG")
    ellipse_mismatch = safe_bool("ELL_MISMATCH")

    filt = safe_float("FILTER_CORRECTION")

    r_maskfrac = safe_float("R_PROFILE_MASKFRAC_MAX")
    h_maskfrac = safe_float("H_PROFILE_MASKFRAC_MAX")

    h_npix = safe_float("H_HAPY_NPIX")
    h_ngood = safe_float("H_PROFILE_NGOOD")
    h_snr = safe_float("H_HAPY_SNP_DET")

    r50 = safe_float("R50_ARCSEC")
    h50 = safe_float("H50_ARCSEC")
    hmax = safe_float("H_MAXDET_ARCSEC")

    r_ngood = safe_float("R_PROFILE_NGOOD")

    phot_ok = safe_bool("PHOT_OK")
    rprof_ok = safe_bool("R_PROFILE_OK")
    hprof_ok = safe_bool("H_PROFILE_OK")

    morph_ok = safe_bool("HAPY_MORPH_OK")
    r_sm_ok = safe_bool("R_SM_OK")
    h_sm_ok = safe_bool("H_SM_OK")

    gal_nc_ok = safe_bool("GAL_NC_OK")
    gal_cv_ok = safe_bool("GAL_CV_OK")

    # -----------------------------
    # warnings
    # -----------------------------
    tab["WARN_MASK"] = (
        mask_warn |
        (np.isfinite(r_maskfrac) & (r_maskfrac > 0.3)) |
        (np.isfinite(h_maskfrac) & (h_maskfrac > 0.3))
    )

    tab["WARN_BRIGHT_STAR"] = bright_star

    tab["WARN_FILTER"] = np.isfinite(filt) & (filt > 1.2)

    tab["WARN_ELLIPSE"] = ellipse_mismatch

    tab["WARN_WEAK_HA"] = (
        (np.isfinite(h_npix) & (h_npix < 50)) |
        (np.isfinite(h_ngood) & (h_ngood < 8)) |
        (np.isfinite(h_snr) & (h_snr < 3))
    )

    # -----------------------------
    # usability flags (0/1/2)
    # -----------------------------
    USE_R = np.zeros(n, dtype=int)
    USE_HA = np.zeros(n, dtype=int)
    USE_HM = np.zeros(n, dtype=int)
    USE_GF = np.zeros(n, dtype=int)

    # R structure
    good_r = (
        phot_ok & rprof_ok &
        np.isfinite(r50) & (r50 > 0) &
        np.isfinite(r_ngood) & (r_ngood >= 20)
    )

    USE_R[good_r] = 2
    USE_R[good_r & tab["WARN_MASK"]] = 1

    # Hα extent
    good_ha = (
        phot_ok & hprof_ok &
        np.isfinite(h50) & (h50 > 0) &
        np.isfinite(hmax) & (hmax > 0) &
        np.isfinite(h_npix) & (h_npix >= 50) &
        np.isfinite(h_ngood) & (h_ngood >= 8) &
        (~tab["WARN_FILTER"])
    )

    USE_HA[good_ha] = 2
    USE_HA[good_ha & tab["WARN_WEAK_HA"]] = 1

    # Hα morphology
    good_hm = (
        morph_ok & h_sm_ok &
        np.isfinite(h_npix) & (h_npix >= 100) &
        (~tab["WARN_FILTER"])
    )

    USE_HM[good_hm] = 2
    USE_HM[good_hm & tab["WARN_WEAK_HA"]] = 1

    # GALFIT
    good_gf = gal_nc_ok | gal_cv_ok
    USE_GF[good_gf] = 2

    # -----------------------------
    # attach to table
    # -----------------------------
    tab["USE_R_STRUCTURE"] = USE_R
    tab["USE_HA_EXTENT"] = USE_HA
    tab["USE_HA_MORPH"] = USE_HM
    tab["USE_GALFIT"] = USE_GF

    # -----------------------------
    # QC tier
    # -----------------------------
    tier = np.full(n, "F", dtype="U1")

    technical_ok = phot_ok

    for i in range(n):
        if not technical_ok[i]:
            tier[i] = "F"
        elif USE_R[i] == 0 or USE_HA[i] == 0:
            tier[i] = "D"
        elif (
            tab["WARN_MASK"][i] or
            tab["WARN_BRIGHT_STAR"][i] or
            tab["WARN_FILTER"][i] or
            tab["WARN_ELLIPSE"][i]
        ):
            tier[i] = "C"
        elif USE_HM[i] < 2:
            tier[i] = "B"
        else:
            tier[i] = "A"

    tab["QC_TIER"] = tier

    return tab
# ----------------------------------------------------------------------
# flag discovery and masks
# ----------------------------------------------------------------------

def find_ok_columns(tab: Table) -> list[str]:
    return [c for c in tab.colnames if c.endswith("_OK") ]

def find_flag_columns(tab: Table) -> list[str]:
    return [c for c in tab.colnames if c.endswith("_FLAG")]

def build_qc_masks(tab: Table, max_ha_filter_correction: float = 1.2) -> dict[str, np.ndarray]:
    flags = build_row_qc_flags(tab, max_ha_filter_correction=max_ha_filter_correction)

    masks = {}
    masks["MASK_OK"] = flags["MASK_OK"]
    masks["PHOT_OK"] = flags["PHOT_OK"]
    masks["PSF_OK"] = flags["PSF_OK"]
    masks["R_PROFILE_OK"] = flags["R_PROFILE_OK"]
    masks["H_PROFILE_OK"] = flags["H_PROFILE_OK"]
    masks["R_SM_OK"] = flags["R_SM_OK"]
    masks["H_SM_OK"] = flags["H_SM_OK"]
    masks["GAL_NC_OK"] = flags["GAL_NC_OK"]
    masks["GAL_CV_OK"] = flags["GAL_CV_OK"]
    masks["FILTER_WARNING"] = flags["FILTER_WARNING"]

    masks["mask_phot_ok"] = flags["MASK_PHOT_OK"]
    masks["profile_ok"] = flags["PROFILE_OK"]
    masks["statmorph_ok"] = flags["R_SM_OK"] & flags["H_SM_OK"] & (~flags["FILTER_WARNING"])
    masks["galfit_any_ok"] = flags["GALFIT_ANY_OK"]
    masks["galfit_both_ok"] = flags["GALFIT_BOTH_OK"]
    masks["science_ready"] = flags["SCIENCE_READY"]
    masks["problem"] = flags["TECHNICAL_PROBLEM"] | flags["SCIENCE_PROBLEM"]

    return masks

# def build_qc_masks(tab: Table, max_ha_filter_correction: float = 1.2) -> dict[str, np.ndarray]:
#     masks = {}

#     # core flags
#     masks["PSF_OK"] = safe_bool_array(tab, "PSF_OK")
#     masks["MASK_OK"] = safe_bool_array(tab, "MASK_OK")
#     masks["PHOT_OK"] = safe_bool_array(tab, "PHOT_OK")
#     masks["HAPY_MORPH_OK"] = safe_bool_array(tab, "HAPY_MORPH_OK")    
#     masks["R_PROFILE_OK"] = safe_bool_array(tab, "R_PROFILE_OK")
#     masks["H_PROFILE_OK"] = safe_bool_array(tab, "H_PROFILE_OK")
#     masks["R_SM_OK"] = safe_bool_array(tab, "R_SM_OK")
#     masks["H_SM_OK"] = safe_bool_array(tab, "H_SM_OK")
#     masks["GAL_NC_OK"] = safe_bool_array(tab, "GAL_NC_OK")
#     masks["GAL_CV_OK"] = safe_bool_array(tab, "GAL_CV_OK")

#     filt_col = first_existing_col(tab, ["FILTER_CORRECTION", "FILT_COR"])
#     if filt_col is not None:
#         filtcor = safe_float_array(tab, filt_col)
#     else:
#         filtcor = np.full(len(tab), np.nan)

#     masks["FILTER_WARNING"] = np.isfinite(filtcor) & (filtcor > max_ha_filter_correction)

#     # science-oriented subsets
#     masks["mask_phot_ok"] = masks["MASK_OK"] & masks["PHOT_OK"]
#     masks["technical_problem"] = ~(masks["MASK_OK"] & masks["PHOT_OK"])
#     masks["science_problem"] = (
#         ~masks["profile_ok"] |
#         masks["FILTER_WARNING"] |
#         safe_bool_array(tab, "ELL0_MASK_WARN") |
#         safe_bool_array(tab, "BRIGHT_STAR_FLAG") |
#         safe_bool_array(tab, "ELL_MISMATCH")
#         )
#     masks["profile_ok"] = (
#         masks["MASK_OK"] &
#         masks["PHOT_OK"] &
#         masks["R_PROFILE_OK"] &
#         masks["H_PROFILE_OK"]
#     )

#     masks["statmorph_ok"] = (
#         masks["profile_ok"] &
#         masks["R_SM_FLAG"] &
#         masks["H_SM_FLAG"] &
#         (~masks["FILTER_WARNING"])
#     )

#     masks["galfit_any_ok"] = (
#         masks["MASK_OK"] &
#         masks["PHOT_OK"] &
#         (masks["GAL_NC_OK"] | masks["GAL_CV_OK"])
#     )

#     masks["galfit_both_ok"] = (
#         masks["MASK_OK"] &
#         masks["PHOT_OK"] &
#         masks["GAL_NC_OK"] &
#         masks["GAL_CV_OK"]
#     )

#     r50 = safe_float_array(tab, "R50_ARCSEC")
#     h50 = safe_float_array(tab, "H50_ARCSEC")
#     hmax = safe_float_array(tab, "H_MAXDET_ARCSEC")
#     hnpix = safe_float_array(tab, "H_HAPY_NPIX")
#     hngood = safe_float_array(tab, "H_PROFILE_NGOOD")
#     rngood = safe_float_array(tab, "R_PROFILE_NGOOD")

#     masks["r_structure_good"] = (
#         masks["PHOT_OK"] &
#         masks["R_PROFILE_OK"] &
#         np.isfinite(r50) & (r50 > 0) &
#         np.isfinite(rngood) & (rngood >= 20)
#         )

    
#     masks["ha_extent_good"] = (
#         masks["PHOT_OK"] &
#         masks["H_PROFILE_OK"] &
#         np.isfinite(h50) & (h50 > 0) &
#         np.isfinite(hmax) & (hmax > 0) &
#         np.isfinite(hnpix) & (hnpix >= 50) &
#         np.isfinite(hngood) & (hngood >= 8) &
#         (~masks["FILTER_WARNING"])
#         )
    
#     masks["science_ready"] = (
#         masks["r_structure_good"] &
#         masks["ha_extent_good"]
#         )
    
#     masks["morph_good"] = (
#         masks["HAPY_MORPH_OK"] &
#         masks["R_SM_OK"] &
#         masks["H_SM_OK"] &
#         (~masks["FILTER_WARNING"])
#     )

#     # TODO - may want to define this differently
#     masks["problem"] = masks["technical_problem"]

#     return masks


# ----------------------------------------------------------------------
# text summaries
# ----------------------------------------------------------------------

def write_text_summary(tab: Table, masks: dict[str, np.ndarray], outpath: Path, scheme) -> None:
    n = len(tab)
    with open(outpath, "w") as fh:
        fh.write("HAPY QC SUMMARY\n")
        fh.write("================\n\n")
        fh.write(f"Total rows: {n}\n")
        if scheme == "virgo":
            if "VFID" in tab.colnames:
                fh.write(f"Unique galaxies: {len(set(tab['VFID']))}")
        elif scheme == "agc":
            if "OBJID" in tab.colnames:
                fh.write(f"Unique galaxies: {len(set(tab['OBJID']))}")

        
        #if "VFID" in tab.colnames:
        #    fh.write(f"Unique VFIDs: {len(set(tab['VFID']))}\n")
        fh.write("\n")

        fh.write("Pipeline flags\n")
        fh.write("--------------\n")
        for f in find_ok_columns(tab):
            col = safe_bool_array(tab, f)
            ntrue = np.sum(col)
            nfalse = np.sum(~col)
            pct = 100.0 * ntrue / n if n > 0 else np.nan
            fh.write(f"{f:20s}: {ntrue:4d} OK  | {nfalse:4d} FAIL  ({pct:5.1f}%)\n")

        fh.write("Warning flags\n")
        fh.write("--------------\n")
        for f in find_flag_columns(tab):
            #col = safe_bool_array(tab, f)
            nzero = np.sum(col == 0)
            nfalse = np.sum(~col)
            pct = 100.0 * ntrue / n if n > 0 else np.nan
            fh.write(f"{f:20s}: {ntrue:4d} Clean  | {nfalse:4d} Flagged  ({pct:5.1f}%)\n")
            
        fh.write("\nQC subsets\n")
        fh.write("----------\n")
        for name in ["mask_phot_ok", "profile_ok", "statmorph_ok", "galfit_any_ok", "galfit_both_ok", "science_ready", "problem"]:
            m = masks[name]
            count = np.sum(m)
            pct = 100.0 * count / n if n > 0 else np.nan
            fh.write(f"{name:20s}: {count:4d} ({pct:5.1f}%)\n")

        if "STATUS" in tab.colnames:
            fh.write("\nSTATUS counts\n")
            fh.write("-------------\n")
            vals, counts = np.unique(np.array(tab["STATUS"]).astype(str), return_counts=True)
            for v, c in zip(vals, counts):
                fh.write(f"{v:20s}: {c}\n")

        if "STAGE" in tab.colnames:
            fh.write("\nSTAGE counts\n")
            fh.write("------------\n")
            vals, counts = np.unique(np.array(tab["STAGE"]).astype(str), return_counts=True)
            for v, c in zip(vals, counts):
                fh.write(f"{v:20s}: {c}\n")


# ----------------------------------------------------------------------
# plotting
# ----------------------------------------------------------------------
def plot_qc_tier_scatter(tab, xcol, ycol, outpath):
    x = safe_float_array(tab, xcol)
    y = safe_float_array(tab, ycol)
    tier = np.array(tab["QC_TIER"]).astype(str)

    good = np.isfinite(x) & np.isfinite(y)

    fig = plt.figure(figsize=(6, 5))
    ax = plt.gca()

    for t in QC_TIER_ORDER:
        m = good & (tier == t)
        if np.sum(m) == 0:
            continue
        ax.scatter(x[m], y[m], s=16, alpha=0.7, color=QC_TIER_PALETTE[t], label=t)

    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.legend(title="QC Tier")
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    

def plot_flag_completion(tab: Table, outpath: Path) -> None:
    flags = find_ok_columns(tab)
    if len(flags) == 0:
        return

    counts = []
    pcts = []
    n = len(tab)

    for f in flags:
        col = safe_bool_array(tab, f)
        counts.append(np.sum(col))
        pcts.append(100.0 * np.sum(col) / n if n > 0 else np.nan)

    fig = plt.figure(figsize=(12, 5))
    ax = plt.gca()
    ax.bar(range(len(flags)), pcts)
    ax.set_xticks(range(len(flags)))
    ax.set_xticklabels(flags, rotation=45, ha="right")
    ax.set_ylabel("Success fraction (%)")
    ax.set_title("Pipeline completion fractions")
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_hist(tab: Table, col: str, outpath: Path, title: str | None = None, logx: bool = False, mask: np.ndarray | None = None) -> None:
    if col not in tab.colnames:
        print(f"WARNING: missing column {col}, skipping {outpath.name}")
        return

    x = safe_float_array(tab, col)
    good = np.isfinite(x)
    if mask is not None:
        good &= mask

    if logx:
        good &= (x > 0)
        x = np.log10(x)

    if np.sum(good) == 0:
        print(f"WARNING: no finite values for {col}, skipping {outpath.name}")
        return

    fig = plt.figure(figsize=(6, 4))
    plt.hist(x[good], bins=30)
    plt.xlabel(f"log10({col})" if logx else col)
    plt.ylabel("N")
    plt.title(title if title is not None else col)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_compare(tab: Table, xcol: str, ycol: str, outpath: Path, title: str | None = None,
                 logx: bool = False, logy: bool = False, mask: np.ndarray | None = None) -> None:
    if xcol not in tab.colnames or ycol not in tab.colnames:
        print(f"WARNING: missing {xcol} or {ycol}, skipping {outpath.name}")
        return

    x = safe_float_array(tab, xcol)
    y = safe_float_array(tab, ycol)

    good = np.isfinite(x) & np.isfinite(y)
    if mask is not None:
        good &= mask
    if logx:
        good &= (x > 0)
    if logy:
        good &= (y > 0)

    if np.sum(good) == 0:
        print(f"WARNING: no finite comparison values for {xcol}, {ycol}")
        return

    fig = plt.figure(figsize=(5.5, 5.5))
    ax = plt.gca()
    ax.scatter(x[good], y[good], s=12, alpha=0.7)

    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    if not logx and not logy:
        lo = np.nanmin([np.nanmin(x[good]), np.nanmin(y[good])])
        hi = np.nanmax([np.nanmax(x[good]), np.nanmax(y[good])])
        ax.plot([lo, hi], [lo, hi], "k--", lw=1)

    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(title if title is not None else f"{ycol} vs {xcol}")
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_filter_warning_vs_ha(tab: Table, masks: dict[str, np.ndarray], outpath: Path) -> None:
    col = first_existing_col(tab, ["FILTER_CORRECTION", "FILT_COR"])
    if col is None or "H_TOT_FLUX_CGS" not in tab.colnames:
        return

    x = safe_float_array(tab, col)
    y = safe_float_array(tab, "H_TOT_FLUX_CGS")

    good = np.isfinite(x) & np.isfinite(y) & (y > 0)
    if np.sum(good) == 0:
        return

    warn = masks["FILTER_WARNING"]

    fig = plt.figure(figsize=(6, 5))
    ax = plt.gca()
    ax.scatter(x[good & ~warn], y[good & ~warn], s=12, alpha=0.7, label="clean")
    ax.scatter(x[good & warn], y[good & warn], s=12, alpha=0.7, label="filter warning")
    ax.set_yscale("log")
    ax.set_xlabel(col)
    ax.set_ylabel("H_TOT_FLUX_CGS")
    ax.set_title("Hα flux vs filter correction")
    ax.legend()
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_raincloud_by_telescope(tab, col, outpath, logx=False):
    if col not in tab.colnames or "TELESCOPE" not in tab.colnames:
        print(f"WARNING: missing {col} or TELESCOPE")
        return

    telescopes = sorted(set(np.array(tab["TELESCOPE"]).astype(str)))
    values = []

    x = np.array(tab[col], dtype=float)
    tel = np.array(tab["TELESCOPE"]).astype(str)

    for t in telescopes:
        vals = x[tel == t]
        vals = vals[np.isfinite(vals)]
        if logx:
            vals = vals[vals > 0]
            vals = np.log10(vals)
        values.append(vals)

    xlabel = f"log10({col})" if logx else col
    fig, ax = raincloud_by_group(values, telescopes, xlabel=xlabel, title=f"{col} by telescope")
    if fig is not None:
        fig.savefig(outpath, dpi=150)
        plt.close(fig)

def plot_failure_fraction_vs_bright_star_distance(tab, outpath):
    """
    Plot failure fraction vs normalized distance to nearest bright star.
    """

    required = [
        "BRIGHT_STAR_DIST_ARCSEC",
        "BRIGHT_STAR_NEAREST_MASKRAD_ARCSEC",
        "PHOT_OK",
        "R_PROFILE_OK",
        "H_PROFILE_OK",
        "R_SM_FLAG",
        "H_SM_FLAG",
        "GAL_NC_OK",
        "GAL_CV_OK",
    ]
    for col in required:
        if col not in tab.colnames:
            print(f"WARNING: missing {col}, skipping bright-star failure plot")
            return

    dist = np.array(tab["BRIGHT_STAR_DIST_ARCSEC"], dtype=float)
    rad = np.array(tab["BRIGHT_STAR_NEAREST_MASKRAD_ARCSEC"], dtype=float)

    good = np.isfinite(dist) & np.isfinite(rad) & (rad > 0)
    if np.sum(good) == 0:
        print("WARNING: no valid bright-star distance data")
        return

    x = dist[good] / rad[good]

    failure_defs = {
        "phot fail": ~np.array(tab["PHOT_OK"][good], dtype=bool),
        "r profile fail": ~np.array(tab["R_PROFILE_OK"][good], dtype=bool),
        "ha profile fail": ~np.array(tab["H_PROFILE_OK"][good], dtype=bool),
        "r statmorph fail": ~np.array(tab["R_SM_FLAG"][good], dtype=bool),
        "ha statmorph fail": ~np.array(tab["H_SM_FLAG"][good], dtype=bool),
        "galfit nc fail": ~np.array(tab["GAL_NC_OK"][good], dtype=bool),
        "galfit cv fail": ~np.array(tab["GAL_CV_OK"][good], dtype=bool),
    }

    bins = np.array([0, 0.5, 1, 1.5, 2, 3, 5, 10])
    xcen = 0.5 * (bins[:-1] + bins[1:])

    fig = plt.figure(figsize=(8, 6))
    ax = plt.gca()

    for label, failflag in failure_defs.items():
        frac = np.full(len(xcen), np.nan)

        for i in range(len(bins) - 1):
            m = (x >= bins[i]) & (x < bins[i+1])
            if np.sum(m) > 0:
                frac[i] = np.mean(failflag[m])

        ax.plot(xcen, frac, marker="o", label=label)

    ax.axvline(1.0, color="k", ls="--", lw=1)
    ax.set_xlabel("Nearest bright-star distance / mask radius")
    ax.set_ylabel("Failure fraction")
    ax.set_title("Pipeline failure rate vs bright-star proximity")
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_qc_dashboard_v1(tab, outpath):
    """
    Make a compact QC dashboard for HAPY merged results.

    Requires:
      - add_science_columns(tab)
      - add_qc_flags(tab)

    Output:
      - one PNG summary figure
    """
    import numpy as np
    import matplotlib.pyplot as plt

    n = len(tab)

    def safe_float(colname):
        if colname not in tab.colnames:
            return np.full(n, np.nan)
        out = np.full(n, np.nan)
        for i, v in enumerate(tab[colname]):
            try:
                out[i] = float(v)
            except Exception:
                pass
        return out

    def safe_str(colname, default=""):
        if colname not in tab.colnames:
            return np.full(n, default, dtype="U8")
        out = np.full(n, default, dtype="U8")
        for i, v in enumerate(tab[colname]):
            out[i] = str(v)
        return out

    tier = safe_str("QC_TIER", default="?")
    h50_r50 = safe_float("H50_R50_RATIO")
    hmax_r25 = safe_float("H_MAXDET_R25_RATIO")
    r50 = safe_float("R50_ARCSEC")
    h50 = safe_float("H50_ARCSEC")
    h_npix = safe_float("H_HAPY_NPIX")
    h_snr = safe_float("H_HAPY_SNP_DET")
    h_gini = safe_float("H_HAPY_GINI")
    h_m20 = safe_float("H_HAPY_M20")

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.ravel()

    # --------------------------------------------------
    # 1. QC tier counts
    # --------------------------------------------------
    ax = axes[0]
    #tiers = ["A", "B", "C", "D", "F"]
    #counts = [np.sum(tier == t) for t in tiers]

    tiers = [t for t in QC_TIER_ORDER if np.any(tier == t)]
    counts = [np.sum(tier == t) for t in tiers]
    colors = [QC_TIER_PALETTE[t] for t in tiers]
    ax.bar(range(len(tiers)), counts,color=colors)
    ax.set_xticks(range(len(tiers)))
    ax.set_xticklabels(tiers)
    ax.set_ylabel("N")
    ax.set_title("QC tier counts")

    # --------------------------------------------------
    # 2. H50 / R50
    # --------------------------------------------------
    ax = axes[1]
    good = np.isfinite(h50_r50) & (h50_r50 > 0)
    if np.sum(good) > 0:
        ax.hist(h50_r50[good], bins=30)
    ax.set_xlabel("H50_R50_RATIO")
    ax.set_ylabel("N")
    ax.set_title("Hα / stellar half-light ratio")

    # --------------------------------------------------
    # 3. Hmax / R25
    # --------------------------------------------------
    ax = axes[2]
    good = np.isfinite(hmax_r25) & (hmax_r25 > 0)
    if np.sum(good) > 0:
        ax.hist(hmax_r25[good], bins=30)
    ax.set_xlabel("H_MAXDET_R25_RATIO")
    ax.set_ylabel("N")
    ax.set_title("Hα max extent / R25")

    # --------------------------------------------------
    # 4. H50 vs R50
    # --------------------------------------------------
    ax = axes[3]
    good = np.isfinite(r50) & np.isfinite(h50) & (r50 > 0) & (h50 > 0)
    if np.sum(good) > 0:
        ax.scatter(r50[good], h50[good], s=12, alpha=0.7)
        lo = np.nanmin([np.nanmin(r50[good]), np.nanmin(h50[good])])
        hi = np.nanmax([np.nanmax(r50[good]), np.nanmax(h50[good])])
        ax.plot([lo, hi], [lo, hi], "k--", lw=1)
    ax.set_xlabel("R50_ARCSEC")
    ax.set_ylabel("H50_ARCSEC")
    ax.set_title("Hα vs stellar half-light radius")

    # --------------------------------------------------
    # 5. Halpha robustness
    # --------------------------------------------------
    ax = axes[4]
    good = np.isfinite(h_npix) & np.isfinite(h_snr) & (h_npix > 0)
    if np.sum(good) > 0:
        ax.scatter(h_npix[good], h_snr[good], s=12, alpha=0.7)
        ax.set_xscale("log")
    ax.set_xlabel("H_HAPY_NPIX")
    ax.set_ylabel("H_HAPY_SNP_DET")
    ax.set_title("Hα detection robustness")

    # --------------------------------------------------
    # 6. Halpha morphology sanity
    # --------------------------------------------------
    ax = axes[5]
    good = np.isfinite(h_gini) & np.isfinite(h_m20)
    if np.sum(good) > 0:
        ax.scatter(h_gini[good], h_m20[good], s=12, alpha=0.7)
    ax.set_xlabel("H_HAPY_GINI")
    ax.set_ylabel("H_HAPY_M20")
    ax.set_title("Hα morphology plane")

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_qc_dashboard(tab, outpath):
    """
    Make a compact QC dashboard for HAPY merged results.

    Requires:
      - add_science_columns(tab)
      - add_qc_flags(tab)

    Output:
      - one PNG summary figure
    """
    import numpy as np
    import matplotlib.pyplot as plt

    n = len(tab)

    def safe_float(colname):
        if colname not in tab.colnames:
            return np.full(n, np.nan)
        out = np.full(n, np.nan)
        for i, v in enumerate(tab[colname]):
            try:
                out[i] = float(v)
            except Exception:
                pass
        return out

    def safe_str(colname, default=""):
        if colname not in tab.colnames:
            return np.full(n, default, dtype="U8")
        out = np.full(n, default, dtype="U8")
        for i, v in enumerate(tab[colname]):
            out[i] = str(v)
        return out

    def scatter_by_tier(ax, x, y, tier, xlabel, ylabel, title, one_to_one=False, logx=False, logy=False):
        good = np.isfinite(x) & np.isfinite(y)
        if logx:
            good &= (x > 0)
        if logy:
            good &= (y > 0)

        if np.sum(good) == 0:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            return

        # Matplotlib default tab colors

        for t in QC_TIER_ORDER:
            m = good & (tier == t)
            if np.sum(m) == 0:
                continue
            ax.scatter(x[m], y[m], s=14, alpha=0.75, label=t, color=QC_TIER_PALETTE[t])

        if logx:
            ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")

        if one_to_one and np.sum(good) > 0:
            if logx or logy:
                lo = np.nanmin([np.nanmin(x[good]), np.nanmin(y[good])])
                hi = np.nanmax([np.nanmax(x[good]), np.nanmax(y[good])])
                if np.isfinite(lo) and np.isfinite(hi) and lo > 0 and hi > 0:
                    ax.plot([lo, hi], [lo, hi], "k--", lw=1)
            else:
                lo = np.nanmin([np.nanmin(x[good]), np.nanmin(y[good])])
                hi = np.nanmax([np.nanmax(x[good]), np.nanmax(y[good])])
                ax.plot([lo, hi], [lo, hi], "k--", lw=1)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    tier = safe_str("QC_TIER", default="?")
    r50 = safe_float("R50_ARCSEC")
    h50 = safe_float("H50_ARCSEC")
    h_npix = safe_float("H_HAPY_NPIX")
    h_snr = safe_float("H_HAPY_SNP_DET")
    h_gini = safe_float("H_HAPY_GINI")
    r_gini = safe_float("R_HAPY_GINI")
    h_m20 = safe_float("H_HAPY_M20")
    r_m20 = safe_float("R_HAPY_M20")

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.ravel()

    # --------------------------------------------------
    # 1. QC tier counts
    # --------------------------------------------------
    ax = axes[0]

    tiers = [t for t in QC_TIER_ORDER if np.any(tier == t)]
    counts = [np.sum(tier == t) for t in tiers]
    colors = [QC_TIER_PALETTE[t] for t in tiers]
    ax.bar(range(len(tiers)), counts, color=colors)
    ax.set_xticks(range(len(tiers)))
    ax.set_xticklabels(tiers)
    ax.set_ylabel("N")
    ax.set_title("QC tier counts")

    # --------------------------------------------------
    # 2. H50 vs R50
    # --------------------------------------------------
    scatter_by_tier(
        axes[1], r50, h50, tier,
        xlabel="R50_ARCSEC",
        ylabel="H50_ARCSEC",
        title="Hα vs stellar half-light radius",
        one_to_one=True,
    )

    # --------------------------------------------------
    # 3. Halpha robustness
    # --------------------------------------------------
    scatter_by_tier(
        axes[2], h_npix, h_snr, tier,
        xlabel="H_HAPY_NPIX",
        ylabel="H_HAPY_SNP_DET",
        title="Hα detection robustness",
        logx=True,
    )

    # --------------------------------------------------
    # 4. Halpha vs stellar Gini
    # --------------------------------------------------
    scatter_by_tier(
        axes[3], r_gini, h_gini, tier,
        xlabel="R_HAPY_GINI",
        ylabel="H_HAPY_GINI",
        title="Hα vs stellar Gini",
        one_to_one=True,
    )

    # --------------------------------------------------
    # 5. Halpha vs stellar M20
    # --------------------------------------------------
    scatter_by_tier(
        axes[4], r_m20, h_m20, tier,
        xlabel="R_HAPY_M20",
        ylabel="H_HAPY_M20",
        title="Hα vs stellar M20",
        one_to_one=True,
    )

    # --------------------------------------------------
    # 6. Halpha morphology plane
    # --------------------------------------------------
    scatter_by_tier(
        axes[5], h_gini, h_m20, tier,
        xlabel="H_HAPY_GINI",
        ylabel="H_HAPY_M20",
        title="Hα morphology plane",
    )

    # one legend for all scatter panels
    handles, labels = axes[5].get_legend_handles_labels()
    if len(handles) > 0:
        fig.legend(handles, labels, title="QC_TIER", loc="upper right")

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
# ----------------------------------------------------------------------
# subset writing
# ----------------------------------------------------------------------

def write_subsets(tab: Table, masks: dict[str, np.ndarray], outdir: Path) -> None:
    subset_names = [
        "mask_phot_ok",
        "profile_ok",
        "statmorph_ok",
        "galfit_any_ok",
        "galfit_both_ok",
        "science_ready",
        "problem",
    ]
    for name in subset_names:
        sub = tab[masks[name]]
        sub.write(outdir / f"{name}.fits", format="fits", overwrite=True)

# ----------------------------------------------------------------------
# output tables
# ----------------------------------------------------------------------
def write_qc_summary_table(tab, outpath):
    """
    Write a compact all-galaxy QC summary table.

    Assumes add_science_columns(tab) and add_qc_flags(tab) were already run.
    """
    from astropy.table import Table
    import numpy as np

    wanted = [
        # identity
        "VFID", "GALNAME", "OBJID", "TELESCOPE", "DATEOBS", "POINTING", "TAG",

        # core stage flags
        "STATUS", "PSF_OK", "MASK_OK", "PHOT_OK",
        "R_PROFILE_OK", "H_PROFILE_OK",
        "R_SM_OK", "H_SM_OK", "HAPY_MORPH_OK",
        "GAL_NC_OK", "GAL_CV_OK",

        # warnings / QC
        "WARN_MASK", "WARN_BRIGHT_STAR", "WARN_FILTER",
        "WARN_ELLIPSE", "WARN_WEAK_HA",
        "USE_R_STRUCTURE", "USE_HA_EXTENT", "USE_HA_MORPH", "USE_GALFIT",
        "QC_TIER",

        # core sizes
        "R25_ARCSEC", "R50_ARCSEC", "R75_ARCSEC",
        "H25_ARCSEC", "H50_ARCSEC", "H75_ARCSEC",
        "H_MAXDET_ARCSEC",
        "R_PETRO_R50_ARCSEC", "H_PETRO_R50_ARCSEC",

        # science ratios
        "H50_R50_RATIO", "H75_R75_RATIO",
        "H_MAXDET_R25_RATIO", "H_PETRO_R50_RATIO",

        # fluxes / magnitudes
        "R24_MAG", "R25_ISO_MAG", "GAL_MAG",
        "H_TOT_FLUX_CGS", "H_R24_FLUX_CGS",
        "R_C30", "H_C30_R24",

        # morphology
        "R_HAPY_GINI", "H_HAPY_GINI",
        "R_HAPY_M20", "H_HAPY_M20",
        "DELTA_GINI", "DELTA_M20",

        # robustness
        "R_PROFILE_NGOOD", "H_PROFILE_NGOOD",
        "R_PROFILE_MASKFRAC_MAX", "H_PROFILE_MASKFRAC_MAX",
        "H_HAPY_NPIX", "H_HAPY_SNP_DET",
        "ELL0_MASKFRAC",

        # galfit
        "GAL_RE_ARCSEC", "GAL_N", "GAL_BA", "GAL_PA", "GAL_CHISQ",
    ]

    keep = [c for c in wanted if c in tab.colnames]
    out = tab[keep].copy()
    out.write(outpath, format="fits", overwrite=True)

def write_outlier_tables(tab, outdir, n_outliers=25):
    """
    Write simple outlier tables for extent, flux/magnitude, and morphology.
    Assumes add_science_columns(tab) and add_qc_flags(tab) already ran.
    """
    import numpy as np

    def top_bottom_indices(x, n=25, positive_only=False):
        x = np.array(x, dtype=float)
        good = np.isfinite(x)
        if positive_only:
            good &= (x > 0)

        idx = np.where(good)[0]
        if len(idx) == 0:
            return np.array([], dtype=int)

        vals = x[idx]
        order = np.argsort(vals)

        lo = idx[order[:min(n, len(order))]]
        hi = idx[order[-min(n, len(order)):]]
        return np.unique(np.concatenate([lo, hi]))

    idcols = [c for c in ["VFID", "GALNAME", "OBJID", "TAG", "QC_TIER"] if c in tab.colnames]

    groups = {
        "outliers_extent_H50_R50.fits": ("H50_R50_RATIO", True),
        "outliers_extent_HMAX_R25.fits": ("H_MAXDET_R25_RATIO", True),
        "outliers_flux_Htot.fits": ("H_TOT_FLUX_CGS", True),
        "outliers_flux_HR24.fits": ("H_R24_FLUX_CGS", True),
        "outliers_mag_R24.fits": ("R24_MAG", False),
        "outliers_mag_GAL.fits": ("GAL_MAG", False),
        "outliers_concentration_R.fits": ("R_C30", False),
        "outliers_concentration_H.fits": ("H_C30_R24", False),
        "outliers_morph_DGini.fits": ("DELTA_GINI", False),
        "outliers_morph_DM20.fits": ("DELTA_M20", False),
    }

    extra = [c for c in [
        "R25_ARCSEC", "R50_ARCSEC", "H50_ARCSEC", "H_MAXDET_ARCSEC",
        "R24_MAG", "GAL_MAG", "H_TOT_FLUX_CGS", "H_R24_FLUX_CGS",
        "R_C30", "H_C30_R24",
        "R_HAPY_GINI", "H_HAPY_GINI", "R_HAPY_M20", "H_HAPY_M20",
        "DELTA_GINI", "DELTA_M20",
        "WARN_MASK", "WARN_BRIGHT_STAR", "WARN_FILTER", "WARN_ELLIPSE", "WARN_WEAK_HA",
        "R_PROFILE_NGOOD", "H_PROFILE_NGOOD", "H_HAPY_NPIX", "H_HAPY_SNP_DET"
    ] if c in tab.colnames]

    for fname, (col, positive_only) in groups.items():
        if col not in tab.colnames:
            continue
        idx = top_bottom_indices(tab[col], n=n_outliers, positive_only=positive_only)
        if len(idx) == 0:
            continue
        keep = idcols + [col] + [c for c in extra if c != col]
        sub = tab[idx][keep]
        sub.write(outdir / fname, format="fits", overwrite=True)
        
# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Survey-level QC for merged HAPY results.")
    parser.add_argument("table", help="Merged HAPY results table, e.g. merged_results.fits")
    parser.add_argument("--outdir", default="qc", help="Output directory")
    parser.add_argument(
        "--scheme",
        choices=["virgo", "agc"],
        required=True,
        help="Pipeline stage whose results should be merged."
    )
    parser.add_argument(
        "--max-ha-filter-correction",
        type=float,
        default=1.2,
        help="Threshold above which FILTER_CORRECTION is flagged as a warning for Halpha science",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    tab = Table.read(args.table)
    print(f"Read {len(tab)} rows from {args.table}")

    # -- add columns
    tab = add_science_columns(tab)
    tab = add_qc_flags(tab)

    # -- write output tables
    write_qc_summary_table(tab, outdir / "qc_summary_table.fits")
    write_outlier_tables(tab, outdir, n_outliers=25)

    masks = build_qc_masks(tab, max_ha_filter_correction=args.max_ha_filter_correction)

    # write text summary
    write_text_summary(tab, masks, outdir / "qc_summary.txt", args.scheme)

    # write subsets
    write_subsets(tab, masks, outdir)

    # plots
    plot_flag_completion(tab, outdir / "flag_completion.png")

    r_fwhm_col = first_populated_col(tab, ["R_FWHM_PSF", "R_FHWM_PSF"])
    h_fwhm_col = first_populated_col(tab, ["H_FWHM_PSF", "H_FHWM_PSF"])

    if r_fwhm_col is not None:
        plot_hist(tab, r_fwhm_col, outdir / "r_fwhm_hist.png", title="R-band FWHM")
    if h_fwhm_col is not None:
        plot_hist(tab, h_fwhm_col, outdir / "ha_fwhm_hist.png", title="Halpha FWHM")

    plot_hist(tab, "R24_MAG", outdir / "R24_MAG_hist.png", title="R24 magnitude", mask=masks["mask_phot_ok"])
    plot_hist(tab, "R25_ISO_MAG", outdir / "R25_ISO_MAG_hist.png", title="R25 isophotal magnitude", mask=masks["mask_phot_ok"])
    plot_hist(tab, "GAL_MAG", outdir / "GAL_MAG_hist.png", title="GALFIT magnitude", mask=masks["galfit_any_ok"])

    plot_hist(tab, "H_TOT_FLUX_CGS", outdir / "H_TOT_FLUX_CGS_hist.png",
              title="Halpha total flux", logx=True, mask=masks["mask_phot_ok"])
    plot_hist(tab, "H_R24_FLUX_CGS", outdir / "H_R24_FLUX_CGS_hist.png",
              title="Halpha R24 flux", logx=True, mask=masks["mask_phot_ok"])

    plot_hist(tab, "H_MAXDET_ARCSEC", outdir / "H_MAXDET_ARCSEC_hist.png",
              title="Halpha max detection radius", mask=masks["mask_phot_ok"])
    plot_hist(tab, "GAL_RE_ARCSEC", outdir / "GAL_RE_hist.png",
              title="GALFIT effective radius", mask=masks["galfit_any_ok"])

    plot_compare(tab, "R24_MAG", "GAL_MAG", outdir / "R24_vs_GAL_MAG.png",
                 title="R24 magnitude vs GALFIT magnitude",
                 mask=masks["galfit_any_ok"])

    plot_compare(tab, "H_R24_FLUX_CGS", "H_TOT_FLUX_CGS", outdir / "H_R24_vs_H_TOT.png",
                 title="Halpha R24 flux vs total flux",
                 logx=True, logy=True, mask=masks["mask_phot_ok"])

    plot_compare(tab, "R24_ARCSEC", "GAL_RE_ARCSEC", outdir / "R24_ARCSEC_vs_GAL_RE.png",
                 title="R24 radius vs GALFIT Re",
                 mask=masks["galfit_any_ok"])

    plot_failure_fraction_vs_bright_star_distance(tab, outdir / "FAILURES_VS_BRIGHT_STAR_DIST.png")
    plot_raincloud_by_telescope(tab, "R_FWHM_PSF", outdir / "raincloud_R_FWHM_by_telescope.png")
    plot_raincloud_by_telescope(tab, "H_FWHM_PSF", outdir / "raincloud_H_FWHM_by_telescope.png")
    plot_raincloud_by_telescope(tab, "FILTER_CORRECTION", outdir / "raincloud_FILTER_CORRECTION_by_telescope.png")
    plot_raincloud_by_telescope(tab, "H_TOT_FLUX_CGS", outdir / "raincloud_H_TOT_FLUX_by_telescope.png", logx=True)

    plot_filter_warning_vs_ha(tab, masks, outdir / "filter_correction_vs_ha_flux.png")

    plot_qc_dashboard(tab, outdir / "qc_dashboard.png")
    
    print(f"Wrote QC products to {outdir}")
    
    tab[tab["QC_TIER"] == "A"].write(outdir / "qc_tier_A.fits", format="fits", overwrite=True)
    tab[np.isin(tab["QC_TIER"], ["A", "B"])].write(outdir / "qc_tier_AB.fits", format="fits", overwrite=True)

if __name__ == "__main__":
    main()
