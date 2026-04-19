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
from hapy.utils.results_table import safe_bool_array, safe_float_array, safe_str_array
from hapy.utils.results_table import first_existing_col, first_populated_col#, build_row_qc_flags, 
from hapy.utils.results_table import ensure_dir, median_and_mad
from hapy.utils.results_table import prepare_analysis_table
from hapy.utils.plotting import QC_TIER_ORDER, QC_TIER_PALETTE

# ----------------------------------------------------------------------
# some science
# ----------------------------------------------------------------------
def get_ratio(a,b):
    return a/b




# ----------------------------------------------------------------------
# flag discovery and masks
# ----------------------------------------------------------------------

def find_ok_columns(tab: Table) -> list[str]:
    return [c for c in tab.colnames if c.endswith("_OK") ]

def find_flag_columns(tab: Table) -> list[str]:
    return [c for c in tab.colnames if c.endswith("_FLAG")]


# ----------------------------------------------------------------------
# text summaries
# ----------------------------------------------------------------------

def write_text_summary(tab: Table, outpath: Path, scheme) -> None:
    """
    Write a plain-text QC summary from a prepared analysis table.

    Assumes prepare_analysis_table(tab) has already been run, so derived
    QC/science columns like PROFILE_OK, SCIENCE_READY, etc. are present.
    """
    n = len(tab)

    with open(outpath, "w") as fh:
        fh.write("HAPY QC SUMMARY\n")
        fh.write("================\n\n")
        fh.write(f"Total rows: {n}\n")

        if scheme == "virgo":
            if "VFID" in tab.colnames:
                vals = safe_str_array(tab, "VFID", default="")
                vals = [v for v in vals if v != ""]
                fh.write(f"Unique galaxies: {len(set(vals))}\n")
        elif scheme == "agc":
            if "OBJID" in tab.colnames:
                vals = safe_str_array(tab, "OBJID", default="")
                vals = [v for v in vals if v != ""]
                fh.write(f"Unique galaxies: {len(set(vals))}\n")

        fh.write("\n")

        fh.write("Pipeline flags\n")
        fh.write("--------------\n")
        for f in find_ok_columns(tab):
            col = safe_bool_array(tab, f)
            ntrue = np.sum(col)
            nfalse = np.sum(~col)
            pct = 100.0 * ntrue / n if n > 0 else np.nan
            fh.write(f"{f:20s}: {ntrue:4d} OK  | {nfalse:4d} FAIL  ({pct:5.1f}%)\n")

        fh.write("\nWarning flags\n")
        fh.write("-------------\n")
        for f in find_flag_columns(tab):
            col = safe_bool_array(tab, f)
            nflag = np.sum(col)
            nclean = np.sum(~col)
            pct = 100.0 * nflag / n if n > 0 else np.nan
            fh.write(f"{f:20s}: {nclean:4d} Clean  | {nflag:4d} Flagged  ({pct:5.1f}%)\n")

        fh.write("\nQC subsets\n")
        fh.write("----------\n")
        subset_cols = [
            "MASK_PHOT_OK",
            "PROFILE_OK",
            "GALFIT_ANY_OK",
            "GALFIT_BOTH_OK",
            "R_STRUCTURE_GOOD",
            "HA_EXTENT_GOOD",
            "HA_MORPH_GOOD",
            "SCIENCE_READY",
            "TECHNICAL_PROBLEM",
            "SCIENCE_PROBLEM",
        ]

        # optional legacy/secondary summary
        if "R_SM_OK" in tab.colnames and "H_SM_OK" in tab.colnames and "FILTER_WARNING" in tab.colnames:
            statmorph_ok = (
                safe_bool_array(tab, "R_SM_OK") &
                safe_bool_array(tab, "H_SM_OK") &
                (~safe_bool_array(tab, "FILTER_WARNING"))
            )
            count = np.sum(statmorph_ok)
            pct = 100.0 * count / n if n > 0 else np.nan
            fh.write(f"{'STATMORPH_OK':20s}: {count:4d} ({pct:5.1f}%)\n")

        for name in subset_cols:
            if name not in tab.colnames:
                continue
            m = safe_bool_array(tab, name)
            count = np.sum(m)
            pct = 100.0 * count / n if n > 0 else np.nan
            fh.write(f"{name:20s}: {count:4d} ({pct:5.1f}%)\n")

        # convenience combined problem count
        if "TECHNICAL_PROBLEM" in tab.colnames and "SCIENCE_PROBLEM" in tab.colnames:
            problem = safe_bool_array(tab, "TECHNICAL_PROBLEM") | safe_bool_array(tab, "SCIENCE_PROBLEM")
            count = np.sum(problem)
            pct = 100.0 * count / n if n > 0 else np.nan
            fh.write(f"{'PROBLEM':20s}: {count:4d} ({pct:5.1f}%)\n")

        if "QC_TIER" in tab.colnames:
            fh.write("\nQC_TIER counts\n")
            fh.write("--------------\n")
            vals, counts = np.unique(safe_str_array(tab, "QC_TIER", default=""), return_counts=True)
            for v, c in zip(vals, counts):
                if v == "":
                    continue
                fh.write(f"{v:20s}: {c}\n")

        if "STATUS" in tab.colnames:
            fh.write("\nSTATUS counts\n")
            fh.write("-------------\n")
            vals, counts = np.unique(safe_str_array(tab, "STATUS", default=""), return_counts=True)
            for v, c in zip(vals, counts):
                if v == "":
                    continue
                fh.write(f"{v:20s}: {c}\n")

        if "STAGE" in tab.colnames:
            fh.write("\nSTAGE counts\n")
            fh.write("------------\n")
            vals, counts = np.unique(safe_str_array(tab, "STAGE", default=""), return_counts=True)
            for v, c in zip(vals, counts):
                if v == "":
                    continue
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
        "BRIGHT_STAR_MASKRAD_ARCSEC",
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
    rad = np.array(tab["BRIGHT_STAR_MASKRAD_ARCSEC"], dtype=float)

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

def write_subsets(tab: Table, outdir: Path) -> None:
    """
    Write standard QC subsets from a prepared analysis table.
    Assumes prepare_analysis_table(tab) has already been run.
    """

    subset_map = {
        "mask_phot_ok": "MASK_PHOT_OK",
        "profile_ok": "PROFILE_OK",
        "galfit_any_ok": "GALFIT_ANY_OK",
        "galfit_both_ok": "GALFIT_BOTH_OK",
        "r_structure_good": "R_STRUCTURE_GOOD",
        "ha_extent_good": "HA_EXTENT_GOOD",
        "ha_morph_good": "HA_MORPH_GOOD",
        "science_ready": "SCIENCE_READY",
    }

    for outname, colname in subset_map.items():
        if colname not in tab.colnames:
            continue

        m = safe_bool_array(tab, colname)
        sub = tab[m]
        sub.write(outdir / "tables" / "subsets" / f"{outname}.fits", format="fits", overwrite=True)

    # --------------------------------------------------
    # Combined "problem" subset
    # --------------------------------------------------
    if "TECHNICAL_PROBLEM" in tab.colnames and "SCIENCE_PROBLEM" in tab.colnames:
        problem = (
            safe_bool_array(tab, "TECHNICAL_PROBLEM") |
            safe_bool_array(tab, "SCIENCE_PROBLEM")
        )
        tab[problem].write(outdir / "tables" / "subsets" / "problem.fits", format="fits", overwrite=True)
        tab[problem].write(outdir / "tables" / "subsets" / "problem.csv", format="ascii.csv", overwrite=True)

    # --------------------------------------------------
    # Optional: legacy statmorph subset (keep for now)
    # --------------------------------------------------
    if all(c in tab.colnames for c in ["R_SM_OK", "H_SM_OK", "FILTER_WARNING"]):
        statmorph_ok = (
            safe_bool_array(tab, "R_SM_OK") &
            safe_bool_array(tab, "H_SM_OK") &
            (~safe_bool_array(tab, "FILTER_WARNING"))
        )
        tab[statmorph_ok].write(outdir / "tables" / "subsets" / "statmorph_ok.fits", format="fits", overwrite=True)

        
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
    #out.write(outpath.replace(".fits",".csv"), format="ascii.csv", overwrite=True)


def write_outlier_tables(tab: Table, outdir: Path, n_outliers: int = 25) -> None:
    """
    Write outlier tables for extent, flux/magnitude, morphology, and QC diagnostics.

    Assumes prepare_analysis_table(tab) has already been run.
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
        "outliers_extent_HPETRO_R50.fits": ("H_PETRO_R50_RATIO", True),

        "outliers_flux_Htot.fits": ("H_TOT_FLUX_CGS", True),
        "outliers_flux_HR24.fits": ("H_R24_FLUX_CGS", True),

        "outliers_mag_R24.fits": ("R24_MAG", False),
        "outliers_mag_GAL.fits": ("GAL_MAG", False),

        "outliers_concentration_R.fits": ("R_C30", False),
        "outliers_concentration_H.fits": ("H_C30_R24", False),

        "outliers_morph_DGini.fits": ("DELTA_GINI", False),
        "outliers_morph_DM20.fits": ("DELTA_M20", False),
        "outliers_morph_DAsym.fits": ("DELTA_ASYM", False),

        "outliers_halpha_fillfrac.fits": ("H_HAPY_FILLFRAC", False),
        "outliers_halpha_npix.fits": ("H_HAPY_NPIX", True),
        "outliers_halpha_snp_det.fits": ("H_HAPY_SNP_DET", False),
    }

    extra = [c for c in [
        # sizes / extents
        "R25_ARCSEC", "R50_ARCSEC", "R75_ARCSEC",
        "H50_ARCSEC", "H75_ARCSEC", "H_MAXDET_ARCSEC",
        "R_PETRO_R50_ARCSEC", "H_PETRO_R50_ARCSEC",
        "H50_R50_RATIO", "H75_R75_RATIO", "H_MAXDET_R25_RATIO", "H_PETRO_R50_RATIO",

        # flux / mags
        "R24_MAG", "GAL_MAG", "H_TOT_FLUX_CGS", "H_R24_FLUX_CGS",

        # concentration / morphology
        "R_C30", "H_C30_R24",
        "R_HAPY_GINI", "H_HAPY_GINI",
        "R_HAPY_M20", "H_HAPY_M20",
        "R_HAPY_ASYM", "H_HAPY_ASYM",
        "DELTA_GINI", "DELTA_M20", "DELTA_ASYM",
        "H_HAPY_FILLFRAC",

        # QC / warnings
        "MASK_PHOT_OK", "PROFILE_OK",
        "R_STRUCTURE_GOOD", "HA_EXTENT_GOOD", "HA_MORPH_GOOD",
        "SCIENCE_READY", "TECHNICAL_PROBLEM", "SCIENCE_PROBLEM",
        "WARN_MASK", "WARN_WEAK_HA",
        "BRIGHT_STAR_FLAG", "FILTER_WARNING", "ELL_MISMATCH",

        # diagnostics
        "R_PROFILE_NGOOD", "H_PROFILE_NGOOD",
        "R_PROFILE_MASKFRAC_MAX", "H_PROFILE_MASKFRAC_MAX",
        "H_HAPY_NPIX", "H_HAPY_SNP_DET",
    ] if c in tab.colnames]

    for fname, (col, positive_only) in groups.items():
        if col not in tab.colnames:
            continue

        idx = top_bottom_indices(tab[col], n=n_outliers, positive_only=positive_only)
        if len(idx) == 0:
            continue

        keep = idcols + [col] + [c for c in extra if c != col]
        sub = tab[idx][keep]
        sub.write(outdir / "tables" / "outliers" /  fname, format="fits", overwrite=True)

        
def write_review_table_old(tab: Table, outdir: Path, scheme: str) -> None:
    """
    Build a compact CSV table for Google Sheets QC review.
    Assumes prepare_analysis_table(tab) has already run.
    """

    import numpy as np

    # -------------------------
    # Build WEBPAGE column
    # -------------------------

    if "TAG" in tab.colnames:
        base = f"http://199.223.247.130/fits/{scheme}/cutouts"
        tags = safe_str_array(tab, "TAG", default="")
        webpage_list = [
            f"{base}/{tag}/{tag}.html" if tag not in ("", "nan", "None") else ""
            for tag in tags
        ]
    else:
        webpage_list = [""] * len(tab)

    maxlen = max([len(s) for s in webpage_list], default=1)

    

    # -------------------------
    # Columns to include
    # -------------------------
    cols = [
        # identity
        "VFID", "GALNAME", "OBJID", "TAG",
        # QC
        "QC_TIER", "SCIENCE_READY",
        "R_STRUCTURE_GOOD", "HA_EXTENT_GOOD", "HA_MORPH_GOOD",
        # warnings
        "WARN_MASK", "WARN_WEAK_HA",
        "BRIGHT_STAR_FLAG", "FILTER_WARNING", "ELL_MISMATCH",
        "SEVERE_CEN_ANY", "WARN_CEN_ANY",
        # metrics
        "H50_R50_RATIO", "H_MAXDET_R25_RATIO",
        "DELTA_GINI", "DELTA_M20", "DELTA_ASYM",
        "H_HAPY_FILLFRAC", "H_HAPY_NPIX",
        # review
        "REVIEW_PRIORITY",
    ]

    cols = [c for c in cols if c in tab.colnames]

    review = tab[cols].copy()
    
    #review["WEBPAGE"] = webpage

    # -------------------------
    # Add human-review columns
    # -------------------------
    n = len(review)

    review["REVIEWED"] = np.full(n, "", dtype=object)
    #review["REVIEW_PRIORITY"] = np.full(n, "", dtype=object)
    review["VIS_CLASS"] = np.full(n, "", dtype=object)
    review["CATALOG_USE"] = np.full(n, "", dtype=object)
    review["VIS_NOTE"] = np.full(n, "", dtype=object)

    # -- mask columns
    review["MASK_FIX_NEEDED"] = np.full(n, "", dtype=object)
    review["MASK_FIXED"] = np.full(n, "", dtype=object)
    review["MASK_ISSUE"] = np.full(n, "", dtype=object)
    review["MASK_NOTE"] = np.full(n, "", dtype=object)    
    #review["VIS_NOTE"] = np.full(n, "", dtype=object)

    # -------------------------
    # Set default priority
    # -------------------------

 

    # add webpage last so url doesn't block other columns
    review["WEBPAGE"] = np.array(webpage_list, dtype=f"<U{maxlen}")    
    # -------------------------
    # Write outputs
    # -------------------------
    outpath = outdir / "tables" / "review"
    ensure_dir(outpath)
    Nhigh = np.sum(review["REVIEW_PRIORITY"] == "high")

    print(f"Number of high priority in {outpath} = {Nhigh}")
    #review.write(outpath / "review_master.fits", format="fits", overwrite=True)
    print("writing ",outpath / "review_sample.csv")
    review.write(outpath / "review_sample.csv", format="ascii.csv", overwrite=True)


def write_review_table(tab: Table, outdir: Path, scheme: str) -> None:
    """
    Build a compact CSV table for Google Sheets QC review.
    Assumes prepare_analysis_table(tab) has already run.
    """

    import numpy as np

    # -------------------------
    # Build WEBPAGE column
    # -------------------------
    if "TAG" in tab.colnames:
        base = f"http://199.223.247.130/fits/{scheme}/cutouts"
        tags = safe_str_array(tab, "TAG", default="")
        webpage_list = [
            f"{base}/{tag}/{tag}.html" if tag not in ("", "nan", "None") else ""
            for tag in tags
        ]
    else:
        webpage_list = [""] * len(tab)

    maxlen = max([len(s) for s in webpage_list], default=1)

    # -------------------------
    # Columns to include
    # Match review-priority drivers + useful diagnostics
    # -------------------------
    cols = [
        # identity
        "VFID", "GALNAME", "OBJID", "TAG",

        # summary / overall
        "QC_TIER", "SCIENCE_READY", "REVIEW_PRIORITY",

        # core success/failure
        "PHOT_OK", "HAPY_MORPH_OK",
        "R_PROFILE_OK", "H_PROFILE_OK",
        "GAL_NC_OK", "GAL_CV_OK",

        # trigger columns used in review priority
        "BRIGHT_STAR_FLAG",
        "WARN_MASK",
        "SEVERE_CEN_ANY",
        "ELL_MISMATCH",
        "WARN_WEAK_HA",
        "FILTER_WARNING",
        "WARN_CEN_ANY",
        "WARN_R_PROFILE_PEAK",
        "WARN_CUTOUT_MISSING",
        "WARN_CUTOUT_MISSING_SHAPE",

        # center diagnostics
        "DOFF_IN_PHOT_ARCSEC",
        "DOFF_IN_GAL_ARCSEC",
        "DOFF_IN_GALC_ARCSEC",
        "DOFF_PHOT_GAL_ARCSEC",
        "DOFF_PHOT_GALC_ARCSEC",
        "DOFF_GAL_GALC_ARCSEC",

        # cutout coverage diagnostics
        "CUTOUT_ELL0_MISSING_FRAC_R",
        "CUTOUT_ELL0_MISSING_FRAC_H",
        "CUTOUT_ELL0_MISSING_FRAC_MAX",
        "CUTOUT_ELL0_NPIX_TOTAL_R",
        "CUTOUT_ELL0_NPIX_TOTAL_H",
        "CUTOUT_ELL0_NPIX_ONIMAGE_R",
        "CUTOUT_ELL0_NPIX_ONIMAGE_H",
        "CUTOUT_ELL0_NPIX_GOOD_R",
        "CUTOUT_ELL0_NPIX_GOOD_H",

        # profile / shape diagnostics
        "R_PROFILE_NONCENTRAL_PEAK",
        "R_PROFILE_PEAK_BIN",
        "R_PROFILE_PEAK_SMA",

        # science/QC context
        "R_STRUCTURE_GOOD", "HA_EXTENT_GOOD", "HA_MORPH_GOOD",
        "H50_R50_RATIO", "H_MAXDET_R25_RATIO",
        "DELTA_GINI", "DELTA_M20", "DELTA_ASYM",
        "H_HAPY_FILLFRAC", "H_HAPY_NPIX",
    ]

    cols = [c for c in cols if c in tab.colnames]
    review = tab[cols].copy()

    # -------------------------
    # Add human-review columns
    # -------------------------
    n = len(review)

    review["REVIEWED"] = np.full(n, "", dtype=object)
    review["VIS_CLASS"] = np.full(n, "", dtype=object)
    review["CATALOG_USE"] = np.full(n, "", dtype=object)
    review["VIS_NOTE"] = np.full(n, "", dtype=object)

    # mask / coverage follow-up
    review["MASK_FIX_NEEDED"] = np.full(n, "", dtype=object)
    review["MASK_FIXED"] = np.full(n, "", dtype=object)
    review["MASK_ISSUE"] = np.full(n, "", dtype=object)
    review["MASK_NOTE"] = np.full(n, "", dtype=object)

    # -------------------------
    # Add webpage last so URL doesn't block other columns
    # -------------------------
    review["WEBPAGE"] = np.array(webpage_list, dtype=f"<U{maxlen}")

    # -------------------------
    # Write outputs
    # -------------------------
    outpath = outdir / "tables" / "review"
    ensure_dir(outpath)

    if "REVIEW_PRIORITY" in review.colnames:
        Nhigh = np.sum(review["REVIEW_PRIORITY"] == "high")
        print(f"Number of high priority in {outpath} = {Nhigh}")

    print("writing ", outpath / "review_sample.csv")
    review.write(outpath / "review_sample.csv", format="ascii.csv", overwrite=True)
    
def print_review_priority_drivers(tab):
    """
    Print counts for the individual flags that drive review priority.

    Assumes `tab` already contains all derived columns from
    `prepare_analysis_table()` / `get_review_priority()`.

    Parameters
    ----------
    tab : astropy.table.Table
        Analysis table with QC / warning columns.
    """
    import numpy as np

    n = len(tab)
    if n == 0:
        print("Table is empty.")
        return

    def _safe_bool(name):
        if name not in tab.colnames:
            return np.zeros(n, dtype=bool)
        arr = np.asarray(tab[name])
        if arr.dtype == bool:
            return arr
        out = np.zeros(n, dtype=bool)
        good = arr == arr
        out[good] = arr[good].astype(bool)
        return out

    high_terms = {
        "NOT_PHOT_OK": ~_safe_bool("PHOT_OK"),
        "NOT_HAPY_MORPH_OK": ~_safe_bool("HAPY_MORPH_OK"),
        "BRIGHT_STAR_FLAG": _safe_bool("BRIGHT_STAR_FLAG"),
        "WARN_MASK": _safe_bool("WARN_MASK"),
        "SEVERE_CEN_ANY": _safe_bool("SEVERE_CEN_ANY"),
        "WARN_CUTOUT_MISSING_SHAPE": _safe_bool("WARN_CUTOUT_MISSING_SHAPE"),
    }

    medium_terms = {
        "ELL_MISMATCH": _safe_bool("ELL_MISMATCH"),
        "WARN_WEAK_HA": _safe_bool("WARN_WEAK_HA"),
        "FILTER_WARNING": _safe_bool("FILTER_WARNING"),
        "WARN_CEN_ANY": _safe_bool("WARN_CEN_ANY"),
        "WARN_R_PROFILE_PEAK": _safe_bool("WARN_R_PROFILE_PEAK"),
        "WARN_CUTOUT_MISSING": _safe_bool("WARN_CUTOUT_MISSING"),
    }

    high = np.zeros(n, dtype=bool)
    for flag in high_terms.values():
        high |= flag

    medium = np.zeros(n, dtype=bool)
    for flag in medium_terms.values():
        medium |= flag

    priority = np.full(n, "low", dtype="U16")
    priority[high] = "high"
    priority[medium & (~high)] = "medium"

    unique, counts = np.unique(priority, return_counts=True)
    summary = dict(zip(unique, counts))

    print("\nREVIEW PRIORITY SUMMARY")
    print(summary)

    print("\nHIGH PRIORITY DRIVERS")
    for name, flag in high_terms.items():
        print(f"{name:28s}: total={np.sum(flag):4d}  in_high={np.sum(flag & high):4d}")

    print("\nMEDIUM PRIORITY DRIVERS")
    for name, flag in medium_terms.items():
        print(f"{name:28s}: total={np.sum(flag):4d}  in_medium={np.sum(flag & (medium & ~high)):4d}")

    print("\nOVERLAP AMONG HIGH DRIVERS")
    high_names = list(high_terms.keys())
    for i in range(len(high_names)):
        for j in range(i + 1, len(high_names)):
            n_ij = np.sum(high_terms[high_names[i]] & high_terms[high_names[j]])
            if n_ij > 0:
                print(f"{high_names[i]:28s} & {high_names[j]:28s}: {n_ij:4d}")

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




    tables_dir = outdir / "tables"
    plots_dir = outdir / "plots"

    ensure_dir(tables_dir)
    ensure_dir(plots_dir)

    # optional subdirs
    ensure_dir(tables_dir / "subsets")
    ensure_dir(tables_dir / "outliers")
    ensure_dir(tables_dir / "summary")


    tab = Table.read(args.table)
    print(f"Read {len(tab)} rows from {args.table}")

    # -- add columns for qc and science
    tab = prepare_analysis_table(tab)
    
    # -- write output tables
    write_qc_summary_table(tab, outdir / "tables" / "summary" / "qc_summary_table.fits")
    write_outlier_tables(tab, outdir, n_outliers=25)


    # write text summary
    write_text_summary(tab, outdir / "qc_summary.txt", args.scheme)

    # write subsets
    write_subsets(tab, outdir)

    # plots
    plot_flag_completion(tab, outdir / "flag_completion.png")

    r_fwhm_col = first_populated_col(tab, ["R_FWHM_PSF", "R_FHWM_PSF"])
    h_fwhm_col = first_populated_col(tab, ["H_FWHM_PSF", "H_FHWM_PSF"])

    #if r_fwhm_col is not None:
    #    plot_hist(tab, r_fwhm_col, outdir / "r_fwhm_hist.png", title="R-band FWHM")
    #if h_fwhm_col is not None:
    #    plot_hist(tab, h_fwhm_col, outdir / "ha_fwhm_hist.png", title="Halpha FWHM")

    #plot_hist(tab, "R24_MAG", outdir / "R24_MAG_hist.png", title="R24 magnitude", mask=masks["mask_phot_ok"])
    #plot_hist(tab, "R25_ISO_MAG", outdir / "R25_ISO_MAG_hist.png", title="R25 isophotal magnitude", mask=masks["mask_phot_ok"])
    #plot_hist(tab, "GAL_MAG", outdir / "GAL_MAG_hist.png", title="GALFIT magnitude", mask=masks["galfit_any_ok"])

    #plot_hist(tab, "H_TOT_FLUX_CGS", outdir / "H_TOT_FLUX_CGS_hist.png",
    #          title="Halpha total flux", logx=True, mask=masks["mask_phot_ok"])
    #plot_hist(tab, "H_R24_FLUX_CGS", outdir / "H_R24_FLUX_CGS_hist.png",
    #          title="Halpha R24 flux", logx=True, mask=masks["mask_phot_ok"])

    #plot_hist(tab, "H_MAXDET_ARCSEC", outdir / "H_MAXDET_ARCSEC_hist.png",
    #          title="Halpha max detection radius", mask=masks["mask_phot_ok"])
    #plot_hist(tab, "GAL_RE_ARCSEC", outdir / "GAL_RE_hist.png",
    #          title="GALFIT effective radius", mask=masks["galfit_any_ok"])

    #plot_compare(tab, "R24_MAG", "GAL_MAG", outdir / "R24_vs_GAL_MAG.png",
    #             title="R24 magnitude vs GALFIT magnitude",
    #             mask=masks["galfit_any_ok"])

    #plot_compare(tab, "H_R24_FLUX_CGS", "H_TOT_FLUX_CGS", outdir / "H_R24_vs_H_TOT.png",
    #             title="Halpha R24 flux vs total flux",
    #             logx=True, logy=True, mask=masks["mask_phot_ok"])

    #plot_compare(tab, "R24_ARCSEC", "GAL_RE_ARCSEC", outdir / "R24_ARCSEC_vs_GAL_RE.png",
    #             title="R24 radius vs GALFIT Re",
    #             mask=masks["galfit_any_ok"])

    #plot_failure_fraction_vs_bright_star_distance(tab, plots_dir / "FAILURES_VS_BRIGHT_STAR_DIST.png")
    plot_raincloud_by_telescope(tab, "R_FWHM_PSF", plots_dir / "raincloud_R_FWHM_by_telescope.png")
    plot_raincloud_by_telescope(tab, "H_FWHM_PSF", plots_dir / "raincloud_H_FWHM_by_telescope.png")
    #plot_raincloud_by_telescope(tab, "FILTER_CORRECTION", outdir / "raincloud_FILTER_CORRECTION_by_telescope.png")
    #plot_raincloud_by_telescope(tab, "H_TOT_FLUX_CGS", outdir / "raincloud_H_TOT_FLUX_by_telescope.png", logx=True)

    #plot_filter_warning_vs_ha(tab, masks, outdir / "filter_correction_vs_ha_flux.png")

    plot_qc_dashboard(tab, plots_dir / "qc_dashboard.png")
    
    print(f"Wrote QC products to {outdir}")
    write_review_table(tab, outdir, args.scheme)
    
    tab[tab["QC_TIER"] == "A"].write(outdir / "tables" / "subsets" / "qc_tier_A.fits", format="fits", overwrite=True)
    tab[np.isin(tab["QC_TIER"], ["A", "B"])].write(outdir / "tables" / "subsets" /"qc_tier_AB.fits", format="fits", overwrite=True)

    print()
    print_review_priority_drivers(tab)

if __name__ == "__main__":
    main()
