#!/usr/bin/env python

"""
science_firstlook.py

First-look science plots for merged HAPY results.

This script is intentionally lightweight and exploratory. It is meant to:
- read the merged HAPY results table
- add a few derived science/QC columns if needed
- apply a simple QC-based sample selection
- generate a compact set of first-science plots

Recommended use
---------------
python science_firstlook.py merged_results.fits --outdir science_firstlook --sample AB

Sample options
--------------
A     : QC_TIER == "A"
AB    : QC_TIER in ["A", "B"]   [default]
ABC   : QC_TIER in ["A", "B", "C"]
ALL   : all rows

Notes
-----
- This script assumes the merged table already contains the raw HAPY outputs.
- If derived columns like H50_R50_RATIO or QC_TIER are missing, it computes them.
- The goal is fast exploratory science plots, not final publication figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from hapy.utils.plotting import QC_TIER_ORDER, QC_TIER_PALETTE

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def safe_float_array(tab: Table, colname: str, default=np.nan) -> np.ndarray:
    if colname not in tab.colnames:
        return np.full(len(tab), default, dtype=float)

    col = tab[colname]
    try:
        if hasattr(col, "filled"):
            col = col.filled(default)
    except Exception:
        pass

    out = np.full(len(tab), default, dtype=float)
    for i, v in enumerate(col):
        try:
            out[i] = float(v)
        except Exception:
            out[i] = default
    return out


def safe_bool_array(tab: Table, colname: str, default: bool = False) -> np.ndarray:
    if colname not in tab.colnames:
        return np.full(len(tab), default, dtype=bool)

    col = tab[colname]
    try:
        if hasattr(col, "filled"):
            col = col.filled(default)
    except Exception:
        pass

    out = np.zeros(len(tab), dtype=bool)
    for i, v in enumerate(col):
        if v is None:
            out[i] = default
        elif isinstance(v, (bool, np.bool_)):
            out[i] = bool(v)
        else:
            s = str(v).strip().lower()
            if s in ("true", "t", "1", "yes", "y"):
                out[i] = True
            elif s in ("false", "f", "0", "no", "n", "", "none", "nan"):
                out[i] = False
            else:
                out[i] = default
    return out


def safe_str_array(tab: Table, colname: str, default: str = "") -> np.ndarray:
    if colname not in tab.colnames:
        return np.full(len(tab), default, dtype="U16")

    out = np.full(len(tab), default, dtype="U16")
    for i, v in enumerate(tab[colname]):
        if v is None:
            out[i] = default
        else:
            out[i] = str(v)
    return out


def add_science_columns(tab: Table) -> Table:
    """
    Add core derived science columns if missing.
    """

    def safe_ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
        out = np.full(len(num), np.nan, dtype=float)
        good = np.isfinite(num) & np.isfinite(den) & (den > 0)
        out[good] = num[good] / den[good]
        return out

    # raw columns
    r50 = safe_float_array(tab, "R50_ARCSEC")
    h50 = safe_float_array(tab, "H50_ARCSEC")

    r75 = safe_float_array(tab, "R75_ARCSEC")
    h75 = safe_float_array(tab, "H75_ARCSEC")

    r25 = safe_float_array(tab, "R25_ARCSEC")
    hmax = safe_float_array(tab, "H_MAXDET_ARCSEC")

    r_p50 = safe_float_array(tab, "R_PETRO_R50_ARCSEC")
    h_p50 = safe_float_array(tab, "H_PETRO_R50_ARCSEC")

    r_gini = safe_float_array(tab, "R_HAPY_GINI")
    h_gini = safe_float_array(tab, "H_HAPY_GINI")

    r_m20 = safe_float_array(tab, "R_HAPY_M20")
    h_m20 = safe_float_array(tab, "H_HAPY_M20")

    # derived
    if "H50_R50_RATIO" not in tab.colnames:
        tab["H50_R50_RATIO"] = safe_ratio(h50, r50)
    if "H75_R75_RATIO" not in tab.colnames:
        tab["H75_R75_RATIO"] = safe_ratio(h75, r75)
    if "H_MAXDET_R25_RATIO" not in tab.colnames:
        tab["H_MAXDET_R25_RATIO"] = safe_ratio(hmax, r25)
    if "H_PETRO_R50_RATIO" not in tab.colnames:
        tab["H_PETRO_R50_RATIO"] = safe_ratio(h_p50, r_p50)

    if "DELTA_GINI" not in tab.colnames:
        tab["DELTA_GINI"] = h_gini - r_gini
    if "DELTA_M20" not in tab.colnames:
        tab["DELTA_M20"] = h_m20 - r_m20

    return tab


def add_qc_flags(tab: Table) -> Table:
    """
    Add a minimal QC tier if missing.
    If QC_TIER already exists, leave it alone.
    """

    if "QC_TIER" in tab.colnames:
        return tab

    n = len(tab)

    phot_ok = safe_bool_array(tab, "PHOT_OK")
    rprof_ok = safe_bool_array(tab, "R_PROFILE_OK")
    hprof_ok = safe_bool_array(tab, "H_PROFILE_OK")
    morph_ok = safe_bool_array(tab, "HAPY_MORPH_OK")
    h_sm_ok = safe_bool_array(tab, "H_SM_OK")
    gal_nc_ok = safe_bool_array(tab, "GAL_NC_OK")
    gal_cv_ok = safe_bool_array(tab, "GAL_CV_OK")

    bright_star = safe_bool_array(tab, "BRIGHT_STAR_FLAG")
    mask_warn = safe_bool_array(tab, "ELL0_MASK_WARN")
    ell_warn = safe_bool_array(tab, "ELL_MISMATCH")

    filt = safe_float_array(tab, "FILTER_CORRECTION")
    warn_filter = np.isfinite(filt) & (filt > 1.2)

    r50 = safe_float_array(tab, "R50_ARCSEC")
    h50 = safe_float_array(tab, "H50_ARCSEC")
    hmax = safe_float_array(tab, "H_MAXDET_ARCSEC")
    h_npix = safe_float_array(tab, "H_HAPY_NPIX")
    h_ngood = safe_float_array(tab, "H_PROFILE_NGOOD")
    h_snr = safe_float_array(tab, "H_HAPY_SNP_DET")
    r_ngood = safe_float_array(tab, "R_PROFILE_NGOOD")

    use_r = np.zeros(n, dtype=int)
    use_ha = np.zeros(n, dtype=int)
    use_hm = np.zeros(n, dtype=int)
    use_gf = np.zeros(n, dtype=int)

    good_r = (
        phot_ok & rprof_ok &
        np.isfinite(r50) & (r50 > 0) &
        np.isfinite(r_ngood) & (r_ngood >= 20)
    )
    use_r[good_r] = 2
    use_r[good_r & mask_warn] = 1

    good_ha = (
        phot_ok & hprof_ok &
        np.isfinite(h50) & (h50 > 0) &
        np.isfinite(hmax) & (hmax > 0) &
        np.isfinite(h_npix) & (h_npix >= 50) &
        np.isfinite(h_ngood) & (h_ngood >= 8) &
        (~warn_filter)
    )
    use_ha[good_ha] = 2

    weak_ha = (
        (np.isfinite(h_npix) & (h_npix < 50)) |
        (np.isfinite(h_ngood) & (h_ngood < 8)) |
        (np.isfinite(h_snr) & (h_snr < 3))
    )
    use_ha[good_ha & weak_ha] = 1

    good_hm = (
        morph_ok & h_sm_ok &
        np.isfinite(h_npix) & (h_npix >= 100) &
        (~warn_filter)
    )
    use_hm[good_hm] = 2
    use_hm[good_hm & weak_ha] = 1

    use_gf[gal_nc_ok | gal_cv_ok] = 2

    tab["USE_R_STRUCTURE"] = use_r
    tab["USE_HA_EXTENT"] = use_ha
    tab["USE_HA_MORPH"] = use_hm
    tab["USE_GALFIT"] = use_gf

    tier = np.full(n, "F", dtype="U1")
    for i in range(n):
        if not phot_ok[i]:
            tier[i] = "F"
        elif use_r[i] == 0 or use_ha[i] == 0:
            tier[i] = "D"
        elif mask_warn[i] or bright_star[i] or warn_filter[i] or ell_warn[i]:
            tier[i] = "C"
        elif use_hm[i] < 2:
            tier[i] = "B"
        else:
            tier[i] = "A"

    tab["QC_TIER"] = tier
    return tab


def select_sample(tab: Table, sample_name: str) -> np.ndarray:
    """
    Return boolean mask for sample selection.
    """
    tier = safe_str_array(tab, "QC_TIER", default="F")

    sample_name = sample_name.upper()
    if sample_name == "A":
        return tier == "A"
    if sample_name == "AB":
        return np.isin(tier, ["A", "B"])
    if sample_name == "ABC":
        return np.isin(tier, ["A", "B", "C"])
    if sample_name == "ALL":
        return np.ones(len(tab), dtype=bool)

    raise ValueError(f"Unknown sample selection: {sample_name}")


def scatter_by_tier(ax, x, y, tier, xlabel, ylabel, title,
                    one_to_one=False, logx=False, logy=False):
    good = np.isfinite(x) & np.isfinite(y)
    if logx:
        good &= (x > 0)
    if logy:
        good &= (y > 0)



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
        lo = np.nanmin([np.nanmin(x[good]), np.nanmin(y[good])])
        hi = np.nanmax([np.nanmax(x[good]), np.nanmax(y[good])])
        if np.isfinite(lo) and np.isfinite(hi):
            ax.plot([lo, hi], [lo, hi], "k--", lw=1)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def plot_firstlook_dashboard(tab: Table, outpath: Path, sample_label: str) -> None:
    tier = safe_str_array(tab, "QC_TIER", default="?")
    r50 = safe_float_array(tab, "R50_ARCSEC")
    h50 = safe_float_array(tab, "H50_ARCSEC")
    h_npix = safe_float_array(tab, "H_HAPY_NPIX")
    h_snr = safe_float_array(tab, "H_HAPY_SNP_DET")
    h_gini = safe_float_array(tab, "H_HAPY_GINI")
    r_gini = safe_float_array(tab, "R_HAPY_GINI")
    h_m20 = safe_float_array(tab, "H_HAPY_M20")
    r_m20 = safe_float_array(tab, "R_HAPY_M20")

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.ravel()

    # 1. tier counts
    ax = axes[0]


    tiers = [t for t in QC_TIER_ORDER if np.any(tier == t)]
    counts = [np.sum(tier == t) for t in tiers]
    colors = [QC_TIER_PALETTE[t] for t in tiers]
    ax.bar(range(len(tiers)), counts, color=colors)
    ax.set_xticks(range(len(tiers)))
    ax.set_xticklabels(tiers)
    ax.set_ylabel("N")
    ax.set_title(f"QC tiers ({sample_label})")

    # 2. H50 vs R50
    scatter_by_tier(
        axes[1], r50, h50, tier,
        xlabel="R50_ARCSEC", ylabel="H50_ARCSEC",
        title="Hα vs stellar half-light radius",
        one_to_one=True
    )

    # 3. Halpha robustness
    scatter_by_tier(
        axes[2], h_npix, h_snr, tier,
        xlabel="H_HAPY_NPIX", ylabel="H_HAPY_SNP_DET",
        title="Hα detection robustness",
        logx=True
    )

    # 4. Halpha vs stellar Gini
    scatter_by_tier(
        axes[3], r_gini, h_gini, tier,
        xlabel="R_HAPY_GINI", ylabel="H_HAPY_GINI",
        title="Hα vs stellar Gini",
        one_to_one=True
    )

    # 5. Halpha vs stellar M20
    scatter_by_tier(
        axes[4], r_m20, h_m20, tier,
        xlabel="R_HAPY_M20", ylabel="H_HAPY_M20",
        title="Hα vs stellar M20",
        one_to_one=True
    )

    # 6. Delta plane
    scatter_by_tier(
        axes[5],
        safe_float_array(tab, "DELTA_GINI"),
        safe_float_array(tab, "DELTA_M20"),
        tier,
        xlabel="DELTA_GINI",
        ylabel="DELTA_M20",
        title="Differential morphology plane"
    )
    plt.sca(axes[5])
    plt.axhline(y=0,ls='--',color='k')
    plt.axvline(x=0,ls='--',color='k')    
    handles, labels = axes[5].get_legend_handles_labels()
    if len(handles) > 0:
        fig.legend(handles, labels, title="QC_TIER", loc="upper right")

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_scatter(tab: Table, xcol: str, ycol: str, outpath: Path,
                 title: str | None = None,
                 one_to_one: bool = False,
                 logx: bool = False,
                 logy: bool = False) -> None:
    x = safe_float_array(tab, xcol)
    y = safe_float_array(tab, ycol)
    tier = safe_str_array(tab, "QC_TIER", default="?")

    fig = plt.figure(figsize=(6, 5))
    ax = plt.gca()
    scatter_by_tier(ax, x, y, tier, xcol, ycol,
                    title or f"{ycol} vs {xcol}",
                    one_to_one=one_to_one, logx=logx, logy=logy)
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 0:
        ax.legend(title="QC_TIER")
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_hist(tab: Table, col: str, outpath: Path, title: str | None = None, logx: bool = False) -> None:
    x = safe_float_array(tab, col)
    good = np.isfinite(x)
    if logx:
        good &= (x > 0)
        x = np.log10(x)

    if np.sum(good) == 0:
        return

    fig = plt.figure(figsize=(6, 4))
    plt.hist(x[good], bins=30)
    plt.xlabel(f"log10({col})" if logx else col)
    plt.ylabel("N")
    plt.title(title or col)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def write_summary_subset(tab: Table, outpath: Path) -> None:
    wanted = [
        "VFID", "GALNAME", "OBJID", "TELESCOPE", "DATEOBS", "TAG",
        "QC_TIER",
        "R50_ARCSEC", "H50_ARCSEC", "H50_R50_RATIO",
        "H_MAXDET_ARCSEC", "H_MAXDET_R25_RATIO",
        "R_HAPY_GINI", "H_HAPY_GINI", "DELTA_GINI",
        "R_HAPY_M20", "H_HAPY_M20", "DELTA_M20",
        "H_HAPY_NPIX", "H_HAPY_SNP_DET",
        "USE_R_STRUCTURE", "USE_HA_EXTENT", "USE_HA_MORPH",
    ]
    keep = [c for c in wanted if c in tab.colnames]
    tab[keep].write(outpath, format="fits", overwrite=True)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="First-look science plots for merged HAPY results.")
    parser.add_argument("table", help="Merged HAPY results table, e.g. merged_results.fits")
    parser.add_argument("--outdir", default="science_firstlook", help="Output directory")
    parser.add_argument("--sample", default="AB", choices=["A", "AB", "ABC", "ALL"],
                        help="QC sample to use")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    plotdir = outdir / "plots"
    tabledir = outdir / "tables"
    ensure_dir(outdir)
    ensure_dir(plotdir)
    ensure_dir(tabledir)

    tab = Table.read(args.table)
    print(f"Read {len(tab)} rows from {args.table}")

    tab = add_science_columns(tab)
    tab = add_qc_flags(tab)

    sample_mask = select_sample(tab, args.sample)
    sub = tab[sample_mask]
    print(f"Selected {len(sub)} rows for sample {args.sample}")

    # save selected sample table
    write_summary_subset(sub, tabledir / f"science_firstlook_{args.sample}.fits")

    # dashboard
    plot_firstlook_dashboard(sub, plotdir / f"dashboard_{args.sample}.png", args.sample)

    # core science plots
    plot_scatter(
        sub, "R_HAPY_GINI", "H_HAPY_GINI",
        plotdir / f"HGini_vs_RGini_{args.sample}.png",
        title="Hα vs stellar Gini",
        one_to_one=True
    )

    plot_scatter(
        sub, "R_HAPY_M20", "H_HAPY_M20",
        plotdir / f"HM20_vs_RM20_{args.sample}.png",
        title="Hα vs stellar M20",
        one_to_one=True
    )

    plot_scatter(
        sub, "DELTA_GINI", "DELTA_M20",
        plotdir / f"DeltaGini_vs_DeltaM20_{args.sample}.png",
        title="Differential morphology plane"
    )

    plot_scatter(
        sub, "H50_R50_RATIO", "DELTA_GINI",
        plotdir / f"H50R50_vs_DeltaGini_{args.sample}.png",
        title="Extent ratio vs differential Gini"
    )

    plot_scatter(
        sub, "H50_R50_RATIO", "DELTA_M20",
        plotdir / f"H50R50_vs_DeltaM20_{args.sample}.png",
        title="Extent ratio vs differential M20"
    )

    plot_scatter(
        sub, "R50_ARCSEC", "H50_ARCSEC",
        plotdir / f"H50_vs_R50_{args.sample}.png",
        title="Hα vs stellar half-light radius",
        one_to_one=True
    )

    plot_hist(
        sub, "H50_R50_RATIO",
        plotdir / f"H50_R50_RATIO_hist_{args.sample}.png",
        title="Hα / stellar half-light ratio"
    )

    plot_hist(
        sub, "H_MAXDET_R25_RATIO",
        plotdir / f"HMAX_R25_RATIO_hist_{args.sample}.png",
        title="Hα max extent / R25"
    )

    print(f"Wrote first-look science products to {outdir}")


if __name__ == "__main__":
    main()
