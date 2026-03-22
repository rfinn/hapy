#!/usr/bin/env python

"""
validate_dashboards.py

Category-level validation dashboards for merged HAPY results.

This script complements validate_measurements.py by making family-specific
dashboards for:
- fluxes / magnitudes
- profile fitting (r and Halpha)
- statmorph / morphology (r and Halpha)
- GALFIT single-component r-band fits
- image quality (jointplots)

Examples
--------
python validate_dashboards.py merged_results.fits --sample AB --hue QC_TIER
python validate_dashboards.py merged_results.fits --sample ALL --hue TELESCOPE
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from astropy.table import Table

from hapy.utils.results_table import prepare_results_table, select_sample
from hapy.utils.plotting import (
    QC_TIER_ORDER,
    QC_TIER_PALETTE,
    PAIRPLOT_LABELS,
    style_pairplot,
    enforce_qc_tier,
)
from validate_measurements import (
    ensure_dir,
    safe_float_array,
    make_dataframe,
    clean_pairplot_df,
    _robust_limits,
    _annotate_pairgrid,
    pairplot_family,
    
    )



def plot_joint(
    tab: Table,
    xcol: str,
    ycol: str,
    outpath: Path,
    hue: str = "QC_TIER",
    logx: bool = False,
    logy: bool = False,
    title: str | None = None,
):
    if xcol not in tab.colnames or ycol not in tab.colnames:
        print(f"WARNING: missing {xcol} or {ycol}")
        return

    x = safe_float_array(tab, xcol)
    y = safe_float_array(tab, ycol)

    good = np.isfinite(x) & np.isfinite(y)
    if logx:
        good &= x > 0
    if logy:
        good &= y > 0

    df = pd.DataFrame({
        xcol: x[good],
        ycol: y[good],
    })

    if hue in tab.colnames:
        df[hue] = np.array(tab[hue])[good]

    if hue == "QC_TIER" and hue in df.columns:
        df = enforce_qc_tier(df)

    joint_kwargs = dict(
        data=df,
        x=xcol,
        y=ycol,
        hue=hue if hue in df.columns else None,
        kind="scatter",
        height=6,
    )

    if hue == "QC_TIER" and hue in df.columns:
        joint_kwargs["hue_order"] = QC_TIER_ORDER
        joint_kwargs["palette"] = QC_TIER_PALETTE

    g = sns.jointplot(**joint_kwargs)

    if logx:
        g.ax_joint.set_xscale("log")
    if logy:
        g.ax_joint.set_yscale("log")

    g.ax_joint.set_title(title or f"{ycol} vs {xcol}")
    plt.tight_layout()
    g.figure.savefig(outpath, dpi=150)
    plt.close(g.figure)


# ----------------------------------------------------------------------
# dataframe builders
# ----------------------------------------------------------------------

def build_fluxmag_df(tab: Table):
    cols = [
        "R24_MAG",
        "GAL_MAG",
        "R24_FLUX_CGS",
        "R_PETRO_FLUX_CGS",
        "H_TOT_FLUX_CGS",
        "H_R24_FLUX_CGS",
        "H30R24_FLUX_CGS",
        "QC_TIER",
        "TELESCOPE",
    ]

    log_cols = [
        "R24_FLUX_CGS",
        "R_PETRO_FLUX_CGS",
        "H_TOT_FLUX_CGS",
        "H_R24_FLUX_CGS",
        "H30R24_FLUX_CGS",
    ]

    df = make_dataframe(tab, cols)
    df = clean_pairplot_df(df, positive_cols=log_cols, log_cols=log_cols)
    return df, log_cols


def build_r_profile_df(tab: Table):
    cols = [
        "R50_ARCSEC",
        "R_PETRO_R50_ARCSEC",
        "R_EXPFIT_RE_ARCSEC",
        "R_LOGFIT_RE_ARCSEC",
        "R_PETRO_CON",
        "R_C30",
        "QC_TIER",
        "TELESCOPE",
    ]
    size_cols = [c for c in cols if c not in ["QC_TIER", "TELESCOPE", "R_PETRO_CON", "R_C30"]]
    log_cols = size_cols + ["R_PETRO_CON", "R_C30"]

    df = make_dataframe(tab, cols)
    df = clean_pairplot_df(df, positive_cols=log_cols, log_cols=log_cols)
    return df, log_cols


def build_h_profile_df(tab: Table):
    cols = [
        "H50_ARCSEC",
        "H_PETRO_R50_ARCSEC",
        "H_EXPFIT_RE_ARCSEC",
        "H_LOGFIT_RE_ARCSEC",
        "H_MAXDET_ARCSEC",
        "H_R95_R24_ARCSEC",
        "H_PETRO_CON",
        "QC_TIER",
        "TELESCOPE",
    ]
    log_cols = [c for c in cols if c not in ["QC_TIER", "TELESCOPE"]]

    df = make_dataframe(tab, cols)
    df = clean_pairplot_df(df, positive_cols=log_cols, log_cols=log_cols)
    return df, log_cols


def build_r_statmorph_df(tab: Table):
    cols = [
        "R_HAPY_GINI",
        "R_HAPY_M20",
        "R_SM_GINI",
        "R_SM_M20",
        "R_SM_R50",
        "R_SM_RHALF_ELLIP",
        "R_SM_C",
        "QC_TIER",
        "TELESCOPE",
    ]
    log_cols = ["R_SM_R50", "R_SM_RHALF_ELLIP", "R_SM_C"]

    df = make_dataframe(tab, cols)
    df = clean_pairplot_df(
        df,
        positive_cols=["R_SM_R50", "R_SM_RHALF_ELLIP", "R_SM_C"],
        nonzero_cols=["R_SM_C"],
        log_cols=log_cols,
    )
    return df, log_cols


def build_h_statmorph_df(tab: Table):
    cols = [
        "H_HAPY_GINI",
        "H_HAPY_M20",
        "H_SM_GINI",
        "H_SM_M20",
        "H_SM_R50",
        "H_SM_RHALF_ELLIP",
        "H_SM_C",
        "QC_TIER",
        "TELESCOPE",
    ]
    log_cols = ["H_SM_R50", "H_SM_RHALF_ELLIP", "H_SM_C"]

    df = make_dataframe(tab, cols)
    df = clean_pairplot_df(
        df,
        positive_cols=["H_SM_R50", "H_SM_RHALF_ELLIP", "H_SM_C"],
        nonzero_cols=["H_SM_C"],
        log_cols=log_cols,
    )
    return df, log_cols


def build_galfit_df(tab: Table):
    cols = [
        "R24_MAG",
        "GAL_MAG",
        "R50_ARCSEC",
        "GAL_RE",
        "GAL_N",
        "ELL0_BA",
        "GAL_BA",
        "QC_TIER",
        "TELESCOPE",
    ]
    log_cols = ["R50_ARCSEC", "GAL_RE"]

    df = make_dataframe(tab, cols)
    df = clean_pairplot_df(df, positive_cols=log_cols, log_cols=log_cols)
    return df, log_cols

def build_cgalfit_df(tab: Table):
    cols = [
        "R24_MAG",
        "GAL_CMAG",
        "R50_ARCSEC",
        "GAL_CRE",
        "GAL_CN",
        "ELL0_BA",
        "GAL_CBA",
        "QC_TIER",
        "TELESCOPE",
    ]
    log_cols = ["R50_ARCSEC", "GAL_RE"]

    df = make_dataframe(tab, cols)
    df = clean_pairplot_df(df, positive_cols=log_cols, log_cols=log_cols)
    return df, log_cols

# ----------------------------------------------------------------------
# main runner
# ----------------------------------------------------------------------

def run_dashboards(tab: Table, outdir: Path, hue: str):
    # flux / magnitude
    df, log_cols = build_fluxmag_df(tab)
    pairplot_family(
        df, hue,
        outdir / f"dashboard_fluxmag_{hue}.png",
        "Photometry / flux / concentration",
        log_cols=log_cols,
        use_robust_limits=False,
    )

    # profile r
    df, log_cols = build_r_profile_df(tab)
    pairplot_family(
        df, hue,
        outdir / f"dashboard_r_profile_{hue}.png",
        "r-band profile fitting",
        log_cols=log_cols,
        use_robust_limits=False,
    )

    # profile Halpha
    df, log_cols = build_h_profile_df(tab)
    pairplot_family(
        df, hue,
        outdir / f"dashboard_h_profile_{hue}.png",
        "Halpha profile fitting",
        log_cols=log_cols,
        use_robust_limits=False,
    )

    # statmorph r
    df, log_cols = build_r_statmorph_df(tab)
    pairplot_family(
        df, hue,
        outdir / f"dashboard_r_statmorph_{hue}.png",
        "r-band statmorph / morphology",
        log_cols=log_cols,
        use_robust_limits=False,
    )

    # statmorph Halpha
    df, log_cols = build_h_statmorph_df(tab)
    pairplot_family(
        df, hue,
        outdir / f"dashboard_h_statmorph_{hue}.png",
        "Halpha statmorph / morphology",
        log_cols=log_cols,
        use_robust_limits=False,
    )

    # GALFIT
    df, log_cols = build_galfit_df(tab)
    pairplot_family(
        df, hue,
        outdir / f"dashboard_galfit_{hue}.png",
        "GALFIT single-component Sérsic validation",
        log_cols=log_cols,
        use_robust_limits=False,
    )

    # GALFIT
    df, log_cols = build_cgalfit_df(tab)
    pairplot_family(
        df, hue,
        outdir / f"dashboard_galfit_{hue}.png",
        "GALFIT+Convolution single-component Sérsic validation",
        log_cols=log_cols,
        use_robust_limits=False,
    )
    
    # image quality
    plot_joint(
        tab,
        "R_FWHM_PSF",
        "H_FWHM_PSF",
        outdir / f"joint_fwhm_{hue}.png",
        hue=hue,
        title="R vs Halpha PSF FWHM",
    )

    plot_joint(
        tab,
        "R_SKYSTD_PHYS",
        "H_SKYSTD_PHYS",
        outdir / f"joint_skynoise_{hue}.png",
        hue=hue,
        logx=True,
        logy=True,
        title="R vs Halpha sky noise",
    )


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Validation dashboards for merged HAPY results.")
    parser.add_argument("table", help="Merged HAPY results table")
    parser.add_argument("--outdir", default="validate_dashboards", help="Output directory")
    parser.add_argument("--sample", default="AB", choices=["A", "AB", "ABC", "ALL"])
    parser.add_argument("--hue", default="QC_TIER", choices=["QC_TIER", "TELESCOPE"])
    args = parser.parse_args()

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    tab = Table.read(args.table)
    print(f"Read {len(tab)} rows from {args.table}")

    tab = prepare_results_table(tab)
    sub = tab[select_sample(tab, args.sample)]
    print(f"Selected {len(sub)} rows for sample {args.sample}")

    run_dashboards(sub, outdir, args.hue)

    print(f"Wrote validation dashboards to {outdir}")


if __name__ == "__main__":
    main()
