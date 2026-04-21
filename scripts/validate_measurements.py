#!/usr/bin/env python

"""
validate_measurements.py

Measurement-validation plots for merged HAPY results.

This script focuses on internal consistency among:
- size metrics
- magnitudes / fluxes
- morphology measurements

It is intended for validation, not operational QC.

Example
-------
python validate_measurements.py merged_results.fits --outdir validation --sample AB
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
import seaborn as sns
import matplotlib.pyplot as plt

from hapy.utils.results_table import (
    ensure_dir,
    safe_float_array,
    safe_bool_array,
    safe_str_array,
    prepare_analysis_table,
    select_sample,
)
from hapy.utils.plotting import (
    QC_TIER_ORDER, QC_TIER_PALETTE, PAIRPLOT_LABELS,
    style_pairplot, enforce_qc_tier,
    make_dataframe, clean_pairplot_df,
    pairplot_family, plot_difference_hist,
    _robust_limits, _annotate_pairgrid
    )

  
    
#######################################
# Plots
#######################################

def plot_wrapped_angle_difference_hist(
    tab,
    col1,
    col2,
    outpath,
    title=None,
    xlabel=None,
    wrap_deg=180.0,
    convert1=None,
    convert2=None,
    bad_sentinels=(-99, -999, 99, 999),
    qclip=None,
):
    """
    Plot histogram of wrapped angle differences between two matched columns.

    Parameters
    ----------
    tab : astropy.table.Table
        Input table.

    col1, col2 : str
        Angle columns in degrees.
        Histogram is for wrapped(col2 - col1).

    outpath : str or Path
        Output figure path.

    title : str or None
        Plot title.

    xlabel : str or None
        X-axis label.

    wrap_deg : float
        Wrapping period in degrees.
        Use 180 for position angles with 180-degree symmetry.
        Use 360 for fully directed angles.

    convert1, convert2 : callable or None
        Optional functions applied to col1 / col2 values before differencing.
        Example:
            convert2 = hapy.geometry.adapters.pa_ccw_north_to_photutils_theta

    bad_sentinels : tuple
        Sentinel values to treat as invalid.

    qclip : tuple or None
        Optional display clipping range in quantiles, e.g. (0.01, 0.99).

    Notes
    -----
    The wrapped difference is mapped to:
        [-wrap_deg/2, +wrap_deg/2)
    """
    import numpy as np
    import matplotlib.pyplot as plt

    def _apply_scalar_converter(arr, func):
        good = np.isfinite(arr)
        out = arr.copy()
        out[good] = np.array([func(v) for v in arr[good]], dtype=float)
        return out

    def _safe_float_array(tab, colname, default=np.nan):
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

    def _wrap_diff_deg(diff, wrap_deg=180.0):
        return ((diff + 0.5 * wrap_deg) % wrap_deg) - 0.5 * wrap_deg

    a1 = _safe_float_array(tab, col1)
    a2 = _safe_float_array(tab, col2)

    for bad in bad_sentinels:
        a1[np.isclose(a1, bad, equal_nan=False)] = np.nan
        a2[np.isclose(a2, bad, equal_nan=False)] = np.nan

    if convert1 is not None:
        a1 = _apply_scalar_converter(a1, convert1)

    if convert2 is not None:
        a2 = _apply_scalar_converter(a2, convert2)
        
        

    good = np.isfinite(a1) & np.isfinite(a2)
    if np.sum(good) == 0:
        print(f"WARNING: no finite matched values for {col1} and {col2}")
        return

    diff = _wrap_diff_deg(a2[good] - a1[good], wrap_deg=wrap_deg)
    diff = diff[np.isfinite(diff)]

    if len(diff) == 0:
        print(f"WARNING: no finite wrapped differences for {col1} and {col2}")
        return

    plot_diff = diff
    if qclip is not None and len(diff) > 5:
        lo, hi = np.nanquantile(diff, qclip)
        plot_diff = diff[(diff >= lo) & (diff <= hi)]

    mean = np.nanmean(diff)
    med = np.nanmedian(diff)
    std = np.nanstd(diff)
    n = len(diff)

    fig = plt.figure(figsize=(6, 4.5))
    ax = plt.gca()

    bins = np.linspace(-0.5 * wrap_deg, 0.5 * wrap_deg, 31)
    ax.hist(plot_diff, bins=bins)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.axvline(med, color="k", ls=":", lw=1)

    if xlabel is None:
        xlabel = f"wrapped({col2} - {col1}) [deg]"

    ax.set_xlabel(xlabel)
    ax.set_ylabel("N")
    ax.set_title(title if title is not None else f"Wrapped angle difference: {col2} - {col1}")

    ax.text(
        0.97, 0.97,
        f"N={n}\nmean={mean:.3g}\nmedian={med:.3g}\nstd={std:.3g}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.8, edgecolor="none"),
    )

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_hist_robust(
    tab,
    col,
    outpath,
    title=None,
    logx=False,
    positive_only=False,
    bad_sentinels=(-99, -999, 99, 999),
    qrange=(0.01, 0.99),
    bins=30,
):
    """
    Plot histogram using robust quantile-based x limits.

    This prevents a few extreme outliers from forcing useless binning.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    def _safe_float_array(tab, colname, default=np.nan):
        if colname not in tab.colnames:
            return np.full(len(tab), default, dtype=float)

        coldata = tab[colname]
        try:
            if hasattr(coldata, "filled"):
                coldata = coldata.filled(default)
        except Exception:
            pass

        out = np.full(len(tab), default, dtype=float)
        for i, v in enumerate(coldata):
            try:
                out[i] = float(v)
            except Exception:
                out[i] = default
        return out

    x = _safe_float_array(tab, col)

    # remove sentinel values
    for bad in bad_sentinels:
        x[np.isclose(x, bad, equal_nan=False)] = np.nan

    good = np.isfinite(x)

    if positive_only:
        good &= (x > 0)

    x = x[good]
    if len(x) == 0:
        print(f"WARNING: no finite values for {col}")
        return

    if logx:
        x = x[x > 0]
        if len(x) == 0:
            print(f"WARNING: no positive values for log histogram of {col}")
            return
        x = np.log10(x)

    # robust display range
    lo, hi = np.nanquantile(x, qrange)

    if not np.isfinite(lo) or not np.isfinite(hi):
        print(f"WARNING: could not determine plotting range for {col}")
        return

    if lo == hi:
        delta = 0.5 if lo == 0 else 0.05 * abs(lo)
        lo -= delta
        hi += delta

    fig = plt.figure(figsize=(6, 4.5))
    ax = plt.gca()
    ax.hist(x, bins=bins, range=(lo, hi))
    ax.set_xlabel(f"log10({col})" if logx else col)
    ax.set_ylabel("N")
    ax.set_title(title if title is not None else col)

    ax.text(
        0.97, 0.97,
        f"N={len(x)}\nmed={np.nanmedian(x):.3g}\nstd={np.nanstd(x):.3g}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.8, edgecolor="none"),
    )

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
# main validation products
# ----------------------------------------------------------------------

def build_r_full_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "R24_ARCSEC",
        "R25_ARCSEC",
        "R75_ARCSEC",
        "R_SM_R80_ARCSEC",
        "R_SM_RMAX_CIRCLE_ARCSEC",
        "R_SM_RMAX_ELLIP_ARCSEC",        
        "QC_TIER",
        "TELESCOPE",        
    ]
    df = make_dataframe(tab, cols)
    size_cols = [c for c in cols if c != "QC_TIER"]
    log_cols = size_cols
    df = clean_pairplot_df(df,positive_cols=size_cols,
                                 log_cols=log_cols,
                                 )
    
    return  df, log_cols
        
 
def build_r_half_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "R50_ARCSEC",
        "R_PETRO_R50_ARCSEC",
        "R_EXPFIT_RE_ARCSEC",
        "GAL_RE_ARCSEC",
        "R_SM_R50_ARCSEC",
        "R_SM_RHALF_ELLIP_ARCSEC",        
        "QC_TIER",
        "TELESCOPE",
    ]
    df = make_dataframe(tab, cols)
    size_cols = [c for c in cols if c != "QC_TIER"]
    log_cols = size_cols
    df = clean_pairplot_df(df,positive_cols=size_cols,
                                 log_cols=log_cols,
                                 )
    
    return  df, log_cols




def build_h_full_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "H75_ARCSEC",
        "H_MAXDET_ARCSEC",
        "H_R95_R24_ARCSEC",
        "H_SM_R80_ARCSEC",
        "H_ISO5E17_ARCSEC",
        "H_ISO17E18_ARCSEC",        
        "H_SM_RMAX_CIRCLE_ARCSEC",
        "H_SM_RMAX_ELLIP_ARCSEC",
        "QC_TIER",
        "TELESCOPE",

    ]
    df = make_dataframe(tab, cols)
    size_cols = [c for c in cols if c != "QC_TIER"]
    log_cols = size_cols
    df = clean_pairplot_df(df,positive_cols=size_cols,
                                 log_cols=log_cols,
                                 )
    
    return  df, log_cols

def build_h_half_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "H25_ARCSEC",
        "H50_ARCSEC",
        "H_PETRO_R50_ARCSEC",
        "H_EXPFIT_RE_ARCSEC",
        "H_SM_R50_ARCSEC",
        "H_SM_RHALF_ELLIP_ARCSEC",        
        "QC_TIER",
        "TELESCOPE",
    ]
    df = make_dataframe(tab, cols)
    size_cols = [c for c in cols if c != "QC_TIER"]
    log_cols = size_cols
    df = clean_pairplot_df(df,positive_cols=size_cols,
                                 log_cols=log_cols,
                                 )
    
    return  df, log_cols
    


def build_r_fluxmag_df(tab: Table) -> pd.DataFrame:
    cols = [
        "R24_MAG",
        "R25_ISO_MAG",
        "R25P5_MAG",
        "R_PETRO_MAG",
        "GAL_MAG",
        "R24_FLUX_CGS",
        "R_PETRO_FLUX_CGS",
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    positive = ["R24_FLUX_CGS", "R_PETRO_FLUX_CGS"]
    log_cols = positive
    df = clean_pairplot_df(df,positive_cols=positive,
                                 log_cols=log_cols,
                                 )
    
    return  df, log_cols
    

def build_h_fluxmag_df(tab: Table) -> pd.DataFrame:
    cols = [
        "H_TOT_FLUX_CGS",
        "H_R24_FLUX_CGS",
        "H30R24_FLUX_CGS",
        "H_ISO5E17_FLUX_CGS",
        "H_ISO17E18_FLUX_CGS",
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    positive = [c for c in cols if c not in ("QC_TIER",)]
    log_cols = positive
    df = clean_pairplot_df(df,positive_cols=positive,
                                 log_cols=log_cols,
                                 )
    
    return  df, log_cols
    


def build_r_morph_df(tab: Table) -> pd.DataFrame:
    cols = [
        "R_HAPY_GINI",
        "R_HAPY_M20",
        "R_HAPY_ASYM",
        "R_SM_GINI",
        "R_SM_M20",
        "R_SM_C",
        "R_SM_A",
        "R_C30",
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    log_cols=None
    return clean_pairplot_df(df, positive_cols=["R_SM_C"],), log_cols


def build_h_morph_df(tab: Table) -> pd.DataFrame:
    cols = [
        "H_HAPY_GINI",
        "H_HAPY_M20",
        "R_HAPY_ASYM",
        "H_SM_GINI",
        "H_SM_M20",
        "H_SM_C",
        "H_SM_A",
        "H_C30_R24",
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    log_cols=None
    return clean_pairplot_df(df, positive_cols=["H_SM_C"],), log_cols



# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Validation plots for HAPY measurements.")
    parser.add_argument("table", help="Merged HAPY results table")
    parser.add_argument("--outdir", default="validation", help="Output directory")
    parser.add_argument("--sample", default="AB", choices=["A", "AB", "ABC", "ALL"],
                        help="QC sample to use")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    tab = Table.read(args.table)
    print(f"Read {len(tab)} rows from {args.table}")

    #tab = add_science_columns(tab)
    #tab = add_qc_flags(tab)
    tab = prepare_analysis_table(tab)
    sample_mask = select_sample(tab, args.sample)
    sub = tab[sample_mask]
    print(f"Selected {len(sub)} rows for sample {args.sample}")

    # Pairplots

    ##################################################################
    ## SIZE METRICS
    ##################################################################    
    r_full_size_df, log_cols = build_r_full_size_df(sub)
    pairplot_family(r_full_size_df, "QC_TIER", outdir / f"pairplot_r_full_sizes_{args.sample}.png",
    #pairplot_family(r_full_size_df, "TELESCOPE", outdir / f"pairplot_r_full_sizes_{args.sample}.png",
                    f"r-band full size validation ({args.sample})", annotate_ratio=True, log_cols=log_cols)
    r_half_size_df, log_cols = build_r_half_size_df(sub)
    pairplot_family(r_half_size_df, "QC_TIER", outdir / f"pairplot_r_half_sizes_{args.sample}.png",
    #pairplot_family(r_half_size_df, "TELESCOPE", outdir / f"pairplot_r_half_sizes_{args.sample}.png",    
                    f"r-band half size validation ({args.sample})", annotate_ratio=True, log_cols=log_cols)
    
    h_full_size_df, log_cols = build_h_full_size_df(sub)
    pairplot_family(h_full_size_df, "QC_TIER", outdir / f"pairplot_h_full_sizes_{args.sample}.png",
                    f"Hα size full validation ({args.sample})", annotate_ratio=True, log_cols=log_cols)
    
    h_half_size_df, log_cols = build_h_half_size_df(sub)
    pairplot_family(h_half_size_df, "QC_TIER", outdir / f"pairplot_h_half_sizes_{args.sample}.png",
                    f"Hα size half validation ({args.sample})", annotate_ratio=True, log_cols=log_cols)

    r_fluxmag_df, log_cols = build_r_fluxmag_df(sub)    
    pairplot_family(r_fluxmag_df, "QC_TIER", outdir / f"pairplot_r_fluxmag_{args.sample}.png",
                    f"r-band magnitude/flux validation ({args.sample})", annotate_ratio=True, log_cols=log_cols)

    h_fluxmag_df, log_cols = build_h_fluxmag_df(sub)
    pairplot_family(h_fluxmag_df, "QC_TIER", outdir / f"pairplot_h_flux_{args.sample}.png",
                    f"Hα flux validation ({args.sample})", annotate_ratio=True, log_cols=log_cols)

    r_morph_df, log_cols = build_r_morph_df(sub)
    pairplot_family(r_morph_df, "QC_TIER", outdir / f"pairplot_r_morph_{args.sample}.png",
                    f"r-band morphology validation ({args.sample})", annotate_ratio=False, log_cols=log_cols)

    h_morph_df, log_cols = build_h_morph_df(sub)
    pairplot_family(h_morph_df, "QC_TIER", outdir / f"pairplot_h_morph_{args.sample}.png",
                    f"Hα morphology validation ({args.sample})", annotate_ratio=False, log_cols=log_cols)


    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    ax = plt.gca()
    plot_difference_hist(
        sub,
        "R50_ARCSEC",
        "R_PETRO_R50_ARCSEC",
        outdir / "diff_R50_vs_RPETRO50.png",
        title="r-band half-light radius differences",
        plotsingle=False,
        ax=ax,
        )
    plt.subplot(1,2,2)
    ax = plt.gca()
    plot_difference_hist(
        sub,
        "R50_ARCSEC",
        "R_EXPFIT_RE_ARCSEC",
        outdir / "fracdiff_R50_vs_REXPFIT.png",
        title="r-band effective radius fractional differences",
        relative=True,
        positive_only=True,
        plotsingle=False,
        ax=ax
        )
    plt.savefig(outdir / "diff_RHALF_measurements.png")
    
    plot_difference_hist(
        sub,
        "ELL0_BA",
        "ELLIP_BA",
        outdir / "diff_ELL0BA_vs_ELLIPBA.png",
        title="Ellipse axis-ratio differences",
        )

    plot_difference_hist(
        sub,
        "ELL0_PA_DEG",
        "ELLIP_PA_DEG",
        outdir / "diff_ELL0PA_vs_ELLIPPA.png",
        title="Ellipse PA differences",
        qclip=(0.01, 0.99),
        )

    from hapy.geometry.adapters import pa_ccw_north_to_photutils_theta

    ellip_pa_deg = np.degrees(safe_float_array(tab, "ELLIP_THETA_RAD"))
    tab["ELLIP_PA_DEG"] = ellip_pa_deg
    plot_wrapped_angle_difference_hist(
        sub,
        "ELL0_PA_DEG",
        "GAL_PA",
        outdir / "diff_ELL0PA_vs_GALPA.png",
        title="Ellipse vs GALFIT PA differences",
        xlabel="Wrapped PA difference [deg]",
        wrap_deg=90.0,
        convert2=pa_ccw_north_to_photutils_theta,
        )
    print(f"Wrote validation products to {outdir}")


if __name__ == "__main__":
    main()
