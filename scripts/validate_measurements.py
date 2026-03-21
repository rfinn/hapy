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
    safe_float_array,
    safe_bool_array,
    safe_str_array,
    add_science_columns,
    add_qc_flags,
    select_sample,
)

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def make_dataframe(tab: Table, columns: list[str]) -> pd.DataFrame:
    data = {}
    for col in columns:
        if col in tab.colnames:
            if tab[col].dtype.kind in ("b",):
                data[col] = np.array(tab[col], dtype=bool)
            else:
                try:
                    data[col] = np.array(tab[col], dtype=float)
                except Exception:
                    data[col] = np.array(tab[col]).astype(str)
    return pd.DataFrame(data)


def clean_pairplot_df(
    df: pd.DataFrame,
    positive_cols: list[str] | None = None,
    bad_sentinels: tuple[float, ...] = (-99, -999, 99, 999),
) -> pd.DataFrame:
    """
    Clean dataframe for pairplots.

    - replaces common sentinel values with NaN
    - enforces positivity for selected columns
    - drops rows with NaN in plotted columns
    """
    out = df.copy()
    positive_cols = positive_cols or []

    # replace common sentinel values
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            for bad in bad_sentinels:
                out.loc[np.isclose(out[col], bad, equal_nan=False), col] = np.nan

    # require positive values where appropriate
    for col in positive_cols:
        if col in out.columns:
            out.loc[out[col] <= 0, col] = np.nan

    return out.dropna()

def _robust_limits(x, qlo=0.01, qhi=0.99, pad_frac=0.05):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return None

    lo, hi = np.nanquantile(x, [qlo, qhi])
    if not np.isfinite(lo) or not np.isfinite(hi):
        return None
    if lo == hi:
        delta = 0.5 if lo == 0 else 0.05 * abs(lo)
        return lo - delta, hi + delta

    pad = pad_frac * (hi - lo)
    return lo - pad, hi + pad

def _annotate_pairgrid(g, df: pd.DataFrame, hue: str | None = None,
                       qlo=0.01, qhi=0.99):
    """
    Add robust limits and ratio annotations to off-diagonal panels
    in a seaborn PairGrid / pairplot.
    """
    vars_ = list(g.x_vars)

    for i, yvar in enumerate(vars_):
        for j, xvar in enumerate(vars_):
            ax = g.axes[i, j]
            if ax is None:
                continue

            # robust limits for all panels
            if xvar in df.columns:
                xlim = _robust_limits(df[xvar].values, qlo=qlo, qhi=qhi)
                if xlim is not None:
                    ax.set_xlim(xlim)

            if yvar in df.columns:
                ylim = _robust_limits(df[yvar].values, qlo=qlo, qhi=qhi)
                if ylim is not None:
                    ax.set_ylim(ylim)

            # only annotate off-diagonal panels
            if i == j:
                continue

            x = np.asarray(df[xvar], dtype=float)
            y = np.asarray(df[yvar], dtype=float)
            good = np.isfinite(x) & np.isfinite(y) & (x != 0)

            if np.sum(good) == 0:
                continue

            ratio = y[good] / x[good]
            ratio = ratio[np.isfinite(ratio)]

            if len(ratio) == 0:
                continue

            mean_ratio = np.nanmean(ratio)
            std_ratio = np.nanstd(ratio)
            n = len(ratio)

            ax.text(
                0.04, 0.96,
                f"N={n}\n⟨y/x⟩={mean_ratio:.2f}\nσ={std_ratio:.2f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"),
            )

def pairplot_family(
    df: pd.DataFrame,
    hue: str,
    outpath: Path,
    title: str,
    corner: bool = True,
    annotate_ratio: bool = True,
):
    if len(df) < 3 or len(df.columns) < 2:
        print(f"WARNING: insufficient data for {outpath.name}")
        return

    sns.set_context("talk")

    g = sns.pairplot(
        df,
        hue=hue if hue in df.columns else None,
        corner=corner,
        diag_kind="hist",
        plot_kws=dict(s=22, alpha=0.7),
        diag_kws=dict(bins=20),
    )

    if annotate_ratio:
        _annotate_pairgrid(g, df, hue=hue if hue in df.columns else None)
    else:
        # still apply robust limits
        _annotate_pairgrid(g, df, hue=hue if hue in df.columns else None)
        # remove text labels if desired
        for axrow in g.axes:
            for ax in axrow:
                if ax is None:
                    continue
                for txt in list(ax.texts):
                    txt.remove()

    g.figure.suptitle(title, y=1.02)
    g.figure.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(g.figure)
    
def plot_difference_hist(
    tab,
    col1,
    col2,
    outpath,
    title=None,
    xlabel=None,
    relative=False,
    positive_only=False,
    bad_sentinels=(-99, -999, 99, 999),
    qclip=None,
):
    """
    Plot histogram of differences between two matched columns.

    Parameters
    ----------
    tab : astropy.table.Table
        Input table.

    col1, col2 : str
        Columns to compare. Histogram is for:
            diff = col2 - col1
        unless relative=True, in which case:
            diff = (col2 - col1) / col1

    outpath : str or Path
        Output figure path.

    title : str or None
        Plot title.

    xlabel : str or None
        X-axis label.

    relative : bool
        If True, plot fractional difference.

    positive_only : bool
        If True, require both columns > 0.

    bad_sentinels : tuple
        Sentinel values to treat as invalid.

    qclip : tuple or None
        Optional quantile clipping range, e.g. (0.01, 0.99)
        to suppress extreme outliers from setting the x-range.

    Notes
    -----
    Annotates the panel with:
      - N
      - median difference
      - standard deviation
    """
    import numpy as np
    import matplotlib.pyplot as plt

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

    x1 = _safe_float_array(tab, col1)
    x2 = _safe_float_array(tab, col2)

    # remove sentinel values
    for bad in bad_sentinels:
        x1[np.isclose(x1, bad, equal_nan=False)] = np.nan
        x2[np.isclose(x2, bad, equal_nan=False)] = np.nan

    good = np.isfinite(x1) & np.isfinite(x2)

    if positive_only:
        good &= (x1 > 0) & (x2 > 0)

    if np.sum(good) == 0:
        print(f"WARNING: no finite matched values for {col1} and {col2}")
        return

    if relative:
        good &= (x1 != 0)
        diff = (x2[good] - x1[good]) / x1[good]
        if xlabel is None:
            xlabel = f"({col2} - {col1}) / {col1}"
    else:
        diff = x2[good] - x1[good]
        if xlabel is None:
            xlabel = f"{col2} - {col1}"

    diff = diff[np.isfinite(diff)]
    if len(diff) == 0:
        print(f"WARNING: no finite differences for {col1} and {col2}")
        return

    # optional clipping for display only
    plot_diff = diff
    if qclip is not None and len(diff) > 5:
        lo, hi = np.nanquantile(diff, qclip)
        plot_diff = diff[(diff >= lo) & (diff <= hi)]

    med = np.nanmedian(diff)
    std = np.nanstd(diff)
    mean = np.nanmean(diff)
    n = len(diff)

    fig = plt.figure(figsize=(6, 4.5))
    ax = plt.gca()

    ax.hist(plot_diff, bins=30)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.axvline(med, color="k", ls=":", lw=1)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("N")
    ax.set_title(title if title is not None else f"{col2} - {col1}")

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
# ----------------------------------------------------------------------
# main validation products
# ----------------------------------------------------------------------

def build_r_full_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "R24_ARCSEC",
        "R25_ARCSEC",
        "R75_ARCSEC",
        "R_SM_R80",
        "R_SM_RMAX_CIRCLE",
        "R_SM_RMAX_ELLIP",        
        "QC_TIER",
    ]
    if "GAL_RE" in tab.colnames:
        cols.insert(-1, "GAL_RE")
    df = make_dataframe(tab, cols)
    return clean_pairplot_df(df, positive_cols=[c for c in cols if c != "QC_TIER"])

def build_r_half_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "R50_ARCSEC",
        "R_PETRO_R50_ARCSEC",
        "R_EXPFIT_RE_ARCSEC",
        "GAL_RE",
        "R_SM_R50",
        "R_SM_RHALF_ELLIP",        
        "QC_TIER",
    ]
    if "GAL_RE" in tab.colnames:
        cols.insert(-1, "GAL_RE")
    df = make_dataframe(tab, cols)
    return clean_pairplot_df(df, positive_cols=[c for c in cols if c != "QC_TIER"])


def build_h_full_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "H75_ARCSEC",
        "H_MAXDET_ARCSEC",
        "H_R95_R24_ARCSEC",
        "H_SM_R80",
        "H_SM_RMAX_CIRCLE",
        "H_SM_RMAX_ELLIP",                
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    return clean_pairplot_df(df, positive_cols=[c for c in cols if c != "QC_TIER"])

def build_h_half_size_df(tab: Table) -> pd.DataFrame:
    cols = [
        "H25_ARCSEC",
        "H50_ARCSEC",
        "H_PETRO_R50_ARCSEC",
        "H_EXPFIT_RE_ARCSEC",
        "H_SM_R50",
        "H_SM_RHALF_ELLIP",        
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    return clean_pairplot_df(df, positive_cols=[c for c in cols if c != "QC_TIER"])


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
    positive = ["R24_FLUX_CGS", "R_PETRO_FLUX_CGS", "R_C30"]
    return clean_pairplot_df(df, positive_cols=positive)


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
    return clean_pairplot_df(df, positive_cols=positive)


def build_r_morph_df(tab: Table) -> pd.DataFrame:
    cols = [
        "R_HAPY_GINI",
        "R_HAPY_M20",
        "R_SM_GINI",
        "R_SM_M20",
        #"R_ASYM",
        "R_C30",
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    return clean_pairplot_df(df)


def build_h_morph_df(tab: Table) -> pd.DataFrame:
    cols = [
        "H_HAPY_GINI",
        "H_HAPY_M20",
        "H_SM_GINI",
        "H_SM_M20",
        #"H_ASYM",
        "H_C30_R24",
        "QC_TIER",
    ]
    df = make_dataframe(tab, cols)
    return clean_pairplot_df(df)


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

    tab = add_science_columns(tab)
    tab = add_qc_flags(tab)

    sample_mask = select_sample(tab, args.sample)
    sub = tab[sample_mask]
    print(f"Selected {len(sub)} rows for sample {args.sample}")

    # Pairplots

    ##################################################################
    ## SIZE METRICS
    ##################################################################    
    r_full_size_df = build_r_full_size_df(sub)
    pairplot_family(r_full_size_df, "QC_TIER", outdir / f"pairplot_r_full_sizes_{args.sample}.png",
                    f"r-band size validation ({args.sample})", annotate_ratio=True)
    r_half_size_df = build_r_half_size_df(sub)
    pairplot_family(r_half_size_df, "QC_TIER", outdir / f"pairplot_r_half_sizes_{args.sample}.png",
                    f"r-band size validation ({args.sample})", annotate_ratio=True)
    
    h_full_size_df = build_h_full_size_df(sub)
    pairplot_family(h_full_size_df, "QC_TIER", outdir / f"pairplot_h_full_sizes_{args.sample}.png",
                    f"Hα size validation ({args.sample})", annotate_ratio=True)
    
    h_half_size_df = build_h_half_size_df(sub)
    pairplot_family(h_half_size_df, "QC_TIER", outdir / f"pairplot_h_half_sizes_{args.sample}.png",
                    f"Hα size validation ({args.sample})", annotate_ratio=True)

    r_fluxmag_df = build_r_fluxmag_df(sub)    
    pairplot_family(r_fluxmag_df, "QC_TIER", outdir / f"pairplot_r_fluxmag_{args.sample}.png",
                    f"r-band magnitude/flux validation ({args.sample})", annotate_ratio=True)

    h_fluxmag_df = build_h_fluxmag_df(sub)
    pairplot_family(h_fluxmag_df, "QC_TIER", outdir / f"pairplot_h_flux_{args.sample}.png",
                    f"Hα flux validation ({args.sample})", annotate_ratio=True)

    r_morph_df = build_r_morph_df(sub)
    pairplot_family(r_morph_df, "QC_TIER", outdir / f"pairplot_r_morph_{args.sample}.png",
                    f"r-band morphology validation ({args.sample})", annotate_ratio=False)

    h_morph_df = build_h_morph_df(sub)
    pairplot_family(h_morph_df, "QC_TIER", outdir / f"pairplot_h_morph_{args.sample}.png",
                    f"Hα morphology validation ({args.sample})", annotate_ratio=False)


    plot_difference_hist(
        sub,
        "R50_ARCSEC",
        "R_PETRO_R50_ARCSEC",
        outdir / "diff_R50_vs_RPETRO50.png",
        title="r-band half-light radius differences",
        )

    plot_difference_hist(
        sub,
        "R50_ARCSEC",
        "R_EXPFIT_RE_ARCSEC",
        outdir / "fracdiff_R50_vs_REXPFIT.png",
        title="r-band effective radius fractional differences",
        relative=True,
        positive_only=True,
        )
    
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
        wrap_deg=180.0,
        convert2=pa_ccw_north_to_photutils_theta,
        )
    print(f"Wrote validation products to {outdir}")


if __name__ == "__main__":
    main()
