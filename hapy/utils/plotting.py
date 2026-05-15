import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from astropy.table import Table
from pathlib import Path
#from hapy.utils.results_table import safe_ratio
# -------------------------
# Global plotting constants
# -------------------------

QC_TIER_ORDER = ["A", "B", "C", "D", "E", "F"]

QC_TIER_PALETTE = {
    "A": "#1b9e77",
    "B": "#b2df8a",
    "C": "#e6ab02",
    "D": "#d95f02",
    "E": "#7570b3",
    "F": "#e7298a",
}

PAIRPLOT_LABELS = {

    # -------------------------
    # R-band sizes
    # -------------------------
    "R50_ARCSEC": r"$R_{50}^{R}$ (arcsec)",
    "R25_ARCSEC": r"$R_{25}^{R}$ (arcsec)",
    "R75_ARCSEC": r"$R_{75}^{R}$ (arcsec)",
    "R_PETRO_R50_ARCSEC": r"$R_{50,\mathrm{Petro}}^{R}$",
    "R_EXPFIT_RE_ARCSEC": r"$R_e^{R}$ (exp fit)",
    "R_LOGFIT_RE_ARCSEC": r"$R_{e}^{R}$ (log fit)",    
    "R_SM_R50_ARCSEC": r"$R_{50}^{R,\mathrm{SM}}$",
    "R_SM_R80_ARCSEC": r"$R_{80}^{R,\mathrm{SM}}$",    
    "R_SM_RHALF_ELLIP_ARCSEC": r"$R_{1/2}^{R,\mathrm{SM,ellip}}$",
    "R_SM_RHALF_CIRCLE_ARCSEC": r"$R_{1/2}^{R,\mathrm{SM,circle}}$",
    "R_SM_RMAX_ELLIP_ARCSEC": r"$R_{max}^{R,\mathrm{SM,ellip}}$",
    "R_SM_RMAX_CIRCLE_ARCSEC": r"$R_{max}^{R,\mathrm{SM,circle}}$",

    # -------------------------
    # Hα sizes
    # -------------------------
    "H50_ARCSEC": r"$R_{50}^{\mathrm{H\alpha}}$ (arcsec)",
    "H25_ARCSEC": r"$R_{25}^{\mathrm{H\alpha}}$",
    "H75_ARCSEC": r"$R_{75}^{\mathrm{H\alpha}}$",
    "H_R95_R24_ARCSEC": r"$R_{R95,R24}^{\mathrm{H\alpha}}$",    
    "H_MAXDET_ARCSEC": r"$R_{\max}^{\mathrm{H\alpha}}$",
    "H_PETRO_R50_ARCSEC": r"$R_{50,\mathrm{Petro}}^{\mathrm{H\alpha}}$",
    "H_LOGFIT_RE_ARCSEC": r"$R_{e}^{\mathrm{H\alpha}}$ (log fit)",
    "H_EXPFIT_RE_ARCSEC": r"$R_e^{\mathrm{H\alpha}}$ (exp fit)",
    "H_SM_R50_ARCSEC": r"$R_{50}^{\mathrm{H\alpha},\mathrm{SM}}$", 
    "H_SM_R80_ARCSEC": r"$R_{80}^{\mathrm{H\alpha},\mathrm{SM}}$",   
    "H_SM_RHALF_ELLIP_ARCSEC": r"$R_{1/2}^{\mathrm{H\alpha},\mathrm{SM,ellip}}$",
    "H_SM_RMAX_ELLIP_ARCSEC": r"$R_{max}^{\mathrm{H\alpha},\mathrm{SM,ellip}}$",    
    "H_SM_RMAX_CIRCLE_ARCSEC": r"$R_{max}^{\mathrm{H\alpha},\mathrm{SM,circle}}$",
    "H_ISO17E18_ARCSEC": r"$R_{ISO17E18}^{\mathrm{H\alpha}}$",
    "H_ISO5E17_ARCSEC": r"$R_{ISO5E17}^{\mathrm{H\alpha}}$",    
    
    # -------------------------
    # Morphology (HAPY + statmorph)
    # -------------------------
    "R_HAPY_GINI": r"$G^{R}$",
    "H_HAPY_GINI": r"$G^{\mathrm{H\alpha}}$",
    "R_HAPY_M20": r"$M_{20}^{R}$",
    "H_HAPY_M20": r"$M_{20}^{\mathrm{H\alpha}}$",

    "R_PETRO_CON": r"$C_{Petro}^{R}$",
    "H_PETRO_CON": r"$C_{Petro}^{\mathrm{H\alpha}}$",
    
    "R_SM_GINI": r"$G^{R}_{\mathrm{SM}}$",
    "H_SM_GINI": r"$G^{\mathrm{H\alpha}}_{\mathrm{SM}}$",
    "R_SM_M20": r"$M_{20,\mathrm{SM}}^{R}$",
    "H_SM_M20": r"$M_{20,\mathrm{SM}}^{\mathrm{H\alpha}}$",
    "R_SM_C": r"$C_{\mathrm{SM}}^{R}$",
    "H_SM_C": r"$C_{\mathrm{SM}}^{\mathrm{H\alpha}}$",

    "R_ASYM": r"$A^{R}$",
    "H_ASYM": r"$A^{\mathrm{H\alpha}}$",

    # -------------------------
    # Fluxes
    # -------------------------
    "R24_FLUX_CGS": r"$F_{R24}^{R}$",
    "R_PETRO_FLUX_CGS": r"$F_{\mathrm{Petro}}^{R}$",

    "H_TOT_FLUX_CGS": r"$F_{\mathrm{tot}}^{\mathrm{H\alpha}}$",
    "H_R24_FLUX_CGS": r"$F_{R24}^{\mathrm{H\alpha}}$",
    "H_ISO17E18_FLUX_CGS": r"$F_{ISO17E18}^{\mathrm{H\alpha}}$",
    "H_ISO5E17_FLUX_CGS": r"$F_{ISO5E17}^{\mathrm{H\alpha}}$",    
    "H30R24_FLUX_CGS": r"$F_{0.3R_{24}}^{\mathrm{H\alpha}}$",

    # -------------------------
    # Concentration
    # -------------------------
    "R_C30": r"$C_{30}^{R}$",
    "H_C30_R24": r"$C_{30}^{\mathrm{H\alpha}}$",

    # -------------------------
    # Magnitudes
    # -------------------------
    "R24_MAG": r"$m_{R24}^{R}$",
    "R25_ISO_MAG": r"$m_{25}^{R}$",
    "GAL_MAG": r"$m_{\mathrm{GALFIT}}$",
    "R_PETRO_MAG": r"$m_{\mathrm{Petro}}$",
    "R25P5_MAG": r"$m_{\mathrm{R25,P5}}$",        

    # -------------------------
    # Structural (GALFIT)
    # -------------------------
    "GAL_RE_ARCSEC": r"$R_e^{\mathrm{GALFIT}}$",
    "GAL_N": r"$n$ (Sérsic)",
    "GAL_BA": r"$b/a$",
    "GAL_PA": r"$\mathrm{PA}$ (deg)",

    "GAL_CRE_ARCSEC": r"$R_e^{\mathrm{GALFIT}}$",
    "GAL_CN": r"$n$ (Sérsic)",
    "GAL_CBA": r"$b/a$",
    "GAL_CPA": r"$\mathrm{PA}$ (deg)",
    
    # -------------------------
    # QC / meta
    # -------------------------
    "QC_TIER": "QC Tier",
    "TELESCOPE": "Telescope",
    "R_FWHM_PSF": r"$R \ FWHM\ (arcsec)$",
    "H_FWHM_PSF": r"$H\alpha \  FWHM \ (arcsec)$",    
    "R_SKYSTD_PHYS": r"$\sigma_{\mathrm{sky}}^{R} (erg~s^{-1}~cm^{-2})$",
    "H_SKYSTD_PHYS": r"$\sigma_{\mathrm{sky}}^{\mathrm{H\alpha}} (erg~s^{-1}~cm^{-2})$",    
}

def enforce_qc_tier(df):
    if "QC_TIER" in df.columns:
        df["QC_TIER"] = pd.Categorical(
            df["QC_TIER"],
            categories=QC_TIER_ORDER,
            ordered=True,
        )
    return df
# =========================
# Styling helpers
# =========================

def pretty_label(name: str) -> str:
    return PAIRPLOT_LABELS.get(name, name)

def style_pairplot(
    g,
    label_map=None,
    base_labelsize=22,
    base_ticksize=22,
    min_size=16,
    panel_size=2.5,
):
    """
    Style a seaborn PairGrid / pairplot.

    Features
    --------
    - automatically scales fonts with number of variables
    - automatically resizes the figure
    - optionally replaces long variable names with shorter labels

    Parameters
    ----------
    g : seaborn PairGrid
        Output of sns.pairplot(...)

    label_map : dict or None
        Optional mapping from column name -> display label

    base_labelsize : int
        Label size for small grids

    base_ticksize : int
        Tick size for small grids

    min_size : int
        Minimum font size

    panel_size : float
        Size in inches per panel
    """
    import numpy as np

    if label_map is None:
        label_map = {}

    # number of variables in the grid
    nvars = len(g.x_vars)

    # resize figure automatically
    g.fig.set_size_inches(panel_size * nvars, panel_size * nvars)

    # continuous scaling
    scale = min(1.0, 3.5 / max(nvars, 1))
    labelsize = max(min_size, base_labelsize * scale)
    ticksize = max(min_size - 1, base_ticksize * scale)

    def pretty_label(text):
        return label_map.get(text, text)

    # style all axes
    for i, row in enumerate(g.axes):
        for j, ax in enumerate(row):
            if ax is None:
                continue

            # only set outer labels in a pairplot
            xlabel = ax.get_xlabel()
            ylabel = ax.get_ylabel()

            if xlabel:
                ax.set_xlabel(pretty_label(xlabel), fontsize=labelsize)
            if ylabel:
                ax.set_ylabel(pretty_label(ylabel), fontsize=labelsize)

            ax.tick_params(axis="both", labelsize=ticksize)

    # legend handling
    if g._legend is not None:
        title = g._legend.get_title().get_text()
        g._legend.set_title(pretty_label(title))
        g._legend.get_title().set_fontsize(labelsize + 5)
        for text in g._legend.texts:
            text.set_fontsize(labelsize)

    # slightly improve spacing
    g.fig.tight_layout()

def style_jointplot(g, label_fs=14, tick_fs=11):
    try:
        xlab = PAIRPLOT_LABELS[g.ax_joint.get_xlabel()]
        ylab = PAIRPLOT_LABELS[g.ax_joint.get_ylabel()]
    except KeyError:
        xlab = g.ax_joint.get_xlabel()
        ylab = g.ax_joint.get_ylabel()
        
    g.ax_joint.set_xlabel(xlab, fontsize=label_fs)
    g.ax_joint.set_ylabel(ylab, fontsize=label_fs)
    g.ax_joint.tick_params(axis="both", labelsize=tick_fs)

    g.ax_marg_x.tick_params(axis="both", labelsize=tick_fs)
    g.ax_marg_y.tick_params(axis="both", labelsize=tick_fs)
    
def style_jointplot_legend(g, title_fs=13, text_fs=11, marker_size=40):
    leg = g.ax_joint.get_legend()
    if leg is None:
        return
    #leg.set_bbox_to_anchor((1.05, 1))
    leg._loc = 4  # lower left

    leg.set_title(leg.get_title().get_text(), prop={'size': title_fs})

    for text in leg.get_texts():
        text.set_fontsize(text_fs)

    for handle in leg.legend_handles:
        try:
            handle.set_sizes([marker_size])
        except Exception:
            pass

def raincloud_by_group(
    values_by_group,
    group_labels,
    xlabel,
    title=None,
    figsize=(10, 6),
    jitter=0.05,
    alpha_points=0.15,
    alpha_violin=0.35,
    alpha_box=0.35,
):
    """
    Make a horizontal raincloud plot.

    Parameters
    ----------
    values_by_group : list of 1D arrays
        One array of values per group.
    group_labels : list of str
        Labels for each group.
    xlabel : str
        Label for x-axis.
    title : str, optional
        Plot title.
    """
    #colors = [QC_TIER_PALETTE.get(t, "gray") for t in group_labels]

    fig, ax = plt.subplots(figsize=figsize)

    # drop non-finite values
    clean_values = []
    clean_labels = []
    for vals, lab in zip(values_by_group, group_labels):
        vals = np.asarray(vals)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        clean_values.append(vals)
        clean_labels.append(lab)

    if len(clean_values) == 0:
        print(f"WARNING: no valid values for raincloud plot: {title}")
        plt.close(fig)
        return None, None

    # boxplot
    bp = ax.boxplot(
        clean_values,
        patch_artist=True,
        vert=False,
        showfliers=False,
        medianprops=dict(color="k", linewidth=1.5),
        widths=0.12,
    )

    for patch in bp["boxes"]:
        patch.set_alpha(alpha_box)

    # violin plot
    vp = ax.violinplot(
        clean_values,
        points=300,
        showmeans=False,
        showextrema=False,
        showmedians=False,
        vert=False,
    )

    for i, body in enumerate(vp["bodies"]):
        yvals = body.get_paths()[0].vertices[:, 1]
        body.get_paths()[0].vertices[:, 1] = np.clip(yvals, i + 1, i + 1.55)
        body.set_alpha(alpha_violin)

    # jittered points
    rng = np.random.default_rng(12345)
    for i, vals in enumerate(clean_values):
        y = np.full(len(vals), i + 0.78)
        y += rng.uniform(low=-jitter, high=jitter, size=len(y))
        ax.scatter(vals, y, s=10, alpha=alpha_points)

    ax.set_yticks(np.arange(1, len(clean_labels) + 1))
    ax.set_yticklabels(clean_labels)
    ax.set_xlabel(xlabel)

    if title is not None:
        ax.set_title(title)

    plt.tight_layout()
    return fig, ax


def jointplot_with_hue(
    x,
    y,
    category,
    xname="x",
    yname="y",
    catname="category",
    kind="scatter",
    height=6,
    alpha=0.6,
    s=20,
):
    """
    Create a seaborn joint plot with hue from x, y, and a categorical array.

    Parameters
    ----------
    x, y : array-like
        Data arrays
    category : array-like
        Categorical labels (e.g., QC tier, telescope, etc.)
    xname, yname : str
        Axis labels
    catname : str
        Name of categorical column
    kind : str
        "scatter", "kde", or "hist"
    height : float
        Size of plot
    alpha : float
        Transparency for scatter
    s : float
        Marker size
    """

    x = np.asarray(x)
    y = np.asarray(y)
    category = np.asarray(category)

    # -----------------------------
    # remove NaNs
    # -----------------------------
    good = np.isfinite(x) & np.isfinite(y)

    ptab = Table([x[good],y[good],category[good]],names=[xname,yname,catname])
    df = ptab.to_pandas()
    # df = pd.DataFrame({
    #     xname: x[good],
    #     yname: y[good],
    #     catname: category[good],
    # })

    # -----------------------------
    # create jointplot
    # -----------------------------
    g = sns.jointplot(
        data=df,
        x=xname,
        y=yname,
        hue=catname,
        kind=kind,
        height=height,
        joint_kws=dict(alpha=alpha, s=s) if kind == "scatter" else None,
    )

    # -----------------------------
    # nicer styling
    # -----------------------------
    g.ax_joint.set_xlabel(xname)
    g.ax_joint.set_ylabel(yname)

    return g
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
    log_cols: list[str] | None = None,
    nonzero_cols: list[str] | None = None,
    bad_sentinels: tuple[float, ...] = (-99, -999, 99, 999),
) -> pd.DataFrame:
    """
    Clean dataframe for pairplots.

    - replaces common sentinel values with NaN
    - enforces positivity for selected numeric columns
    - optionally log10-transforms selected numeric columns
    - drops rows with NaN in plotted columns
    """
    out = df.copy()
    positive_cols = positive_cols or []
    log_cols = log_cols or []
    nonzero_cols = nonzero_cols or []

    # replace sentinel values only in numeric columns
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            for bad in bad_sentinels:
                out.loc[np.isclose(out[col], bad, equal_nan=False), col] = np.nan

    # require > 0 only for numeric columns
    for col in positive_cols:
        if col in out.columns and pd.api.types.is_numeric_dtype(out[col]):
            out.loc[out[col] <= 0, col] = np.nan

    # require != 0 only for numeric columns
    for col in nonzero_cols:
        if col in out.columns and pd.api.types.is_numeric_dtype(out[col]):
            out.loc[out[col] == 0, col] = np.nan

    # log-transform only numeric columns
    for col in log_cols:
        if col in out.columns and pd.api.types.is_numeric_dtype(out[col]):
            out.loc[out[col] <= 0, col] = np.nan
            out[col] = np.log10(out[col])

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

def _annotate_pairgrid(
    g,
    df: pd.DataFrame,
    hue: str | None = None,
    qlo=0.01,
    qhi=0.99,
    log_cols: list[str] | None = None,
    use_robust_limits: bool = True,
    use_robust_diag_limits: bool = False,
):
    """
    Add robust limits and panel annotations to off-diagonal panels
    in a seaborn PairGrid / pairplot.

    For linear quantities:
        annotate mean(y/x) and std(y/x)

    For logged quantities:
        annotate median(y-x) and std(y-x)
        where y-x = log10(y/x)
    """
    vars_ = list(g.x_vars)
    log_cols = log_cols or []

    for i, yvar in enumerate(vars_):
        for j, xvar in enumerate(vars_):
            ax = g.axes[i, j]
            if ax is None:
                continue

            
            # robust limits for off-diagonal scatter panels
            if use_robust_limits and i != j:
                if xvar in df.columns:
                    xlim = _robust_limits(df[xvar].values, qlo=qlo, qhi=qhi)
                    if xlim is not None:
                        ax.set_xlim(xlim)

                if yvar in df.columns:
                    ylim = _robust_limits(df[yvar].values, qlo=qlo, qhi=qhi)
                    if ylim is not None:
                        ax.set_ylim(ylim)

            # robust limits for diagonal histogram panels only
            if use_robust_diag_limits and i == j:
                if xvar in df.columns:
                    xlim = _robust_limits(df[xvar].values, qlo=qlo, qhi=qhi)
                    if xlim is not None:
                        ax.set_xlim(xlim)

            
            if i == j:
                continue

            x = np.asarray(df[xvar], dtype=float)
            y = np.asarray(df[yvar], dtype=float)
            good = np.isfinite(x) & np.isfinite(y)

            if np.sum(good) == 0:
                continue

            x = x[good]
            y = y[good]
            n = len(x)

            x_logged = xvar in log_cols
            y_logged = yvar in log_cols

            if x_logged and y_logged:
                d = y - x
                med_d = np.nanmedian(d)
                std_d = np.nanstd(d)
                text = (
                    f"N={n}\n"
                    f"med Δlog={med_d:.2f}\n"
                    f"σΔlog={std_d:.2f}"
                    )
 
            else:
                good_ratio = x != 0
                if np.sum(good_ratio) == 0:
                    continue
                ratio = y[good_ratio] / x[good_ratio]
                ratio = ratio[np.isfinite(ratio)]

                if len(ratio) == 0:
                    continue

                mean_ratio = np.nanmean(ratio)
                std_ratio = np.nanstd(ratio)

                text = (
                    f"N={n}\n"
                    f"⟨y/x⟩={mean_ratio:.2f}\n"
                    f"σ(y/x)={std_ratio:.2f}"
                )

            ax.text(
                0.04, 0.96,
                text,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    alpha=0.7,
                    edgecolor="none",
                ),
            )



    
#######################################
# Plots
#######################################
def pairplot_family(
    df: pd.DataFrame,
    hue: str,
    outpath: Path,
    title: str,
    corner: bool = True,
    annotate_ratio: bool = True,
    log_cols: list[str] | None = None,
    use_robust_limits: bool = True,   # NEW
    use_robust_diag_limits: bool = False,
):
    if len(df) < 3 or len(df.columns) < 2:
        print(f"WARNING: insufficient data for {outpath.name}")
        return

    pairplot_kwargs = dict(
        data=df,
        hue=hue if hue in df.columns else None,
        corner=corner,
        diag_kind="hist",
        plot_kws=dict(s=22, alpha=0.7),
        diag_kws=dict(bins=20),
    )

    if hue == "QC_TIER":
        from hapy.utils.plotting import QC_TIER_ORDER, QC_TIER_PALETTE, enforce_qc_tier
        df = enforce_qc_tier(df.copy())
        pairplot_kwargs["data"] = df
        pairplot_kwargs["hue_order"] = QC_TIER_ORDER
        pairplot_kwargs["palette"] = QC_TIER_PALETTE
        
    sns.set_context("talk")

    g = sns.pairplot(**pairplot_kwargs)

    
    if annotate_ratio:
        _annotate_pairgrid(
            g,
            pairplot_kwargs["data"],
            hue=hue if hue in df.columns else None,
            log_cols=log_cols,
            use_robust_limits=use_robust_limits,
            use_robust_diag_limits=use_robust_diag_limits,
            )
        #_annotate_pairgrid(g, pairplot_kwargs["data"], hue=hue if hue in df.columns else None, log_cols=log_cols)
    else:
        # still apply robust limits
        _annotate_pairgrid(
            g,
            pairplot_kwargs["data"],
            hue=hue if hue in df.columns else None,
            use_robust_limits=use_robust_limits,
            use_robust_diag_limits=use_robust_diag_limits,
            )
        #_annotate_pairgrid(g, pairplot_kwargs["data"], hue=hue if hue in df.columns else None)
        # remove text labels if desired
        for axrow in g.axes:
            for ax in axrow:
                if ax is None:
                    continue
                for txt in list(ax.texts):
                    txt.remove()
    style_pairplot(g, label_map=PAIRPLOT_LABELS)
    
    g.figure.suptitle(title, y=.97)
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
    plotsingle=True,
    ax=None,
    idx1=None,
    idx2=None,
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

    idx1 : array
        Optional array of indices; used for plotting duplicates

    idx2 : array
        Optional array of indices; used for plotting duplicates

    Notes
    -----
    Annotates the panel with:
      - N
      - median difference
      - standard deviation
    """
    import numpy as np
    import matplotlib.pyplot as plt

        
    if not plotsingle and ax is None:
        print("WARNING: must provide an axis if plotsingle=True")
        return
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

    if idx1 is not None and idx2 is not None:
        x1 = x1[idx1]
        x2 = x2[idx2]        

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

    if plotsingle:
        fig = plt.figure(figsize=(6, 4.5))
        ax = plt.gca()

    ax.hist(plot_diff, bins=30)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.axvline(med, color="k", ls=":", lw=1)

    ax.set_xlabel(xlabel,fontsize=10)
    ax.set_ylabel("N",fontsize=10)
    ax.set_title(title if title is not None else f"{col2} - {col1}",fontsize=10)

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
    if plotsingle:
        fig.savefig(outpath, dpi=150)
        #plt.close(fig)



import numpy as np
import matplotlib.pyplot as plt


def plot_with_residuals(
    x1, y1,
    x2=None, y2=None,
    xlabel1="x",
    ylabel1="y",
    xlabel2=None,
    ylabel2=None,
    title1=None,
    title2=None,
    label1="data",
    label2=None,
    residual="diff",   # "diff", "frac", or "logratio"
    one_to_one=True,
    figsize=None,
    s=20,
    alpha=0.7,
    outfile=None,
):
    """
    Plot y vs x in a top panel and residuals in a smaller bottom panel.

    If x2/y2 are supplied, make a 2-column version.

    Parameters
    ----------
    residual : {"diff", "frac", "logratio"}
        diff     : y - x
        frac     : (y - x) / x
        logratio : log10(y / x)
    """

    def clean_arrays(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        good = np.isfinite(x) & np.isfinite(y)
        if residual in ["frac", "logratio"]:
            good &= x != 0
        if residual == "logratio":
            good &= (x > 0) & (y > 0)
        return x[good], y[good]

    def calc_resid(x, y):
        if residual == "diff":
            return y - x, r"$y-x$"
        if residual == "frac":
            return (y - x) / x, r"$(y-x)/x$"
        if residual == "logratio":
            return np.log10(y / x), r"$\log_{10}(y/x)$"
        raise ValueError("residual must be 'diff', 'frac', or 'logratio'")

    def draw_column(ax_top, ax_res, x, y, xlabel, ylabel, title=None, label=None):
        x, y = clean_arrays(x, y)
        res, res_label = calc_resid(x, y)

        ax_top.scatter(x, y, s=s, alpha=alpha, label=label)

        if one_to_one and len(x) > 0:
            lo = np.nanmin([np.nanmin(x), np.nanmin(y)])
            hi = np.nanmax([np.nanmax(x), np.nanmax(y)])
            ax_top.plot([lo, hi], [lo, hi], "k--", lw=1)
            ax_top.set_xlim(lo, hi)
            ax_top.set_ylim(lo, hi)

        ax_res.axhline(0, color="k", lw=1, ls="--")
        ax_res.scatter(x, res, s=s, alpha=alpha)

        ax_top.set_ylabel(ylabel)
        ax_res.set_xlabel(xlabel)
        ax_res.set_ylabel(res_label)

        if title is not None:
            ax_top.set_title(title)

        if label is not None:
            ax_top.legend()

        if len(res) > 0:
            med = np.nanmedian(res)
            mad = np.nanmedian(np.abs(res - med))
            std = np.nanstd(res)

            ax_res.text(
                0.04, 0.96,
                (
                    f"$\\mathrm{{med}}={med:.3g}$\n"
                    f"$\\mathrm{{MAD}}={mad:.3g}$\n"
                    f"$\\sigma={std:.3g}$\n"
                    f"$N={len(res)}$"
                ),
                transform=ax_res.transAxes,
                ha="left",
                va="top",
                fontsize=9,
            )

    ncols = 2 if (x2 is not None and y2 is not None) else 1

    if figsize is None:
        figsize = (12, 6) if ncols == 2 else (6, 6)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        2, ncols,
        height_ratios=[3, 1],
        hspace=0.05,
        wspace=0.3,
    )

    ax1_top = fig.add_subplot(gs[0, 0])
    ax1_res = fig.add_subplot(gs[1, 0], sharex=ax1_top)

    draw_column(
        ax1_top, ax1_res,
        x1, y1,
        xlabel=xlabel1,
        ylabel=ylabel1,
        title=title1,
        label=label1,
    )

    plt.setp(ax1_top.get_xticklabels(), visible=False)

    if ncols == 2:
        ax2_top = fig.add_subplot(gs[0, 1])
        ax2_res = fig.add_subplot(gs[1, 1], sharex=ax2_top)

        draw_column(
            ax2_top, ax2_res,
            x2, y2,
            xlabel=xlabel2 or xlabel1,
            ylabel=ylabel2 or ylabel1,
            title=title2,
            label=label2,
        )

        plt.setp(ax2_top.get_xticklabels(), visible=False)

    if outfile is not None:
        fig.savefig(outfile, dpi=150, bbox_inches="tight")

    return fig

def maybe_log(arr, do_log=False):
    arr = np.asarray(arr, dtype=float)
    if not do_log:
        return arr, np.isfinite(arr)

    out = np.full(len(arr), np.nan)
    good = np.isfinite(arr) & (arr > 0)
    out[good] = np.log10(arr[good])
    return out, good

def plot_with_residuals(
    x1, y1,
    x2=None, y2=None,
    xlabel1="x",
    ylabel1="y",
    xlabel2=None,
    ylabel2=None,
    title1=None,
    title2=None,
    residual="diff",      # "diff", "frac", "logratio"
    logx=False,
    logy=False,
    one_to_one=True,
    figsize=None,
    s=20,
    alpha=0.7,
    outfile=None,
    pids=None,
):
    """
    Plot y vs x in a top panel and residuals in a smaller bottom panel.

    If x2/y2 are supplied, make a 2-column version.

    Notes
    -----
    If logx=True and logy=True with residual="diff", the residual is:
        log10(y_original) - log10(x_original) = log10(y_original / x_original)
    """
    import numpy as np
    import matplotlib.pyplot as plt

    def clean_xy(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        good = np.isfinite(x) & np.isfinite(y)

        if logx:
            good &= x > 0
        if logy:
            good &= y > 0

        x = x[good]
        y = y[good]

        if logx:
            x = np.log10(x)
        if logy:
            y = np.log10(y)

        return x, y

    def get_residual(x, y):
        if residual == "diff":
            return y - x, r"$y - x$"
        elif residual == "frac":
            good = x != 0
            out = np.full(len(x), np.nan)
            out[good] = (y[good] - x[good]) / x[good]
            return out, r"$(y - x)/x$"
        elif residual == "logratio":
            good = (x > 0) & (y > 0)
            out = np.full(len(x), np.nan)
            out[good] = np.log10(y[good] / x[good])
            return out, r"$\log_{10}(y/x)$"
        else:
            raise ValueError("residual must be 'diff', 'frac', or 'logratio'")

    def draw_panel(ax_top, ax_res, x, y, xlabel, ylabel, title=None):
        x, y = clean_xy(x, y)
        res, res_label = get_residual(x, y)

        good_res = np.isfinite(res)

        ax_top.scatter(x, y, s=s, alpha=alpha)

        if one_to_one and len(x) > 0:
            lo = np.nanmin([np.nanmin(x), np.nanmin(y)])
            hi = np.nanmax([np.nanmax(x), np.nanmax(y)])
            if np.isfinite(lo) and np.isfinite(hi):
                pad = 0.03 * (hi - lo) if hi > lo else 0.1
                ax_top.plot([lo, hi], [lo, hi], "k--", lw=1)
                ax_top.set_xlim(lo - pad, hi + pad)
                ax_top.set_ylim(lo - pad, hi + pad)

        ax_top.set_ylabel(ylabel)
        if title is not None:
            ax_top.set_title(title)

        ax_res.axhline(0, color="k", lw=1, ls="--")
        ax_res.scatter(x[good_res], res[good_res], s=s, alpha=alpha)

        ax_res.set_xlabel(xlabel)
        ax_res.set_ylabel(res_label)

        if np.any(good_res):
            med = np.nanmedian(res[good_res])
            mad = np.nanmedian(np.abs(res[good_res] - med))
            std = np.nanstd(res[good_res])


            nmad = 1.4826 * mad

            ax_res.text(
                0.04, 0.96,
                (
                    f"$\\mathrm{{med}}={med:.3g}$\n"
                    f"$\\mathrm{{MAD}}={mad:.3g}$\n"
                    f"$\\mathrm{{NMAD}}={nmad:.3g}$\n"
                    f"$\\sigma={std:.3g}$\n"
                    f"$N={np.sum(good_res)}$"
                ),
                transform=ax_res.transAxes,
                ha="left",
                va="top",
                fontsize=9,
            )

            # print outliers
            #outlier = np.abs(residual) > 5 * std
            #print(np.arange(len(

        plt.setp(ax_top.get_xticklabels(), visible=False)

    two_panel = (x2 is not None) and (y2 is not None)
    ncols = 2 if two_panel else 1

    if figsize is None:
        figsize = (12, 6) if two_panel else (6, 6)

    fig, axes = plt.subplots(
        2,
        ncols,
        figsize=figsize,
        gridspec_kw={"height_ratios": [3, 1]},
        sharex="col",
        squeeze=False,
    )

    draw_panel(
        axes[0, 0],
        axes[1, 0],
        x1,
        y1,
        xlabel=xlabel1,
        ylabel=ylabel1,
        title=title1,
    )

    if two_panel:
        draw_panel(
            axes[0, 1],
            axes[1, 1],
            x2,
            y2,
            xlabel=xlabel2 or xlabel1,
            ylabel=ylabel2 or ylabel1,
            title=title2,
        )

    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile, dpi=150, bbox_inches="tight")

    return fig, axes
