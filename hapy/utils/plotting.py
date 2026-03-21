import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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
    "R_SM_R50": r"$R_{50}^{R,\mathrm{SM}}$",
    "R_SM_R80": r"$R_{80}^{R,\mathrm{SM}}$",    
    "R_SM_RHALF_ELLIP": r"$R_{1/2}^{R,\mathrm{SM,ellip}}$",
    "R_SM_RHALF_CIRCLE": r"$R_{1/2}^{R,\mathrm{SM,circle}}$",
    "R_SM_RMAX_ELLIP": r"$R_{max}^{R,\mathrm{SM,ellip}}$",
    "R_SM_RMAX_CIRCLE": r"$R_{max}^{R,\mathrm{SM,circle}}$",

    # -------------------------
    # Hα sizes
    # -------------------------
    "H50_ARCSEC": r"$R_{50}^{\mathrm{H\alpha}}$ (arcsec)",
    "H25_ARCSEC": r"$R_{25}^{\mathrm{H\alpha}}$",
    "H75_ARCSEC": r"$R_{75}^{\mathrm{H\alpha}}$",
    "H_R95_R24_ARCSEC": r"$R_{R95,R24}^{\mathrm{H\alpha}}$",    
    "H_MAXDET_ARCSEC": r"$R_{\max}^{\mathrm{H\alpha}}$",
    "H_PETRO_R50_ARCSEC": r"$R_{50,\mathrm{Petro}}^{\mathrm{H\alpha}}$",
    "H_EXPFIT_RE_ARCSEC": r"$R_e^{\mathrm{H\alpha}}$ (exp fit)",
    "H_SM_R50": r"$R_{50}^{\mathrm{H\alpha},\mathrm{SM}}$", 
    "H_SM_R80": r"$R_{80}^{\mathrm{H\alpha},\mathrm{SM}}$",   
    "H_SM_RHALF_ELLIP": r"$R_{1/2}^{\mathrm{H\alpha},\mathrm{SM,ellip}}$",
    "H_SM_RMAX_ELLIP": r"$R_{max}^{\mathrm{H\alpha},\mathrm{SM,ellip}}$",    
    "H_SM_RMAX_CIRCLE": r"$R_{max}^{\mathrm{H\alpha},\mathrm{SM,circle}}$",
    # -------------------------
    # Morphology (HAPY + statmorph)
    # -------------------------
    "R_HAPY_GINI": r"$G^{R}$",
    "H_HAPY_GINI": r"$G^{\mathrm{H\alpha}}$",
    "R_HAPY_M20": r"$M_{20}^{R}$",
    "H_HAPY_M20": r"$M_{20}^{\mathrm{H\alpha}}$",

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
    "GAL_RE": r"$R_e^{\mathrm{GALFIT}}$",
    "GAL_N": r"$n$ (Sérsic)",
    "GAL_BA": r"$b/a$",
    "GAL_PA": r"$\mathrm{PA}$ (deg)",

    # -------------------------
    # QC / meta
    # -------------------------
    "QC_TIER": "QC Tier",
    "TELESCOPE": "Telescope",
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
