#!/usr/bin/env python

"""
qc_duplicates.py

Quality-control analysis for duplicate HAPY observations.

This script:
- reads a merged HAPY results table
- identifies galaxies with duplicate observations
- builds all pairwise duplicate comparisons
- separates same-telescope and cross-telescope pairs
- applies parameter-family-specific quality cuts
- writes summary tables
- generates pair and residual plots

Example
-------
python qc_duplicates.py merged_results.fits --outdir qc_duplicates

Notes
-----
- R-band duplicate analyses do NOT exclude rows with large Halpha filter correction.
- Halpha duplicate analyses DO exclude rows with FILTER_CORRECTION > threshold
  via FILTER_WARNING.
- GALFIT NC and CV are analyzed separately.
- This script is robust to the current typo bug where FWHM values may be stored in
  R_FHWM / H_FHWM instead of R_FWHM / H_FWHM.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
from hapy.utils.plotting import raincloud_by_group


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def safe_bool_array(tab: Table, colname: str, default: bool = False) -> np.ndarray:
    """Return a boolean numpy array for a table column, or default if missing."""
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


def safe_float_array(tab: Table, colname: str, default=np.nan) -> np.ndarray:
    """Return a float numpy array for a table column, or default if missing."""
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


def first_existing_col(tab: Table, names: list[str]) -> str | None:
    """Return the first existing column from a list of candidate names."""
    for name in names:
        if name in tab.colnames:
            return name
    return None


def first_populated_col(tab: Table, names: list[str]) -> str | None:
    """
    Return the first existing column that has at least one finite value.
    Useful for handling old typo columns like R_FHWM / H_FHWM.
    """
    for name in names:
        if name not in tab.colnames:
            continue
        vals = safe_float_array(tab, name)
        if np.any(np.isfinite(vals)):
            return name
    return None


def finite_pair_mask(x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
    return np.isfinite(x1) & np.isfinite(x2)


def median_and_mad(dx: np.ndarray) -> tuple[float, float]:
    """Return median and MAD-like robust scatter."""
    good = np.isfinite(dx)
    if np.sum(good) == 0:
        return np.nan, np.nan
    med = np.nanmedian(dx[good])
    mad = np.nanmedian(np.abs(dx[good] - med))
    return med, mad


# ----------------------------------------------------------------------
# duplicate pair construction
# ----------------------------------------------------------------------

def build_duplicate_pairs(tab: Table, id_col: str = "VFID") -> Table:
    """
    Build all pairwise duplicate observation pairs for rows sharing the same VFID.

    Returns a table with columns:
      idx1, idx2, VFID, TEL1, TEL2, same_tel
    """
    if id_col not in tab.colnames:
        raise ValueError(f"Missing required column '{id_col}'")

    tel_col = first_existing_col(tab, ["TELESCOPE", "TEL"])

    groups = {}
    for i, vfid in enumerate(tab[id_col]):
        groups.setdefault(str(vfid), []).append(i)

    rows = []
    for vfid, inds in groups.items():
        if len(inds) < 2:
            continue
        for i1, i2 in itertools.combinations(inds, 2):
            if tel_col is not None:
                tel1 = str(tab[tel_col][i1])
                tel2 = str(tab[tel_col][i2])
                same_tel = (tel1 == tel2)
            else:
                tel1 = ""
                tel2 = ""
                same_tel = False

            rows.append((i1, i2, vfid, tel1, tel2, same_tel))

    if len(rows) == 0:
        return Table(
            names=("idx1", "idx2", "VFID", "TEL1", "TEL2", "same_tel"),
            dtype=("i4", "i4", "U32", "U16", "U16", "bool"),
        )

    return Table(
        rows=rows,
        names=("idx1", "idx2", "VFID", "TEL1", "TEL2", "same_tel"),
    )


# ----------------------------------------------------------------------
# row and pair masks
# ----------------------------------------------------------------------

def build_row_flags(tab: Table, max_ha_filter_correction: float = 1.2) -> dict[str, np.ndarray]:
    """
    Construct row-level QC flags.

    R-band analyses remain usable even when Halpha filter correction is bad.
    Halpha analyses are more restrictive.
    """
    flags = {}

    flags["MASK_OK"] = safe_bool_array(tab, "MASK_OK")
    flags["PHOT_OK"] = safe_bool_array(tab, "PHOT_OK")
    flags["PSF_OK"] = safe_bool_array(tab, "PSF_OK")
    flags["R_PROFILE_OK"] = safe_bool_array(tab, "R_PROFILE_OK")
    flags["HA_PROFILE_OK"] = safe_bool_array(tab, "HA_PROFILE_OK")
    flags["R_SM_FLAG"] = safe_bool_array(tab, "R_SM_FLAG")
    flags["H_SM_FLAG"] = safe_bool_array(tab, "H_SM_FLAG")
    flags["GAL_NC_OK"] = safe_bool_array(tab, "GAL_NC_OK")
    flags["GAL_CV_OK"] = safe_bool_array(tab, "GAL_CV_OK")

    filt_col = first_existing_col(tab, ["FILTER_CORRECTION", "FILT_COR"])
    if filt_col is not None:
        filtcor = safe_float_array(tab, filt_col)
    else:
        filtcor = np.full(len(tab), np.nan)

    flags["FILTER_CORR"] = filtcor
    flags["FILTER_WARNING"] = np.isfinite(filtcor) & (filtcor > max_ha_filter_correction)

    # broader row families
    flags["R_DUP_OK"] = flags["MASK_OK"] & flags["PHOT_OK"]
    flags["HA_DUP_OK"] = flags["MASK_OK"] & flags["PHOT_OK"] & (~flags["FILTER_WARNING"])

    flags["R_SM_OK"] = flags["R_SM_FLAG"]
    flags["H_SM_OK"] = flags["H_SM_FLAG"] & (~flags["FILTER_WARNING"])

    flags["GALFIT_NC_OK"] = flags["GAL_NC_OK"]
    flags["GALFIT_CV_OK"] = flags["GAL_CV_OK"]
    flags["GALFIT_ANY_OK"] = flags["GAL_NC_OK"] | flags["GAL_CV_OK"]

    return flags


def pair_mask_from_row_flag(pairtab: Table, rowflag: np.ndarray) -> np.ndarray:
    """Require both rows in a duplicate pair to pass a row-level flag."""
    return rowflag[pairtab["idx1"]] & rowflag[pairtab["idx2"]]


# ----------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------

def summarize_pair_subset(name: str, pairmask: np.ndarray, pairtab: Table) -> dict:
    n_pairs = int(np.sum(pairmask))
    same_tel = int(np.sum(pairmask & pairtab["same_tel"]))
    cross_tel = int(np.sum(pairmask & (~pairtab["same_tel"])))
    return dict(
        subset=name,
        n_pairs=n_pairs,
        n_same_tel=same_tel,
        n_cross_tel=cross_tel,
    )


def write_subset_summary(outdir: Path, pairtab: Table, pairmasks: dict[str, np.ndarray]) -> None:
    rows = [summarize_pair_subset(name, mask, pairtab) for name, mask in pairmasks.items()]
    Table(rows=rows).write(outdir / "duplicate_pair_summary.ecsv",
                           format="ascii.ecsv", overwrite=True)


def write_duplicate_pairs(outdir: Path, pairtab: Table) -> None:
    pairtab.write(outdir / "duplicate_pairs.ecsv", format="ascii.ecsv", overwrite=True)


# ----------------------------------------------------------------------
# plotting
# ----------------------------------------------------------------------

def plot_pair_grid(
    tab: Table,
    pairtab: Table,
    pairmask: np.ndarray,
    cols: list[str],
    outpath: Path,
    color_by: str | None = None,
    title: str = "",
    residual: bool = False,
) -> None:
    """
    Make pair plots or residual plots for a family of columns.
    """
    use = np.where(pairmask)[0]
    if len(use) == 0:
        print(f"WARNING: no pairs for plot {outpath.name}")
        return

    ncol = len(cols)
    ncols = 4
    nrows = int(np.ceil(ncol / ncols))

    fig = plt.figure(figsize=(4 * ncols, 3.5 * nrows))
    plt.subplots_adjust(hspace=0.45, wspace=0.35)

    cb_handle = None

    if color_by is not None and color_by in tab.colnames:
        cval = safe_float_array(tab, color_by)
        cpair = 0.5 * (cval[pairtab["idx1"][use]] + cval[pairtab["idx2"][use]])
    else:
        cpair = None

    for i, col in enumerate(cols, start=1):
        ax = plt.subplot(nrows, ncols, i)

        if col not in tab.colnames:
            ax.text(0.5, 0.5, f"Missing:\n{col}", ha="center", va="center")
            ax.set_axis_off()
            continue

        x1 = safe_float_array(tab, col)[pairtab["idx1"][use]]
        x2 = safe_float_array(tab, col)[pairtab["idx2"][use]]
        good = finite_pair_mask(x1, x2)

        if np.sum(good) == 0:
            ax.text(0.5, 0.5, f"No finite pairs:\n{col}", ha="center", va="center")
            ax.set_axis_off()
            continue

        x1 = x1[good]
        x2 = x2[good]

        cc = cpair[good] if cpair is not None else None

        if residual:
            dx = x2 - x1
            if cc is not None:
                cb_handle = ax.scatter(x1, dx, c=cc, s=18, alpha=0.75)
            else:
                ax.scatter(x1, dx, s=18, alpha=0.75)

            med, mad = median_and_mad(dx)
            ax.axhline(0, color="k", ls="--", lw=1)
            if np.isfinite(mad):
                ax.axhline(+mad, color="k", ls=":", lw=1)
                ax.axhline(-mad, color="k", ls=":", lw=1)

            ax.text(
                0.04, 0.96,
                f"med={med:.3g}\nmad={mad:.3g}\nN={len(dx)}",
                transform=ax.transAxes,
                ha="left", va="top", fontsize=9,
            )
            ax.set_xlabel(col)
            ax.set_ylabel("obs2 - obs1")

        else:
            if cc is not None:
                cb_handle = ax.scatter(x1, x2, c=cc, s=18, alpha=0.75)
            else:
                ax.scatter(x1, x2, s=18, alpha=0.75)

            lo = np.nanmin([np.nanmin(x1), np.nanmin(x2)])
            hi = np.nanmax([np.nanmax(x1), np.nanmax(x2)])
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo, hi], [lo, hi], "k--", lw=1)

            dx = x2 - x1
            med, mad = median_and_mad(dx)
            ax.text(
                0.04, 0.96,
                f"med={med:.3g}\nmad={mad:.3g}\nN={len(dx)}",
                transform=ax.transAxes,
                ha="left", va="top", fontsize=9,
            )
            ax.set_xlabel("obs1")
            ax.set_ylabel("obs2")

        ax.set_title(col)

    if title:
        fig.suptitle(title, fontsize=14)

    if cb_handle is not None and color_by is not None:
        cbar = fig.colorbar(cb_handle, ax=fig.axes, fraction=0.025, pad=0.02)
        cbar.set_label(color_by)

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_subset_histogram(pairtab: Table, pairmasks: dict[str, np.ndarray], outpath: Path) -> None:
    names = list(pairmasks.keys())
    counts = [int(np.sum(pairmasks[k])) for k in names]

    fig = plt.figure(figsize=(10, 5))
    plt.bar(range(len(names)), counts)
    plt.xticks(range(len(names)), names, rotation=45, ha="right")
    plt.ylabel("Number of duplicate pairs")
    plt.title("Duplicate pair counts by subset")
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
# parameter summaries
# ----------------------------------------------------------------------

def summarize_parameters(
    tab: Table,
    pairtab: Table,
    pairmask: np.ndarray,
    cols: list[str],
    subset_name: str,
) -> Table:
    rows = []
    use = np.where(pairmask)[0]

    for col in cols:
        if col not in tab.colnames:
            rows.append((subset_name, col, 0, np.nan, np.nan))
            continue

        x1 = safe_float_array(tab, col)[pairtab["idx1"][use]]
        x2 = safe_float_array(tab, col)[pairtab["idx2"][use]]
        good = finite_pair_mask(x1, x2)

        if np.sum(good) == 0:
            rows.append((subset_name, col, 0, np.nan, np.nan))
            continue

        dx = x2[good] - x1[good]
        med, mad = median_and_mad(dx)
        rows.append((subset_name, col, int(np.sum(good)), med, mad))

    return Table(
        rows=rows,
        names=("subset", "column", "n_pairs", "median_offset", "mad_offset"),
    )

def plot_telescope_pair_residuals(
    tab,
    pairtab,
    pairmask,
    col,
    outpath,
    ylabel=None,
):
    """
    Plot duplicate residuals grouped by telescope pair.

    Only cross-telescope pairs are used.
    Residual is defined as obs2 - obs1.
    """

    if col not in tab.colnames:
        print(f"WARNING: missing column {col}, skipping {outpath}")
        return

    use = pairmask & (~pairtab["same_tel"])
    if np.sum(use) == 0:
        print(f"WARNING: no cross-telescope pairs for {col}")
        return

    idx1 = pairtab["idx1"][use]
    idx2 = pairtab["idx2"][use]

    x1 = safe_float_array(tab, col)[idx1]
    x2 = safe_float_array(tab, col)[idx2]
    tel1 = np.array(pairtab["TEL1"][use]).astype(str)
    tel2 = np.array(pairtab["TEL2"][use]).astype(str)

    good = np.isfinite(x1) & np.isfinite(x2)
    if np.sum(good) == 0:
        print(f"WARNING: no finite cross-telescope pairs for {col}")
        return

    x1 = x1[good]
    x2 = x2[good]
    tel1 = tel1[good]
    tel2 = tel2[good]

    residual = x2 - x1
    pair_names = np.array(["-".join(sorted([a, b])) for a, b in zip(tel1, tel2)])

    unique_pairs = sorted(set(pair_names))
    xpos = {p: i for i, p in enumerate(unique_pairs)}

    fig = plt.figure(figsize=(1.6 * max(len(unique_pairs), 4), 6))
    ax = plt.gca()

    for p in unique_pairs:
        m = pair_names == p
        y = residual[m]
        x = np.full(np.sum(m), xpos[p], dtype=float)

        # small horizontal jitter
        jitter = np.random.uniform(-0.12, 0.12, size=len(y))
        ax.scatter(x + jitter, y, alpha=0.7, s=24)

        med = np.nanmedian(y)
        mad = np.nanmedian(np.abs(y - med))

        ax.plot([xpos[p] - 0.2, xpos[p] + 0.2], [med, med], "k-", lw=2)
        ax.plot([xpos[p], xpos[p]], [med - mad, med + mad], "k-", lw=2)

        ax.text(
            xpos[p], med,
            f"  N={len(y)}",
            fontsize=8,
            va="bottom",
            ha="left",
        )

    ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_xticks(range(len(unique_pairs)))
    ax.set_xticklabels(unique_pairs, rotation=45, ha="right")
    ax.set_ylabel(ylabel if ylabel is not None else f"{col}: obs2 - obs1")
    ax.set_title(f"Cross-telescope duplicate residuals: {col}")
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_telescope_group_residuals(
    tab,
    pairtab,
    pairmask,
    col,
    outpath,
    mode="cross",
    ylabel=None,
    logdiff=False,
):
    """
    Plot duplicate residuals grouped by telescope pairing.

    Parameters
    ----------
    mode : {"cross", "same"}
        cross -> use only cross-telescope pairs, grouped as BOK-HDI, etc.
        same  -> use only same-telescope pairs, grouped as BOK, HDI, etc.

    logdiff : bool
        If True, plot log10(obs2) - log10(obs1) instead of obs2 - obs1.
        Useful for strictly positive flux-like quantities.
    """

    if col not in tab.colnames:
        print(f"WARNING: missing column {col}, skipping {outpath}")
        return

    if mode == "cross":
        use = pairmask & (~pairtab["same_tel"])
    elif mode == "same":
        use = pairmask & pairtab["same_tel"]
    else:
        raise ValueError("mode must be 'cross' or 'same'")

    if np.sum(use) == 0:
        print(f"WARNING: no {mode}-telescope pairs for {col}")
        return

    idx1 = pairtab["idx1"][use]
    idx2 = pairtab["idx2"][use]

    x1 = safe_float_array(tab, col)[idx1]
    x2 = safe_float_array(tab, col)[idx2]
    tel1 = np.array(pairtab["TEL1"][use]).astype(str)
    tel2 = np.array(pairtab["TEL2"][use]).astype(str)

    if logdiff:
        good = np.isfinite(x1) & np.isfinite(x2) & (x1 > 0) & (x2 > 0)
    else:
        good = np.isfinite(x1) & np.isfinite(x2)

    if np.sum(good) == 0:
        print(f"WARNING: no finite {mode}-telescope pairs for {col}")
        return

    x1 = x1[good]
    x2 = x2[good]
    tel1 = tel1[good]
    tel2 = tel2[good]

    if logdiff:
        residual = np.log10(x2) - np.log10(x1)
    else:
        residual = x2 - x1

    if mode == "cross":
        groups = np.array(["-".join(sorted([a, b])) for a, b in zip(tel1, tel2)])
        title = f"Cross-telescope duplicate residuals: {col}"
    else:
        groups = tel1
        title = f"Same-telescope duplicate residuals: {col}"

    unique_groups = sorted(set(groups))
    xpos = {g: i for i, g in enumerate(unique_groups)}

    fig = plt.figure(figsize=(1.6 * max(len(unique_groups), 4), 6))
    ax = plt.gca()

    rng = np.random.default_rng(12345)

    for g in unique_groups:
        m = groups == g
        y = residual[m]
        x = np.full(np.sum(m), xpos[g], dtype=float)

        jitter = rng.uniform(-0.12, 0.12, size=len(y))
        ax.scatter(x + jitter, y, alpha=0.7, s=24)

        med = np.nanmedian(y)
        mad = np.nanmedian(np.abs(y - med))

        ax.plot([xpos[g] - 0.2, xpos[g] + 0.2], [med, med], "k-", lw=2)
        ax.plot([xpos[g], xpos[g]], [med - mad, med + mad], "k-", lw=2)

        ax.text(
            xpos[g], med,
            f"  N={len(y)}",
            fontsize=8,
            va="bottom",
            ha="left",
        )

    ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_xticks(range(len(unique_groups)))
    ax.set_xticklabels(unique_groups, rotation=45, ha="right")

    if ylabel is None:
        if logdiff:
            ylabel = f"log10(obs2) - log10(obs1): {col}"
        else:
            ylabel = f"obs2 - obs1: {col}"

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_rainclouds(unique_pairs, pair_names, residual):
    values_by_group = []
    labels = []

    for pair_name in unique_pairs:
        mask = pair_names == pair_name
        values_by_group.append(residual[mask])
        labels.append(pair_name)

    fig, ax = raincloud_by_group(
        values_by_group,
        labels,
        xlabel="obs2 - obs1",
        title="Cross-telescope residuals: R24_MAG"
        )
# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="QC analysis for duplicate HAPY observations.")
    parser.add_argument("table", help="Merged HAPY results table (e.g. merged_results.fits)")
    parser.add_argument("--outdir", default="qc_duplicates", help="Output directory")
    parser.add_argument(
        "--max-ha-filter-correction",
        type=float,
        default=1.2,
        help="Maximum Halpha filter correction for 'good' Halpha duplicate comparisons",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    tab = Table.read(args.table)
    print(f"Read {len(tab)} rows from {args.table}")

    if "VFID" not in tab.colnames:
        raise RuntimeError("Merged results table must contain VFID for duplicate analysis.")

    pairtab = build_duplicate_pairs(tab, id_col="VFID")
    print(f"Found {len(pairtab)} duplicate pairs")

    write_duplicate_pairs(outdir, pairtab)

    flags = build_row_flags(tab, max_ha_filter_correction=args.max_ha_filter_correction)

    # resolve FWHM columns, including typo fallback
    r_fwhm_col = first_populated_col(tab, ["R_FWHM", "R_FHWM"])
    h_fwhm_col = first_populated_col(tab, ["H_FWHM", "H_FHWM"])

    if r_fwhm_col is not None:
        print(f"Using {r_fwhm_col} for R-band FWHM coloring")
    if h_fwhm_col is not None:
        print(f"Using {h_fwhm_col} for Halpha FWHM coloring")

    # pair-level masks
    all_pairs = np.ones(len(pairtab), dtype=bool)

    pairmasks = {
        "all_pairs": all_pairs,
        "same_tel_pairs": pairtab["same_tel"],
        "cross_tel_pairs": ~pairtab["same_tel"],
        "r_pairs": pair_mask_from_row_flag(pairtab, flags["R_DUP_OK"]),
        "ha_pairs": pair_mask_from_row_flag(pairtab, flags["HA_DUP_OK"]),
        "r_sm_pairs": pair_mask_from_row_flag(pairtab, flags["R_SM_OK"]),
        "h_sm_pairs": pair_mask_from_row_flag(pairtab, flags["H_SM_OK"]),
        "galfit_nc_pairs": pair_mask_from_row_flag(pairtab, flags["GALFIT_NC_OK"]),
        "galfit_cv_pairs": pair_mask_from_row_flag(pairtab, flags["GALFIT_CV_OK"]),
        "galfit_any_pairs": pair_mask_from_row_flag(pairtab, flags["GALFIT_ANY_OK"]),
    }

    # Halpha warning subset
    filter_warning = flags["FILTER_WARNING"]
    pairmasks["ha_warning_pairs"] = (
        filter_warning[pairtab["idx1"]] | filter_warning[pairtab["idx2"]]
    )

    write_subset_summary(outdir, pairtab, pairmasks)
    plot_subset_histogram(pairtab, pairmasks, outdir / "duplicate_pair_counts.png")

    # parameter families using real schema
    r_cols = [
        "R24_MAG", "R25_ISO_MAG", "R25P5_MAG", "R_PETRO_MAG",
        "R24_ARCSEC", "R25_ISO_ARCSEC", "R25P5_ARCSEC",
        "R25_ARCSEC", "R50_ARCSEC", "R75_ARCSEC",
        "R24_FLUX_CGS", "R_PETRO_FLUX_CGS",
        "R_C30", "R_PETRO_CON",
        "ELLIP_GINI_DET", "R_M20", "R_ASYM",
    ]

    ha_cols = [
        "HA_TOT_FLUX_CGS", "HA_R24_FLUX_CGS", "HA30R24_FLUX_CGS",
        "HA_ISO5E17_FLUX_CGS", "HA_ISO17E18_FLUX_CGS",
        "HA_MAXDET_ARCSEC", "HA25_ARCSEC", "HA50_ARCSEC", "HA75_ARCSEC",
        "HA_ISO5E17_ARCSEC", "HA_ISO17E18_ARCSEC",
        "HA_C30_R24", "HA_PETRO_CON",
        "H_M20", "H_ASYM",
    ]

    r_sm_cols = [
        "R_SM_GINI", "R_SM_M20", "R_SM_C", "R_SM_A", "R_SM_S",
        "R_SM_RPETRO_ELLIP", "R_SM_RHALF_ELLIP",
        "R_SM_R20", "R_SM_R50", "R_SM_R80",
    ]

    h_sm_cols = [
        "H_SM_GINI", "H_SM_M20", "H_SM_C", "H_SM_A", "H_SM_S",
        "H_SM_RPETRO_ELLIP", "H_SM_RHALF_ELLIP",
        "H_SM_R20", "H_SM_R50", "H_SM_R80",
    ]

    galfit_nc_cols = [
        "GAL_MAG", "GAL_RE", "GAL_N", "GAL_BA", "GAL_PA", "GAL_CHISQ",
    ]

    galfit_cv_cols = [
        "GAL_CMAG", "GAL_CRE", "GAL_CN", "GAL_CBA", "GAL_CPA", "GAL_CCHISQ",
    ]

    # plots
    plot_pair_grid(
        tab, pairtab, pairmasks["r_pairs"], r_cols,
        outdir / "r_pairs.png",
        color_by=r_fwhm_col,
        title="R-band duplicate comparisons",
        residual=False,
    )
    plot_pair_grid(
        tab, pairtab, pairmasks["r_pairs"], r_cols,
        outdir / "r_residuals.png",
        color_by=r_fwhm_col,
        title="R-band duplicate residuals",
        residual=True,
    )

    plot_pair_grid(
        tab, pairtab, pairmasks["ha_pairs"], ha_cols,
        outdir / "ha_pairs.png",
        color_by=h_fwhm_col,
        title="Halpha duplicate comparisons (clean sample)",
        residual=False,
    )
    plot_pair_grid(
        tab, pairtab, pairmasks["ha_pairs"], ha_cols,
        outdir / "ha_residuals.png",
        color_by=h_fwhm_col,
        title="Halpha duplicate residuals (clean sample)",
        residual=True,
    )

    plot_pair_grid(
        tab, pairtab, pairmasks["r_sm_pairs"], r_sm_cols,
        outdir / "r_statmorph_pairs.png",
        color_by=r_fwhm_col,
        title="R-band statmorph duplicate comparisons",
        residual=False,
    )
    plot_pair_grid(
        tab, pairtab, pairmasks["r_sm_pairs"], r_sm_cols,
        outdir / "r_statmorph_residuals.png",
        color_by=r_fwhm_col,
        title="R-band statmorph duplicate residuals",
        residual=True,
    )

    plot_pair_grid(
        tab, pairtab, pairmasks["h_sm_pairs"], h_sm_cols,
        outdir / "ha_statmorph_pairs.png",
        color_by=h_fwhm_col,
        title="Halpha statmorph duplicate comparisons",
        residual=False,
    )
    plot_pair_grid(
        tab, pairtab, pairmasks["h_sm_pairs"], h_sm_cols,
        outdir / "ha_statmorph_residuals.png",
        color_by=h_fwhm_col,
        title="Halpha statmorph duplicate residuals",
        residual=True,
    )

    plot_pair_grid(
        tab, pairtab, pairmasks["galfit_nc_pairs"], galfit_nc_cols,
        outdir / "galfit_nc_pairs.png",
        color_by=r_fwhm_col,
        title="GALFIT NC duplicate comparisons",
        residual=False,
    )
    plot_pair_grid(
        tab, pairtab, pairmasks["galfit_nc_pairs"], galfit_nc_cols,
        outdir / "galfit_nc_residuals.png",
        color_by=r_fwhm_col,
        title="GALFIT NC duplicate residuals",
        residual=True,
    )

    plot_pair_grid(
        tab, pairtab, pairmasks["galfit_cv_pairs"], galfit_cv_cols,
        outdir / "galfit_cv_pairs.png",
        color_by=r_fwhm_col,
        title="GALFIT CV duplicate comparisons",
        residual=False,
    )
    plot_pair_grid(
        tab, pairtab, pairmasks["galfit_cv_pairs"], galfit_cv_cols,
        outdir / "galfit_cv_residuals.png",
        color_by=r_fwhm_col,
        title="GALFIT CV duplicate residuals",
        residual=True,
    )

    plot_telescope_pair_residuals(
        tab, pairtab, pairmasks["r_pairs"],
        "R24_MAG",
        outdir / "telpair_R24_MAG.png",
        ylabel="Δ R24_MAG (obs2 - obs1)"
        )

    plot_telescope_pair_residuals(
        tab, pairtab, pairmasks["r_pairs"],
        "GAL_MAG",
        outdir / "telpair_GAL_MAG.png",
        ylabel="Δ GAL_MAG (obs2 - obs1)"
        )

    plot_telescope_pair_residuals(
        tab, pairtab, pairmasks["galfit_nc_pairs"],
        "GAL_RE",
        outdir / "telpair_GAL_RE.png",
        ylabel="Δ GAL_RE (obs2 - obs1)"
        )
    plot_telescope_pair_residuals(
        tab, pairtab, pairmasks["ha_pairs"],
        "HA_TOT_FLUX_CGS",
        outdir / "telpair_HA_TOT_FLUX_CGS.png",
        ylabel="Δ HA_TOT_FLUX_CGS (obs2 - obs1)"
        )

    plot_telescope_pair_residuals(
        tab, pairtab, pairmasks["ha_pairs"],
        "HA_R24_FLUX_CGS",
        outdir / "telpair_HA_R24_FLUX_CGS.png",
        ylabel="Δ HA_R24_FLUX_CGS (obs2 - obs1)"
        )

    plot_telescope_group_residuals(
        tab, pairtab, pairmasks["r_pairs"],
        "R24_MAG",
        outdir / "telpair_cross_R24_MAG.png",
        mode="cross",
        ylabel="Δ R24_MAG (obs2 - obs1)"
        )

    plot_telescope_group_residuals(
        tab, pairtab, pairmasks["galfit_nc_pairs"],
        "GAL_RE",
        outdir / "telpair_cross_GAL_RE.png",
        mode="cross",
        ylabel="Δ GAL_RE (obs2 - obs1)"
        )

    plot_telescope_group_residuals(
        tab, pairtab, pairmasks["ha_pairs"],
        "HA_TOT_FLUX_CGS",
        outdir / "telpair_cross_HA_TOT_FLUX_CGS.png",
        mode="cross",
        logdiff=True,
        ylabel="Δ log10(HA_TOT_FLUX_CGS)"
        )

    plot_telescope_group_residuals(
        tab, pairtab, pairmasks["r_pairs"],
        "R24_MAG",
        outdir / "telpair_same_R24_MAG.png",
        mode="same",
        ylabel="Δ R24_MAG (obs2 - obs1)"
        )

    plot_telescope_group_residuals(
        tab, pairtab, pairmasks["galfit_nc_pairs"],
        "GAL_RE",
        outdir / "telpair_same_GAL_RE.png",
        mode="same",
        ylabel="Δ GAL_RE (obs2 - obs1)"
        )

    plot_telescope_group_residuals(
        tab, pairtab, pairmasks["ha_pairs"],
        "HA_TOT_FLUX_CGS",
        outdir / "telpair_same_HA_TOT_FLUX_CGS.png",
        mode="same",
        logdiff=True,
        ylabel="Δ log10(HA_TOT_FLUX_CGS)"
        )
    # summary tables
    summary_tabs = [
        summarize_parameters(tab, pairtab, pairmasks["r_pairs"], r_cols, "r_pairs"),
        summarize_parameters(tab, pairtab, pairmasks["ha_pairs"], ha_cols, "ha_pairs"),
        summarize_parameters(tab, pairtab, pairmasks["r_sm_pairs"], r_sm_cols, "r_sm_pairs"),
        summarize_parameters(tab, pairtab, pairmasks["h_sm_pairs"], h_sm_cols, "h_sm_pairs"),
        summarize_parameters(tab, pairtab, pairmasks["galfit_nc_pairs"], galfit_nc_cols, "galfit_nc_pairs"),
        summarize_parameters(tab, pairtab, pairmasks["galfit_cv_pairs"], galfit_cv_cols, "galfit_cv_pairs"),
    ]

    summary = vstack(summary_tabs, metadata_conflicts="silent")
    summary.write(outdir / "duplicate_parameter_summary.ecsv",
                  format="ascii.ecsv", overwrite=True)

    print(f"Wrote duplicate QC products to {outdir}")


if __name__ == "__main__":
    main()
