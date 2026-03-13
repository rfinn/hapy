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

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


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


def first_existing_col(tab: Table, names: list[str]) -> str | None:
    for name in names:
        if name in tab.colnames:
            return name
    return None


def first_populated_col(tab: Table, names: list[str]) -> str | None:
    for name in names:
        if name not in tab.colnames:
            continue
        vals = safe_float_array(tab, name)
        if np.any(np.isfinite(vals)):
            return name
    return None


# ----------------------------------------------------------------------
# flag discovery and masks
# ----------------------------------------------------------------------

def find_status_columns(tab: Table) -> list[str]:
    return sorted([
        c for c in tab.colnames
        if c.endswith("_OK") or c.endswith("_FLAG")
    ])


def build_qc_masks(tab: Table, max_ha_filter_correction: float = 1.2) -> dict[str, np.ndarray]:
    masks = {}

    # core flags
    masks["MASK_OK"] = safe_bool_array(tab, "MASK_OK")
    masks["PHOT_OK"] = safe_bool_array(tab, "PHOT_OK")
    masks["PSF_OK"] = safe_bool_array(tab, "PSF_OK")
    masks["R_PROFILE_OK"] = safe_bool_array(tab, "R_PROFILE_OK")
    masks["HA_PROFILE_OK"] = safe_bool_array(tab, "HA_PROFILE_OK")
    masks["R_SM_FLAG"] = safe_bool_array(tab, "R_SM_FLAG")
    masks["H_SM_FLAG"] = safe_bool_array(tab, "H_SM_FLAG")
    masks["GAL_NC_OK"] = safe_bool_array(tab, "GAL_NC_OK")
    masks["GAL_CV_OK"] = safe_bool_array(tab, "GAL_CV_OK")

    filt_col = first_existing_col(tab, ["FILTER_CORRECTION", "FILT_COR"])
    if filt_col is not None:
        filtcor = safe_float_array(tab, filt_col)
    else:
        filtcor = np.full(len(tab), np.nan)

    masks["FILTER_WARNING"] = np.isfinite(filtcor) & (filtcor > max_ha_filter_correction)

    # science-oriented subsets
    masks["mask_phot_ok"] = masks["MASK_OK"] & masks["PHOT_OK"]

    masks["profile_ok"] = (
        masks["MASK_OK"] &
        masks["PHOT_OK"] &
        masks["R_PROFILE_OK"] &
        masks["HA_PROFILE_OK"]
    )

    masks["statmorph_ok"] = (
        masks["profile_ok"] &
        masks["R_SM_FLAG"] &
        masks["H_SM_FLAG"] &
        (~masks["FILTER_WARNING"])
    )

    masks["galfit_any_ok"] = (
        masks["MASK_OK"] &
        masks["PHOT_OK"] &
        (masks["GAL_NC_OK"] | masks["GAL_CV_OK"])
    )

    masks["galfit_both_ok"] = (
        masks["MASK_OK"] &
        masks["PHOT_OK"] &
        masks["GAL_NC_OK"] &
        masks["GAL_CV_OK"]
    )

    masks["science_ready"] = (
        masks["profile_ok"] &
        masks["R_SM_FLAG"] &
        masks["H_SM_FLAG"] &
        (masks["GAL_NC_OK"] | masks["GAL_CV_OK"]) &
        (~masks["FILTER_WARNING"])
    )

    masks["problem"] = ~masks["mask_phot_ok"]

    return masks


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
        for f in find_status_columns(tab):
            col = safe_bool_array(tab, f)
            ntrue = np.sum(col)
            nfalse = np.sum(~col)
            pct = 100.0 * ntrue / n if n > 0 else np.nan
            fh.write(f"{f:20s}: {ntrue:4d} OK  | {nfalse:4d} FAIL  ({pct:5.1f}%)\n")

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

def plot_flag_completion(tab: Table, outpath: Path) -> None:
    flags = find_status_columns(tab)
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
    if col is None or "HA_TOT_FLUX_CGS" not in tab.colnames:
        return

    x = safe_float_array(tab, col)
    y = safe_float_array(tab, "HA_TOT_FLUX_CGS")

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
    ax.set_ylabel("HA_TOT_FLUX_CGS")
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
        "HA_PROFILE_OK",
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
        "ha profile fail": ~np.array(tab["HA_PROFILE_OK"][good], dtype=bool),
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

    masks = build_qc_masks(tab, max_ha_filter_correction=args.max_ha_filter_correction)

    # write text summary
    write_text_summary(tab, masks, outdir / "qc_summary.txt", args.scheme)

    # write subsets
    write_subsets(tab, masks, outdir)

    # plots
    plot_flag_completion(tab, outdir / "flag_completion.png")

    r_fwhm_col = first_populated_col(tab, ["R_FWHM", "R_FHWM"])
    h_fwhm_col = first_populated_col(tab, ["H_FWHM", "H_FHWM"])

    if r_fwhm_col is not None:
        plot_hist(tab, r_fwhm_col, outdir / "r_fwhm_hist.png", title="R-band FWHM")
    if h_fwhm_col is not None:
        plot_hist(tab, h_fwhm_col, outdir / "ha_fwhm_hist.png", title="Halpha FWHM")

    plot_hist(tab, "R24_MAG", outdir / "R24_MAG_hist.png", title="R24 magnitude", mask=masks["mask_phot_ok"])
    plot_hist(tab, "R25_ISO_MAG", outdir / "R25_ISO_MAG_hist.png", title="R25 isophotal magnitude", mask=masks["mask_phot_ok"])
    plot_hist(tab, "GAL_MAG", outdir / "GAL_MAG_hist.png", title="GALFIT magnitude", mask=masks["galfit_any_ok"])

    plot_hist(tab, "HA_TOT_FLUX_CGS", outdir / "HA_TOT_FLUX_CGS_hist.png",
              title="Halpha total flux", logx=True, mask=masks["mask_phot_ok"])
    plot_hist(tab, "HA_R24_FLUX_CGS", outdir / "HA_R24_FLUX_CGS_hist.png",
              title="Halpha R24 flux", logx=True, mask=masks["mask_phot_ok"])

    plot_hist(tab, "HA_MAXDET_ARCSEC", outdir / "HA_MAXDET_ARCSEC_hist.png",
              title="Halpha max detection radius", mask=masks["mask_phot_ok"])
    plot_hist(tab, "GAL_RE", outdir / "GAL_RE_hist.png",
              title="GALFIT effective radius", mask=masks["galfit_any_ok"])

    plot_compare(tab, "R24_MAG", "GAL_MAG", outdir / "R24_vs_GAL_MAG.png",
                 title="R24 magnitude vs GALFIT magnitude",
                 mask=masks["galfit_any_ok"])

    plot_compare(tab, "HA_R24_FLUX_CGS", "HA_TOT_FLUX_CGS", outdir / "HA_R24_vs_HA_TOT.png",
                 title="Halpha R24 flux vs total flux",
                 logx=True, logy=True, mask=masks["mask_phot_ok"])

    plot_compare(tab, "R24_ARCSEC", "GAL_RE", outdir / "R24_ARCSEC_vs_GAL_RE.png",
                 title="R24 radius vs GALFIT Re",
                 mask=masks["galfit_any_ok"])

    plot_failure_fraction_vs_bright_star_distance(tab, outdir / "FAILURES_VS_BRIGHT_STAR_DIST.png")
    plot_raincloud_by_telescope(tab, "R_FWHM", outdir / "raincloud_R_FWHM_by_telescope.png")
    plot_raincloud_by_telescope(tab, "H_FWHM", outdir / "raincloud_H_FWHM_by_telescope.png")
    plot_raincloud_by_telescope(tab, "FILTER_CORRECTION", outdir / "raincloud_FILTER_CORRECTION_by_telescope.png")
    plot_raincloud_by_telescope(tab, "HA_TOT_FLUX_CGS", outdir / "raincloud_HA_TOT_FLUX_by_telescope.png", logx=True)

    plot_filter_warning_vs_ha(tab, masks, outdir / "filter_correction_vs_ha_flux.png")

    print(f"Wrote QC products to {outdir}")


if __name__ == "__main__":
    main()
