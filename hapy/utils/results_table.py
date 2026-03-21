#!/usr/bin/env python

"""
results_table.py

Utilities for working with merged HAPY results tables.

Provides:
- safe column extraction helpers
- derived science quantities
- QC flags and tiers
- sample selection

This module is intended to be the single source of truth for
table-level logic used across QC, validation, and science scripts.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table
from hapy.utils.astro import KE_SFR

# ----------------------------------------------------------------------
# QC defaults
# ----------------------------------------------------------------------

QC_DEFAULTS = {
    "r_profile_ngood_min": 20,
    "ha_profile_ngood_min": 8,
    "ha_npix_min_extent": 50,
    "ha_npix_min_morph": 100,
    "ha_snr_det_min": 3.0,
    "filter_correction_warn": 1.2,
}


# ----------------------------------------------------------------------
# basic helpers
# ----------------------------------------------------------------------

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


# ----------------------------------------------------------------------
# science columns
# ----------------------------------------------------------------------

def _safe_ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    out = np.full(len(num), np.nan, dtype=float)
    good = np.isfinite(num) & np.isfinite(den) & (den > 0)
    out[good] = num[good] / den[good]
    return out


def add_science_columns(tab: Table) -> Table:
    """
    Add derived science columns used in analysis.
    """

    r50 = safe_float_array(tab, "R50_ARCSEC")
    h50 = safe_float_array(tab, "H50_ARCSEC")
    r25 = safe_float_array(tab, "R25_ARCSEC")
    hmax = safe_float_array(tab, "H_MAXDET_ARCSEC")

    r_p50 = safe_float_array(tab, "R_PETRO_R50_ARCSEC")
    h_p50 = safe_float_array(tab, "H_PETRO_R50_ARCSEC")

    r_gini = safe_float_array(tab, "R_HAPY_GINI")
    h_gini = safe_float_array(tab, "H_HAPY_GINI")

    r_m20 = safe_float_array(tab, "R_HAPY_M20")
    h_m20 = safe_float_array(tab, "H_HAPY_M20")

    if "H50_R50_RATIO" not in tab.colnames:
        tab["H50_R50_RATIO"] = _safe_ratio(h50, r50)

    if "H_MAXDET_R25_RATIO" not in tab.colnames:
        tab["H_MAXDET_R25_RATIO"] = _safe_ratio(hmax, r25)

    if "H_PETRO_R50_RATIO" not in tab.colnames:
        tab["H_PETRO_R50_RATIO"] = _safe_ratio(h_p50, r_p50)

    if "DELTA_GINI" not in tab.colnames:
        tab["DELTA_GINI"] = h_gini - r_gini

    if "DELTA_M20" not in tab.colnames:
        tab["DELTA_M20"] = h_m20 - r_m20

    # Optional Halpha SFR column
    redshift_col = None
    for c in ["REDSHIFT", "Z", "ZDIST", "ZDIST", "vr"]:
        if c in tab.colnames:
            redshift_col = c
            break

    if "LOG_SFR_HA" not in tab.colnames and redshift_col is not None:
        haflux = safe_float_array(tab, "H_TOT_FLUX_CGS")
        z = safe_float_array(tab, redshift_col)

        # If only velocity exists, convert to redshift approximately
        if redshift_col.lower() == "vr":
            z = z / 3.0e5

        tab["LOG_SFR_HA"] = KE_SFR(haflux, z)
        
    return tab


# ----------------------------------------------------------------------
# QC flags
# ----------------------------------------------------------------------

def add_qc_flags(tab: Table, qc: dict | None = None) -> Table:
    """
    Add QC tiering and usability flags.

    Parameters
    ----------
    tab : astropy.table.Table
        Input merged results table.

    qc : dict or None
        Optional override dictionary for QC thresholds.
        Keys may include:
          - r_profile_ngood_min
          - ha_profile_ngood_min
          - ha_npix_min_extent
          - ha_npix_min_morph
          - ha_snr_det_min
          - filter_correction_warn
    """

    if "QC_TIER" in tab.colnames:
        return tab

    cfg = QC_DEFAULTS.copy()
    if qc is not None:
        cfg.update(qc)

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
    warn_filter = np.isfinite(filt) & (filt > cfg["filter_correction_warn"])

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
        phot_ok &
        rprof_ok &
        np.isfinite(r50) & (r50 > 0) &
        np.isfinite(r_ngood) & (r_ngood >= cfg["r_profile_ngood_min"])
    )
    use_r[good_r] = 2
    use_r[good_r & mask_warn] = 1

    good_ha = (
        phot_ok &
        hprof_ok &
        np.isfinite(h50) & (h50 > 0) &
        np.isfinite(hmax) & (hmax > 0) &
        np.isfinite(h_npix) & (h_npix >= cfg["ha_npix_min_extent"]) &
        np.isfinite(h_ngood) & (h_ngood >= cfg["ha_profile_ngood_min"]) &
        (~warn_filter)
    )
    use_ha[good_ha] = 2

    weak_ha = (
        (np.isfinite(h_npix) & (h_npix < cfg["ha_npix_min_extent"])) |
        (np.isfinite(h_ngood) & (h_ngood < cfg["ha_profile_ngood_min"])) |
        (np.isfinite(h_snr) & (h_snr < cfg["ha_snr_det_min"]))
    )
    use_ha[good_ha & weak_ha] = 1

    good_hm = (
        morph_ok &
        h_sm_ok &
        np.isfinite(h_npix) & (h_npix >= cfg["ha_npix_min_morph"]) &
        (~warn_filter)
    )
    use_hm[good_hm] = 2
    use_hm[good_hm & weak_ha] = 1

    use_gf[gal_nc_ok | gal_cv_ok] = 2

    tab["USE_R_STRUCTURE"] = use_r
    tab["USE_HA_EXTENT"] = use_ha
    tab["USE_HA_MORPH"] = use_hm
    tab["USE_GALFIT"] = use_gf

    if "WARN_FILTER" not in tab.colnames:
        tab["WARN_FILTER"] = warn_filter
    if "WARN_MASK" not in tab.colnames:
        tab["WARN_MASK"] = mask_warn
    if "WARN_BRIGHT_STAR" not in tab.colnames:
        tab["WARN_BRIGHT_STAR"] = bright_star
    if "WARN_ELLIPSE" not in tab.colnames:
        tab["WARN_ELLIPSE"] = ell_warn
    if "WARN_WEAK_HA" not in tab.colnames:
        tab["WARN_WEAK_HA"] = weak_ha

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


# ----------------------------------------------------------------------
# pipeline interface
# ----------------------------------------------------------------------

def prepare_results_table(tab: Table, qc: dict | None = None) -> Table:
    """
    Apply standard derived science and QC columns.

    Parameters
    ----------
    tab : astropy.table.Table
        Input merged results table.

    qc : dict or None
        Optional QC-threshold overrides passed to add_qc_flags().
    """
    tab = add_science_columns(tab)
    tab = add_qc_flags(tab, qc=qc)
    return tab


def select_sample(tab: Table, sample: str = "AB") -> np.ndarray:
    """
    Select QC subset.

    Parameters
    ----------
    sample : str
        One of: A, AB, ABC, ALL
    """
    tier = safe_str_array(tab, "QC_TIER", default="F")
    sample = sample.upper()

    if sample == "A":
        return tier == "A"
    elif sample == "AB":
        return np.isin(tier, ["A", "B"])
    elif sample == "ABC":
        return np.isin(tier, ["A", "B", "C"])
    elif sample == "ALL":
        return np.ones(len(tab), dtype=bool)

    raise ValueError(f"Unknown sample: {sample}")
