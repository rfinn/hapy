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
from astropy.table import Column


#from hapy.utils.astro import KE_SFR_from_redshift

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

def safe_str_array(tab: Table, colname: str, default: str = "") -> np.ndarray:
    """Return a clean string numpy array for a table column."""
    if colname not in tab.colnames:
        return np.full(len(tab), default, dtype=object)

    col = tab[colname]

    try:
        if hasattr(col, "filled"):
            col = col.filled(default)
    except Exception:
        pass

    out = np.full(len(tab), default, dtype=object)

    for i, v in enumerate(col):
        if v is None:
            continue

        s = str(v).strip()

        if s.lower() in ("none", "nan", "null", ""):
            continue

        out[i] = s

    return out

def safe_ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """Safely divide two float arrays."""
    out = np.full(len(num), np.nan, dtype=float)
    good = np.isfinite(num) & np.isfinite(den) & (den > 0)
    out[good] = num[good] / den[good]
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

def median_and_mad(dx: np.ndarray) -> tuple[float, float]:
    """Return median and MAD-like robust scatter."""
    good = np.isfinite(dx)
    if np.sum(good) == 0:
        return np.nan, np.nan
    med = np.nanmedian(dx[good])
    mad = np.nanmedian(np.abs(dx[good] - med))
    return med, mad

def get_std(dx: np.ndarray) -> tuple[float, float]:
    """Return median and MAD-like robust scatter."""
    good = np.isfinite(dx)
    if np.sum(good) == 0:
        return np.nan, np.nan
    std = np.nanstd(dx[good])
    return std


# ----------------------------------------------------------------------
# science columns
# ----------------------------------------------------------------------

def add_science_columns(tab: Table) -> Table:
    """
    Add core derived science columns if missing.
    Safe to call multiple times.
    """
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

    # optional, but useful for consistency with validation/science plots
    r_asym = safe_float_array(tab, "R_HAPY_ASYM")
    h_asym = safe_float_array(tab, "H_HAPY_ASYM")

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
    if "DELTA_ASYM" not in tab.colnames:
        tab["DELTA_ASYM"] = h_asym - r_asym

    return tab




# ----------------------------------------------------------------------
# QC flags
# ----------------------------------------------------------------------

def add_duplicate_metadata(tab, id_col="VFID"):
    ids = np.array(tab[id_col]).astype(str)

    unique, counts = np.unique(ids, return_counts=True)
    count_map = dict(zip(unique, counts))

    ndup = np.array([count_map[i] for i in ids])

    tab["N_DUP"] = ndup
    tab["IS_DUPLICATE"] = ndup > 1

    return tab




def add_center_offset_columns(tab):
    """
    Add center-offset columns comparing:
      - input coords vs photutils center
      - input coords vs GALFIT center
      - input coords vs GALFIT+conv center
      - photutils center vs GALFIT center
      - photutils center vs GALFIT+conv center
      - GALFIT center vs GALFIT+conv center

    Offsets are added in both pixels and arcsec.

    Expected input columns
    ----------------------
    ELL0_XC, ELL0_YC
    ELLIP_XCENTROID, ELLIP_YCENTROID
    GAL_XC, GAL_YC
    GAL_CXC, GAL_CYC
    PIXSCALE
    """

    def _get(name):
        return np.asarray(tab[name], dtype=float)

    def _offsets(x1, y1, x2, y2, pixscale):
        dx = x2 - x1
        dy = y2 - y1
        dr_pix = np.hypot(dx, dy)

        bad = (
            ~np.isfinite(x1) | ~np.isfinite(y1) |
            ~np.isfinite(x2) | ~np.isfinite(y2)
        )
        dr_pix[bad] = np.nan

        dr_arcsec = dr_pix * pixscale
        bad_arc = bad | ~np.isfinite(pixscale) | (pixscale <= 0)
        dr_arcsec[bad_arc] = np.nan

        dx[bad] = np.nan
        dy[bad] = np.nan
        return dx, dy, dr_pix, dr_arcsec

    xin = _get("ELL0_XC")
    yin = _get("ELL0_YC")

    xph = _get("ELLIP_XCENTROID")
    yph = _get("ELLIP_YCENTROID")

    xg = _get("GAL_XC")
    yg = _get("GAL_YC")

    xgc = _get("GAL_CXC")
    ygc = _get("GAL_CYC")

    pixscale = _get("PIXSCALE")

    pairs = [
        ("IN_PHOT", xin, yin, xph, yph),
        ("IN_GAL", xin, yin, xg, yg),
        ("IN_GALC", xin, yin, xgc, ygc),
        ("PHOT_GAL", xph, yph, xg, yg),
        ("PHOT_GALC", xph, yph, xgc, ygc),
        ("GAL_GALC", xg, yg, xgc, ygc),
    ]

    for tag, x1, y1, x2, y2 in pairs:
        dx, dy, dr_pix, dr_arcsec = _offsets(x1, y1, x2, y2, pixscale)

        cols = {
            f"DX_{tag}_PIX": dx,
            f"DY_{tag}_PIX": dy,
            f"DOFF_{tag}_PIX": dr_pix,
            f"DOFF_{tag}_ARCSEC": dr_arcsec,
        }

        for name, values in cols.items():
            if name in tab.colnames:
                tab[name] = values
            else:
                tab.add_column(Column(values, name=name))

    return tab


def add_center_offset_flags(tab, warn_arcsec=2.0, severe_arcsec=5.0):
    """
    Add QC flags based on center offsets in arcsec.

    Uses these scalar offset columns (created by add_center_offset_columns):
      DOFF_IN_PHOT_ARCSEC
      DOFF_IN_GAL_ARCSEC
      DOFF_IN_GALC_ARCSEC
      DOFF_PHOT_GAL_ARCSEC
      DOFF_PHOT_GALC_ARCSEC
      DOFF_GAL_GALC_ARCSEC

    Parameters
    ----------
    warn_arcsec : float
        Threshold for warning flags.
    severe_arcsec : float
        Threshold for severe flags.
    """

    def _flag_from_offset(colname, thresh):
        vals = np.asarray(tab[colname], dtype=float)
        return np.isfinite(vals) & (vals > thresh)

    flag_map_warn = {
        "WARN_CEN_IN_PHOT": "DOFF_IN_PHOT_ARCSEC",
        "WARN_CEN_IN_GAL": "DOFF_IN_GAL_ARCSEC",
        "WARN_CEN_IN_GALC": "DOFF_IN_GALC_ARCSEC",
        "WARN_CEN_PHOT_GAL": "DOFF_PHOT_GAL_ARCSEC",
        "WARN_CEN_PHOT_GALC": "DOFF_PHOT_GALC_ARCSEC",
        "WARN_CEN_GAL_GALC": "DOFF_GAL_GALC_ARCSEC",
    }

    flag_map_severe = {
        "SEVERE_CEN_IN_PHOT": "DOFF_IN_PHOT_ARCSEC",
        "SEVERE_CEN_IN_GAL": "DOFF_IN_GAL_ARCSEC",
        "SEVERE_CEN_IN_GALC": "DOFF_IN_GALC_ARCSEC",
        "SEVERE_CEN_PHOT_GAL": "DOFF_PHOT_GAL_ARCSEC",
        "SEVERE_CEN_PHOT_GALC": "DOFF_PHOT_GALC_ARCSEC",
        "SEVERE_CEN_GAL_GALC": "DOFF_GAL_GALC_ARCSEC",
    }

    for flagname, offcol in flag_map_warn.items():
        vals = _flag_from_offset(offcol, warn_arcsec)
        if flagname in tab.colnames:
            tab[flagname] = vals
        else:
            tab.add_column(Column(vals, name=flagname))

    for flagname, offcol in flag_map_severe.items():
        vals = _flag_from_offset(offcol, severe_arcsec)
        if flagname in tab.colnames:
            tab[flagname] = vals
        else:
            tab.add_column(Column(vals, name=flagname))

    # combined summary flags
    warn_any = np.zeros(len(tab), dtype=bool)
    severe_any = np.zeros(len(tab), dtype=bool)

    for name in flag_map_warn:
        warn_any |= np.asarray(tab[name], dtype=bool)

    for name in flag_map_severe:
        severe_any |= np.asarray(tab[name], dtype=bool)

    if "WARN_CEN_ANY" in tab.colnames:
        tab["WARN_CEN_ANY"] = warn_any
    else:
        tab.add_column(Column(warn_any, name="WARN_CEN_ANY"))

    if "SEVERE_CEN_ANY" in tab.colnames:
        tab["SEVERE_CEN_ANY"] = severe_any
    else:
        tab.add_column(Column(severe_any, name="SEVERE_CEN_ANY"))

    return tab

def build_row_qc_flags(tab, max_ha_filter_correction: float = 1.2) -> dict[str, np.ndarray]:
    """
    Unified row-level QC flags for merged HAPY results.

    Returns a dict of boolean and float arrays that both qc_results.py
    and qc_duplicates.py can use.
    """
    flags = {}

    # -----------------------------
    # Core stage booleans
    # -----------------------------
    flags["PSF_OK"] = safe_bool_array(tab, "PSF_OK")
    flags["MASK_OK"] = safe_bool_array(tab, "MASK_OK")
    flags["PHOT_OK"] = safe_bool_array(tab, "PHOT_OK")
    flags["R_PROFILE_OK"] = safe_bool_array(tab, "R_PROFILE_OK")
    flags["H_PROFILE_OK"] = safe_bool_array(tab, "H_PROFILE_OK")
    flags["WARN_R_PROFILE_PEAK"] = safe_bool_array(tab, "R_PROFILE_NONCENTRAL_PEAK")
    flags["R_SM_OK"] = safe_bool_array(tab, "R_SM_OK")
    flags["H_SM_OK"] = safe_bool_array(tab, "H_SM_OK")
    flags["GAL_NC_OK"] = safe_bool_array(tab, "GAL_NC_OK")
    flags["GAL_CV_OK"] = safe_bool_array(tab, "GAL_CV_OK")
    flags["HAPY_MORPH_OK"] = safe_bool_array(tab, "HAPY_MORPH_OK")

    # -----------------------------
    # Warning booleans
    # -----------------------------
    flags["BRIGHT_STAR_FLAG"] = safe_bool_array(tab, "BRIGHT_STAR_FLAG")
    flags["ELL0_MASK_WARN"] = safe_bool_array(tab, "ELL0_MASK_WARN")
    flags["ELL_MISMATCH"] = safe_bool_array(tab, "ELL_MISMATCH")

    # -----------------------------
    # Filter correction
    # -----------------------------
    filt_col = first_existing_col(tab, ["FILTER_CORRECTION", "FILT_COR"])
    if filt_col is not None:
        filtcor = safe_float_array(tab, filt_col)
    else:
        filtcor = np.full(len(tab), np.nan)

    flags["FILTER_CORR"] = filtcor
    flags["FILTER_WARNING"] = np.isfinite(filtcor) & (filtcor > max_ha_filter_correction)

    # -----------------------------
    # Numeric QC helpers
    # -----------------------------
    flags["R_PROFILE_NGOOD"] = safe_float_array(tab, "R_PROFILE_NGOOD")
    flags["H_PROFILE_NGOOD"] = safe_float_array(tab, "H_PROFILE_NGOOD")
    flags["MASKFRAC_GUESS_ELLIPSE"] = safe_float_array(tab, "MASKFRAC_GUESS_ELLIPSE")
    flags["R_PROFILE_MASKFRAC_MAX"] = safe_float_array(tab, "R_PROFILE_MASKFRAC_MAX")    
    flags["H_PROFILE_MASKFRAC_MAX"] = safe_float_array(tab, "H_PROFILE_MASKFRAC_MAX")
    flags["H_HAPY_NPIX"] = safe_float_array(tab, "H_HAPY_NPIX")
    flags["H_HAPY_SNP_DET"] = safe_float_array(tab, "H_HAPY_SNP_DET")
    flags["R50_ARCSEC"] = safe_float_array(tab, "R50_ARCSEC")
    flags["H50_ARCSEC"] = safe_float_array(tab, "H50_ARCSEC")
    flags["H_MAXDET_ARCSEC"] = safe_float_array(tab, "H_MAXDET_ARCSEC")

    # -----------------------------
    # Science-oriented row families
    # -----------------------------
    flags["MASK_PHOT_OK"] = flags["MASK_OK"] & flags["PHOT_OK"]

    flags["PROFILE_OK"] = (
        flags["MASK_OK"] &
        flags["PHOT_OK"] &
        flags["R_PROFILE_OK"] &
        flags["H_PROFILE_OK"]
    )

    flags["R_STRUCTURE_GOOD"] = (
        flags["PHOT_OK"] &
        flags["R_PROFILE_OK"] &
        np.isfinite(flags["R50_ARCSEC"]) & (flags["R50_ARCSEC"] > 0) &
        np.isfinite(flags["R_PROFILE_NGOOD"]) & (flags["R_PROFILE_NGOOD"] >= 20) &
        (
            ~np.isfinite(flags["R_PROFILE_MASKFRAC_MAX"]) |
            (flags["R_PROFILE_MASKFRAC_MAX"] < 0.5)
        ) &
        (~flags["ELL0_MASK_WARN"])
    )

    flags["HA_EXTENT_GOOD"] = (
        flags["PHOT_OK"] &
        flags["H_PROFILE_OK"] &
        np.isfinite(flags["H50_ARCSEC"]) & (flags["H50_ARCSEC"] > 0) &
        np.isfinite(flags["H_MAXDET_ARCSEC"]) & (flags["H_MAXDET_ARCSEC"] > 0) &
        np.isfinite(flags["H_HAPY_NPIX"]) & (flags["H_HAPY_NPIX"] >= 50) &
        np.isfinite(flags["H_PROFILE_NGOOD"]) & (flags["H_PROFILE_NGOOD"] >= 8) &
        (
            ~np.isfinite(flags["H_HAPY_SNP_DET"]) |
            (flags["H_HAPY_SNP_DET"] > 3)
        ) &
        (
            ~np.isfinite(flags["H_PROFILE_MASKFRAC_MAX"]) |
            (flags["H_PROFILE_MASKFRAC_MAX"] < 0.5)
        ) &
        (~flags["FILTER_WARNING"])
    )

    flags["HA_MORPH_GOOD"] = (
        flags["HAPY_MORPH_OK"] &
        flags["H_SM_OK"] &
        np.isfinite(flags["H_HAPY_NPIX"]) & (flags["H_HAPY_NPIX"] >= 100) &
        (
            ~np.isfinite(flags["H_HAPY_SNP_DET"]) |
            (flags["H_HAPY_SNP_DET"] > 5)
        ) &
        (~flags["FILTER_WARNING"])
    )

    flags["GALFIT_ANY_OK"] = flags["GAL_NC_OK"] | flags["GAL_CV_OK"]
    flags["GALFIT_BOTH_OK"] = flags["GAL_NC_OK"] & flags["GAL_CV_OK"]

    flags["WARN_MASK"] = (
        flags["ELL0_MASK_WARN"] |
        (np.isfinite(flags["MASKFRAC_GUESS_ELLIPSE"]) & (flags["MASKFRAC_GUESS_ELLIPSE"] > 0.3))
        #(np.isfinite(flags["R_PROFILE_MASKFRAC_MAX"]) & (flags["R_PROFILE_MASKFRAC_MAX"] > 0.3)) |
        #(np.isfinite(flags["H_PROFILE_MASKFRAC_MAX"]) & (flags["H_PROFILE_MASKFRAC_MAX"] > 0.3))
    )

    flags["WARN_WEAK_HA"] = (
        (np.isfinite(flags["H_HAPY_NPIX"]) & (flags["H_HAPY_NPIX"] < 50)) |
        (np.isfinite(flags["H_PROFILE_NGOOD"]) & (flags["H_PROFILE_NGOOD"] < 8)) |
        (np.isfinite(flags["H_HAPY_SNP_DET"]) & (flags["H_HAPY_SNP_DET"] < 3))
    )

    flags["SCIENCE_READY"] = (
        flags["R_STRUCTURE_GOOD"] &
        flags["HA_EXTENT_GOOD"] &
        (~flags["BRIGHT_STAR_FLAG"]) &
        (~flags["ELL_MISMATCH"])
    )

    flags["TECHNICAL_PROBLEM"] = ~flags["MASK_PHOT_OK"]

    flags["SCIENCE_PROBLEM"] = (
        ~flags["PROFILE_OK"] |
        flags["FILTER_WARNING"] |
        flags["WARN_MASK"] |
        flags["WARN_WEAK_HA"] |
        flags["BRIGHT_STAR_FLAG"] |
        flags["ELL_MISMATCH"]
    )

    return flags


# -- this is a duplicate function from when I was building this in scripts
# -- this is currently NOT USED!!!
# -- keeping in case some part is useful down the road...
# def add_qc_flags(tab: Table, qc: dict | None = None) -> Table:
#     """
#     Add QC tiering and usability flags.

#     Parameters
#     ----------
#     tab : astropy.table.Table
#         Input merged results table.

#     qc : dict or None
#         Optional override dictionary for QC thresholds.
#         Keys may include:
#           - r_profile_ngood_min
#           - ha_profile_ngood_min
#           - ha_npix_min_extent
#           - ha_npix_min_morph
#           - ha_snr_det_min
#           - filter_correction_warn
#     """

#     if "QC_TIER" in tab.colnames:
#         return tab

#     cfg = QC_DEFAULTS.copy()
#     if qc is not None:
#         cfg.update(qc)

#     n = len(tab)

#     phot_ok = safe_bool_array(tab, "PHOT_OK")
#     rprof_ok = safe_bool_array(tab, "R_PROFILE_OK")
#     hprof_ok = safe_bool_array(tab, "H_PROFILE_OK")
#     morph_ok = safe_bool_array(tab, "HAPY_MORPH_OK")
#     h_sm_ok = safe_bool_array(tab, "H_SM_OK")
#     gal_nc_ok = safe_bool_array(tab, "GAL_NC_OK")
#     gal_cv_ok = safe_bool_array(tab, "GAL_CV_OK")

#     bright_star = safe_bool_array(tab, "BRIGHT_STAR_FLAG")
#     mask_warn = safe_bool_array(tab, "ELL0_MASK_WARN")
#     ell_warn = safe_bool_array(tab, "ELL_MISMATCH")

#     filt = safe_float_array(tab, "FILTER_CORRECTION")
#     warn_filter = np.isfinite(filt) & (filt > cfg["filter_correction_warn"])

#     r50 = safe_float_array(tab, "R50_ARCSEC")
#     h50 = safe_float_array(tab, "H50_ARCSEC")
#     hmax = safe_float_array(tab, "H_MAXDET_ARCSEC")

#     h_npix = safe_float_array(tab, "H_HAPY_NPIX")
#     h_ngood = safe_float_array(tab, "H_PROFILE_NGOOD")
#     h_snr = safe_float_array(tab, "H_HAPY_SNP_DET")
#     r_ngood = safe_float_array(tab, "R_PROFILE_NGOOD")

#     use_r = np.zeros(n, dtype=int)
#     use_ha = np.zeros(n, dtype=int)
#     use_hm = np.zeros(n, dtype=int)
#     use_gf = np.zeros(n, dtype=int)

#     good_r = (
#         phot_ok &
#         rprof_ok &
#         np.isfinite(r50) & (r50 > 0) &
#         np.isfinite(r_ngood) & (r_ngood >= cfg["r_profile_ngood_min"])
#     )
#     use_r[good_r] = 2
#     use_r[good_r & mask_warn] = 1

#     good_ha = (
#         phot_ok &
#         hprof_ok &
#         np.isfinite(h50) & (h50 > 0) &
#         np.isfinite(hmax) & (hmax > 0) &
#         np.isfinite(h_npix) & (h_npix >= cfg["ha_npix_min_extent"]) &
#         np.isfinite(h_ngood) & (h_ngood >= cfg["ha_profile_ngood_min"]) &
#         (~warn_filter)
#     )
#     use_ha[good_ha] = 2

#     weak_ha = (
#         (np.isfinite(h_npix) & (h_npix < cfg["ha_npix_min_extent"])) |
#         (np.isfinite(h_ngood) & (h_ngood < cfg["ha_profile_ngood_min"])) |
#         (np.isfinite(h_snr) & (h_snr < cfg["ha_snr_det_min"]))
#     )
#     use_ha[good_ha & weak_ha] = 1

#     good_hm = (
#         morph_ok &
#         h_sm_ok &
#         np.isfinite(h_npix) & (h_npix >= cfg["ha_npix_min_morph"]) &
#         (~warn_filter)
#     )
#     use_hm[good_hm] = 2
#     use_hm[good_hm & weak_ha] = 1

#     use_gf[gal_nc_ok | gal_cv_ok] = 2

#     tab["USE_R_STRUCTURE"] = use_r
#     tab["USE_HA_EXTENT"] = use_ha
#     tab["USE_HA_MORPH"] = use_hm
#     tab["USE_GALFIT"] = use_gf

#     if "WARN_FILTER" not in tab.colnames:
#         tab["WARN_FILTER"] = warn_filter
#     if "WARN_MASK" not in tab.colnames:
#         tab["WARN_MASK"] = mask_warn
#     if "WARN_BRIGHT_STAR" not in tab.colnames:
#         tab["WARN_BRIGHT_STAR"] = bright_star
#     if "WARN_ELLIPSE" not in tab.colnames:
#         tab["WARN_ELLIPSE"] = ell_warn
#     if "WARN_WEAK_HA" not in tab.colnames:
#         tab["WARN_WEAK_HA"] = weak_ha

#     center_warn = tab["WARN_CEN_ANY"]
    
#     tier = np.full(n, "F", dtype="U1")

#     for i in range(n):
#         if not phot_ok[i]:
#             tier[i] = "F"
#         elif use_r[i] == 0 or use_ha[i] == 0:
#             tier[i] = "D"
#         elif mask_warn[i] or bright_star[i] or warn_filter[i] or ell_warn[i] or center_warn[i]:
#         tier[i] = "C"
#         elif use_hm[i] < 2:
#             tier[i] = "B"
#         else:
#             tier[i] = "A"

#     tab["QC_TIER"] = tier

#     return tab


# ----------------------------------------------------------------------
# pipeline interface
# ----------------------------------------------------------------------

# def prepare_results_table(tab: Table, qc: dict | None = None) -> Table:
#     """
#     Apply standard derived science and QC columns.

#     Parameters
#     ----------
#     tab : astropy.table.Table
#         Input merged results table.

#     qc : dict or None
#         Optional QC-threshold overrides passed to add_qc_flags().
#     """
#     tab = add_science_columns(tab)
#     tab = add_qc_flags(tab, qc=qc)
#     return tab


def add_qc_columns(tab: Table, max_ha_filter_correction: float = 1.2) -> Table:
    """
    Attach QC columns from build_row_qc_flags() to the table.

    This replaces the older science_firstlook-specific add_qc_flags()
    implementation and uses the canonical QC logic.
    """
    flags = build_row_qc_flags(tab, max_ha_filter_correction=max_ha_filter_correction)
    for key, val in flags.items():
        if key not in tab.colnames:
            tab[key] = val
    return tab


def add_qc_tier(tab: Table) -> Table:
    """
    Add a minimal QC tier if missing.

    Uses canonical QC columns if available; otherwise computes them first.
    Safe to call multiple times.
    """
    if "QC_TIER" in tab.colnames:
        return tab

    if "R_STRUCTURE_GOOD" not in tab.colnames:
        tab = add_qc_columns(tab)

    phot_ok = safe_bool_array(tab, "PHOT_OK")
    use_r = safe_bool_array(tab, "R_STRUCTURE_GOOD")
    use_ha = safe_bool_array(tab, "HA_EXTENT_GOOD")
    use_hm = safe_bool_array(tab, "HA_MORPH_GOOD")

    bright_star = safe_bool_array(tab, "BRIGHT_STAR_FLAG")
    mask_warn = safe_bool_array(tab, "WARN_MASK")
    ell_warn = safe_bool_array(tab, "ELL_MISMATCH")
    warn_filter = safe_bool_array(tab, "FILTER_WARNING")
    center_warn = tab["WARN_CEN_ANY"]
    r_profile_offcenter = tab["WARN_R_PROFILE_PEAK"]
    n = len(tab)
    tier = np.full(n, "F", dtype="U1")

    for i in range(n):
        if not phot_ok[i]:
            tier[i] = "F"
        elif (not use_r[i]) or (not use_ha[i]):
            tier[i] = "D"
        elif mask_warn[i] or bright_star[i] or warn_filter[i] or ell_warn[i] or center_warn[i] or r_profile_offcenter[i]:
            tier[i] = "C"
        elif not use_hm[i]:
            tier[i] = "B"
        else:
            tier[i] = "A"

    tab["QC_TIER"] = tier
    return tab

def add_vfindex(tab):
    vfindex = np.zeros(len(tab))

    for i in range(len(tab)):
        vfindex[i] = int(tab['VFID'][i].replace('VFID',''))
    vfindex = np.array(vfindex,'i')
    tab["VFINDEX"] = vfindex
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


def get_review_priority(tab: Table) -> np.ndarray:
    """
    Review priority focused on core science reliability.

    High = fundamental measurement failure
    Medium = caution / common issues
    Low = not urgent
    """
    if "QC_TIER" not in tab.colnames:
        tab = prepare_analysis_table(tab, copy=False)

    n = len(tab)
    priority = np.full(n, "low", dtype="U16")

    # -----------------------------
    # HIGH: core failures
    # -----------------------------
    high = (
        ~safe_bool_array(tab, "PHOT_OK") |
        ~safe_bool_array(tab, "HAPY_MORPH_OK") |
        safe_bool_array(tab, "BRIGHT_STAR_FLAG") |
        safe_bool_array(tab, "WARN_MASK") |
        safe_bool_array(tab, "SEVERE_CEN_ANY")
    )

    priority[high] = "high"

    # -----------------------------
    # MEDIUM: cautionary
    # -----------------------------
    medium = (
        safe_bool_array(tab, "ELL_MISMATCH") |
        safe_bool_array(tab, "WARN_WEAK_HA") |
        safe_bool_array(tab, "FILTER_WARNING") |
        safe_bool_array(tab, "WARN_CEN_ANY") |
        safe_bool_array(tab, "WARN_R_PROFILE_PEAK")
    )

    priority[medium & (~high)] = "medium"

    return priority


def prepare_analysis_table(
    tab: Table,
    add_qc: bool = True,
    add_tier: bool = True,    
    #add_derived: bool = True,
    add_duplicates=True,    
    add_science: bool = True,    
    copy: bool = True,
) -> Table:
    """
    Standard table preparation for QC, validation, and first-look science.

    Safe to call from qc_results.py, validate_measurements.py,
    validate_duplicates.py, validate_dashboards.py, and science_firstlook.py.
    """
    if copy:
        tab = tab.copy()
        
    if add_qc:
        tab = add_center_offset_columns(tab)
        tab = add_center_offset_flags(tab)
        tab = add_qc_columns(tab)

    if add_tier:
        tab = add_qc_tier(tab)

    #if add_derived:
    #    tab = add_derived_columns(tab)

    if add_duplicates:
        tab = add_duplicate_metadata(tab)

    if add_science:
        tab = add_science_columns(tab)

    # add vfindex
    if ("VFID" in tab.colnames) and ("VFINDEX" not in tab.colnames):
        tab = add_vfindex(tab)
        
    if "REVIEW_PRIORITY" not in tab.colnames:
        tab["REVIEW_PRIORITY"] = get_review_priority(tab)
        priority = tab["REVIEW_PRIORITY"]
        vals, counts = np.unique(priority, return_counts=True)
        print("REVIEW_PRIORITY SUMMARY")
        print(dict(zip(vals, counts)))

        for name in ["ELL_MISMATCH", "FILTER_WARNING", "WARN_MASK", "BRIGHT_STAR_FLAG", "WARN_WEAK_HA"]:
            arr = safe_bool_array(tab, name)
            print(name, np.sum(arr))

    return tab
