from pathlib import Path
import numpy as np
from astropy.table import Table

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

def median_and_mad(dx: np.ndarray) -> tuple[float, float]:
    """Return median and MAD-like robust scatter."""
    good = np.isfinite(dx)
    if np.sum(good) == 0:
        return np.nan, np.nan
    med = np.nanmedian(dx[good])
    mad = np.nanmedian(np.abs(dx[good] - med))
    return med, mad

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
        (np.isfinite(flags["R_PROFILE_MASKFRAC_MAX"]) & (flags["R_PROFILE_MASKFRAC_MAX"] > 0.3)) |
        (np.isfinite(flags["H_PROFILE_MASKFRAC_MAX"]) & (flags["H_PROFILE_MASKFRAC_MAX"] > 0.3))
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
