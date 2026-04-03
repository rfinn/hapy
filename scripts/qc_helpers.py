from pathlib import Path
import numpy as np
from astropy.table import Table

# ----------------------------------------------------------------------
# helper functions
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

    out = np.full(len(tab), default, dtype=bool)

    for i, v in enumerate(col):
        if v is None:
            continue

        if isinstance(v, (bool, np.bool_)):
            out[i] = bool(v)
            continue

        s = str(v).strip().lower()

        if s in ("true", "t", "1", "yes", "y"):
            out[i] = True
        elif s in ("false", "f", "0", "no", "n", "", "none", "nan"):
            out[i] = False
        # else: leave default

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
        if v is None:
            continue

        try:
            out[i] = float(v)
        except Exception:
            continue

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
# table manipulations
# ----------------------------------------------------------------------


def add_derived_columns(tab):
    def f(name):
        return safe_float_array(tab, name)

    # morphology deltas
    tab["DELTA_GINI"] = f("H_HAPY_GINI") - f("R_HAPY_GINI")
    tab["DELTA_M20"]  = f("H_HAPY_M20")  - f("R_HAPY_M20")
    tab["DELTA_ASYM"] = f("H_HAPY_ASYM") - f("R_HAPY_ASYM")

    # size ratio (safe divide)
    r50 = f("R50_ARCSEC")
    h50 = f("H50_ARCSEC")

    ratio = np.full(len(tab), np.nan)
    good = np.isfinite(r50) & np.isfinite(h50) & (r50 > 0)
    ratio[good] = h50[good] / r50[good]

    tab["H_R50_RATIO"] = ratio



    return tab

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
def add_duplicate_metadata(tab, id_col="VFID"):
    ids = np.array(tab[id_col]).astype(str)

    unique, counts = np.unique(ids, return_counts=True)
    count_map = dict(zip(unique, counts))

    ndup = np.array([count_map[i] for i in ids])

    tab["N_DUP"] = ndup
    tab["IS_DUPLICATE"] = ndup > 1

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

    n = len(tab)
    tier = np.full(n, "F", dtype="U1")

    for i in range(n):
        if not phot_ok[i]:
            tier[i] = "F"
        elif (not use_r[i]) or (not use_ha[i]):
            tier[i] = "D"
        elif mask_warn[i] or bright_star[i] or warn_filter[i] or ell_warn[i]:
            tier[i] = "C"
        elif not use_hm[i]:
            tier[i] = "B"
        else:
            tier[i] = "A"

    tab["QC_TIER"] = tier
    return tab


def prepare_analysis_table(
    tab: Table,
    add_qc: bool = True,
    add_tier: bool = True,    
    add_derived: bool = True,
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
        tab = add_qc_columns(tab)

    if add_tier:
        tab = add_qc_tier(tab)

    if add_derived:
        tab = add_derived_columns(tab)

    if add_duplicates:
        tab = add_duplicate_metadata(tab)

    if add_science:
        tab = add_science_columns(tab)

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



def select_sample(tab: Table, sample_name: str) -> np.ndarray:
    """
    Return boolean mask for sample selection.
    """
    tier = safe_str_array(tab, "QC_TIER", default="F")
    tier = np.array([str(t).upper() for t in tier], dtype=object)

    sample_name = sample_name.upper()
    if sample_name == "A":
        return tier == "A"
    if sample_name == "AB":
        return np.isin(tier, ["A", "B"])
    if sample_name == "ABC":
        return np.isin(tier, ["A", "B", "C"])
    if sample_name == "ALL":
        return np.ones(len(tab), dtype=bool)

    raise ValueError(f"Unknown sample selection: {sample_name}")


# def get_review_priority(tab: Table) -> np.ndarray:
#     """
#     Return review priority labels: high / medium / low.

#     High = likely bug / failure / severe contamination / extreme outlier
#     Medium = warning or cautionary case
#     Low = likely fine
#     """
#     if "QC_TIER" not in tab.colnames:
#         tab = prepare_analysis_table(tab, copy=False)

#     n = len(tab)
#     priority = np.full(n, "low", dtype="U16")

#     tier = safe_str_array(tab, "QC_TIER", default="")

#     # -----------------------------
#     # Base priority from QC tier
#     # -----------------------------
#     priority[np.isin(tier, ["D", "F"])] = "high"
#     priority[np.isin(tier, ["C"])] = "medium"

#     # -----------------------------
#     # Promote only severe warning cases
#     # -----------------------------
#     # severe_flags = (
#     #     safe_bool_array(tab, "ELL_MISMATCH") #| safe_bool_array(tab, "FILTER_WARNING")
#     # )
#     # priority[severe_flags] = "high"

#     # -----------------------------
#     # Moderate warnings stay medium
#     # -----------------------------
#     moderate_flags = (
#         safe_bool_array(tab, "WARN_MASK") |
#         safe_bool_array(tab, "BRIGHT_STAR_FLAG") |
#         safe_bool_array(tab, "WARN_WEAK_HA")
#     )
#     promote_medium = moderate_flags & (priority != "high")
#     priority[promote_medium] = "medium"

#     # -----------------------------
#     # Promote extreme outliers to high
#     # -----------------------------
#     if "H50_R50_RATIO" in tab.colnames:
#         ratio = safe_float_array(tab, "H50_R50_RATIO")
#         extreme_ratio = np.isfinite(ratio) & ((ratio > 3.0) | (ratio < 0.2))
#         priority[extreme_ratio] = "high"

#     if "DELTA_GINI" in tab.colnames:
#         dg = safe_float_array(tab, "DELTA_GINI")
#         extreme_dg = np.isfinite(dg) & (np.abs(dg) > 0.5)
#         priority[extreme_dg] = "high"

#     if "DELTA_M20" in tab.colnames:
#         dm20 = safe_float_array(tab, "DELTA_M20")
#         extreme_dm20 = np.isfinite(dm20) & (np.abs(dm20) > 2.0)
#         priority[extreme_dm20] = "high"

#     priorities = ["low", "medium", "high"]
#     print("REVIEW_PRIORITY SUMMARY")
#     for p in priorities:
#         Np = np.sum(priority == p)
#         print(f"\t{p} = {Np} ({Np/len(priority)*100:.1f} )%")
        
#     return priority


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
        safe_bool_array(tab, "WARN_MASK") 
    )

    priority[high] = "high"

    # -----------------------------
    # MEDIUM: cautionary
    # -----------------------------
    medium = (
        safe_bool_array(tab, "ELL_MISMATCH") |
        safe_bool_array(tab, "WARN_MASK") |
        safe_bool_array(tab, "WARN_WEAK_HA") |
        safe_bool_array(tab, "FILTER_WARNING")
    )

    priority[medium & (~high)] = "medium"

    return priority


