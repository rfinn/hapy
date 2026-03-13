#!/usr/bin/env python

# from chatgpt's adaptation of my fit_profile.py

"""
profile_summary.py

Summary measurements derived from HAPY elliptical-aperture photometry tables.

This module is intended to replace the old script-style logic in
hatools/fit_profile.py with importable, array-based utilities that can be
called directly from run_analysis.py.

Assumptions
-----------
Each input photometry table has one row per elliptical aperture, ordered by
increasing semi-major axis, and follows the newer HAPY phot-table schema:

    sma_arcsec
    sma_pix
    flux_cum
    flux_cum_err
    sb_avg
    sb_avg_err
    sb_avg_snr
    flux_cgs
    flux_cgs_err
    mag_cum
    mag_cum_err
    sb_cgs_arcsec2
    sb_cgs_arcsec2_err
    sb_mag_arcsec2
    sb_mag_arcsec2_err
    masked_fraction

Notes
-----
- Missing/unmeasurable quantities are returned as np.nan.
- Errors are first-order estimates only.
- This module does not read or write files; it operates on in-memory tables.
"""

from __future__ import annotations

from typing import Any

import warnings
import numpy as np
from astropy.table import Table
from scipy import optimize


# ----------------------------------------------------------------------
# generic helpers
# ----------------------------------------------------------------------

def _as_float_array(x) -> np.ndarray:
    return np.asarray(x, dtype=float)


def _safe_ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den == 0:
        return np.nan
    return num / den


def ratio_error(a: float, b: float, erra: float, errb: float) -> float:
    """
    Propagate uncertainty for ratio a/b.
    """
    if not all(np.isfinite([a, b, erra, errb])) or b == 0:
        return np.nan
    return np.sqrt((erra / b) ** 2 + (a * errb / b**2) ** 2)


def safe_mag_from_flux(flux: float, magzp: float) -> float:
    """
    Convert scalar flux to magnitude, returning np.nan for non-positive flux.
    """
    if not np.isfinite(flux) or not np.isfinite(magzp) or flux <= 0:
        return np.nan
    return magzp - 2.5 * np.log10(flux)


def _combine_two_mags(m1: float, m2: float) -> float:
    """
    Combine two cumulative magnitudes by averaging in flux space.

    Preserves the spirit of the legacy fit_profile.py behavior.
    """
    if not np.isfinite(m1) and not np.isfinite(m2):
        return np.nan
    if np.isfinite(m1) and not np.isfinite(m2):
        return m1
    if np.isfinite(m2) and not np.isfinite(m1):
        return m2

    f1 = 10.0 ** (-0.4 * m1)
    f2 = 10.0 ** (-0.4 * m2)
    f = 0.5 * (f1 + f2)
    if f <= 0:
        return np.nan
    return -2.5 * np.log10(f)


def _interp_monotonic(x_target: float, x, y) -> float:
    """
    1D interpolation with monotonic sorting and NaN handling.
    """
    x = _as_float_array(x)
    y = _as_float_array(y)

    good = np.isfinite(x) & np.isfinite(y)
    if np.sum(good) < 2:
        return np.nan

    x = x[good]
    y = y[good]

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    xu, idx = np.unique(x, return_index=True)
    yu = y[idx]

    if len(xu) < 2:
        return np.nan

    if x_target < xu[0] or x_target > xu[-1]:
        return np.nan

    return float(np.interp(x_target, xu, yu))


def _crossing_radius(radius, profile, threshold, decreasing=True, mode="first") -> tuple[float, float]:
    """
    Estimate radius where a profile crosses a threshold.

    Parameters
    ----------
    radius : array-like
    profile : array-like
    threshold : float
    decreasing : bool
        True for profiles that generally decrease with radius.
    mode : {"first", "last"}
        Which crossing to use if there are multiple crossings.

    Returns
    -------
    value, err : float, float
        Radius estimate and simple half-bracket uncertainty.
    """
    
    radius = _as_float_array(radius)
    profile = _as_float_array(profile)

    good = np.isfinite(radius) & np.isfinite(profile)
    if np.sum(good) < 2:
        return np.nan, np.nan

    r = radius[good]
    p = profile[good]

    order = np.argsort(r)
    r = r[order]
    p = p[order]

    if decreasing:
        sign = p - threshold
    else:
        sign = threshold - p

    finite = np.isfinite(sign)
    r = r[finite]
    sign = sign[finite]
    if len(sign) < 2:
        return np.nan, np.nan

    # exact hits
    exact = np.where(sign == 0)[0]
    if len(exact) > 0:
        idx = exact[0] if mode == "first" else exact[-1]
        return float(r[idx]), 0.0

    # sign changes
    cross = np.where(sign[:-1] * sign[1:] < 0)[0]
    if len(cross) == 0:
        return np.nan, np.nan

    i = cross[0] if mode == "first" else cross[-1]
    val = 0.5 * (r[i] + r[i + 1])
    err = 0.5 * abs(r[i + 1] - r[i])
    return float(val), float(err)


def _crossing_bracket_indices(profile, threshold, decreasing=True, mode="first") -> tuple[int | None, int | None]:
    """
    Return bracketing indices around a threshold crossing.
    """
    p = _as_float_array(profile)
    good = np.isfinite(p)
    if np.sum(good) < 2:
        return None, None

    idx = np.where(good)[0]
    pg = p[good]

    if decreasing:
        sign = pg - threshold
    else:
        sign = threshold - pg

    exact = np.where(sign == 0)[0]
    if len(exact) > 0:
        j = exact[0] if mode == "first" else exact[-1]
        ii = idx[j]
        return ii, ii

    cross = np.where(sign[:-1] * sign[1:] < 0)[0]
    if len(cross) == 0:
        return None, None

    j = cross[0] if mode == "first" else cross[-1]
    return idx[j], idx[j + 1]


def _flux_radius(radius, flux_cum, frac) -> tuple[float, float]:
    """
    Radius enclosing a fraction of the maximum cumulative flux.
    """
    radius = _as_float_array(radius)
    flux_cum = _as_float_array(flux_cum)

    good = np.isfinite(radius) & np.isfinite(flux_cum)
    if np.sum(good) < 2:
        return np.nan, np.nan

    r = radius[good]
    f = flux_cum[good]

    order = np.argsort(r)
    r = r[order]
    f = f[order]

    total_flux = np.nanmax(f)
    if not np.isfinite(total_flux) or total_flux <= 0:
        return np.nan, np.nan

    target = frac * total_flux
    val = _interp_monotonic(target, f, r)
    return val, np.nan


def _longest_true_run(mask: np.ndarray) -> int:
    """
    Length of longest contiguous True run.
    """
    mask = np.asarray(mask, dtype=bool)
    best = 0
    run = 0
    for val in mask:
        if val:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return int(best)


def _true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """
    Return inclusive index bounds for contiguous True runs.
    """
    mask = np.asarray(mask, dtype=bool)
    runs = []
    start = None

    for i, val in enumerate(mask):
        if val and start is None:
            start = i
        elif not val and start is not None:
            runs.append((start, i - 1))
            start = None

    if start is not None:
        runs.append((start, len(mask) - 1))

    return runs


def _significant_runs(
    sb_snr,
    sb_cgs,
    min_snr=2.0,
    min_run=2,
) -> list[tuple[int, int]]:
    """
    Find contiguous significant runs in an Halpha profile.

    Significant means:
    - finite sb_snr
    - finite sb_cgs
    - sb_snr > min_snr
    - sb_cgs > 0

    Only runs with length >= min_run are returned.
    """
    sb_snr = _as_float_array(sb_snr)
    sb_cgs = _as_float_array(sb_cgs)

    detect = np.isfinite(sb_snr) & np.isfinite(sb_cgs) & (sb_snr > min_snr) & (sb_cgs > 0)
    runs = _true_runs(detect)
    runs = [(i0, i1) for i0, i1 in runs if (i1 - i0 + 1) >= min_run]
    return runs


# ----------------------------------------------------------------------
# simple model fits
# ----------------------------------------------------------------------

def _exp_profile(r, i0, k):
    return i0 * np.exp(-k * r)


def fit_exponential_profile(radius, sb, sb_snr=None, min_snr=2.0) -> dict[str, float]:
    """
    Fit I(r) = I0 exp(-k r) to positive high-SNR points.
    """
    radius = _as_float_array(radius)
    sb = _as_float_array(sb)

    good = np.isfinite(radius) & np.isfinite(sb) & (sb > 0)
    if sb_snr is not None:
        sb_snr = _as_float_array(sb_snr)
        good &= np.isfinite(sb_snr) & (sb_snr > min_snr)

    if np.sum(good) < 4:
        return {
            "EXPFIT_I0": np.nan,
            "EXPFIT_K": np.nan,
            "EXPFIT_RE_ARCSEC": np.nan,
            "EXPFIT_OK": False,
        }

    try:
        popt, _pcov = optimize.curve_fit(_exp_profile, radius[good], sb[good], maxfev=10000)
        i0, k = popt
        re = np.nan if (not np.isfinite(k) or k == 0) else 1.0 / k
        return {
            "EXPFIT_I0": float(i0),
            "EXPFIT_K": float(k),
            "EXPFIT_RE_ARCSEC": float(re) if np.isfinite(re) else np.nan,
            "EXPFIT_OK": True,
        }
    except Exception:
        return {
            "EXPFIT_I0": np.nan,
            "EXPFIT_K": np.nan,
            "EXPFIT_RE_ARCSEC": np.nan,
            "EXPFIT_OK": False,
        }


def fit_log_linear_profile(radius, sb, sb_snr=None, min_snr=2.0) -> dict[str, float]:
    """
    Fit log10(I) as a linear function of radius.
    """
    radius = _as_float_array(radius)
    sb = _as_float_array(sb)

    good = np.isfinite(radius) & np.isfinite(sb) & (sb > 0)
    if sb_snr is not None:
        sb_snr = _as_float_array(sb_snr)
        good &= np.isfinite(sb_snr) & (sb_snr > min_snr)

    if np.sum(good) < 4:
        return {
            "LOGFIT_SLOPE": np.nan,
            "LOGFIT_INTERCEPT": np.nan,
            "LOGFIT_RE_ARCSEC": np.nan,
            "LOGFIT_OK": False,
        }

    try:
        coeff = np.polyfit(radius[good], np.log10(sb[good]), 1)
        slope, intercept = coeff
        re = np.nan if (not np.isfinite(slope) or slope == 0) else -1.0 / slope
        return {
            "LOGFIT_SLOPE": float(slope),
            "LOGFIT_INTERCEPT": float(intercept),
            "LOGFIT_RE_ARCSEC": float(re) if np.isfinite(re) else np.nan,
            "LOGFIT_OK": True,
        }
    except Exception:
        return {
            "LOGFIT_SLOPE": np.nan,
            "LOGFIT_INTERCEPT": np.nan,
            "LOGFIT_RE_ARCSEC": np.nan,
            "LOGFIT_OK": False,
        }


# ----------------------------------------------------------------------
# Petrosian measurements
# ----------------------------------------------------------------------

def petrosian_summary(
    radius_arcsec,
    flux_cum,
    flux_cgs=None,
    flux_cgs_err=None,
    magzp=None,
) -> dict[str, float]:
    """
    Compute Petrosian summary quantities from cumulative flux profile.
    """
    r = _as_float_array(radius_arcsec)
    f = _as_float_array(flux_cum)

    good = np.isfinite(r) & np.isfinite(f) & (r > 0)
    if np.sum(good) < 4:
        return {
            "PETRO_RAD_ARCSEC": np.nan,
            "PETRO_FLUX": np.nan,
            "PETRO_FLUX_CGS": np.nan,
            "PETRO_FLUX_CGS_ERR": np.nan,
            "PETRO_MAG": np.nan,
            "PETRO_R50_ARCSEC": np.nan,
            "PETRO_R90_ARCSEC": np.nan,
            "PETRO_CON": np.nan,
            "PETRO_OK": False,
        }

    # Sort the main working arrays by radius
    r = r[good]
    f = f[good]

    r_order = np.argsort(r)
    r = r[r_order]
    f = f[r_order]

    mean_sb_to_r = f / (np.pi * r**2)

    fl = np.array([_interp_monotonic(0.8 * rr, r, f) for rr in r], dtype=float)
    fu = np.array([_interp_monotonic(1.25 * rr, r, f) for rr in r], dtype=float)

    ann_area = np.pi * (1.25**2 - 0.8**2) * r**2
    ann_sb = (fu - fl) / ann_area
    petro_ratio = ann_sb / mean_sb_to_r

    valid = np.isfinite(petro_ratio)
    if np.sum(valid) < 2 or np.nanmax(petro_ratio[valid]) < 0.2:
        return {
            "PETRO_RAD_ARCSEC": np.nan,
            "PETRO_FLUX": np.nan,
            "PETRO_FLUX_CGS": np.nan,
            "PETRO_FLUX_CGS_ERR": np.nan,
            "PETRO_MAG": np.nan,
            "PETRO_R50_ARCSEC": np.nan,
            "PETRO_R90_ARCSEC": np.nan,
            "PETRO_CON": np.nan,
            "PETRO_OK": False,
        }

    rr = r[valid]
    pr = petro_ratio[valid]

    pr_order = np.argsort(pr)
    pr = pr[pr_order]
    rr = rr[pr_order]

    pr_u, idx = np.unique(pr, return_index=True)
    rr_u = rr[idx]

    if len(pr_u) < 2 or 0.2 < pr_u[0] or 0.2 > pr_u[-1]:
        petro_rad = np.nan
    else:
        petro_rad = float(np.interp(0.2, pr_u, rr_u))

    if not np.isfinite(petro_rad):
        return {
            "PETRO_RAD_ARCSEC": np.nan,
            "PETRO_FLUX": np.nan,
            "PETRO_FLUX_CGS": np.nan,
            "PETRO_FLUX_CGS_ERR": np.nan,
            "PETRO_MAG": np.nan,
            "PETRO_R50_ARCSEC": np.nan,
            "PETRO_R90_ARCSEC": np.nan,
            "PETRO_CON": np.nan,
            "PETRO_OK": False,
        }

    two_rp = 2.0 * petro_rad

    petro_flux = _interp_monotonic(two_rp, r, f)
    if not np.isfinite(petro_flux):
        petro_flux = np.nanmax(f)

    petro_flux_cgs = np.nan
    petro_flux_cgs_err = np.nan

    if flux_cgs is not None:
        fcgs = _as_float_array(flux_cgs)[good]
        fcgs = fcgs[r_order]

        petro_flux_cgs = _interp_monotonic(two_rp, r, fcgs)
        if not np.isfinite(petro_flux_cgs):
            petro_flux_cgs = np.nanmax(fcgs)

        if flux_cgs_err is not None:
            fcgs_err = _as_float_array(flux_cgs_err)[good]
            fcgs_err = fcgs_err[r_order]
            petro_flux_cgs_err = _interp_monotonic(two_rp, r, fcgs_err)

    petro_mag = safe_mag_from_flux(petro_flux, magzp) if magzp is not None else np.nan

    totalfrac = f / petro_flux if np.isfinite(petro_flux) and petro_flux > 0 else np.full_like(f, np.nan)
    petro_r50 = _interp_monotonic(0.5, totalfrac, r)
    petro_r90 = _interp_monotonic(0.9, totalfrac, r)
    petro_con = _safe_ratio(petro_r90, petro_r50)

    return {
        "PETRO_RAD_ARCSEC": float(petro_rad),
        "PETRO_FLUX": float(petro_flux) if np.isfinite(petro_flux) else np.nan,
        "PETRO_FLUX_CGS": float(petro_flux_cgs) if np.isfinite(petro_flux_cgs) else np.nan,
        "PETRO_FLUX_CGS_ERR": float(petro_flux_cgs_err) if np.isfinite(petro_flux_cgs_err) else np.nan,
        "PETRO_MAG": float(petro_mag) if np.isfinite(petro_mag) else np.nan,
        "PETRO_R50_ARCSEC": float(petro_r50) if np.isfinite(petro_r50) else np.nan,
        "PETRO_R90_ARCSEC": float(petro_r90) if np.isfinite(petro_r90) else np.nan,
        "PETRO_CON": float(petro_con) if np.isfinite(petro_con) else np.nan,
        "PETRO_OK": True,
    }


# ----------------------------------------------------------------------
# R-band profile summary
# ----------------------------------------------------------------------

def summarize_r_profile(tab: Table, magzp: float | None = None) -> dict[str, Any]:
    """
    Derive one-row R-band summary quantities from a photometry table.
    """
    results: dict[str, Any] = {}

    sma_arcsec = _as_float_array(tab["sma_arcsec"])
    sma_pix = _as_float_array(tab["sma_pix"])
    flux_cum = _as_float_array(tab["flux_cum"])
    flux_cum_err = _as_float_array(tab["flux_cum_err"])
    flux_cgs = _as_float_array(tab["flux_cgs"])
    flux_cgs_err = _as_float_array(tab["flux_cgs_err"])
    mag_cum = _as_float_array(tab["mag_cum"])
    mag_cum_err = _as_float_array(tab["mag_cum_err"])
    sb_mag = _as_float_array(tab["sb_mag_arcsec2"])
    sb_snr = _as_float_array(tab["sb_avg_snr"])
    sb_cgs = _as_float_array(tab["sb_cgs_arcsec2"])
    masked_fraction = _as_float_array(tab["masked_fraction"]) if "masked_fraction" in tab.colnames else None

    good_profile = np.isfinite(sb_snr) & (sb_snr > 2)
    results["R_PROFILE_OK"] = bool(np.sum(good_profile) >= 4)
    results["R_PROFILE_NGOOD"] = int(np.sum(good_profile))

    if masked_fraction is not None and np.any(np.isfinite(masked_fraction)):
        results["R_PROFILE_MASKFRAC_MAX"] = float(np.nanmax(masked_fraction))
    else:
        results["R_PROFILE_MASKFRAC_MAX"] = np.nan

    # flux radii
    for frac, label in [(0.25, "R25"), (0.50, "R50"), (0.75, "R75")]:
        r_arcsec, _ = _flux_radius(sma_arcsec, flux_cum, frac)
        r_pix, _ = _flux_radius(sma_pix, flux_cum, frac)
        results[f"{label}_ARCSEC"] = float(r_arcsec) if np.isfinite(r_arcsec) else np.nan
        results[f"{label}_PIX"] = float(r_pix) if np.isfinite(r_pix) else np.nan

    # isophotal radii/mags
    for iso, label in [
        (24.0, "R24"),
        (25.0, "R25_ISO"),
        (25.5, "R25P5"),
        (24.16, "R24_VEGA"),
        (25.16, "R25_VEGA"),
    ]:
        rval, rerr = _crossing_radius(sma_arcsec, sb_mag, iso, decreasing=False, mode="first")
        results[f"{label}_ARCSEC"] = float(rval) if np.isfinite(rval) else np.nan
        results[f"{label}_ARCSEC_ERR"] = float(rerr) if np.isfinite(rerr) else np.nan

        a, b = _crossing_bracket_indices(sb_mag, iso, decreasing=False, mode="first")
        if a is None or b is None:
            results[f"{label}_MAG"] = np.nan
            results[f"{label}_MAG_ERR"] = np.nan
        else:
            results[f"{label}_MAG"] = _combine_two_mags(mag_cum[a], mag_cum[b])
            errs = [mag_cum_err[a], mag_cum_err[b]]
            errs = [e for e in errs if np.isfinite(e)]
            results[f"{label}_MAG_ERR"] = float(np.nanmax(errs)) if errs else np.nan

    # C30 based on R24
    r24 = results.get("R24_ARCSEC", np.nan)
    if np.isfinite(r24):
        inner = _interp_monotonic(0.3 * r24, sma_arcsec, flux_cgs)
        outer = _interp_monotonic(r24, sma_arcsec, flux_cgs)
        inner_err = _interp_monotonic(0.3 * r24, sma_arcsec, flux_cgs_err)
        outer_err = _interp_monotonic(r24, sma_arcsec, flux_cgs_err)

        results["R30R24_FLUX_CGS"] = float(inner) if np.isfinite(inner) else np.nan
        results["R30R24_FLUX_CGS_ERR"] = float(inner_err) if np.isfinite(inner_err) else np.nan
        results["R24_FLUX_CGS"] = float(outer) if np.isfinite(outer) else np.nan
        results["R24_FLUX_CGS_ERR"] = float(outer_err) if np.isfinite(outer_err) else np.nan

        c30 = _safe_ratio(inner, outer)
        c30_err = ratio_error(inner, outer, inner_err, outer_err)
        results["R_C30"] = float(c30) if np.isfinite(c30) else np.nan
        results["R_C30_ERR"] = float(c30_err) if np.isfinite(c30_err) else np.nan
    else:
        results["R30R24_FLUX_CGS"] = np.nan
        results["R30R24_FLUX_CGS_ERR"] = np.nan
        results["R24_FLUX_CGS"] = np.nan
        results["R24_FLUX_CGS_ERR"] = np.nan
        results["R_C30"] = np.nan
        results["R_C30_ERR"] = np.nan

    petro = petrosian_summary(
        sma_arcsec,
        flux_cum,
        flux_cgs=flux_cgs,
        flux_cgs_err=flux_cgs_err,
        magzp=magzp,
    )
    for key, val in petro.items():
        results[f"R_{key}"] = val

    expfit = fit_exponential_profile(sma_arcsec, sb_cgs, sb_snr=sb_snr, min_snr=2.0)
    logfit = fit_log_linear_profile(sma_arcsec, sb_cgs, sb_snr=sb_snr, min_snr=2.0)
    for key, val in expfit.items():
        results[f"R_{key}"] = val
    for key, val in logfit.items():
        results[f"R_{key}"] = val

    if np.sum(sb_snr > 1.5) > 2:
        idx = int(np.max(np.where(sb_snr > 1.5)))
        results["R_TOT_MAG_SNR"] = float(mag_cum[idx]) if np.isfinite(mag_cum[idx]) else np.nan
        results["R_TOT_FLUX_CGS_SNR"] = float(flux_cgs[idx]) if np.isfinite(flux_cgs[idx]) else np.nan
        results["R_TOT_FLUX_CGS_SNR_ERR"] = float(flux_cgs_err[idx]) if np.isfinite(flux_cgs_err[idx]) else np.nan
        results["R_SNR_TRUNC_ARCSEC"] = float(sma_arcsec[idx]) if np.isfinite(sma_arcsec[idx]) else np.nan
    else:
        results["R_TOT_MAG_SNR"] = np.nan
        results["R_TOT_FLUX_CGS_SNR"] = np.nan
        results["R_TOT_FLUX_CGS_SNR_ERR"] = np.nan
        results["R_SNR_TRUNC_ARCSEC"] = np.nan

    return results


# ----------------------------------------------------------------------
# Halpha profile summary
# ----------------------------------------------------------------------

def summarize_ha_profile(
    tab: Table,
    r24_arcsec: float | None = None,
    min_snr: float = 2.0,
    min_run: int = 2,
) -> dict[str, Any]:
    """
    Derive one-row Halpha summary quantities from a photometry table.

    Halpha profiles are allowed to be non-monotonic and ring-like. Significant
    radial structure is identified using contiguous runs of annuli with:
        sb_avg_snr > min_snr
        sb_cgs_arcsec2 > 0

    Outer Halpha extent is defined from the outermost significant run.
    """
    results: dict[str, Any] = {}

    sma_arcsec = _as_float_array(tab["sma_arcsec"])
    sma_pix = _as_float_array(tab["sma_pix"])
    flux_cum = _as_float_array(tab["flux_cum"])
    flux_cum_err = _as_float_array(tab["flux_cum_err"])
    flux_cgs = _as_float_array(tab["flux_cgs"])
    flux_cgs_err = _as_float_array(tab["flux_cgs_err"])
    sb_cgs = _as_float_array(tab["sb_cgs_arcsec2"])
    sb_cgs_err = _as_float_array(tab["sb_cgs_arcsec2_err"])
    sb_snr = _as_float_array(tab["sb_avg_snr"])
    masked_fraction = _as_float_array(tab["masked_fraction"]) if "masked_fraction" in tab.colnames else None

    raw_detect = np.isfinite(sb_snr) & np.isfinite(sb_cgs) & (sb_snr > min_snr) & (sb_cgs > 0)
    runs = _significant_runs(sb_snr, sb_cgs, min_snr=min_snr, min_run=min_run)

    results["H_PROFILE_NGOOD"] = int(np.sum(raw_detect))
    results["H_PROFILE_LONGRUN"] = _longest_true_run(raw_detect)
    results["H_NDET_RUNS"] = int(len(runs))
    results["H_PROFILE_OK"] = bool(len(runs) > 0)

    if masked_fraction is not None and np.any(np.isfinite(masked_fraction)):
        results["H_PROFILE_MASKFRAC_MAX"] = float(np.nanmax(masked_fraction))
    else:
        results["H_PROFILE_MASKFRAC_MAX"] = np.nan

    if len(runs) == 0:
        results["H_MAXDET_ARCSEC"] = np.nan
        results["H_MAXDET_PIX"] = np.nan
        results["H_TOT_FLUX_CGS"] = np.nan
        results["H_TOT_FLUX_CGS_ERR"] = np.nan
        results["H_SNR_TRUNC_ARCSEC"] = np.nan

        # no trustworthy profile structure
        results["H25_ARCSEC"] = np.nan
        results["H25_PIX"] = np.nan
        results["H50_ARCSEC"] = np.nan
        results["H50_PIX"] = np.nan
        results["H75_ARCSEC"] = np.nan
        results["H75_PIX"] = np.nan

        for prefix in ["H_ISO5E17", "H_ISO17E18"]:
            results[f"{prefix}_ARCSEC"] = np.nan
            results[f"{prefix}_ARCSEC_ERR"] = np.nan
            results[f"{prefix}_FLUX_CGS"] = np.nan
            results[f"{prefix}_FLUX_CGS_ERR"] = np.nan

        results["H30R24_FLUX_CGS"] = np.nan
        results["H30R24_FLUX_CGS_ERR"] = np.nan
        results["H_R24_FLUX_CGS"] = np.nan
        results["H_R24_FLUX_CGS_ERR"] = np.nan
        results["H_C30_R24"] = np.nan
        results["H_C30_R24_ERR"] = np.nan
        results["H_R95_R24_ARCSEC"] = np.nan

        # fits/petrosian also not trustworthy
        for prefix in [
            "H_PETRO_RAD_ARCSEC", "H_PETRO_FLUX", "H_PETRO_FLUX_CGS",
            "H_PETRO_FLUX_CGS_ERR", "H_PETRO_MAG", "H_PETRO_R50_ARCSEC",
            "H_PETRO_R90_ARCSEC", "H_PETRO_CON"
        ]:
            results[prefix] = np.nan
        results["H_PETRO_OK"] = False

        for prefix in [
            "H_EXPFIT_I0", "H_EXPFIT_K", "H_EXPFIT_RE_ARCSEC",
            "H_LOGFIT_SLOPE", "H_LOGFIT_INTERCEPT", "H_LOGFIT_RE_ARCSEC"
        ]:
            results[prefix] = np.nan
        results["H_EXPFIT_OK"] = False
        results["H_LOGFIT_OK"] = False

        return results

    # Outer edge of last significant run
    last_i0, last_i1 = runs[-1]
    results["H_MAXDET_ARCSEC"] = float(sma_arcsec[last_i1]) if np.isfinite(sma_arcsec[last_i1]) else np.nan
    results["H_MAXDET_PIX"] = float(sma_pix[last_i1]) if np.isfinite(sma_pix[last_i1]) else np.nan

    # Integrated Halpha quantities: use cumulative profile evaluated at outermost significant radius
    maxdet_r = results["H_MAXDET_ARCSEC"]
    results["H_TOT_FLUX_CGS"] = _interp_monotonic(maxdet_r, sma_arcsec, flux_cgs)
    results["H_TOT_FLUX_CGS_ERR"] = _interp_monotonic(maxdet_r, sma_arcsec, flux_cgs_err)
    results["H_SNR_TRUNC_ARCSEC"] = maxdet_r

    # Flux radii relative to cumulative profile out to outermost detected extent
    tot_flux = results["H_TOT_FLUX_CGS"]
    if np.isfinite(tot_flux) and tot_flux > 0:
        for frac, label in [(0.25, "H25"), (0.50, "H50"), (0.75, "H75")]:
            target = frac * tot_flux
            r_arcsec = _interp_monotonic(target, flux_cgs, sma_arcsec)
            r_pix = _interp_monotonic(target, flux_cgs, sma_pix)
            results[f"{label}_ARCSEC"] = float(r_arcsec) if np.isfinite(r_arcsec) else np.nan
            results[f"{label}_PIX"] = float(r_pix) if np.isfinite(r_pix) else np.nan
    else:
        for label in ["H25", "H50", "H75"]:
            results[f"{label}_ARCSEC"] = np.nan
            results[f"{label}_PIX"] = np.nan

    # Isophotal Halpha radii: use outermost crossing
    iso_map = {
        "H_ISO5E17": 5.0e-17,
        "H_ISO17E18": 1.7e-17,
    }

    for prefix, level in iso_map.items():
        rval, rerr = _crossing_radius(sma_arcsec, sb_cgs, level, decreasing=True, mode="last")
        results[f"{prefix}_ARCSEC"] = float(rval) if np.isfinite(rval) else np.nan
        results[f"{prefix}_ARCSEC_ERR"] = float(rerr) if np.isfinite(rerr) else np.nan

        if np.isfinite(rval):
            fval = _interp_monotonic(rval, sma_arcsec, flux_cgs)
            ferr = _interp_monotonic(rval, sma_arcsec, flux_cgs_err)
            results[f"{prefix}_FLUX_CGS"] = float(fval) if np.isfinite(fval) else np.nan
            results[f"{prefix}_FLUX_CGS_ERR"] = float(ferr) if np.isfinite(ferr) else np.nan
        else:
            results[f"{prefix}_FLUX_CGS"] = np.nan
            results[f"{prefix}_FLUX_CGS_ERR"] = np.nan

    # R24-based Halpha quantities
    if r24_arcsec is not None and np.isfinite(r24_arcsec):
        inner = _interp_monotonic(0.3 * r24_arcsec, sma_arcsec, flux_cgs)
        outer = _interp_monotonic(r24_arcsec, sma_arcsec, flux_cgs)
        inner_err = _interp_monotonic(0.3 * r24_arcsec, sma_arcsec, flux_cgs_err)
        outer_err = _interp_monotonic(r24_arcsec, sma_arcsec, flux_cgs_err)

        results["H30R24_FLUX_CGS"] = float(inner) if np.isfinite(inner) else np.nan
        results["H30R24_FLUX_CGS_ERR"] = float(inner_err) if np.isfinite(inner_err) else np.nan
        results["H_R24_FLUX_CGS"] = float(outer) if np.isfinite(outer) else np.nan
        results["H_R24_FLUX_CGS_ERR"] = float(outer_err) if np.isfinite(outer_err) else np.nan

        c30 = _safe_ratio(inner, outer)
        c30_err = ratio_error(inner, outer, inner_err, outer_err)
        results["H_C30_R24"] = float(c30) if np.isfinite(c30) else np.nan
        results["H_C30_R24_ERR"] = float(c30_err) if np.isfinite(c30_err) else np.nan

        if np.isfinite(outer) and outer > 0:
            r95 = _interp_monotonic(0.95 * outer, flux_cgs, sma_arcsec)
            results["H_R95_R24_ARCSEC"] = float(r95) if np.isfinite(r95) else np.nan
        else:
            results["H_R95_R24_ARCSEC"] = np.nan
    else:
        results["H30R24_FLUX_CGS"] = np.nan
        results["H30R24_FLUX_CGS_ERR"] = np.nan
        results["H_R24_FLUX_CGS"] = np.nan
        results["H_R24_FLUX_CGS_ERR"] = np.nan
        results["H_C30_R24"] = np.nan
        results["H_C30_R24_ERR"] = np.nan
        results["H_R95_R24_ARCSEC"] = np.nan

    # Petrosian quantities are not especially natural for multi-run Halpha, but keep them
    # for simple cases only
    if len(runs) == 1:
        petro = petrosian_summary(
            sma_arcsec,
            flux_cum,
            flux_cgs=flux_cgs,
            flux_cgs_err=flux_cgs_err,
            magzp=None,
        )
        for key, val in petro.items():
            results[f"H_{key}"] = val
    else:
        results["H_PETRO_RAD_ARCSEC"] = np.nan
        results["H_PETRO_FLUX"] = np.nan
        results["H_PETRO_FLUX_CGS"] = np.nan
        results["H_PETRO_FLUX_CGS_ERR"] = np.nan
        results["H_PETRO_MAG"] = np.nan
        results["H_PETRO_R50_ARCSEC"] = np.nan
        results["H_PETRO_R90_ARCSEC"] = np.nan
        results["H_PETRO_CON"] = np.nan
        results["H_PETRO_OK"] = False

    # Fits only for simple single-run profiles
    if len(runs) == 1:
        expfit = fit_exponential_profile(sma_arcsec, sb_cgs, sb_snr=sb_snr, min_snr=min_snr)
        logfit = fit_log_linear_profile(sma_arcsec, sb_cgs, sb_snr=sb_snr, min_snr=min_snr)
        for key, val in expfit.items():
            results[f"H_{key}"] = val
        for key, val in logfit.items():
            results[f"H_{key}"] = val
    else:
        results["H_EXPFIT_I0"] = np.nan
        results["H_EXPFIT_K"] = np.nan
        results["H_EXPFIT_RE_ARCSEC"] = np.nan
        results["H_EXPFIT_OK"] = False

        results["H_LOGFIT_SLOPE"] = np.nan
        results["H_LOGFIT_INTERCEPT"] = np.nan
        results["H_LOGFIT_RE_ARCSEC"] = np.nan
        results["H_LOGFIT_OK"] = False

    return results


# ----------------------------------------------------------------------
# convenience wrapper
# ----------------------------------------------------------------------

def summarize_dual_profiles(
    rtab: Table | None = None,
    hatab: Table | None = None,
    r_magzp: float | None = None,
) -> dict[str, Any]:
    """
    Summarize R and Halpha photometry tables together.
    """
    results: dict[str, Any] = {}

    r24_arcsec = np.nan
    if rtab is not None:
        rres = summarize_r_profile(rtab, magzp=r_magzp)
        results.update(rres)
        r24_arcsec = rres.get("R24_ARCSEC", np.nan)

    if hatab is not None:
        hares = summarize_ha_profile(hatab, r24_arcsec=r24_arcsec)
        results.update(hares)

    return results
