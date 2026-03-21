#!/usr/bin/env python

"""
astro.py

Astrophysical helper functions used by HAPY utilities.
"""

from __future__ import annotations

import numpy as np
from astropy.cosmology import WMAP9 as cosmo


def KE_SFR_from_redshift(haflux, redshift):
    """
    Convert Halpha flux to log10(SFR / Msun yr^-1)
    using the Kennicutt & Evans (2012) calibration.

    Parameters
    ----------
    haflux : float or array-like
        Halpha flux in erg s^-1 cm^-2

    redshift : float or array-like
        Redshift

    Returns
    -------
    logsfr : float or ndarray
        log10(SFR / Msun yr^-1)

    Notes
    -----
    Uses:
        log10(SFR / Msun yr^-1) = log10(L_Ha / erg s^-1) - 41.27

    Returns NaN where inputs are invalid or haflux <= 0.
    """
    haflux = np.asarray(haflux, dtype=float)
    redshift = np.asarray(redshift, dtype=float)

    L = haflux * (4.0 * np.pi * cosmo.luminosity_distance(redshift).cgs.value**2)

    logsfr = np.full(np.shape(L), np.nan, dtype=float)
    good = np.isfinite(L) & (L > 0)
    logsfr[good] = np.log10(L[good]) - 41.27

    if logsfr.ndim == 0:
        return float(logsfr)

    return logsfr


def KE_SFR_from_distance(haflux, dist_mpc):
    import numpy as np

    haflux = np.asarray(haflux, dtype=float)
    dist_mpc = np.asarray(dist_mpc, dtype=float)

    d_cm = dist_mpc * 3.085677581e24
    L = haflux * (4.0 * np.pi * d_cm**2)

    logsfr = np.full(np.shape(L), np.nan, dtype=float)
    good = np.isfinite(L) & (L > 0)
    logsfr[good] = np.log10(L[good]) - 41.27

    if logsfr.ndim == 0:
        return float(logsfr)

    return logsfr
