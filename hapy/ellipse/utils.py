from photutils.segmentation import detect_threshold, detect_sources, SourceCatalog
from astropy.stats import SigmaClip

import numpy as np

def ellipse_geom_missing(params):
    """True if we cannot build an ellipse at all (size/shape missing)."""
    sma = params.get("sma_arcsec")
    ba  = params.get("ba")
    def bad(x): return x is None or (isinstance(x, (int, float)) and (not np.isfinite(x) or x <= 0))
    return bad(sma) or bad(ba)

def pa_is_placeholder(params):
    """True if PA is present but likely a placeholder (AGC default)."""
    scheme = (params.get("scheme") or "").lower()
    pa = params.get("pa_deg")
    if pa is None or not isinstance(pa, (int, float)) or not np.isfinite(pa):
        return True
    return (scheme == "agc" and float(pa) == 0.0)

def infer_ellipse_from_r_cutout(r_data, r_wcs=None, nsigma=2.5, npixels=20):
    """
    Return an EllipseParams guess in pixel coords (theta_deg CCW from +x).
    Heuristic: choose the detection closest to image center.
    """
    ny, nx = r_data.shape
    xc0, yc0 = nx / 2.0, ny / 2.0

    base_mask = ~np.isfinite(r_data)
    sigma_clip = SigmaClip(sigma=3.0, maxiters=5)
    thr = detect_threshold(r_data, nsigma=nsigma, sigma_clip=sigma_clip, mask=base_mask)
    segm = detect_sources(r_data, thr, npixels=npixels, mask=base_mask)
    if segm is None:
        return None

    cat = SourceCatalog(r_data, segm, mask=base_mask)

    # pick object closest to cutout center
    x = np.array(cat.xcentroid)
    y = np.array(cat.ycentroid)
    d2 = (x - xc0)**2 + (y - yc0)**2
    i = int(np.nanargmin(d2))

    obj = cat[i]

    # shape
    ba = float(obj.semiminor_sigma.value / obj.semimajor_sigma.value) if obj.semimajor_sigma.value > 0 else 1.0
    ba = float(np.clip(ba, 0.05, 1.0))

    # size: use a multiple of sigma as a rough SMA (tune factor as needed)
    sma_pix = float(obj.semimajor_sigma.value * 3.0)  # 3-sigma as a starting ellipse
    sma_pix = max(sma_pix, 5.0)

    # orientation: photutils orientation is radians CCW from +x (already your EllipseParams theta)
    theta_deg = float(np.degrees(obj.orientation.to_value(u.rad)) % 180.0)

    return EllipseParams(
        xc=float(obj.xcentroid.value),
        yc=float(obj.ycentroid.value),
        sma_pix=sma_pix,
        ba=ba,
        theta_deg=theta_deg
    )
