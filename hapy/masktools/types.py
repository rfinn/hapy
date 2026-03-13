#from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class EllipseParams:
    """Ellipse in pixel coords.  theta_deg is angle CCW from +x axis"""
    xc: float
    yc: float
    sma_pix: float
    ba: float
    theta_deg: float # angle in deg from +x axis


@dataclass(frozen=True)
class MaskFractionResult:
    frac_masked: float
    n_total: int
    n_masked: int
    n_unmasked: int
    
def build_ell0_from_metadata(params, pixscale):
    # params: dict from metadata.json
    sma_pix = params["sma_arcsec"] / pixscale
    pa_north = params["pa_deg"]  # CCW from North
    theta_deg = pa_ccw_north_to_photutils_theta(pa_north)  # uses (90 + pa) % 180 for your axes
    return EllipseParams(
        xc=params.get("xc", None),
        yc=params.get("yc", None),
        sma_pix=sma_pix,
        ba=params["ba"],
        pa_deg=theta_deg,  # photutils theta deg CCW from +x
    )
