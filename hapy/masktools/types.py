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
