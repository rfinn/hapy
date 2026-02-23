#from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class EllipseParams:
    """Ellipse in pixel coords."""
    xc: float
    yc: float
    sma_pix: float
    ba: float
    pa_deg: float
