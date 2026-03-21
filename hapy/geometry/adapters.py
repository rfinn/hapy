from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import numpy as np

from hapy.galfittools.results import GalfitResult
from hapy.masktools.types import EllipseParams

#@dataclass(frozen=True, slots=True)
#class EllipseGeometry:
#    """Geometry in photutils convention: PA deg CCW from +x axis."""
#    xc: float
#    yc: float
#    ba: float
#    pa_deg: float
#    sky: Optional[float] = None



def geometry_from_galfit(
    res: GalfitResult,
    component: int = 1,
    clamp_ba: bool = True,
    min_ba: float = 0.05,
) -> EllipseParams:
    """
    Convert a GalfitResult into EllipseGeometry (photutils conventions).
    """
    if component == 1:
        c = res.comp1
    elif component == 2:
        if res.comp2 is None:
            raise ValueError("Requested component=2 but res.comp2 is None")
        c = res.comp2
    else:
        raise ValueError("component must be 1 or 2")

    ba = float(c.ba)
    if clamp_ba:
        ba = float(np.clip(ba, min_ba, 1.0))

    return EllipseParams(
        xc=float(c.xc),
        yc=float(c.yc),
        ba=ba,
        theta_deg=pa_ccw_north_to_photutils_theta(float(c.pa)),
        #sky=float(res.sky) if hasattr(res, "sky") else None,
    )


def pa_ccw_north_deg_to_photutils_theta_rad(pa_deg):
    """
    Convert standard internal PA convention to photutils theta.

    Internal convention:
      PA_DEG = CCW from North (+y), degrees, periodic over 180 deg. Assumes pixel +x is West and +y is North.

    Photutils convention:
      theta = radians CCW from +x axis.

    Assumes image axes: +x=East, +y=North.
    """
    if pa_deg is None:
        return None
    try:
        pa = float(pa_deg)
    except Exception:
        return None
    return np.deg2rad(pa_ccw_north_to_photutils_theta(pa))

def pa_ccw_north_to_photutils_theta(pa_deg: float) -> float:
    """
    Convert internal PA_DEG (deg CCW from North/+y) to photutils theta (deg CCW from +x). 

    Assumes pixel +x is West and +y is North.

    Ellipse angles are 180-deg periodic.
    """

    return float((90.0 + pa_deg) % 180.0)


def photutils_theta_to_pa_ccw_north(theta_deg: float) -> float:
    """
    Convert photutils theta (deg CCW from +x) to internal PA_DEG (deg CCW from North/+y).
    """
    return float((90.0 + theta_deg) % 180.0)



