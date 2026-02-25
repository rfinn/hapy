from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import numpy as np

from hapy.galfittools.results import GalfitResult


@dataclass(frozen=True, slots=True)
class EllipseGeometry:
    """Geometry in photutils convention: PA deg CCW from +x axis."""
    xc: float
    yc: float
    ba: float
    pa_deg: float
    sky: Optional[float] = None


def galfit_pa_to_photutils_pa(pa_galfit_deg: float) -> float:
    """
    Convert GALFIT PA (deg CCW from +y) -> photutils PA (deg CCW from +x).
    """
    return float((90.0 - pa_galfit_deg) % 360.0)


def geometry_from_galfit(
    res: GalfitResult,
    component: int = 1,
    clamp_ba: bool = True,
    min_ba: float = 0.05,
) -> EllipseGeometry:
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

    return EllipseGeometry(
        xc=float(c.xc),
        yc=float(c.yc),
        ba=ba,
        pa_deg=galfit_pa_to_photutils_pa(float(c.pa)),
        sky=float(res.sky) if hasattr(res, "sky") else None,
    )


def galfit_pa_to_photutils_pa(pa_galfit_deg: float) -> float:
    # +y -> +x
    return (90.0 - pa_galfit_deg) % 360.0

def photutils_pa_to_galfit_pa(pa_phot_deg: float) -> float:
    # +x -> +y
    return (90.0 - pa_phot_deg) % 360.0
