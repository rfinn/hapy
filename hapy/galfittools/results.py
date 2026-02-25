from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True, slots=True)
class GalfitComponent:
    xc: float
    yc: float
    mag: float
    re: float
    n: float
    ba: float
    pa: float

    xc_err: float = 0.0
    yc_err: float = 0.0
    mag_err: float = 0.0
    re_err: float = 0.0
    n_err: float = 0.0
    ba_err: float = 0.0
    pa_err: float = 0.0

    numerical_error_flag: int = 0


@dataclass(frozen=True, slots=True)
class GalfitResult:
    ncomp: int
    comp1: GalfitComponent
    comp2: Optional[GalfitComponent] = None

    sky: float = 0.0
    sky_err: float = 0.0

    chi2nu: float = 0.0
    error: float = 0.0  # GALFIT's "ERROR" keyword if present/parsed

    # optional asymmetry (Fourier mode)
    f1: Optional[float] = None
    f1_err: Optional[float] = None
    f1_pa: Optional[float] = None
    f1_pa_err: Optional[float] = None
