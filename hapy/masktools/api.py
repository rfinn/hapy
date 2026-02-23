"""
Public API for the headless masking engine.

maskgui (and any future CLI/tests) should import ONLY from here,
not from masktools internals.
"""

from .engine import MaskEngine
from .types import EllipseParams

__all__ = ["MaskEngine", "EllipseParams"]
