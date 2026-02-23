from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Any

import numpy as np

# import statmorph ONLY here
import statmorph
from statmorph.utils.image_diagnostics import make_figure

from astropy.utils import lazyproperty
import scipy.ndimage as ndi


class MyStatmorph(statmorph.SourceMorphology):
    """Statmorph subclass that forces gini segmap behavior."""

    @lazyproperty
    def _segmap_gini(self):
        segmap = np.array(self._segmap.data == 1, "i")
        return segmap[self._slice_stamp]

    def print(self):
        ''' adding a print method to print out the instance variables '''
        for k in self.__dict__.keys():
            if k.startswith('_'):
                continue
            print(f"{k}: {self.__dict__[k]}")
            
        

@dataclass
class MorphologyResult:
    morph_r: Any
    morph_img2: Optional[Any] = None
    fig_r: Optional[Any] = None
    fig_img2: Optional[Any] = None


def make_object_segmap(segmentation_data: np.ndarray, object_label: int, smooth_size: int = 10) -> np.ndarray:
    """Return a 0/1 segmap for the central object only."""
    segmap = segmentation_data == object_label
    segmap_float = ndi.uniform_filter(np.float64(segmap), size=smooth_size)
    return np.array(segmap_float > 0.5, "i")


def run_statmorph_single(
    image: np.ndarray,
    segmap: np.ndarray,
    gain: float = 1.0,
    mask: Optional[np.ndarray] = None,
    psf: Optional[np.ndarray] = None,
    cutout_extent: float = 1.5,
    make_fig: bool = True,
) -> tuple[Any, Optional[Any]]:
    """
    Run statmorph on a single image and segmap (0/1).
    Returns (morph, fig_or_None).
    """
    label = 1  # statmorph expects an integer label

    if psf is None:
        morph = MyStatmorph(image, segmap, label, gain=gain, mask=mask, cutout_extent=cutout_extent)
    else:
        morph = MyStatmorph(image, segmap, label, gain=gain, mask=mask, psf=psf, cutout_extent=cutout_extent)

    fig = make_figure(morph) if make_fig else None
    return morph, fig


def run_statmorph_for_photometry(
    *,
    image: np.ndarray,
    segmentation_data: np.ndarray,
    object_label: int,
    gain: float = 1.0,
    mask: Optional[np.ndarray] = None,
    psf: Optional[np.ndarray] = None,
    image2: Optional[np.ndarray] = None,
    psf2: Optional[np.ndarray] = None,
    make_fig: bool = True,
) -> MorphologyResult:
    """
    High-level helper: build single-object segmap and run statmorph on image and optional image2.
    """
    segmap = make_object_segmap(segmentation_data, object_label)

    morph_r, fig_r = run_statmorph_single(
        image=image, segmap=segmap, gain=gain, mask=mask, psf=psf, make_fig=make_fig
    )

    morph2 = fig2 = None
    if image2 is not None:
        morph2, fig2 = run_statmorph_single(
            image=image2, segmap=segmap, gain=gain, mask=mask, psf=psf2, make_fig=make_fig
        )

    return MorphologyResult(morph_r=morph_r, morph_img2=morph2, fig_r=fig_r, fig_img2=fig2)


def m20_from_sourcecatalog(cat, idx):
    """
    Compute M20 for a single object in a photutils SourceCatalog.

    Uses only pixels in the object's segment, ignoring masked pixels.
    """
    label = cat.label[idx]
    seg = (cat.segment[idx] == label)

    data_ma = cat.data_ma[idx]          # MaskedArray
    data = np.asarray(data_ma.filled(0.0), dtype=float)
    # exclude masked pixels
    good = seg & ~data_ma.mask

    if not np.any(good):
        return np.nan

    # centroid in cutout coordinates
    xc, yc = cat.cutout_centroid[idx]

    yy, xx = np.indices(data.shape)
    r2 = (xx - xc) ** 2 + (yy - yc) ** 2

    flux = data[good]
    r2_good = r2[good]

    # If flux can be negative (sky-subtracted), M20 becomes ill-defined.
    # Common approach: restrict to positive-flux pixels.
    pos = flux > 0
    if not np.any(pos):
        return np.nan

    flux = flux[pos]
    r2_good = r2_good[pos]

    mtot = np.sum(flux * r2_good)
    if mtot <= 0:
        return np.nan

    # Find the subset of pixels containing the brightest 20% of the total flux.
    order = np.argsort(flux)[::-1]  # descending
    flux_sorted = flux[order]
    r2_sorted = r2_good[order]

    cumsum = np.cumsum(flux_sorted)
    frac = cumsum / cumsum[-1]

    top20 = frac <= 0.20
    # ensure at least one pixel
    if not np.any(top20):
        top20[np.argmax(flux_sorted)] = True

    m20 = np.log10(np.sum(flux_sorted[top20] * r2_sorted[top20]) / mtot)
    return float(m20)
