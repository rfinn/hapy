"""
H-alpha clump detection and measurement utilities for HAPY.

This module detects individual H-alpha star-forming clumps within the
R-band segmentation footprint of the central galaxy.

The intended workflow is:

1. Use the existing R-band segmentation image from EllipsePhotometry.
2. Select the central galaxy footprint using object_index.
3. Restrict the H-alpha image to that footprint.
4. Run photutils detect_sources/deblend_sources on the restricted H-alpha image.
5. Measure the resulting H-alpha clumps with SourceCatalog.
6. Save the clump catalog table and clump segmentation map for later analysis.
7. Return summary quantities for storage in the main HAPY catalog.


2026-06-27
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import inspect
import warnings

import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable

from photutils.segmentation import (
    SegmentationImage,
    SourceCatalog,
    detect_sources,
    deblend_sources,
)

from photutils.detection import find_peaks

from hapy.imagetools.imutils import calculate_background_photutils, get_bbox_from_mask, make_limits_square

__all__ = [
    "ClumpDetectionConfig",
    "ClumpSummary",
    "ClumpAnalysisResult",
    "make_object_footprint",
    "mask_halpha_to_footprint",
    "detect_halpha_clumps",
    "make_clump_catalog",
    "summarize_clumps",
    "save_clump_outputs",
    "analyze_halpha_clumps",
]


# ---------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------



@dataclass
class ClumpDetectionConfig:
    """
    Configuration parameters for H-alpha clump detection.
    """

    # Final H-alpha clump detection threshold:
    # threshold = local background median + nsigma * local background std
    nsigma: float = 3.0
    npixels: int = 5

    # Parameters passed to hapy.imagetools.imutils.calculate_background_photutils
    background_grow_radius: float = 10
    background_npixels: int = 10
    background_mask_nsigma: float = 2.0
    background_clip_sigma: float = 3.0
    background_clip_maxiters: int = 5

    # Deblending
    deblend: bool = True
    contrast: float = 0.001

    # Testing showed linear + 64 is better
    # e.g. VFID1802-CGCG218-049-BOK-20210418-VFID1821
    mode: str = "linear"
    nlevels: int = 64

    connectivity: int = 8
    nproc: int = 1
    relabel: bool = True
    progress_bar: bool = False

    # Optional cuts used only for summary statistics.
    # The full clump catalog is still saved.
    min_flux: float | None = None
    min_area: int | None = None

    # Local peak finding
    find_peaks: bool = True
    peak_box_size: int = 5
    peak_min_separation: int | None = None
    peak_nmax: int | None = None

    # Nuclear radius in multiples of the measured H-alpha FWHM.
    # Nuclear clump = clump centroid within this radius of the adopted galaxy center.
    nuclear_radius_fwhm: float = 1.5

    # Optional point-source-like detection
    find_point_sources: bool = False
    point_source_method: str = "dao"
    point_source_fwhm: float | None = None
    point_source_threshold_nsigma: float = 5.0

    # Output control
    save_catalog: bool = True
    save_segmentation: bool = True
    save_peaks: bool = True
    save_point_sources: bool = True
    save_diagnostic: bool = True

    table_format: str = "ascii.ecsv"
    diagnostic_format: str = "png"
    diagnostic_percent: float = 99.5
    plot_kron_apertures: bool = False

    
@dataclass
class ClumpSummary:
    """
    Compact summary of the H-alpha clump population.

    These are the quantities intended for storage in the main HAPY catalog.
    """    
    n_clumps: int = 0
    n_clumps_good: int = 0

    n_peaks: int = 0
    n_peaks_in_clumps: int = 0
    peak_max: float = np.nan
    peak_sum: float = np.nan

    n_point_sources: int = 0
    point_source_flux_sum: float = np.nan

    threshold: float = np.nan
    background_mean: float = np.nan
    background_median: float = np.nan
    background_rms: float = np.nan

    total_halpha_flux_footprint: float = np.nan
    total_clump_flux: float = np.nan
    clump_flux_fraction: float = np.nan

    brightest_clump_flux: float = np.nan
    brightest_clump_fraction: float = np.nan

    top2_clump_fraction: float = np.nan
    top3_clump_fraction: float = np.nan

    has_nuclear: bool = False
    n_nuclear: int = 0
    nuclear_flux_frac: float = np.nan

    bright_to_second_flux: float = np.nan

    clump_area_pix: int = 0
    footprint_area_pix: int = 0
    clump_area_fraction: float = np.nan

    flux_weighted_xcentroid: float = np.nan
    flux_weighted_ycentroid: float = np.nan

    brightest_xcentroid: float = np.nan
    brightest_ycentroid: float = np.nan

    catalog_path: str = ""
    segmentation_path: str = ""
    peaks_path: str = ""
    point_sources_path: str = ""
    diagnostic_path: str = ""

    # clump positions
    rmin_pix: float = np.nan
    rmax_pix: float = np.nan
    rmed_pix: float = np.nan
    rmean_pix: float = np.nan
    rfluxwt_pix: float = np.nan
    rbright_pix: float = np.nan

    dx_bright_pix: float = np.nan
    dy_bright_pix: float = np.nan
    

@dataclass
class ClumpAnalysisResult:
    """
    Container returned by analyze_halpha_clumps.

    Attributes
    ----------
    catalog : SourceCatalog or None
        In-memory photutils SourceCatalog for the H-alpha clumps.

    table : QTable
        Persistent table version of the clump catalog.

    segm : SegmentationImage or None
        Final H-alpha clump segmentation image, after deblending if requested.

    segm_detected : SegmentationImage or None
        Initial detected-source segmentation image before deblending.

    footprint_mask : ndarray
        Boolean mask selecting the central galaxy R-band segmentation footprint.

    detect_mask : ndarray
        Boolean mask used during H-alpha clump detection.

    summary : ClumpSummary
        Compact summary intended for the main HAPY catalog.
    """

    catalog: object | None
    table: QTable

    peak_table: QTable | None = None
    point_source_table: QTable | None = None

    segm: object | None = None
    segm_detected: object | None = None

    footprint_mask: np.ndarray | None = None
    detect_mask: np.ndarray | None = None

    summary: ClumpSummary | None = None
    


    def to_hapy_columns(self, prefix: str = "HCL_") -> dict:
        """
        Return summary quantities as a dictionary for the main HAPY catalog.

        Parameters
        ----------
        prefix : str
            Prefix for output column names.

        Returns
        -------
        dict
            Dictionary mapping HAPY catalog column names to scalar values.
        """

        out = {}

        for key, value in asdict(self.summary).items():
            colname = prefix + key.upper()
            out[colname] = value

        return out


# ---------------------------------------------------------------------
# Small compatibility helpers
# ---------------------------------------------------------------------


def _call_with_supported_kwargs(func, *args, **kwargs):
    """
    Call a function using only keyword arguments supported by its signature.

    This helps keep the module compatible with small photutils API changes,
    e.g., nlevels vs n_levels and nproc vs n_processes.
    """

    sig = inspect.signature(func)
    allowed = set(sig.parameters)

    clean_kwargs = {
        key: value
        for key, value in kwargs.items()
        if key in allowed
    }

    return func(*args, **clean_kwargs)


def _as_array(segm_or_array):
    """
    Return the underlying data array from a SegmentationImage or ndarray.
    """

    if hasattr(segm_or_array, "data"):
        return np.asarray(segm_or_array.data)

    return np.asarray(segm_or_array)


def _as_float_array(values):
    """
    Convert an astropy Column/Quantity/array-like object to a float ndarray.
    """

    if hasattr(values, "value"):
        values = values.value

    return np.asarray(values, dtype=float)


def _empty_clump_table() -> QTable:
    """
    Return an empty clump catalog table with basic expected columns.
    """

    tab = QTable()
    tab["label"] = np.array([], dtype=int)
    tab["is_good"] = np.array([], dtype=bool)
    return tab


# ---------------------------------------------------------------------
# Footprint and masking
# ---------------------------------------------------------------------

def make_footprint_weightimage(
    hdata,
    footprint_mask: np.ndarray,
    mask: np.ndarray | None = None,
) -> np.ndarray:
    """
    Make a simple weight image that is positive only for usable pixels
    inside the central R-band segmentation footprint.

    This is used by calculate_background_photutils so the H-alpha clump
    threshold is estimated only from the central galaxy footprint.

    Parameters
    ----------
    hdata : ndarray
        H-alpha image.

    footprint_mask : ndarray of bool
        True for pixels belonging to the central R-band galaxy footprint.

    mask : ndarray of bool or None
        Optional bad-pixel/object mask. True means masked.

    Returns
    -------
    weightimage : ndarray
        Integer weight image. Good pixels have value 1. Bad pixels have 0.
    """

    hdata = np.asarray(hdata)
    footprint_mask = np.asarray(footprint_mask, dtype=bool)

    good = footprint_mask & np.isfinite(hdata)

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)

        if mask.shape != hdata.shape:
            raise ValueError(
                "mask and hdata must have the same shape. "
                f"Got {mask.shape} and {hdata.shape}."
            )

        good &= ~mask

    weightimage = np.zeros(hdata.shape, dtype=np.int16)
    weightimage[good] = 1

    return weightimage

def make_object_footprint(rsegm, object_index: int) -> np.ndarray:
    """
    Return boolean mask for the central galaxy R-band segmentation footprint.

    Parameters
    ----------
    rsegm : SegmentationImage or ndarray
        R-band segmentation image.

    object_index : int
        Label corresponding to the central galaxy.

    Returns
    -------
    footprint_mask : ndarray of bool
        True for pixels belonging to the central galaxy footprint.
    """

    segdata = _as_array(rsegm)

    if object_index is None:
        raise ValueError("object_index must not be None.")

    footprint_mask = segdata == int(object_index)

    if not np.any(footprint_mask):
        warnings.warn(
            f"No pixels found for object_index={object_index}. "
            "Returning an empty footprint.",
            RuntimeWarning,
        )

    return footprint_mask


def mask_halpha_to_footprint(
    hdata,
    footprint_mask: np.ndarray,
    mask: np.ndarray | None = None,
):
    """
    Build an H-alpha detection image and detection mask.

    Pixels outside the R-band footprint, non-finite pixels, and optional
    existing mask pixels are excluded from detection.

    Parameters
    ----------
    hdata : ndarray
        H-alpha image.

    footprint_mask : ndarray of bool
        Boolean mask selecting the central galaxy R-band footprint.

    mask : ndarray of bool or None
        Optional existing bad-pixel/object mask. True means masked.

    Returns
    -------
    data_detect : ndarray
        H-alpha image for detection. Invalid pixels are replaced with 0
        initially; they will be set below the threshold later.

    detect_mask : ndarray of bool
        Boolean mask where True pixels should be ignored.
    """

    hdata = np.asarray(hdata, dtype=float)
    footprint_mask = np.asarray(footprint_mask, dtype=bool)

    if hdata.shape != footprint_mask.shape:
        raise ValueError(
            "hdata and footprint_mask must have the same shape. "
            f"Got {hdata.shape} and {footprint_mask.shape}."
        )

    detect_mask = ~footprint_mask
    detect_mask |= ~np.isfinite(hdata)

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)

        if mask.shape != hdata.shape:
            raise ValueError(
                "mask and hdata must have the same shape. "
                f"Got {mask.shape} and {hdata.shape}."
            )

        detect_mask |= mask

    data_detect = np.array(hdata, dtype=float, copy=True)
    data_detect[detect_mask] = 0.0

    return data_detect, detect_mask


# ---------------------------------------------------------------------
# Detection and catalog construction
# ---------------------------------------------------------------------
def summarize_nuclear_clumps(summary, table, nuclear_radius_pix=None):
    """
    Add nuclear-clump summary quantities.

    A nuclear clump is defined as a good clump whose centroid lies within
    nuclear_radius_pix of the adopted galaxy center.

    If the nuclear radius is defined and there is positive total good-clump
    flux, then nuclear_flux_frac is zero when no nuclear clumps are present.
    """

    if table is None or len(table) == 0:
        return summary

    if nuclear_radius_pix is None or not np.isfinite(nuclear_radius_pix):
        return summary

    if "r_gal_pix" not in table.colnames:
        return summary

    flux_col = _get_column(table, ["segment_flux", "source_sum", "flux"])
    if flux_col is None:
        return summary

    if "is_good" in table.colnames:
        good = np.asarray(table["is_good"], dtype=bool)
    else:
        good = np.ones(len(table), dtype=bool)

    r = _as_float_array(table["r_gal_pix"])
    flux = _as_float_array(flux_col)

    nuclear = good & np.isfinite(r) & (r <= float(nuclear_radius_pix))

    table["is_nuclear"] = nuclear

    summary.nuclear_radius_pix = float(nuclear_radius_pix)
    summary.has_nuclear = bool(np.any(nuclear))
    summary.n_nuclear = int(np.sum(nuclear))

    good_flux = good & np.isfinite(flux) & (flux > 0)
    total_flux = np.nansum(flux[good_flux])

    if total_flux > 0:
        nuclear_flux = np.nansum(flux[nuclear & np.isfinite(flux) & (flux > 0)])
        summary.nuclear_flux_frac = float(nuclear_flux / total_flux)

    return summary



def add_clump_offset_columns(
    table,
    xcenter: float,
    ycenter: float,
):
    """
    Add clump centroid offsets relative to the adopted galaxy center.

    Parameters
    ----------
    table : QTable
        Clump SourceCatalog table.

    xcenter, ycenter : float
        Adopted galaxy center in pixel coordinates.

    Returns
    -------
    table : QTable
        Input table with added dx/dy/r/theta columns.
    """

    if table is None or len(table) == 0:
        return table

    x_col = None
    y_col = None

    for candidate in ["xcentroid", "x_centroid"]:
        if candidate in table.colnames:
            x_col = candidate
            break

    for candidate in ["ycentroid", "y_centroid"]:
        if candidate in table.colnames:
            y_col = candidate
            break

    if x_col is None or y_col is None:
        return table

    x = _as_float_array(table[x_col])
    y = _as_float_array(table[y_col])

    dx = x - float(xcenter)
    dy = y - float(ycenter)

    r = np.sqrt(dx**2 + dy**2)
    theta = np.degrees(np.arctan2(dy, dx))

    table["dx_gal_pix"] = dx
    table["dy_gal_pix"] = dy
    table["r_gal_pix"] = r
    table["theta_gal_deg"] = theta

    return table

def summarize_clump_offsets(summary, table):
    """
    Add clump centroid offset summary quantities to ClumpSummary.
    """

    if table is None or len(table) == 0:
        return summary

    if "r_gal_pix" not in table.colnames:
        return summary

    if "is_good" in table.colnames:
        good = np.asarray(table["is_good"], dtype=bool)
    else:
        good = np.ones(len(table), dtype=bool)

    if np.sum(good) == 0:
        return summary

    r = _as_float_array(table["r_gal_pix"])[good]
    dx = _as_float_array(table["dx_gal_pix"])[good]
    dy = _as_float_array(table["dy_gal_pix"])[good]

    flux_col = _get_column(table, ["segment_flux", "source_sum", "flux"])
    if flux_col is not None:
        flux = _as_float_array(flux_col)[good]
    else:
        flux = np.ones(np.sum(good), dtype=float)

    finite = np.isfinite(r)

    if not np.any(finite):
        return summary

    summary.rmin_pix = float(np.nanmin(r[finite]))
    summary.rmax_pix = float(np.nanmax(r[finite]))
    summary.rmed_pix = float(np.nanmedian(r[finite]))
    summary.rmean_pix = float(np.nanmean(r[finite]))

    good_weight = finite & np.isfinite(flux) & (flux > 0)

    if np.any(good_weight):
        summary.rfluxwt_pix = float(
            np.average(r[good_weight], weights=flux[good_weight])
        )

        # Brightest clump offset.
        imax = np.argmax(flux[good_weight])

        r_good = r[good_weight]
        dx_good = dx[good_weight]
        dy_good = dy[good_weight]

        summary.rbright_pix = float(r_good[imax])
        summary.dx_bright_pix = float(dx_good[imax])
        summary.dy_bright_pix = float(dy_good[imax])

    return summary

def estimate_clump_threshold(
    hdata: np.ndarray,
    footprint_mask: np.ndarray,
    mask: np.ndarray | None,
    config: ClumpDetectionConfig,
):
    """
    Estimate the H-alpha clump detection threshold using HAPY's
    calculate_background_photutils helper.

    The background is estimated only from valid pixels inside the central
    galaxy R-band footprint by passing a footprint weight image.

    Returns
    -------
    threshold : float
        Detection threshold for H-alpha clumps.

    background_mean : float
        Sigma-clipped mean background.

    background_median : float
        Sigma-clipped median background.

    background_rms : float
        Sigma-clipped background standard deviation.

    weightimage : ndarray
        Simple weight image used for the background calculation.
    """

    hdata = np.asarray(hdata, dtype=float)

    weightimage = make_footprint_weightimage(
        hdata,
        footprint_mask,
        mask=mask,
    )

    if np.sum(weightimage > 0) == 0:
        return np.nan, np.nan, np.nan, np.nan, weightimage

    mean, median, std = calculate_background_photutils(
        hdata,
        grow_radius=config.background_grow_radius,
        npixels=config.background_npixels,
        weightimage=weightimage,
        nsigma=config.background_mask_nsigma,
        clip_sigma=config.background_clip_sigma,
        clip_maxiters=config.background_clip_maxiters,
    )

    threshold = median + config.nsigma * std

    return (
        float(threshold),
        float(mean),
        float(median),
        float(std),
        weightimage,
    )

def detect_halpha_clumps(
    hdata: np.ndarray,
    data_detect: np.ndarray,
    footprint_mask: np.ndarray,
    detect_mask: np.ndarray,
    config: ClumpDetectionConfig,
    mask: np.ndarray | None = None,
):
    """
    Detect and optionally deblend H-alpha clumps.

    Parameters
    ----------
    hdata : ndarray
        Original H-alpha image.

    data_detect : ndarray
        H-alpha image prepared for detection.

    footprint_mask : ndarray of bool
        Central R-band galaxy footprint.

    detect_mask : ndarray of bool
        Boolean detection mask. True means ignore pixel.

    config : ClumpDetectionConfig
        Detection and deblending parameters.

    mask : ndarray or None
        Optional external mask. True means masked.

    Returns
    -------
    segm : SegmentationImage or None
        Final segmentation image after optional deblending.

    segm_detected : SegmentationImage or None
        Initial segmentation image before deblending.

    threshold : float
        Detection threshold.

    background_mean : float
        Sigma-clipped mean background.

    background_median : float
        Sigma-clipped median background.

    background_rms : float
        Sigma-clipped background standard deviation.

    weightimage : ndarray
        Weight image used for background estimation.
    """

    (
        threshold,
        background_mean,
        background_median,
        background_rms,
        weightimage,
    ) = estimate_clump_threshold(
        hdata,
        footprint_mask,
        mask,
        config,
    )

    if not np.isfinite(threshold):
        return (
            None,
            None,
            threshold,
            background_mean,
            background_median,
            background_rms,
            weightimage,
        )

    work = np.array(data_detect, dtype=float, copy=True)

    # Make sure masked pixels are safely below threshold.
    good = ~detect_mask
    scale = np.nanstd(work[good]) if np.any(good) else 1.0

    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0

    fill_value = threshold - 10.0 * scale - 1.0
    work[detect_mask] = fill_value

    segm_detected = _call_with_supported_kwargs(
        detect_sources,
        work,
        threshold,
        config.npixels,
        connectivity=config.connectivity,
        mask=detect_mask,
    )

    if segm_detected is None:
        return (
            None,
            None,
            threshold,
            background_mean,
            background_median,
            background_rms,
            weightimage,
        )

    segm = segm_detected

    if config.deblend:
        segm = _call_with_supported_kwargs(
            deblend_sources,
            work,
            segm_detected,
            config.npixels,
            labels=None,
            nlevels=config.nlevels,
            n_levels=config.nlevels,
            contrast=config.contrast,
            mode=config.mode,
            connectivity=config.connectivity,
            relabel=config.relabel,
            nproc=config.nproc,
            n_processes=config.nproc,
            progress_bar=config.progress_bar,
        )

    return (
        segm,
        segm_detected,
        threshold,
        background_mean,
        background_median,
        background_rms,
        weightimage,
    )




def make_clump_catalog(
    hdata,
    segm,
    error=None,
    mask=None,
    wcs=None,
    progress_bar: bool = False,
):
    """
    Create a photutils SourceCatalog for the H-alpha clumps.

    Parameters
    ----------
    hdata : ndarray
        Original H-alpha image used for measurements.

    segm : SegmentationImage
        Final clump segmentation image.

    error : ndarray or None
        Optional H-alpha uncertainty image.

    mask : ndarray or None
        Optional mask. True means masked.

    wcs : astropy.wcs.WCS or None
        Optional WCS.

    progress_bar : bool
        Whether SourceCatalog should show a progress bar.

    Returns
    -------
    catalog : SourceCatalog
        Photutils source catalog.
    """

    if segm is None:
        return None

    kwargs = {
        "error": error,
        "mask": mask,
        "wcs": wcs,
        "progress_bar": progress_bar,
    }

    return _call_with_supported_kwargs(
        SourceCatalog,
        np.asarray(hdata, dtype=float),
        segm,
        **kwargs,
    )


# ---------------------------------------------------------------------
# Summary calculations
# ---------------------------------------------------------------------


def _get_column(table: QTable, names, default=None):
    """
    Return the first available column from a list of candidate names.
    """

    for name in names:
        if name in table.colnames:
            return table[name]

    return default


def _add_good_clump_flag(
    table: QTable,
    config: ClumpDetectionConfig,
) -> np.ndarray:
    """
    Add an is_good column to the clump table.

    The full catalog is always preserved. This flag only controls which
    clumps enter the compact summary statistics.
    """

    n = len(table)

    if n == 0:
        table["is_good"] = np.array([], dtype=bool)
        return table["is_good"]

    good = np.ones(n, dtype=bool)

    flux = _get_column(table, ["segment_flux", "source_sum", "flux"])
    area = _get_column(table, ["area", "segment_area"])

    if config.min_flux is not None and flux is not None:
        flux_arr = _as_float_array(flux)
        good &= np.isfinite(flux_arr)
        good &= flux_arr >= config.min_flux

    if config.min_area is not None and area is not None:
        area_arr = _as_float_array(area)
        good &= np.isfinite(area_arr)
        good &= area_arr >= config.min_area

    table["is_good"] = good

    return good

def summarize_clump_flux_ratios(summary, table):
    """
    Add clump flux dominance ratios.
    """

    if table is None or len(table) == 0:
        return summary

    flux_col = _get_column(table, ["segment_flux", "source_sum", "flux"])
    if flux_col is None:
        return summary

    if "is_good" in table.colnames:
        good = np.asarray(table["is_good"], dtype=bool)
    else:
        good = np.ones(len(table), dtype=bool)

    flux = _as_float_array(flux_col)
    flux = flux[good & np.isfinite(flux) & (flux > 0)]

    if len(flux) < 2:
        summary.bright_to_second_flux = np.nan
        return summary

    flux_sorted = np.sort(flux)[::-1]
    summary.bright_to_second_flux = float(flux_sorted[0] / flux_sorted[1])

    return summary

def summarize_clumps(
    hdata,
    footprint_mask: np.ndarray,
    detect_mask: np.ndarray,
    segm,
    catalog,
    config: ClumpDetectionConfig,
    threshold: float = np.nan,
    background_mean: float = np.nan,
    background_median: float = np.nan,
    background_rms: float = np.nan,
) -> tuple[ClumpSummary, QTable]:
    

    """
    Summarize the detected H-alpha clumps.

    Parameters
    ----------
    hdata : ndarray
        Original H-alpha image.

    footprint_mask : ndarray of bool
        Central galaxy R-band footprint.

    detect_mask : ndarray of bool
        Detection mask.

    segm : SegmentationImage or None
        H-alpha clump segmentation.

    catalog : SourceCatalog or None
        Photutils SourceCatalog.

    config : ClumpDetectionConfig
        Detection configuration.

    threshold, background, background_rms : float
        Detection statistics.

    Returns
    -------
    summary : ClumpSummary
        Compact summary quantities.

    table : QTable
        SourceCatalog table with an added is_good column.
    """

    hdata = np.asarray(hdata, dtype=float)
    footprint_mask = np.asarray(footprint_mask, dtype=bool)
    detect_mask = np.asarray(detect_mask, dtype=bool)

    valid_footprint = footprint_mask & ~detect_mask & np.isfinite(hdata)

    footprint_area_pix = int(np.sum(valid_footprint))

    if footprint_area_pix > 0:
        total_halpha_flux_footprint = float(np.nansum(hdata[valid_footprint]))
    else:
        total_halpha_flux_footprint = np.nan

    summary = ClumpSummary(
        threshold=threshold,
        background_mean=background_mean,
        background_median=background_median,
        background_rms=background_rms,
        total_halpha_flux_footprint=total_halpha_flux_footprint,
        footprint_area_pix=footprint_area_pix,
    )



    if catalog is None or segm is None:
        table = _empty_clump_table()
        return summary, table

    table = catalog.to_table()

    if len(table) == 0:
        table["is_good"] = np.array([], dtype=bool)
        return summary, table

    good = _add_good_clump_flag(table, config)

    label_col = _get_column(table, ["label", "id"])
    flux_col = _get_column(table, ["segment_flux", "source_sum", "flux"])
    x_col = _get_column(table, ["xcentroid", "x_centroid", "xcen"])
    y_col = _get_column(table, ["ycentroid", "y_centroid", "ycen"])

    if label_col is None or flux_col is None:
        return summary, table

    labels = np.asarray(label_col, dtype=int)
    flux = _as_float_array(flux_col)

    good &= np.isfinite(flux)

    n_clumps = len(table)
    n_good = int(np.sum(good))

    summary.n_clumps = int(n_clumps)
    summary.n_clumps_good = n_good

    if n_good == 0:
        return summary, table

    good_labels = labels[good]
    good_flux = flux[good]

    segdata = _as_array(segm)

    clump_area_pix = int(np.sum(np.isin(segdata, good_labels)))
    total_clump_flux = float(np.nansum(good_flux))

    summary.clump_area_pix = clump_area_pix
    summary.total_clump_flux = total_clump_flux

    if summary.footprint_area_pix > 0:
        summary.clump_area_fraction = clump_area_pix / summary.footprint_area_pix

    if np.isfinite(summary.total_halpha_flux_footprint) and summary.total_halpha_flux_footprint > 0:
        summary.clump_flux_fraction = total_clump_flux / summary.total_halpha_flux_footprint

    # Sort clumps by flux from brightest to faintest.
    order = np.argsort(good_flux)[::-1]
    flux_sorted = good_flux[order]
    labels_sorted = good_labels[order]

    if len(flux_sorted) > 0:
        summary.brightest_clump_flux = float(flux_sorted[0])

        if total_clump_flux > 0:
            summary.brightest_clump_fraction = float(flux_sorted[0] / total_clump_flux)
            summary.top2_clump_fraction = float(np.nansum(flux_sorted[:2]) / total_clump_flux)
            summary.top3_clump_fraction = float(np.nansum(flux_sorted[:3]) / total_clump_flux)

    # Brightest clump centroid.
    if x_col is not None and y_col is not None:
        x = _as_float_array(x_col)[good]
        y = _as_float_array(y_col)[good]

        x_sorted = x[order]
        y_sorted = y[order]

        if len(x_sorted) > 0:
            summary.brightest_xcentroid = float(x_sorted[0])
            summary.brightest_ycentroid = float(y_sorted[0])

        # Flux-weighted centroid of all good clumps.
        pos_good = np.isfinite(good_flux) & (good_flux > 0)
        pos_good &= np.isfinite(x)
        pos_good &= np.isfinite(y)

        if np.any(pos_good):
            weights = good_flux[pos_good]
            summary.flux_weighted_xcentroid = float(
                np.average(x[pos_good], weights=weights)
            )
            summary.flux_weighted_ycentroid = float(
                np.average(y[pos_good], weights=weights)
            )

    # Add a rank column for convenience.
    rank = np.full(len(table), -1, dtype=int)
    for i, lab in enumerate(labels_sorted, start=1):
        rank[labels == lab] = i

    table["flux_rank"] = rank

    return summary, table

def find_halpha_peaks(
    data_detect: np.ndarray,
    detect_mask: np.ndarray,
    threshold: float,
    segm=None,
    config: ClumpDetectionConfig | None = None,
    wcs=None,
) -> QTable:
    """
    Find local H-alpha peaks within the central R-band footprint.

    This is a companion diagnostic to the segmentation-based clump catalog.
    Peaks identify local maxima, while clumps are defined by segmented regions.
    """

    if config is None:
        config = ClumpDetectionConfig()

    if not np.isfinite(threshold):
        return QTable()

    kwargs = dict(
        box_size=config.peak_box_size,
        mask=detect_mask,
        wcs=wcs,
    )

    if config.peak_nmax is not None:
        # photutils has used npeaks in some versions. The compatibility
        # wrapper will ignore unsupported alternatives.
        kwargs["npeaks"] = config.peak_nmax
        kwargs["n_peaks"] = config.peak_nmax

    peaks = _call_with_supported_kwargs(
        find_peaks,
        data_detect,
        threshold,
        **kwargs,
    )

    if peaks is None:
        return QTable()

    peaks = QTable(peaks)

    # Add clump segment label for each peak.
    if segm is not None and len(peaks) > 0:
        segdata = _as_array(segm)

        x_col = "x_peak" if "x_peak" in peaks.colnames else None
        y_col = "y_peak" if "y_peak" in peaks.colnames else None

        if x_col is not None and y_col is not None:
            x = np.asarray(peaks[x_col], dtype=int)
            y = np.asarray(peaks[y_col], dtype=int)

            inside = (
                (x >= 0)
                & (x < segdata.shape[1])
                & (y >= 0)
                & (y < segdata.shape[0])
            )

            labels = np.full(len(peaks), -1, dtype=int)
            labels[inside] = segdata[y[inside], x[inside]]

            peaks["clump_label"] = labels
            peaks["in_clump"] = labels > 0

    return peaks

def summarize_peak_table(
    summary: ClumpSummary,
    peak_table: QTable | None,
) -> ClumpSummary:
    """
    Add local-peak summary quantities to ClumpSummary.
    """

    if peak_table is None or len(peak_table) == 0:
        return summary

    summary.n_peaks = int(len(peak_table))

    if "in_clump" in peak_table.colnames:
        summary.n_peaks_in_clumps = int(np.sum(peak_table["in_clump"]))

    peak_col = None
    for candidate in ["peak_value", "peak", "value"]:
        if candidate in peak_table.colnames:
            peak_col = candidate
            break

    if peak_col is not None:
        peak_values = _as_float_array(peak_table[peak_col])
        finite = np.isfinite(peak_values)

        if np.any(finite):
            summary.peak_max = float(np.nanmax(peak_values[finite]))
            summary.peak_sum = float(np.nansum(peak_values[finite]))

    return summary
# ---------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------


def save_clump_outputs(
    result: ClumpAnalysisResult,
    output_dir,
    basename: str,
    overwrite: bool = True,
    table_format: str = "ascii.ecsv",
) -> ClumpAnalysisResult:
    """
    Save the clump SourceCatalog table and segmentation image.

    Parameters
    ----------
    result : ClumpAnalysisResult
        Result returned by analyze_halpha_clumps.

    output_dir : str or Path
        Directory where outputs should be written.

    basename : str
        Base filename, usually the HAPY cutout tag.

    overwrite : bool
        Whether to overwrite existing files.

    table_format : str
        Astropy table format for the catalog.

    Returns
    -------
    result : ClumpAnalysisResult
        Same result object, with catalog_path and segmentation_path filled in.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    catalog_path = output_dir / f"{basename}-halpha-clumps.ecsv"
    segm_path = output_dir / f"{basename}-halpha-clumps-segm.fits"

    # Save the persistent SourceCatalog table.
    if result.table is not None:
        result.table.meta["HAPY_PRODUCT"] = "halpha_clump_catalog"
        result.table.write(
            catalog_path,
            format=table_format,
            overwrite=overwrite,
        )

        result.summary.catalog_path = str(catalog_path)

    # Save the clump segmentation map.
    if result.segm is not None:
        segm_data = _as_array(result.segm).astype(np.int32)
    else:
        segm_data = np.zeros(result.footprint_mask.shape, dtype=np.int32)

    header = fits.Header()
    header["HAPYPROD"] = "HALPHA_CLUMP_SEGM"
    header["NCLUMP"] = int(result.summary.n_clumps)
    header["NGOOD"] = int(result.summary.n_clumps_good)

    fits.writeto(
        segm_path,
        segm_data,
        header=header,
        overwrite=overwrite,
    )

    result.summary.segmentation_path = str(segm_path)

    if getattr(result.summary, "diagnostic_path", ""):
        header["DIAGPATH"] = result.summary.diagnostic_path



    # Save local peak table, if available.
    if getattr(result, "peak_table", None) is not None:
        if len(result.peak_table) > 0 or True:
            peaks_path = output_dir / f"{basename}-halpha-clump-peaks.ecsv"

            result.peak_table.meta["HAPY_PRODUCT"] = "halpha_clump_peaks"
            result.peak_table.write(
                peaks_path,
                format=table_format,
                overwrite=overwrite,
            )

            result.summary.peaks_path = str(peaks_path)
            
    return result


def save_clump_diagnostic(
    result: ClumpAnalysisResult,
    hdata,
    output_dir,
    basename: str,
    config: ClumpDetectionConfig,
    overwrite: bool = True,
):
    """
    Save a diagnostic image showing the H-alpha image, R-band footprint,
    H-alpha clump segmentation, clump centroids, nuclear clumps, and
    optional peak locations.

    Parameters
    ----------
    result : ClumpAnalysisResult
        Result returned by analyze_halpha_clumps.

    hdata : ndarray
        Original H-alpha image.

    output_dir : str or Path
        Directory where diagnostic image should be written.

    basename : str
        Base filename, usually the HAPY cutout tag.

    config : ClumpDetectionConfig
        Clump detection configuration.

    overwrite : bool
        Whether to overwrite an existing diagnostic image.

    Returns
    -------
    diagnostic_path : pathlib.Path
        Path to saved diagnostic image.
    """

    import warnings
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path

    from hapy.imagetools.imutils import make_masked_display_image_norm

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    diagnostic_path = (
        output_dir
        / f"{basename}-halpha-clumps-diagnostic.{config.diagnostic_format}"
    )

    if diagnostic_path.exists() and not overwrite:
        result.summary.diagnostic_path = str(diagnostic_path)
        return diagnostic_path

    hdata = np.asarray(hdata, dtype=float)

    footprint = np.asarray(result.footprint_mask, dtype=bool)
    detect_mask = np.asarray(result.detect_mask, dtype=bool)

    if result.segm is not None:
        segm_data = _as_array(result.segm)
    else:
        segm_data = np.zeros(hdata.shape, dtype=int)

    # ------------------------------------------------------------
    # Display image and normalization
    # ------------------------------------------------------------
    hshow, norm = make_masked_display_image_norm(
        hdata,
        mask=footprint,
        percent=config.diagnostic_percent,
        stretch="sqrt",
    )

    # Image shown in panel 2: hide pixels excluded from detection.
    detect_show = np.array(hdata, dtype=float, copy=True)
    detect_show[detect_mask] = np.nan
    detect_show[~footprint] = np.nan

    # ------------------------------------------------------------
    # Zoom limits
    # ------------------------------------------------------------
    xlim, ylim = get_bbox_from_mask(
        result.footprint_mask,
        pad=25,
        min_size=50,
    )

    xlim, ylim = make_limits_square(
        xlim,
        ylim,
        image_shape=hdata.shape,
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    ax0, ax1, ax2 = axes

    # ------------------------------------------------------------
    # Panel 1: H-alpha image with R-band footprint
    # ------------------------------------------------------------
    ax0.imshow(hshow, origin="lower", norm=norm, cmap="gray_r")

    try:
        ax0.contour(
            footprint.astype(int),
            levels=[0.5],
            colors="cyan",
            linewidths=1.0,
        )
    except Exception:
        pass

    ax0.set_title("Hα within R-band footprint")

    # ------------------------------------------------------------
    # Panel 2: H-alpha image with clump segmentation boundaries
    # ------------------------------------------------------------
    ax1.imshow(detect_show, origin="lower", norm=norm, cmap="gray")

    if np.nanmax(segm_data) > 0:
        try:
            ax1.contour(
                segm_data > 0,
                levels=[0.5],
                colors="yellow",
                linewidths=1.0,
            )
        except Exception:
            pass

    ax1.set_title("Hα clump segmentation")

    # ------------------------------------------------------------
    # Add clump centroids.
    #
    # Red circles  = good clumps
    # Lime circles = nuclear good clumps
    # ------------------------------------------------------------
    table = result.table

    if table is not None and len(table) > 0:
        x_col = None
        y_col = None

        for candidate in ["xcentroid", "x_centroid"]:
            if candidate in table.colnames:
                x_col = candidate
                break

        for candidate in ["ycentroid", "y_centroid"]:
            if candidate in table.colnames:
                y_col = candidate
                break

        if x_col is not None and y_col is not None:
            x = _as_float_array(table[x_col])
            y = _as_float_array(table[y_col])

            if "is_good" in table.colnames:
                good = np.asarray(table["is_good"], dtype=bool)
            else:
                good = np.ones(len(table), dtype=bool)

            ax1.plot(
                x[good],
                y[good],
                "o",
                ms=5,
                mfc="none",
                mec="red",
                mew=1.2,
                label="good clump",
            )

            if "is_nuclear" in table.colnames:
                is_nuclear = np.asarray(table["is_nuclear"], dtype=bool)
                nuclear_good = good & is_nuclear

                if np.any(nuclear_good):
                    ax1.plot(
                        x[nuclear_good],
                        y[nuclear_good],
                        "o",
                        ms=9,
                        mfc="none",
                        mec="lime",
                        mew=1.8,
                        label="nuclear clump",
                    )

            if "label" in table.colnames:
                labels = np.asarray(table["label"], dtype=int)

                for xi, yi, lab, use in zip(x, y, labels, good):
                    if use and np.isfinite(xi) and np.isfinite(yi):
                        ax1.text(
                            xi + 1,
                            yi + 1,
                            str(lab),
                            color="red",
                            fontsize=7,
                        )

    # Optional Kron apertures
    if config.plot_kron_apertures and result.catalog is not None:
        try:
            result.catalog.plot_kron_apertures(
                ax=ax1,
                color="white",
                lw=1.0,
            )
        except Exception as err:
            warnings.warn(
                f"Could not plot Kron apertures: {err}",
                RuntimeWarning,
            )

    # ------------------------------------------------------------
    # Panel 3: segmentation map, plus peaks if available
    # ------------------------------------------------------------
    if result.segm is not None and hasattr(result.segm, "imshow"):
        result.segm.imshow(ax=ax2)
    else:
        ax2.imshow(segm_data, origin="lower", interpolation="nearest")

    ax2.set_title("Clump labels")

    peak_table = getattr(result, "peak_table", None)

    if peak_table is not None and len(peak_table) > 0:
        x_peak_col = "x_peak" if "x_peak" in peak_table.colnames else None
        y_peak_col = "y_peak" if "y_peak" in peak_table.colnames else None

        if x_peak_col is not None and y_peak_col is not None:
            xpeak = _as_float_array(peak_table[x_peak_col])
            ypeak = _as_float_array(peak_table[y_peak_col])

            ax2.plot(
                xpeak,
                ypeak,
                "+",
                ms=8,
                mew=1.5,
                color="white",
                alpha=0.4,
            )

    # ------------------------------------------------------------
    # Cosmetic cleanup
    # ------------------------------------------------------------
    for ax in axes:
        ax.set_xlabel("x [pix]")
        ax.set_ylabel("y [pix]")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    # ------------------------------------------------------------
    # Figure title
    # ------------------------------------------------------------
    summary = result.summary

    threshold = getattr(summary, "threshold", np.nan)
    if np.isfinite(threshold):
        threshold_text = f"{threshold:.3g}"
    else:
        threshold_text = "NaN"

    has_nuclear = bool(getattr(summary, "has_nuclear", False))
    n_nuclear = int(getattr(summary, "n_nuclear", 0))
    nuclear_flux_frac = getattr(summary, "nuclear_flux_frac", np.nan)

    if np.isfinite(nuclear_flux_frac):
        nuclear_flux_frac_text = f"{nuclear_flux_frac:.2f}"
    else:
        nuclear_flux_frac_text = "NaN"

    nuclear_radius_fwhm = getattr(config, "nuclear_radius_fwhm", np.nan)

    if np.isfinite(nuclear_radius_fwhm):
        nuclear_radius_text = f"{nuclear_radius_fwhm:.1f} FWHM"
    else:
        nuclear_radius_text = "NaN"

    title = (
        f"{basename}\n"
        f"Nclump={summary.n_clumps}, "
        f"Ngood={summary.n_clumps_good}, "
        f"Npeak={getattr(summary, 'n_peaks', 0)}, "
        f"threshold={threshold_text}\n"
        f"nuclear={has_nuclear}, "
        f"Nnuc={n_nuclear}, "
        f"fnuc={nuclear_flux_frac_text}, "
        f"Rnuc={nuclear_radius_text}"
    )

    fig.suptitle(title, fontsize=11)

    fig.savefig(diagnostic_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    result.summary.diagnostic_path = str(diagnostic_path)

    return diagnostic_path





# ---------------------------------------------------------------------
# High-level wrapper
# ---------------------------------------------------------------------

def analyze_halpha_clumps(
    hdata,
    rsegm,
    object_index: int,
    error=None,
    mask=None,
    wcs=None,
    xcenter=None,
    ycenter=None,
    nuclear_radius_pix=None,
    config: ClumpDetectionConfig | None = None,
    output_dir=None,
    basename: str | None = None,
    overwrite: bool = True,
) -> ClumpAnalysisResult:
    

    """
    Detect, deblend, measure, summarize, and optionally save H-alpha clumps.
    """

    if config is None:
        config = ClumpDetectionConfig()

    footprint_mask = make_object_footprint(rsegm, object_index)

    data_detect, detect_mask = mask_halpha_to_footprint(
        hdata,
        footprint_mask,
        mask=mask,
    )

    segdata = _as_array(rsegm)

    # print("DEBUG segmentation")
    # print("  object_index/object_label passed =", object_index)
    # print("  segdata shape =", segdata.shape)
    # print("  segdata dtype =", segdata.dtype)
    # print("  unique labels first 20 =", np.unique(segdata)[:20])
    # print("  unique labels last 20 =", np.unique(segdata)[-20:])
    # print("  footprint pixels =", np.sum(footprint_mask))


    (
        segm,
        segm_detected,
        threshold,
        background_mean,
        background_median,
        background_rms,
        weightimage,
    ) = detect_halpha_clumps(
        hdata,
        data_detect,
        footprint_mask,
        detect_mask,
        config,
        mask=mask,
    )

    # print("DEBUG clumps")
    # print("  footprint pixels:", np.sum(footprint_mask))
    # print("  detect good pixels:", np.sum(~detect_mask))
    # print("  hdata max in footprint:", np.nanmax(np.asarray(hdata)[footprint_mask]))
    # print("  data_detect max good:", np.nanmax(data_detect[~detect_mask]))
    # print("  threshold:", threshold if "threshold" in locals() else "not computed yet")
        

    # For SourceCatalog measurement, include both the user mask and the
    # outside-footprint mask. This ensures clump measurements are restricted
    # to valid central-galaxy pixels.
    source_mask = detect_mask

    catalog = None

    if segm is not None:
        catalog = make_clump_catalog(
            hdata,
            segm,
            error=error,
            mask=source_mask,
            wcs=wcs,
            progress_bar=config.progress_bar,
        )

    summary, table = summarize_clumps(
        hdata,
        footprint_mask,
        detect_mask,
        segm,
        catalog,
        config,
        threshold=threshold,
        background_mean=background_mean,
        background_median=background_median,
        background_rms=background_rms,
    )

    if xcenter is not None and ycenter is not None:
        table = add_clump_offset_columns(
            table,
            xcenter=xcenter,
            ycenter=ycenter,
        )

        summary = summarize_clump_offsets(summary, table)


    summary = summarize_clump_flux_ratios(summary, table)

    if nuclear_radius_pix is not None:
        summary = summarize_nuclear_clumps(
            summary,
            table,
            nuclear_radius_pix=nuclear_radius_pix,
        )
    # ------------------------------------------------------------
    # Optional local peak finding
    # ------------------------------------------------------------
    peak_table = None
    point_source_table = None

    if config.find_peaks:
        peak_table = find_halpha_peaks(
            data_detect,
            detect_mask,
            threshold,
            segm=segm,
            config=config,
            wcs=wcs,
        )
        summary = summarize_peak_table(summary, peak_table)

    # ------------------------------------------------------------
    # Optional point-source-like detection
    # ------------------------------------------------------------
    if config.find_point_sources:
        point_source_table = find_halpha_point_sources(
            data_detect,
            detect_mask,
            background_rms,
            config,
        )
        summary = summarize_point_source_table(summary, point_source_table)

    # ------------------------------------------------------------
    # Store useful metadata in the persistent clump table
    # ------------------------------------------------------------
    table.meta["HAPY_PRODUCT"] = "halpha_clump_catalog"
    table.meta["OBJECT_INDEX"] = int(object_index)

    table.meta["NSIGMA"] = float(config.nsigma)
    table.meta["NPIXELS"] = int(config.npixels)

    table.meta["DEBLEND"] = bool(config.deblend)
    table.meta["NLEVELS"] = int(config.nlevels)
    table.meta["CONTRAST"] = float(config.contrast)

    table.meta["THRESH"] = float(threshold) if np.isfinite(threshold) else np.nan
    table.meta["BKGMEAN"] = float(background_mean) if np.isfinite(background_mean) else np.nan
    table.meta["BKGMED"] = float(background_median) if np.isfinite(background_median) else np.nan
    table.meta["BKGRMS"] = float(background_rms) if np.isfinite(background_rms) else np.nan

    table.meta["FINDPEAK"] = bool(config.find_peaks)
    table.meta["PEAKBOX"] = int(config.peak_box_size)

    if config.peak_min_separation is not None:
        table.meta["PEAKSEP"] = int(config.peak_min_separation)

    table.meta["FINDPSRC"] = bool(config.find_point_sources)

    # ------------------------------------------------------------
    # Package result
    # ------------------------------------------------------------
    result = ClumpAnalysisResult(
        catalog=catalog,
        table=table,
        peak_table=peak_table,
        point_source_table=point_source_table,
        segm=segm,
        segm_detected=segm_detected,
        footprint_mask=footprint_mask,
        detect_mask=detect_mask,
        summary=summary,
    )

    # ------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------
    if output_dir is not None:
        if basename is None:
            raise ValueError("basename must be provided when output_dir is provided.")

        if config.save_catalog or config.save_segmentation:
            result = save_clump_outputs(
                result,
                output_dir=output_dir,
                basename=basename,
                overwrite=overwrite,
                table_format=config.table_format,
            )

        if config.save_diagnostic:
            save_clump_diagnostic(
                result,
                hdata=hdata,
                output_dir=output_dir,
                basename=basename,
                config=config,
                overwrite=overwrite,
            )

    return result

