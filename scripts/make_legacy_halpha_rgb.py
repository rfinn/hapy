#!/usr/bin/env python

# Created with Codex assistance for JWST proposal figure prototyping.
# Adapted from local hapy/havirgo workflows and Legacy Survey RGB stretch logic.

"""
Make a false-color RGB image from aligned Legacy Survey and Halpha FITS images.

Default channel recipe:
    R = Legacy z
    G = Legacy r
    B = Legacy g

The Halpha/SFR image is scaled separately and overlaid in color.  If an HI
moment-zero map is available, HI contours or a filled HI intensity layer can be
drawn using WCSAxes transforms.
"""

import argparse
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib.patheffects as path_effects
import numpy as np

if not hasattr(np, "in1d"):
    np.in1d = np.isin

try:
    import numpy.lib._function_base_impl as np_function_base_impl

    if not hasattr(np_function_base_impl, "_check_interpolation_as_method"):

        def _check_interpolation_as_method(method, interpolation, fname):
            if interpolation is not None:
                if method != "linear":
                    raise TypeError(
                        f"{fname} received both 'method' and 'interpolation'."
                    )
                method = interpolation
            return method

        np_function_base_impl._check_interpolation_as_method = (
            _check_interpolation_as_method
        )
except Exception:
    pass

from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.cosmology import WMAP9 as cosmo
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import AsinhStretch, LinearStretch
from astropy.wcs import WCS
from astropy.wcs import utils as wcs_utils
from scipy import ndimage
from scipy.spatial import ConvexHull


DEFAULT_HI_FILES = {
    "VFID5709": "~/research/Virgo/MeerKAT/final_24Jan2025/J1406_p_0601_lw00_final_1_mom0.fits",
    "VFID6033": "~/research/Virgo/MeerKAT/final_24Jan2025/J1420_p_0336_lw00_final_1_mom0.fits",
    "VFID6091": "~/research/Virgo/MeerKAT/final_24Jan2025/J1420_p_0336_lw00_final_3_mom0.fits",
    "VFID6018": "~/research/Virgo/MeerKAT/final_24Jan2025/J1420_p_0336_lw00_final_6_mom0.fits",
    "VFID6362": "~/research/Virgo/MeerKAT/final_24Jan2025/J1502_p_0150_lw00_final_2_mom0.fits",
    "VFID5889": "~/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/ngc5363_ngc5364_ngc5360.fits",
    "VFID5851": "~/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_2_mom0.fits",
    "VFID5855": "~/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_6_mom0.fits",
    "VFID5842": "~/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_fin_lw05_bpcorr_5_mom0.fits",
    "VFID5859": "~/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/J1355_p_0512_lw05_final_15_mom0.fits",
    "VFID5892": "~/research/Virgo/alma/2023/MeerKAT_ALMA_target_list/ngc5363_ngc5364_ngc5360.fits",
}

DEFAULT_NIRCAM_APERTURES = tuple(
    f"NRC{module}{detector}_FULL" for module in "AB" for detector in range(1, 6)
)


def expand_path(path):
    """
    Expand ~ and environment variables in a user-supplied path.
    """
    return Path(os.path.expandvars(os.path.expanduser(str(path))))


def find_one(directory, pattern, label, preferred_prefix=None):
    """
    Return exactly one file matching a glob pattern.
    """
    matches = sorted(Path(directory).glob(pattern))

    if len(matches) == 0:
        raise FileNotFoundError(f"Could not find {label}: {Path(directory) / pattern}")

    if len(matches) > 1 and preferred_prefix is not None:
        preferred = [m for m in matches if m.name.startswith(preferred_prefix)]
        if len(preferred) == 1:
            return preferred[0]

    if len(matches) > 1:
        names = "\n  ".join(str(m) for m in matches)
        raise ValueError(f"Found multiple {label} matches:\n  {names}")

    return matches[0]


def infer_vfid(*values):
    """
    Infer a VFID string from paths or names.
    """
    for value in values:
        match = re.search(r"VFID\d+", str(value))
        if match:
            return match.group(0)
    return None


def infer_display_name(*values):
    """
    Infer a compact galaxy display name from paths or names.
    """
    for value in values:
        text = str(value)
        match = re.search(r"\b(NGC|IC)[-_ ]?(\d+[A-Za-z]?)\b", text, re.IGNORECASE)
        if match:
            return f"{match.group(1).upper()} {match.group(2)}"
    return None


def parse_csv_floats(value):
    """
    Parse a comma-separated list of floats.
    """
    if value is None:
        return None
    return [float(x) for x in value.split(",") if x.strip()]


def parse_csv_strings(value):
    """
    Parse a comma-separated list of strings.
    """
    if value is None:
        return None

    values = [x.strip() for x in value.split(",") if x.strip()]
    if len(values) == 1 and values[0].lower() in {"nircam", "nircam-full"}:
        return list(DEFAULT_NIRCAM_APERTURES)

    return values


def parse_crop(value):
    """
    Parse crop coordinates as xmin,xmax,ymin,ymax.
    """
    if value is None:
        return None

    vals = [int(float(x)) for x in value.split(",") if x.strip()]
    if len(vals) != 4:
        raise ValueError("--crop must have four values: xmin,xmax,ymin,ymax")

    xmin, xmax, ymin, ymax = vals
    if xmax <= xmin or ymax <= ymin:
        raise ValueError("--crop must satisfy xmax > xmin and ymax > ymin")

    return xmin, xmax, ymin, ymax


def format_param(value):
    """
    Format a parameter value compactly for filenames.
    """
    if value is None:
        return "none"

    if isinstance(value, float):
        text = f"{value:g}"
    else:
        text = str(value)

    text = text.replace(".", "p")
    text = re.sub(r"[^A-Za-z0-9_-]+", "", text)
    return text or "none"


def make_output_name(
    tag,
    recipe,
    crop=None,
    optical_stretch="legacy",
    halpha_color="magenta",
    halpha_alpha=0.55,
    halpha_nsigma=3.0,
    halpha_smooth_fwhm=1.5,
    halpha_min_area=20,
    halpha_detection=True,
    hi_contours=True,
    hi_image=False,
    hi_image_alpha=0.45,
    hi_image_percentiles=(50.0, 99.5),
    hi_image_stretch="asinh",
    hi_image_asinh_a=0.08,
    hi_image_alpha_power=0.7,
    hi_image_cmap="hi-cyan",
    hi_image_interpolation="bilinear",
    hi_image_high_fade_start=None,
    hi_image_high_fade_end=None,
    jwst_footprint=False,
    jwst_dither_footprints=False,
    jwst_pa_v3=None,
    jwst_color="gold",
    jwst_alpha=0.95,
    jwst_linewidth=2.0,
    jwst_stroke_linewidth=None,
    jwst_fill_alpha=0.0,
    jwst_outer_outline=False,
    jwst_outer_color=None,
    jwst_outer_alpha=0.95,
    jwst_outer_linewidth=2.0,
    jwst_outer_stroke_linewidth=0.0,
    jwst_outer_outline_mode="per-dither",
    hide_axes=False,
    corner_label=None,
    corner_label_fontsize=18,
    scale_fontsize=12,
    scale_linewidth=2.0,
):
    """
    Build an output filename that records the main visualization parameters.
    """
    recipe_label = {
        "legacy-grz-halpha-magenta": "legacy-grz-haoverlay",
    }.get(recipe, recipe)
    parts = [tag, recipe_label]
    parts.append(f"opt-{format_param(optical_stretch)}")

    if recipe != "legacy-grz":
        parts.extend(
            [
                f"ha-{format_param(halpha_color)}",
                f"a{format_param(halpha_alpha)}",
            ]
        )

        if halpha_detection:
            parts.extend(
                [
                    f"nsig{format_param(halpha_nsigma)}",
                    f"fwhm{format_param(halpha_smooth_fwhm)}",
                    f"area{format_param(halpha_min_area)}",
                ]
            )
        else:
            parts.append("nohadetect")

    if hi_contours:
        parts.append("hi")
    if hi_image:
        parts.extend(
            [
                "hiimg",
                f"hia{format_param(hi_image_alpha)}",
                (
                    "hip"
                    f"{format_param(hi_image_percentiles[0])}-"
                    f"{format_param(hi_image_percentiles[1])}"
                ),
                f"hist{format_param(hi_image_stretch)}",
                f"hiaa{format_param(hi_image_asinh_a)}",
                f"hiap{format_param(hi_image_alpha_power)}",
                f"hii{format_param(hi_image_interpolation)}",
            ]
        )
        if hi_image_high_fade_start is not None:
            parts.append(f"hifs{format_param(hi_image_high_fade_start)}")
        if hi_image_high_fade_end is not None:
            parts.append(f"hife{format_param(hi_image_high_fade_end)}")
        if hi_image_cmap != "hi-cyan":
            parts.append(f"hicmap{format_param(hi_image_cmap)}")

    if jwst_footprint or jwst_dither_footprints:
        parts.append("nircamdith" if jwst_dither_footprints else "nircam")
        if jwst_pa_v3 is not None:
            parts.append(f"pa{format_param(jwst_pa_v3)}")
        if jwst_color != "gold":
            parts.append(f"jwst{format_param(jwst_color)}")
        if jwst_alpha != 0.95:
            parts.append(f"ja{format_param(jwst_alpha)}")
        if jwst_linewidth != 2.0:
            parts.append(f"jlw{format_param(jwst_linewidth)}")
        if jwst_stroke_linewidth is not None:
            parts.append(f"jsl{format_param(jwst_stroke_linewidth)}")
        if jwst_fill_alpha > 0:
            parts.append(f"jfa{format_param(jwst_fill_alpha)}")
        if jwst_outer_outline:
            parts.append("jouter")
            if jwst_outer_color is not None:
                parts.append(f"joc{format_param(jwst_outer_color)}")
            if jwst_outer_alpha != 0.95:
                parts.append(f"joa{format_param(jwst_outer_alpha)}")
            if jwst_outer_linewidth != 2.0:
                parts.append(f"jolw{format_param(jwst_outer_linewidth)}")
            if jwst_outer_stroke_linewidth != 0.0:
                parts.append(f"josl{format_param(jwst_outer_stroke_linewidth)}")
            if jwst_outer_outline_mode != "per-dither":
                parts.append(f"jom{format_param(jwst_outer_outline_mode)}")

    if hide_axes:
        parts.append("noaxes")

    if corner_label:
        parts.append("label")
        if corner_label_fontsize != 18:
            parts.append(f"lfs{format_param(corner_label_fontsize)}")

    if scale_linewidth != 2.0:
        parts.append(f"slw{format_param(scale_linewidth)}")

    if scale_fontsize != 12:
        parts.append(f"sfs{format_param(scale_fontsize)}")

    if crop is not None:
        parts.append("crop" + "-".join(str(int(x)) for x in crop))
    else:
        parts.append("full")

    return "-".join(parts) + ".png"


def read_image(filename):
    """
    Read a FITS image and return 2D float data plus its header.
    """
    data, header = fits.getdata(filename, header=True)
    data = np.squeeze(np.asarray(data, dtype=float))

    if data.ndim != 2:
        raise ValueError(f"Expected a 2D image in {filename}; found shape {data.shape}")

    return data, header


def read_mask(filename, shape, mode="nonzero"):
    """
    Read a FITS mask and return a boolean array.

    The default mode matches the way masks are used in hafunctions.py: nonzero
    pixels are masked.
    """
    mask_data = np.squeeze(np.asarray(fits.getdata(filename)))

    if mask_data.shape != shape:
        raise ValueError(
            f"Mask shape {mask_data.shape} does not match image shape {shape}: {filename}"
        )

    if mode == "nonzero":
        mask = mask_data != 0
    elif mode == "zero":
        mask = mask_data == 0
    else:
        raise ValueError(f"Unknown mask mode: {mode}")

    return mask


def resolve_mask_file(galaxy_dir, mask_glob=None, mask_file=None, preferred_prefix=None):
    """
    Resolve an optional Halpha mask file.
    """
    if mask_file is not None:
        return expand_path(mask_file)

    if mask_glob is None:
        return None

    return find_one(
        galaxy_dir,
        mask_glob,
        "Halpha mask",
        preferred_prefix=preferred_prefix,
    )


def scale_channel(
    data,
    lower_percent=1.0,
    upper_percent=99.5,
    stretch="asinh",
    asinh_a=0.08,
    clip_negative=False,
):
    """
    Scale one image channel to [0, 1].
    """
    arr = np.asarray(data, dtype=float).copy()

    if clip_negative:
        arr[arr < 0] = 0

    finite = np.isfinite(arr)
    if not np.any(finite):
        return np.zeros_like(arr, dtype=float)

    lo, hi = np.nanpercentile(arr[finite], [lower_percent, upper_percent])

    if not np.isfinite(lo):
        lo = np.nanmin(arr[finite])

    if not np.isfinite(hi) or hi <= lo:
        hi = np.nanmax(arr[finite])

    if not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(arr, dtype=float)

    scaled = np.clip((arr - lo) / (hi - lo), 0, 1)
    scaled[~finite] = 0

    if stretch == "asinh":
        scaled = AsinhStretch(a=asinh_a)(scaled)
    elif stretch == "linear":
        scaled = LinearStretch()(scaled)
    else:
        raise ValueError(f"Unknown stretch: {stretch}")

    return np.asarray(np.clip(scaled, 0, 1), dtype=float)


def legacy_rgb(g_data, r_data, z_data, m=0.03, q=20.0, stretch_factor=1.0):
    """
    Reproduce the Legacy Survey grz RGB stretch used by legacypipe.survey.

    This is adapted from legacypipe.survey.sdss_rgb for bands ['g', 'r', 'z'].
    The images are assumed to be in Legacy Survey flux units.
    """
    bands = ["g", "r", "z"]
    imgs = [g_data, r_data, z_data]
    rgbscales = {
        "g": (2, 6.0 * stretch_factor),
        "r": (1, 3.4 * stretch_factor),
        "z": (0, 2.2 * stretch_factor),
    }

    intensity = np.zeros_like(np.asarray(g_data, dtype=float), dtype=float)
    scaled_imgs = []

    for img, band in zip(imgs, bands):
        _, scale = rgbscales[band]
        scaled = np.maximum(0, np.asarray(img, dtype=float) * scale + m)
        scaled[~np.isfinite(scaled)] = 0
        scaled_imgs.append(scaled)
        intensity += scaled

    intensity /= len(bands)

    if q is not None:
        f_intensity = np.arcsinh(q * intensity) / np.sqrt(q)
        intensity = intensity + (intensity == 0.0) * 1e-6
        intensity = f_intensity / intensity

    h, w = intensity.shape
    rgb = np.zeros((h, w, 3), dtype=np.float32)

    for img, band in zip(imgs, bands):
        plane, scale = rgbscales[band]
        imgplane = (np.asarray(img, dtype=float) * scale + m) * intensity
        imgplane[~np.isfinite(imgplane)] = 0
        rgb[:, :, plane] = np.clip(imgplane, 0, 1)

    return rgb


def remove_small_regions(detection, min_area=0):
    """
    Remove connected detections smaller than min_area pixels.
    """
    if min_area is None or min_area <= 1:
        return detection

    labels, nlabels = ndimage.label(detection)
    if nlabels == 0:
        return detection

    counts = np.bincount(labels.ravel())
    keep = counts >= int(min_area)
    keep[0] = False

    return keep[labels]


def make_halpha_detection_mask(
    halpha_data,
    base_mask=None,
    nsigma=3.0,
    smooth_fwhm=1.5,
    min_area=20,
    verbose=False,
):
    """
    Build a boolean mask for pixels that should not contribute Halpha color.

    The detection threshold is estimated from the unmasked Halpha image with
    sigma-clipped statistics.  Smoothing is used only to define detected
    regions; the original Halpha values are still used for the magenta overlay.
    """
    arr = np.asarray(halpha_data, dtype=float)
    invalid = ~np.isfinite(arr)

    if base_mask is not None:
        invalid = invalid | np.asarray(base_mask, dtype=bool)

    stats_data = arr[~invalid]
    stats_data = stats_data[np.isfinite(stats_data)]

    if stats_data.size == 0:
        return np.ones(arr.shape, dtype=bool)

    mean, median, std = sigma_clipped_stats(stats_data, sigma=3.0, maxiters=5)
    threshold = float(median) + float(nsigma) * float(std)

    detect_image = np.where(invalid, np.nan, arr)
    detect_image = np.where(detect_image > 0, detect_image, 0.0)

    if smooth_fwhm is not None and smooth_fwhm > 0:
        stddev = float(smooth_fwhm) / 2.3548
        kernel = Gaussian2DKernel(stddev)
        detect_image = convolve_fft(
            detect_image,
            kernel,
            nan_treatment="interpolate",
            preserve_nan=True,
            allow_huge=True,
        )

    detection = np.isfinite(detect_image) & (detect_image > threshold) & ~invalid
    detection = remove_small_regions(detection, min_area=min_area)

    if verbose:
        print(
            "Halpha detection: "
            f"mean={mean:.4g}, median={median:.4g}, sigma={std:.4g}, "
            f"threshold={threshold:.4g}, detected={np.mean(detection):.4f}"
        )

    return ~detection


def build_rgb(
    g_data,
    r_data,
    z_data,
    halpha_data,
    recipe="legacy-grz-halpha-magenta",
    optical_stretch="legacy",
    optical_percentiles=(1.0, 99.5),
    legacy_rgb_m=0.03,
    legacy_rgb_q=20.0,
    legacy_rgb_stretch_factor=1.0,
    halpha_percentiles=(1.0, 99.5),
    optical_asinh_a=0.08,
    halpha_asinh_a=0.03,
    halpha_alpha=0.55,
    halpha_color="magenta",
    halpha_mask=None,
    halpha_detection_mask=None,
):
    """
    Construct an RGB array from aligned image channels.
    """
    opt_lo, opt_hi = optical_percentiles
    ha_lo, ha_hi = halpha_percentiles

    if optical_stretch == "legacy":
        optical_rgb = legacy_rgb(
            g_data,
            r_data,
            z_data,
            m=legacy_rgb_m,
            q=legacy_rgb_q,
            stretch_factor=legacy_rgb_stretch_factor,
        )
    elif optical_stretch == "percentile":
        g = scale_channel(
            g_data,
            lower_percent=opt_lo,
            upper_percent=opt_hi,
            asinh_a=optical_asinh_a,
        )
        r = scale_channel(
            r_data,
            lower_percent=opt_lo,
            upper_percent=opt_hi,
            asinh_a=optical_asinh_a,
        )
        z = scale_channel(
            z_data,
            lower_percent=opt_lo,
            upper_percent=opt_hi,
            asinh_a=optical_asinh_a,
        )
        optical_rgb = np.dstack([z, r, g])
    else:
        raise ValueError(f"Unknown optical stretch: {optical_stretch}")

    if halpha_mask is not None:
        halpha_data = np.where(halpha_mask, np.nan, halpha_data)
    if halpha_detection_mask is not None:
        halpha_data = np.where(halpha_detection_mask, np.nan, halpha_data)

    ha = scale_channel(
        halpha_data,
        lower_percent=ha_lo,
        upper_percent=ha_hi,
        asinh_a=halpha_asinh_a,
        clip_negative=True,
    )

    if recipe == "legacy-grz":
        return optical_rgb

    if recipe == "halpha-r-g":
        g = scale_channel(
            g_data,
            lower_percent=opt_lo,
            upper_percent=opt_hi,
            asinh_a=optical_asinh_a,
        )
        r = scale_channel(
            r_data,
            lower_percent=opt_lo,
            upper_percent=opt_hi,
            asinh_a=optical_asinh_a,
        )
        return np.dstack([ha, r, g])

    if recipe != "legacy-grz-halpha-magenta":
        raise ValueError(f"Unknown recipe: {recipe}")

    rgb = optical_rgb
    alpha = np.clip(ha * halpha_alpha, 0, 1)
    halpha_rgb = np.asarray(mcolors.to_rgb(halpha_color), dtype=float)

    return rgb * (1 - alpha[:, :, None]) + halpha_rgb * alpha[:, :, None]


def add_scale_bar(
    ax,
    header,
    redshift=None,
    vr=None,
    barsize=5.0,
    color="white",
    fontsize=12,
    xscale=0.08,
    yscale=0.08,
    dytext=0.06,
    capfrac=0.18,
    linewidth=2.0,
):
    """
    Add a physical scale bar to an image axis.

    This follows the add_scale helper in havirgo/python/hafunctions.py, with the
    recession velocity inferred from REDSHIFT when vr is not supplied.
    """
    if vr is None:
        if redshift is None:
            redshift = header.get("REDSHIFT")
        if redshift is None:
            raise ValueError("Could not add scale bar: no REDSHIFT in header")
        vr = float(redshift) * 3.0e5

    z = float(vr) / 3.0e5
    angular_distance = cosmo.angular_diameter_distance(z)
    kpc_per_arcsec = angular_distance.value * 1000.0 * np.pi / (180.0 * 3600.0)

    if kpc_per_arcsec <= 0:
        raise ValueError(f"Bad physical scale for vr={vr}")

    pscale = (
        wcs_utils.proj_plane_pixel_scales(WCS(header).celestial) * u.deg
    ).to(u.arcsec)
    pscale_arcsec = float(np.abs(pscale[0].value))

    barsize_arcsec = float(barsize) / kpc_per_arcsec
    barsize_pixels = barsize_arcsec / pscale_arcsec

    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    dx = x2 - x1
    dy = y2 - y1

    xline1 = x1 + xscale * dx
    xline2 = xline1 + barsize_pixels
    yline = y1 + yscale * dy

    capsize_pixels = capfrac * barsize_pixels
    ycap1 = yline - 0.5 * capsize_pixels
    ycap2 = yline + 0.5 * capsize_pixels

    xtext = xline1 + 0.5 * barsize_pixels
    ytext = yline - dytext * dy

    ax.plot([xline1, xline2], [yline, yline], "-", lw=linewidth, color=color, zorder=9)
    ax.plot([xline1, xline1], [ycap1, ycap2], "-", lw=linewidth, color=color, zorder=9)
    ax.plot([xline2, xline2], [ycap1, ycap2], "-", lw=linewidth, color=color, zorder=9)
    ax.text(
        xtext,
        ytext,
        f"{barsize:.0f} kpc",
        horizontalalignment="center",
        fontsize=fontsize,
        color=color,
        path_effects=[
            path_effects.withStroke(linewidth=max(1.0, linewidth * 0.6), foreground="black")
        ],
        zorder=10,
    )


def add_corner_label(
    ax,
    label,
    position="upper left",
    color="white",
    fontsize=18,
    pad=0.035,
):
    """
    Add a readable in-image label in axis-fraction coordinates.
    """
    positions = {
        "upper left": (pad, 1.0 - pad, "left", "top"),
        "upper right": (1.0 - pad, 1.0 - pad, "right", "top"),
        "lower left": (pad, pad, "left", "bottom"),
        "lower right": (1.0 - pad, pad, "right", "bottom"),
    }
    if position not in positions:
        raise ValueError(f"Unknown corner-label position: {position}")

    x, y, ha, va = positions[position]
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        color=color,
        fontsize=fontsize,
        fontweight="semibold",
        path_effects=[path_effects.withStroke(linewidth=3.0, foreground="black")],
        zorder=10,
    )


def hide_wcs_axes(ax):
    """
    Hide WCS axes, ticks, labels, and frame for clean presentation figures.
    """
    ax.set_axis_off()
    for coord in ax.coords:
        coord.set_axislabel("")
        coord.set_ticklabel_visible(False)
        coord.set_ticks_visible(False)


def resolve_hi_file(vfid, hi_file=None, use_hi=True, verbose=False):
    """
    Resolve an HI contour file from an explicit path or the default VFID map.
    """
    if not use_hi:
        return None

    if hi_file is not None:
        path = expand_path(hi_file)
    elif vfid in DEFAULT_HI_FILES:
        path = expand_path(DEFAULT_HI_FILES[vfid])
    else:
        if verbose:
            print(f"No default HI file for {vfid}")
        return None

    if not path.exists():
        if verbose:
            print(f"HI file does not exist: {path}")
        return None

    return path


def plot_hi_contours(
    ax,
    hi_file,
    levels=None,
    color="lightskyblue",
    alpha=0.75,
    linewidth=1.0,
    ncontour=5,
):
    """
    Add HI contours to a WCSAxes image.
    """
    hi_data, hi_header = read_image(hi_file)
    hi_wcs = WCS(hi_header).celestial

    if levels is None:
        levels = 3 ** np.arange(1, ncontour + 1) + 1

    ax.contour(
        hi_data,
        levels=levels,
        colors=color,
        alpha=alpha,
        linewidths=linewidth,
        transform=ax.get_transform(hi_wcs),
    )


def positive_hi_scaled_image(
    hi_data,
    lower_percent=50.0,
    upper_percent=99.5,
    stretch="asinh",
    asinh_a=0.08,
):
    """
    Return a [0, 1] image scaled from positive HI emission only.
    """
    arr = np.asarray(hi_data, dtype=float)
    positive = np.isfinite(arr) & (arr > 0)
    scaled_input = np.where(positive, arr, np.nan)

    return scale_channel(
        scaled_input,
        lower_percent=lower_percent,
        upper_percent=upper_percent,
        stretch=stretch,
        asinh_a=asinh_a,
        clip_negative=True,
    )


def hi_cmap_from_name(name):
    """
    Resolve a colormap for filled HI overlays.
    """
    if name == "hi-cyan":
        return mcolors.LinearSegmentedColormap.from_list(
            "hi-cyan",
            ["#002f73", "#006db5", "#00d5ff", "#d8ffff"],
        )

    return plt.get_cmap(name)


def plot_hi_image(
    ax,
    hi_file,
    alpha=0.45,
    percentiles=(50.0, 99.5),
    stretch="asinh",
    asinh_a=0.08,
    alpha_power=0.7,
    cmap="hi-cyan",
    interpolation="bilinear",
    high_fade_start=None,
    high_fade_end=None,
):
    """
    Add a semi-transparent HI moment-zero intensity image to a WCSAxes image.
    """
    hi_data, hi_header = read_image(hi_file)
    hi_wcs = WCS(hi_header).celestial
    hi_lo, hi_hi = percentiles

    scaled = positive_hi_scaled_image(
        hi_data,
        lower_percent=hi_lo,
        upper_percent=hi_hi,
        stretch=stretch,
        asinh_a=asinh_a,
    )
    positive = np.isfinite(hi_data) & (hi_data > 0)
    scaled = np.where(positive, scaled, 0.0)
    alpha_image = np.clip(scaled, 0, 1) ** float(alpha_power)
    if high_fade_start is not None or high_fade_end is not None:
        if high_fade_start is None or high_fade_end is None:
            raise ValueError(
                "Set both --hi-image-high-fade-start and "
                "--hi-image-high-fade-end, or neither."
            )
        fade_start = float(high_fade_start)
        fade_end = float(high_fade_end)
        if fade_end <= fade_start:
            raise ValueError("--hi-image-high-fade-end must be larger than start")
        high_frac = np.clip((scaled - fade_start) / (fade_end - fade_start), 0, 1)
        high_fade = high_frac * high_frac * (3 - 2 * high_frac)
        alpha_image *= 1 - high_fade
    alpha_image = np.where(positive, alpha * alpha_image, 0.0)

    ax.imshow(
        scaled,
        origin="lower",
        cmap=hi_cmap_from_name(cmap),
        alpha=alpha_image,
        interpolation=interpolation,
        transform=ax.get_transform(hi_wcs),
        zorder=2.5,
    )


def resolve_jwst_center(header, data_shape, crop=None, ra=None, dec=None):
    """
    Resolve the JWST footprint center in world coordinates.
    """
    if ra is not None or dec is not None:
        if ra is None or dec is None:
            raise ValueError("Set both --jwst-ra and --jwst-dec, or neither.")
        return float(ra), float(dec)

    if crop is not None:
        xmin, xmax, ymin, ymax = crop
        xcenter = 0.5 * (xmin + xmax)
        ycenter = 0.5 * (ymin + ymax)
    else:
        ny, nx = data_shape
        xcenter = 0.5 * (nx - 1)
        ycenter = 0.5 * (ny - 1)

    return WCS(header).celestial.pixel_to_world_values(xcenter, ycenter)


def build_jwst_footprints(
    ra,
    dec,
    pa_v3=0.0,
    apertures=None,
    reference_aperture="NRCALL_FULL",
    reference_v2=None,
    reference_v3=None,
):
    """
    Build JWST/NIRCam detector footprint polygons using pysiaf.
    """
    try:
        import pysiaf
        from pysiaf.utils import rotations
    except Exception as err:
        raise ImportError(
            "Could not import pysiaf. Install it in the active Python environment "
            "or rerun without --jwst-footprint."
        ) from err

    if apertures is None:
        apertures = DEFAULT_NIRCAM_APERTURES

    siaf = pysiaf.Siaf("NIRCam")
    reference = siaf[reference_aperture]
    if reference_v2 is None:
        reference_v2 = reference.V2Ref
    if reference_v3 is None:
        reference_v3 = reference.V3Ref
    attitude = rotations.attitude_matrix(
        float(reference_v2),
        float(reference_v3),
        float(ra),
        float(dec),
        float(pa_v3),
    )

    footprints = []
    for aperture_name in apertures:
        aperture = siaf[aperture_name]
        aperture.set_attitude_matrix(attitude)
        footprint_ra, footprint_dec = aperture.closed_polygon_points("sky")
        footprints.append((aperture_name, footprint_ra, footprint_dec))

    return footprints


def plot_jwst_footprint_polygons(
    ax,
    footprints,
    color="gold",
    alpha=0.95,
    linewidth=2.0,
    stroke_linewidth=None,
    fill_alpha=0.0,
    outer_outline=False,
    outer_color=None,
    outer_alpha=0.95,
    outer_linewidth=2.0,
    outer_stroke_linewidth=0.0,
):
    """
    Draw precomputed JWST footprint polygons on a WCSAxes image.
    """
    transform = ax.get_transform("world")

    for _, footprint_ra, footprint_dec in footprints:
        if fill_alpha > 0:
            ax.fill(
                footprint_ra,
                footprint_dec,
                color=color,
                alpha=fill_alpha,
                transform=transform,
                zorder=6,
            )
        if stroke_linewidth is None:
            stroke_width = max(1.0, linewidth + 1.5)
        else:
            stroke_width = stroke_linewidth
        stroke_effects = []
        if stroke_width > 0:
            stroke_effects.append(
                path_effects.withStroke(linewidth=stroke_width, foreground="black")
            )

        if linewidth > 0 and alpha > 0:
            ax.plot(
                footprint_ra,
                footprint_dec,
                color=color,
                alpha=alpha,
                linewidth=linewidth,
                transform=transform,
                zorder=7,
                path_effects=stroke_effects,
            )

    if outer_outline:
        outline_ra, outline_dec = jwst_outer_outline(footprints)
        outline_color = outer_color or color
        outer_effects = []
        if outer_stroke_linewidth > 0:
            outer_effects.append(
                path_effects.withStroke(
                    linewidth=outer_stroke_linewidth,
                    foreground="black",
                )
            )
        ax.plot(
            outline_ra,
            outline_dec,
            color=outline_color,
            alpha=outer_alpha,
            linewidth=outer_linewidth,
            transform=transform,
            zorder=8,
            path_effects=outer_effects,
        )


def jwst_outer_outline(footprints):
    """
    Return a convex-hull outline around a set of JWST detector polygons.
    """
    points = []
    for _, footprint_ra, footprint_dec in footprints:
        ra = np.asarray(footprint_ra, dtype=float)
        dec = np.asarray(footprint_dec, dtype=float)
        if ra.size > 1 and dec.size > 1:
            ra = ra[:-1]
            dec = dec[:-1]
        finite = np.isfinite(ra) & np.isfinite(dec)
        points.append(np.column_stack([ra[finite], dec[finite]]))

    if not points:
        raise ValueError("No finite JWST footprint points available for outer outline")

    points = np.vstack(points)
    if points.shape[0] < 3:
        raise ValueError("At least three JWST footprint points are needed for outline")

    dec0 = np.nanmean(points[:, 1])
    scaled = np.column_stack(
        [
            (points[:, 0] - np.nanmean(points[:, 0])) * np.cos(np.deg2rad(dec0)),
            points[:, 1] - dec0,
        ]
    )
    hull = ConvexHull(scaled)
    indices = np.append(hull.vertices, hull.vertices[0])
    return points[indices, 0], points[indices, 1]


def plot_jwst_footprints(
    ax,
    ra,
    dec,
    pa_v3=0.0,
    apertures=None,
    reference_aperture="NRCALL_FULL",
    color="gold",
    alpha=0.95,
    linewidth=2.0,
    stroke_linewidth=None,
    fill_alpha=0.0,
    outer_outline=False,
    outer_color=None,
    outer_alpha=0.95,
    outer_linewidth=2.0,
    outer_stroke_linewidth=0.0,
):
    """
    Overlay JWST/NIRCam footprint polygons on a WCSAxes image.
    """
    footprints = build_jwst_footprints(
        ra,
        dec,
        pa_v3=pa_v3,
        apertures=apertures,
        reference_aperture=reference_aperture,
    )
    plot_jwst_footprint_polygons(
        ax,
        footprints,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
        stroke_linewidth=stroke_linewidth,
        fill_alpha=fill_alpha,
        outer_outline=outer_outline,
        outer_color=outer_color,
        outer_alpha=outer_alpha,
        outer_linewidth=outer_linewidth,
        outer_stroke_linewidth=outer_stroke_linewidth,
    )


def parse_apt_pointing_file(pointing_file, target=None, pointing_types=None):
    """
    Parse unique NIRCam pointing/dither rows from an APT pointing report.

    The APT report repeats the same telescope offsets for each exposure/filter.
    Rows are de-duplicated by target, tile, dither, aperture, and V2/V3 position.
    """
    pointing_file = expand_path(pointing_file)
    if pointing_types is not None:
        pointing_types = {x.upper() for x in pointing_types}

    rows = []
    seen = set()
    current_target = None
    current_visit = None
    section_re = re.compile(r"^\*\s+(\S+)\s+\(Obs\s+([^)]+)\)")
    visit_re = re.compile(r"^\*\*\s+Visit\s+(\S+)")

    with pointing_file.open() as handle:
        for line in handle:
            section_match = section_re.match(line)
            if section_match:
                current_target = section_match.group(1)
                current_visit = None
                continue

            visit_match = visit_re.match(line)
            if visit_match:
                current_visit = visit_match.group(1)
                continue

            tokens = line.split()
            if len(tokens) < 22 or not tokens[0].isdigit():
                continue

            aperture_name = tokens[4]
            if not aperture_name.startswith("NRC"):
                continue

            row_target = tokens[6]
            if target is not None:
                target_text = str(target)
                matches_target = (
                    current_target == target_text
                    or row_target == target_text
                    or (
                        row_target.startswith("VFID")
                        and len(row_target) >= 7
                        and target_text.startswith(row_target)
                    )
                )
                if not matches_target:
                    continue

            try:
                row = {
                    "section_target": current_target,
                    "visit": current_visit,
                    "tar": int(tokens[0]),
                    "tile": int(tokens[1]),
                    "exp": int(tokens[2]),
                    "dither": int(tokens[3]),
                    "aperture": aperture_name,
                    "target_order": tokens[5],
                    "target": row_target,
                    "ra": float(tokens[7]),
                    "dec": float(tokens[8]),
                    "basex": float(tokens[9]),
                    "basey": float(tokens[10]),
                    "dithx": float(tokens[11]),
                    "dithy": float(tokens[12]),
                    "v2": float(tokens[13]),
                    "v3": float(tokens[14]),
                    "idlx": float(tokens[15]),
                    "idly": float(tokens[16]),
                    "level": tokens[17].upper(),
                    "type": tokens[18].upper(),
                }
            except (IndexError, ValueError):
                continue

            if pointing_types is not None and row["level"] not in pointing_types:
                continue

            key = (
                row["section_target"],
                row["target"],
                row["tile"],
                row["dither"],
                row["aperture"],
                round(row["ra"], 8),
                round(row["dec"], 8),
                round(row["v2"], 3),
                round(row["v3"], 3),
            )
            if key in seen:
                continue

            seen.add(key)
            rows.append(row)

    if not rows:
        target_text = f" for {target}" if target is not None else ""
        raise ValueError(f"No NIRCam pointing rows found{target_text} in {pointing_file}")

    return rows


def plot_jwst_dither_footprints(
    ax,
    pointing_file,
    target,
    pa_v3=0.0,
    apertures=None,
    reference_aperture="NRCALL_FULL",
    pointing_types=None,
    color="gold",
    alpha=0.35,
    linewidth=1.2,
    stroke_linewidth=0.0,
    fill_alpha=0.04,
    outer_outline=False,
    outer_color=None,
    outer_alpha=0.95,
    outer_linewidth=2.0,
    outer_stroke_linewidth=0.0,
    outer_outline_mode="per-dither",
    verbose=False,
):
    """
    Overlay NIRCam detector footprints for each unique APT dither/tile position.
    """
    pointings = parse_apt_pointing_file(
        pointing_file,
        target=target,
        pointing_types=pointing_types,
    )
    if verbose:
        print(f"APT NIRCam pointings for {target}: {len(pointings)}")
        for row in pointings:
            print(
                "  "
                f"tile={row['tile']} dith={row['dither']} "
                f"ap={row['aperture']} level={row['level']} "
                f"V2={row['v2']:.3f} V3={row['v3']:.3f}"
            )

    if outer_outline_mode not in {"per-dither", "combined"}:
        raise ValueError("--jwst-outer-outline-mode must be per-dither or combined")

    all_footprints = []
    for row in pointings:
        footprints = build_jwst_footprints(
            row["ra"],
            row["dec"],
            pa_v3=pa_v3,
            apertures=apertures,
            reference_aperture=reference_aperture,
            reference_v2=row["v2"],
            reference_v3=row["v3"],
        )
        all_footprints.extend(footprints)
        plot_jwst_footprint_polygons(
            ax,
            footprints,
            color=color,
            alpha=alpha,
            linewidth=linewidth,
            stroke_linewidth=stroke_linewidth,
            fill_alpha=fill_alpha,
            outer_outline=outer_outline and outer_outline_mode == "per-dither",
            outer_color=outer_color,
            outer_alpha=outer_alpha,
            outer_linewidth=outer_linewidth,
            outer_stroke_linewidth=outer_stroke_linewidth,
        )

    if outer_outline and outer_outline_mode == "combined":
        plot_jwst_footprint_polygons(
            ax,
            all_footprints,
            color=color,
            alpha=0.0,
            linewidth=0.0,
            stroke_linewidth=0.0,
            fill_alpha=0.0,
            outer_outline=True,
            outer_color=outer_color,
            outer_alpha=outer_alpha,
            outer_linewidth=outer_linewidth,
            outer_stroke_linewidth=outer_stroke_linewidth,
        )


def make_legacy_halpha_rgb(
    galaxy_dir,
    halpha_glob="*-sfr-vr.fits",
    legacy_dir="legacy",
    g_glob="*-g-ha.fits",
    r_glob="*-r-ha.fits",
    z_glob="*-z-ha.fits",
    mask_glob="*-R-mask.fits",
    mask_file=None,
    mask_mode="nonzero",
    outfile=None,
    output_dir=None,
    recipe="legacy-grz-halpha-magenta",
    crop=None,
    optical_stretch="legacy",
    optical_percentiles=(1.0, 99.5),
    legacy_rgb_m=0.03,
    legacy_rgb_q=20.0,
    legacy_rgb_stretch_factor=1.0,
    halpha_percentiles=(1.0, 99.5),
    optical_asinh_a=0.08,
    halpha_asinh_a=0.03,
    halpha_alpha=0.55,
    halpha_color="magenta",
    halpha_nsigma=3.0,
    halpha_smooth_fwhm=1.5,
    halpha_min_area=20,
    halpha_detection=True,
    hi_contours=True,
    hi_file=None,
    hi_levels=None,
    hi_color="lightskyblue",
    hi_alpha=0.75,
    hi_linewidth=1.0,
    hi_image=False,
    hi_image_alpha=0.45,
    hi_image_percentiles=(50.0, 99.5),
    hi_image_stretch="asinh",
    hi_image_asinh_a=0.08,
    hi_image_alpha_power=0.7,
    hi_image_cmap="hi-cyan",
    hi_image_interpolation="bilinear",
    hi_image_high_fade_start=None,
    hi_image_high_fade_end=None,
    jwst_footprint=False,
    jwst_ra=None,
    jwst_dec=None,
    jwst_pa_v3=0.0,
    jwst_apertures=None,
    jwst_reference_aperture="NRCALL_FULL",
    jwst_color="gold",
    jwst_alpha=0.95,
    jwst_linewidth=2.0,
    jwst_stroke_linewidth=None,
    jwst_fill_alpha=0.0,
    jwst_outer_outline=False,
    jwst_outer_color=None,
    jwst_outer_alpha=0.95,
    jwst_outer_linewidth=2.0,
    jwst_outer_stroke_linewidth=0.0,
    jwst_outer_outline_mode="per-dither",
    jwst_dither_footprints=False,
    jwst_pointing_file=None,
    jwst_pointing_target=None,
    jwst_pointing_types=None,
    add_scale=True,
    scale_kpc=5.0,
    scale_color="white",
    scale_fontsize=12,
    scale_x=0.08,
    scale_y=0.08,
    scale_dytext=0.06,
    scale_linewidth=2.0,
    title=None,
    show_title=True,
    hide_axes=False,
    corner_label=None,
    corner_label_position="upper left",
    corner_label_color="white",
    corner_label_fontsize=18,
    corner_label_pad=0.035,
    figsize=(8.0, 8.0),
    dpi=200,
    verbose=False,
):
    """
    Make a false-color RGB image for one galaxy directory.
    """
    galaxy_dir = expand_path(galaxy_dir)
    legacy_dir = galaxy_dir / legacy_dir

    halpha_file = find_one(
        galaxy_dir,
        halpha_glob,
        "Halpha image",
        preferred_prefix=galaxy_dir.name,
    )
    g_file = find_one(legacy_dir, g_glob, "Legacy g image")
    r_file = find_one(legacy_dir, r_glob, "Legacy r image")
    z_file = find_one(legacy_dir, z_glob, "Legacy z image")

    vfid = infer_vfid(galaxy_dir.name, halpha_file.name, g_file.name)
    mask_file = resolve_mask_file(
        galaxy_dir,
        mask_glob=mask_glob,
        mask_file=mask_file,
        preferred_prefix=galaxy_dir.name,
    )
    crop = parse_crop(crop) if isinstance(crop, str) else crop

    if outfile is None:
        tag = vfid or galaxy_dir.name
        if output_dir is None:
            output_dir = galaxy_dir
        else:
            output_dir = expand_path(output_dir)
        outfile = output_dir / make_output_name(
            tag,
            recipe,
            crop=crop,
            optical_stretch=optical_stretch,
            halpha_color=halpha_color,
            halpha_alpha=halpha_alpha,
            halpha_nsigma=halpha_nsigma,
            halpha_smooth_fwhm=halpha_smooth_fwhm,
            halpha_min_area=halpha_min_area,
            halpha_detection=halpha_detection,
            hi_contours=hi_contours,
            hi_image=hi_image,
            hi_image_alpha=hi_image_alpha,
            hi_image_percentiles=hi_image_percentiles,
            hi_image_stretch=hi_image_stretch,
            hi_image_asinh_a=hi_image_asinh_a,
            hi_image_alpha_power=hi_image_alpha_power,
            hi_image_cmap=hi_image_cmap,
            hi_image_interpolation=hi_image_interpolation,
            hi_image_high_fade_start=hi_image_high_fade_start,
            hi_image_high_fade_end=hi_image_high_fade_end,
            jwst_footprint=jwst_footprint,
            jwst_dither_footprints=jwst_dither_footprints,
            jwst_pa_v3=jwst_pa_v3,
            jwst_color=jwst_color,
            jwst_alpha=jwst_alpha,
            jwst_linewidth=jwst_linewidth,
            jwst_stroke_linewidth=jwst_stroke_linewidth,
            jwst_fill_alpha=jwst_fill_alpha,
            jwst_outer_outline=jwst_outer_outline,
            jwst_outer_color=jwst_outer_color,
            jwst_outer_alpha=jwst_outer_alpha,
            jwst_outer_linewidth=jwst_outer_linewidth,
            jwst_outer_stroke_linewidth=jwst_outer_stroke_linewidth,
            jwst_outer_outline_mode=jwst_outer_outline_mode,
            hide_axes=hide_axes,
            corner_label=corner_label,
            corner_label_fontsize=corner_label_fontsize,
            scale_fontsize=scale_fontsize,
            scale_linewidth=scale_linewidth,
        )
    else:
        outfile = expand_path(outfile)

    if verbose:
        print(f"galaxy_dir = {galaxy_dir}")
        print(f"vfid = {vfid}")
        print(f"Halpha = {halpha_file}")
        print(f"Legacy g = {g_file}")
        print(f"Legacy r = {r_file}")
        print(f"Legacy z = {z_file}")
        print(f"Halpha mask = {mask_file}")
        print(f"outfile = {outfile}")

    halpha_data, halpha_header = read_image(halpha_file)
    g_data, _ = read_image(g_file)
    r_data, _ = read_image(r_file)
    z_data, _ = read_image(z_file)

    shapes = {arr.shape for arr in [halpha_data, g_data, r_data, z_data]}
    if len(shapes) != 1:
        raise ValueError(f"Input images must already be aligned; found shapes {shapes}")

    halpha_mask = None
    if mask_file is not None:
        halpha_mask = read_mask(mask_file, halpha_data.shape, mode=mask_mode)
        if verbose:
            print(f"masked Halpha pixels = {np.mean(halpha_mask):.3f}")

    halpha_detection_mask = None
    if halpha_detection:
        halpha_detection_mask = make_halpha_detection_mask(
            halpha_data,
            base_mask=halpha_mask,
            nsigma=halpha_nsigma,
            smooth_fwhm=halpha_smooth_fwhm,
            min_area=halpha_min_area,
            verbose=verbose,
        )

    rgb = build_rgb(
        g_data,
        r_data,
        z_data,
        halpha_data,
        recipe=recipe,
        optical_stretch=optical_stretch,
        optical_percentiles=optical_percentiles,
        legacy_rgb_m=legacy_rgb_m,
        legacy_rgb_q=legacy_rgb_q,
        legacy_rgb_stretch_factor=legacy_rgb_stretch_factor,
        halpha_percentiles=halpha_percentiles,
        optical_asinh_a=optical_asinh_a,
        halpha_asinh_a=halpha_asinh_a,
        halpha_alpha=halpha_alpha,
        halpha_color=halpha_color,
        halpha_mask=halpha_mask,
        halpha_detection_mask=halpha_detection_mask,
    )

    image_wcs = WCS(halpha_header).celestial
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(projection=image_wcs)
    ax.imshow(rgb, origin="lower")

    if crop is not None:
        xmin, xmax, ymin, ymax = crop
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

    if hi_image:
        resolved_hi_file = resolve_hi_file(vfid, hi_file=hi_file, use_hi=True, verbose=verbose)
        if resolved_hi_file is not None:
            if verbose:
                print(f"HI image = {resolved_hi_file}")
            plot_hi_image(
                ax,
                resolved_hi_file,
                alpha=hi_image_alpha,
                percentiles=hi_image_percentiles,
                stretch=hi_image_stretch,
                asinh_a=hi_image_asinh_a,
                alpha_power=hi_image_alpha_power,
                cmap=hi_image_cmap,
                interpolation=hi_image_interpolation,
                high_fade_start=hi_image_high_fade_start,
                high_fade_end=hi_image_high_fade_end,
            )

    if hi_contours:
        resolved_hi_file = resolve_hi_file(vfid, hi_file=hi_file, use_hi=True, verbose=verbose)
        if resolved_hi_file is not None:
            if verbose:
                print(f"HI contours = {resolved_hi_file}")
            plot_hi_contours(
                ax,
                resolved_hi_file,
                levels=hi_levels,
                color=hi_color,
                alpha=hi_alpha,
                linewidth=hi_linewidth,
            )

    if jwst_dither_footprints:
        if jwst_pointing_file is None:
            raise ValueError("--jwst-pointing-file is required with --jwst-dither-footprints")
        if jwst_pointing_target is None:
            jwst_pointing_target = vfid
        try:
            if verbose:
                print(
                    "NIRCam dither footprints from "
                    f"{expand_path(jwst_pointing_file)} for {jwst_pointing_target}; "
                    f"PA_V3 = {jwst_pa_v3:g}"
                )
            plot_jwst_dither_footprints(
                ax,
                jwst_pointing_file,
                jwst_pointing_target,
                pa_v3=jwst_pa_v3,
                apertures=jwst_apertures,
                reference_aperture=jwst_reference_aperture,
                pointing_types=jwst_pointing_types,
                color=jwst_color,
                alpha=jwst_alpha,
                linewidth=jwst_linewidth,
                stroke_linewidth=jwst_stroke_linewidth,
                fill_alpha=jwst_fill_alpha,
                outer_outline=jwst_outer_outline,
                outer_color=jwst_outer_color,
                outer_alpha=jwst_outer_alpha,
                outer_linewidth=jwst_outer_linewidth,
                outer_stroke_linewidth=jwst_outer_stroke_linewidth,
                outer_outline_mode=jwst_outer_outline_mode,
                verbose=verbose,
            )
        except Exception as err:
            if verbose:
                print(f"WARNING: could not add JWST dither footprints: {err}")
            else:
                raise

    elif jwst_footprint:
        try:
            jwst_center_ra, jwst_center_dec = resolve_jwst_center(
                halpha_header,
                halpha_data.shape,
                crop=crop,
                ra=jwst_ra,
                dec=jwst_dec,
            )
            if verbose:
                print(
                    "NIRCam footprint center = "
                    f"{jwst_center_ra:.7f}, {jwst_center_dec:.7f}; "
                    f"PA_V3 = {jwst_pa_v3:g}"
                )
            plot_jwst_footprints(
                ax,
                jwst_center_ra,
                jwst_center_dec,
                pa_v3=jwst_pa_v3,
                apertures=jwst_apertures,
                reference_aperture=jwst_reference_aperture,
                color=jwst_color,
                alpha=jwst_alpha,
                linewidth=jwst_linewidth,
                stroke_linewidth=jwst_stroke_linewidth,
                fill_alpha=jwst_fill_alpha,
                outer_outline=jwst_outer_outline,
                outer_color=jwst_outer_color,
                outer_alpha=jwst_outer_alpha,
                outer_linewidth=jwst_outer_linewidth,
                outer_stroke_linewidth=jwst_outer_stroke_linewidth,
            )
        except Exception as err:
            if verbose:
                print(f"WARNING: could not add JWST footprint: {err}")
            else:
                raise

    if add_scale:
        try:
            add_scale_bar(
                ax,
                halpha_header,
                barsize=scale_kpc,
                color=scale_color,
                fontsize=scale_fontsize,
                xscale=scale_x,
                yscale=scale_y,
                dytext=scale_dytext,
                linewidth=scale_linewidth,
            )
        except Exception as err:
            if verbose:
                print(f"WARNING: could not add scale bar: {err}")

    if title is None:
        title = vfid

    if show_title and title:
        ax.set_title(title)

    if corner_label == "auto":
        corner_label = infer_display_name(galaxy_dir.name, halpha_file.name) or vfid

    if corner_label:
        add_corner_label(
            ax,
            corner_label,
            position=corner_label_position,
            color=corner_label_color,
            fontsize=corner_label_fontsize,
            pad=corner_label_pad,
        )

    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_axislabel("RA")
    lat.set_axislabel("Dec")
    lon.set_major_formatter("hh:mm:ss")
    lat.set_major_formatter("dd:mm")

    if hide_axes:
        hide_wcs_axes(ax)

    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    if verbose:
        print(f"wrote {outfile}")

    return outfile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Make a Legacy grz + Halpha false-color RGB image."
    )

    parser.add_argument("galaxy_dir", help="Galaxy cutout directory.")
    parser.add_argument(
        "--halpha-glob",
        default="*-sfr-vr.fits",
        help="Glob for the Halpha/SFR image inside galaxy_dir. Default: *-sfr-vr.fits",
    )
    parser.add_argument(
        "--legacy-dir",
        default="legacy",
        help="Legacy image subdirectory name. Default: legacy",
    )
    parser.add_argument("--g-glob", default="*-g-ha.fits", help="Legacy g-band glob.")
    parser.add_argument("--r-glob", default="*-r-ha.fits", help="Legacy r-band glob.")
    parser.add_argument("--z-glob", default="*-z-ha.fits", help="Legacy z-band glob.")
    parser.add_argument(
        "--mask-glob",
        default="*-R-mask.fits",
        help="Glob for Halpha mask inside galaxy_dir. Default: *-R-mask.fits",
    )
    parser.add_argument(
        "--mask-file",
        default=None,
        help="Optional explicit Halpha mask FITS file.",
    )
    parser.add_argument(
        "--mask-mode",
        default="nonzero",
        choices=["nonzero", "zero"],
        help="Mask convention. Default: nonzero pixels are masked.",
    )
    parser.add_argument(
        "--no-halpha-mask",
        action="store_true",
        help="Do not mask the Halpha channel.",
    )
    parser.add_argument("--outfile", default=None, help="Output PNG filename.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory for auto-named output when --outfile is omitted. "
            "Default: galaxy_dir."
        ),
    )
    parser.add_argument(
        "--recipe",
        default="legacy-grz-halpha-magenta",
        choices=["legacy-grz-halpha-magenta", "legacy-grz", "halpha-r-g"],
        help="RGB recipe. Default: legacy-grz-halpha-magenta",
    )
    parser.add_argument(
        "--crop",
        default=None,
        help="Optional crop as xmin,xmax,ymin,ymax in pixel coordinates.",
    )
    parser.add_argument(
        "--optical-stretch",
        default="legacy",
        choices=["legacy", "percentile"],
        help="Optical RGB stretch. Default: legacy.",
    )
    parser.add_argument(
        "--optical-percentiles",
        default="1,99.5",
        help=(
            "Lower,upper percentiles for Legacy channels when "
            "--optical-stretch percentile. Default: 1,99.5"
        ),
    )
    parser.add_argument(
        "--legacy-rgb-m",
        type=float,
        default=0.03,
        help="Legacy sdss_rgb m parameter. Default: 0.03.",
    )
    parser.add_argument(
        "--legacy-rgb-q",
        type=float,
        default=20.0,
        help="Legacy sdss_rgb Q parameter. Default: 20.",
    )
    parser.add_argument(
        "--legacy-rgb-stretch-factor",
        type=float,
        default=1.0,
        help="Multiplier for Legacy g/r/z RGB scales. Default: 1.",
    )
    parser.add_argument(
        "--halpha-percentiles",
        default="1,99.5",
        help="Lower,upper percentiles for detected Halpha pixels. Default: 1,99.5",
    )
    parser.add_argument(
        "--optical-asinh-a",
        type=float,
        default=0.08,
        help="Asinh a parameter for Legacy channels. Default: 0.08",
    )
    parser.add_argument(
        "--halpha-asinh-a",
        type=float,
        default=0.03,
        help="Asinh a parameter for the Halpha channel. Default: 0.03",
    )
    parser.add_argument(
        "--halpha-alpha",
        type=float,
        default=0.55,
        help="Maximum opacity of the Halpha overlay. Default: 0.55",
    )
    parser.add_argument(
        "--halpha-color",
        default="magenta",
        help="Matplotlib color for the Halpha overlay. Default: magenta.",
    )
    parser.add_argument(
        "--halpha-nsigma",
        type=float,
        default=3.0,
        help="Sigma threshold for Halpha detection mask. Default: 3.",
    )
    parser.add_argument(
        "--halpha-smooth-fwhm",
        type=float,
        default=1.5,
        help="Gaussian FWHM in pixels for Halpha detection only. Default: 1.5.",
    )
    parser.add_argument(
        "--halpha-min-area",
        type=int,
        default=20,
        help="Minimum connected Halpha detection area in pixels. Default: 20.",
    )
    parser.add_argument(
        "--no-halpha-detection",
        action="store_true",
        help="Do not threshold/filter the Halpha overlay.",
    )
    parser.add_argument(
        "--no-hi-contours",
        action="store_true",
        help="Do not draw HI contours.",
    )
    parser.add_argument(
        "--hi-file",
        default=None,
        help="Optional HI moment-zero FITS file. If omitted, use the built-in VFID lookup.",
    )
    parser.add_argument(
        "--hi-levels",
        default=None,
        help="Comma-separated HI contour levels. Default: 3**arange(1,6)+1.",
    )
    parser.add_argument("--hi-color", default="lightskyblue", help="HI contour color.")
    parser.add_argument(
        "--hi-alpha",
        type=float,
        default=0.75,
        help="HI contour alpha. Default: 0.75",
    )
    parser.add_argument(
        "--hi-linewidth",
        type=float,
        default=1.0,
        help="HI contour linewidth. Default: 1.0",
    )
    parser.add_argument(
        "--hi-image",
        action="store_true",
        help="Draw the HI moment-zero map as a filled semi-transparent image layer.",
    )
    parser.add_argument(
        "--hi-image-alpha",
        type=float,
        default=0.45,
        help="Maximum opacity of the filled HI image layer. Default: 0.45.",
    )
    parser.add_argument(
        "--hi-image-percentiles",
        default="50,99.5",
        help=(
            "Lower,upper percentiles for positive HI image scaling. "
            "Default: 50,99.5."
        ),
    )
    parser.add_argument(
        "--hi-image-stretch",
        default="asinh",
        choices=["asinh", "linear"],
        help="Stretch for the filled HI image layer. Default: asinh.",
    )
    parser.add_argument(
        "--hi-image-asinh-a",
        type=float,
        default=0.08,
        help=(
            "Asinh a parameter for the filled HI image layer. Smaller values "
            "bring out fainter HI. Default: 0.08."
        ),
    )
    parser.add_argument(
        "--hi-image-alpha-power",
        type=float,
        default=0.7,
        help=(
            "Power applied to the scaled HI image before alpha mapping. Values "
            "<1 make faint HI more visible. Default: 0.7."
        ),
    )
    parser.add_argument(
        "--hi-image-cmap",
        default="hi-cyan",
        help=(
            "Matplotlib colormap for the filled HI image layer, or hi-cyan for "
            "the built-in blue/cyan map. Default: hi-cyan."
        ),
    )
    parser.add_argument(
        "--hi-image-interpolation",
        default="bilinear",
        help="Matplotlib interpolation for the filled HI image layer. Default: bilinear.",
    )
    parser.add_argument(
        "--hi-image-high-fade-start",
        type=float,
        default=None,
        help=(
            "Scaled HI value where the filled HI layer begins fading out. "
            "Use with --hi-image-high-fade-end."
        ),
    )
    parser.add_argument(
        "--hi-image-high-fade-end",
        type=float,
        default=None,
        help=(
            "Scaled HI value where the filled HI layer is fully suppressed. "
            "Use with --hi-image-high-fade-start."
        ),
    )
    parser.add_argument(
        "--jwst-footprint",
        action="store_true",
        help="Overlay a JWST/NIRCam footprint using pysiaf.",
    )
    parser.add_argument(
        "--jwst-dither-footprints",
        action="store_true",
        help=(
            "Overlay all unique NIRCam dither/tile footprints from an APT "
            "pointing report."
        ),
    )
    parser.add_argument(
        "--jwst-pointing-file",
        default=None,
        help="APT .pointing report used by --jwst-dither-footprints.",
    )
    parser.add_argument(
        "--jwst-pointing-target",
        default=None,
        help=(
            "Target identifier in the APT pointing report. Default: inferred "
            "VFID from galaxy_dir."
        ),
    )
    parser.add_argument(
        "--jwst-pointing-types",
        default=None,
        help=(
            "Optional comma-separated APT pointing row types to include, e.g. "
            "TARGET,DITHER,TILE. Default: include all types, then de-duplicate."
        ),
    )
    parser.add_argument(
        "--jwst-ra",
        type=float,
        default=None,
        help=(
            "JWST footprint center RA in degrees. If omitted, use the crop center "
            "or image center."
        ),
    )
    parser.add_argument(
        "--jwst-dec",
        type=float,
        default=None,
        help=(
            "JWST footprint center Dec in degrees. If omitted, use the crop center "
            "or image center."
        ),
    )
    parser.add_argument(
        "--jwst-pa-v3",
        type=float,
        default=0.0,
        help="JWST PA_V3 angle in degrees east of north. Default: 0.",
    )
    parser.add_argument(
        "--jwst-apertures",
        default="nircam-full",
        help=(
            "Comma-separated NIRCam SIAF apertures, or nircam-full for all ten "
            "detectors. Default: nircam-full."
        ),
    )
    parser.add_argument(
        "--jwst-reference-aperture",
        default="NRCALL_FULL",
        help="SIAF aperture used as the pointing reference. Default: NRCALL_FULL.",
    )
    parser.add_argument(
        "--jwst-color",
        default="gold",
        help="JWST footprint outline color. Default: gold.",
    )
    parser.add_argument(
        "--jwst-alpha",
        type=float,
        default=0.95,
        help="JWST footprint outline alpha. Default: 0.95.",
    )
    parser.add_argument(
        "--jwst-linewidth",
        type=float,
        default=2.0,
        help="JWST footprint outline linewidth. Default: 2.",
    )
    parser.add_argument(
        "--jwst-stroke-linewidth",
        type=float,
        default=None,
        help=(
            "Black stroke width behind JWST footprint outlines. Default: automatic; "
            "use 0 to disable."
        ),
    )
    parser.add_argument(
        "--jwst-fill-alpha",
        type=float,
        default=0.0,
        help="Optional JWST footprint fill alpha. Default: 0.",
    )
    parser.add_argument(
        "--jwst-outer-outline",
        action="store_true",
        help="Draw a convex outer outline around each NIRCam footprint/dither.",
    )
    parser.add_argument(
        "--jwst-outer-color",
        default=None,
        help="Color for the JWST outer footprint outline. Default: use --jwst-color.",
    )
    parser.add_argument(
        "--jwst-outer-alpha",
        type=float,
        default=0.95,
        help="Alpha for the JWST outer footprint outline. Default: 0.95.",
    )
    parser.add_argument(
        "--jwst-outer-linewidth",
        type=float,
        default=2.0,
        help="Linewidth for the JWST outer footprint outline. Default: 2.",
    )
    parser.add_argument(
        "--jwst-outer-stroke-linewidth",
        type=float,
        default=0.0,
        help="Black stroke width behind the JWST outer outline. Default: 0.",
    )
    parser.add_argument(
        "--no-scale",
        action="store_true",
        help="Do not draw a physical scale bar.",
    )
    parser.add_argument(
        "--scale-kpc",
        type=float,
        default=5.0,
        help="Scale-bar size in kpc. Default: 5.",
    )
    parser.add_argument("--scale-color", default="white", help="Scale-bar color.")
    parser.add_argument(
        "--scale-fontsize",
        type=float,
        default=12,
        help="Scale-bar label font size. Default: 12.",
    )
    parser.add_argument(
        "--scale-x",
        type=float,
        default=0.08,
        help="Scale-bar x position as an axis fraction. Default: 0.08.",
    )
    parser.add_argument(
        "--scale-y",
        type=float,
        default=0.08,
        help="Scale-bar y position as an axis fraction. Default: 0.08.",
    )
    parser.add_argument(
        "--scale-dytext",
        type=float,
        default=0.06,
        help="Scale-bar label offset as an axis fraction. Default: 0.06.",
    )
    parser.add_argument(
        "--scale-linewidth",
        type=float,
        default=2.0,
        help="Scale-bar line width. Default: 2.",
    )
    parser.add_argument("--title", default=None, help="Figure title. Default: inferred VFID.")
    parser.add_argument(
        "--no-title",
        action="store_true",
        help="Do not draw the external figure title.",
    )
    parser.add_argument(
        "--hide-axes",
        action="store_true",
        help="Hide WCS axes, ticks, labels, and frame.",
    )
    parser.add_argument(
        "--corner-label",
        nargs="?",
        const="auto",
        default=None,
        help=(
            "Draw an in-image corner label. If no value is provided, infer the "
            "galaxy name from the directory."
        ),
    )
    parser.add_argument(
        "--corner-label-position",
        default="upper left",
        choices=["upper left", "upper right", "lower left", "lower right"],
        help="Corner label position. Default: upper left.",
    )
    parser.add_argument(
        "--corner-label-color",
        default="white",
        help="Corner label text color. Default: white.",
    )
    parser.add_argument(
        "--corner-label-fontsize",
        type=float,
        default=18,
        help="Corner label font size. Default: 18.",
    )
    parser.add_argument(
        "--corner-label-pad",
        type=float,
        default=0.035,
        help="Corner label padding as an axis fraction. Default: 0.035.",
    )
    parser.add_argument(
        "--figsize",
        default="8,8",
        help="Figure size in inches as width,height. Default: 8,8",
    )
    parser.add_argument("--dpi", type=int, default=200, help="Output DPI. Default: 200.")
    parser.add_argument("--verbose", action="store_true", help="Print file details.")

    return parser.parse_args()


def main():
    args = parse_args()

    optical_percentiles = parse_csv_floats(args.optical_percentiles)
    halpha_percentiles = parse_csv_floats(args.halpha_percentiles)
    figsize = parse_csv_floats(args.figsize)
    hi_levels = parse_csv_floats(args.hi_levels)
    hi_image_percentiles = parse_csv_floats(args.hi_image_percentiles)
    jwst_apertures = parse_csv_strings(args.jwst_apertures)
    jwst_pointing_types = parse_csv_strings(args.jwst_pointing_types)

    if len(optical_percentiles) != 2:
        raise SystemExit("--optical-percentiles must have two values")
    if len(halpha_percentiles) != 2:
        raise SystemExit("--halpha-percentiles must have two values")
    if len(figsize) != 2:
        raise SystemExit("--figsize must have two values")
    if len(hi_image_percentiles) != 2:
        raise SystemExit("--hi-image-percentiles must have two values")

    make_legacy_halpha_rgb(
        galaxy_dir=args.galaxy_dir,
        halpha_glob=args.halpha_glob,
        legacy_dir=args.legacy_dir,
        g_glob=args.g_glob,
        r_glob=args.r_glob,
        z_glob=args.z_glob,
        mask_glob=None if args.no_halpha_mask else args.mask_glob,
        mask_file=None if args.no_halpha_mask else args.mask_file,
        mask_mode=args.mask_mode,
        outfile=args.outfile,
        output_dir=args.output_dir,
        recipe=args.recipe,
        crop=args.crop,
        optical_stretch=args.optical_stretch,
        optical_percentiles=tuple(optical_percentiles),
        legacy_rgb_m=args.legacy_rgb_m,
        legacy_rgb_q=args.legacy_rgb_q,
        legacy_rgb_stretch_factor=args.legacy_rgb_stretch_factor,
        halpha_percentiles=tuple(halpha_percentiles),
        optical_asinh_a=args.optical_asinh_a,
        halpha_asinh_a=args.halpha_asinh_a,
        halpha_alpha=args.halpha_alpha,
        halpha_color=args.halpha_color,
        halpha_nsigma=args.halpha_nsigma,
        halpha_smooth_fwhm=args.halpha_smooth_fwhm,
        halpha_min_area=args.halpha_min_area,
        halpha_detection=not args.no_halpha_detection,
        hi_contours=not args.no_hi_contours,
        hi_file=args.hi_file,
        hi_levels=hi_levels,
        hi_color=args.hi_color,
        hi_alpha=args.hi_alpha,
        hi_linewidth=args.hi_linewidth,
        hi_image=args.hi_image,
        hi_image_alpha=args.hi_image_alpha,
        hi_image_percentiles=tuple(hi_image_percentiles),
        hi_image_stretch=args.hi_image_stretch,
        hi_image_asinh_a=args.hi_image_asinh_a,
        hi_image_alpha_power=args.hi_image_alpha_power,
        hi_image_cmap=args.hi_image_cmap,
        hi_image_interpolation=args.hi_image_interpolation,
        hi_image_high_fade_start=args.hi_image_high_fade_start,
        hi_image_high_fade_end=args.hi_image_high_fade_end,
        jwst_footprint=args.jwst_footprint,
        jwst_ra=args.jwst_ra,
        jwst_dec=args.jwst_dec,
        jwst_pa_v3=args.jwst_pa_v3,
        jwst_apertures=jwst_apertures,
        jwst_reference_aperture=args.jwst_reference_aperture,
        jwst_color=args.jwst_color,
        jwst_alpha=args.jwst_alpha,
        jwst_linewidth=args.jwst_linewidth,
        jwst_stroke_linewidth=args.jwst_stroke_linewidth,
        jwst_fill_alpha=args.jwst_fill_alpha,
        jwst_outer_outline=args.jwst_outer_outline,
        jwst_outer_color=args.jwst_outer_color,
        jwst_outer_alpha=args.jwst_outer_alpha,
        jwst_outer_linewidth=args.jwst_outer_linewidth,
        jwst_outer_stroke_linewidth=args.jwst_outer_stroke_linewidth,
        jwst_dither_footprints=args.jwst_dither_footprints,
        jwst_pointing_file=args.jwst_pointing_file,
        jwst_pointing_target=args.jwst_pointing_target,
        jwst_pointing_types=jwst_pointing_types,
        add_scale=not args.no_scale,
        scale_kpc=args.scale_kpc,
        scale_color=args.scale_color,
        scale_fontsize=args.scale_fontsize,
        scale_x=args.scale_x,
        scale_y=args.scale_y,
        scale_dytext=args.scale_dytext,
        scale_linewidth=args.scale_linewidth,
        title=args.title,
        show_title=not args.no_title,
        hide_axes=args.hide_axes,
        corner_label=args.corner_label,
        corner_label_position=args.corner_label_position,
        corner_label_color=args.corner_label_color,
        corner_label_fontsize=args.corner_label_fontsize,
        corner_label_pad=args.corner_label_pad,
        figsize=tuple(figsize),
        dpi=args.dpi,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
