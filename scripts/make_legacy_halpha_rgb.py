#!/usr/bin/env python

"""
Make a false-color RGB image from aligned Legacy Survey and Halpha FITS images.

Default channel recipe:
    R = Legacy z
    G = Legacy r
    B = Legacy g

The Halpha/SFR image is scaled separately and overlaid in color.  If an HI
moment-zero map is available, HI contours are drawn using WCSAxes transforms.
"""

import argparse
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.cosmology import WMAP9 as cosmo
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import AsinhStretch, LinearStretch
from astropy.wcs import WCS
from astropy.wcs import utils as wcs_utils
from scipy import ndimage


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


def parse_csv_floats(value):
    """
    Parse a comma-separated list of floats.
    """
    if value is None:
        return None
    return [float(x) for x in value.split(",") if x.strip()]


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

    ax.plot([xline1, xline2], [yline, yline], "-", lw=linewidth, color=color)
    ax.plot([xline1, xline1], [ycap1, ycap2], "-", lw=linewidth, color=color)
    ax.plot([xline2, xline2], [ycap1, ycap2], "-", lw=linewidth, color=color)
    ax.text(
        xtext,
        ytext,
        f"{barsize:.0f} kpc",
        horizontalalignment="center",
        fontsize=fontsize,
        color=color,
    )


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
    add_scale=True,
    scale_kpc=5.0,
    scale_color="white",
    scale_fontsize=12,
    scale_x=0.08,
    scale_y=0.08,
    scale_dytext=0.06,
    title=None,
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
            )
        except Exception as err:
            if verbose:
                print(f"WARNING: could not add scale bar: {err}")

    if title is None:
        title = vfid

    if title:
        ax.set_title(title)

    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_axislabel("RA")
    lat.set_axislabel("Dec")
    lon.set_major_formatter("hh:mm:ss")
    lat.set_major_formatter("dd:mm")

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
    parser.add_argument("--title", default=None, help="Figure title. Default: inferred VFID.")
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

    if len(optical_percentiles) != 2:
        raise SystemExit("--optical-percentiles must have two values")
    if len(halpha_percentiles) != 2:
        raise SystemExit("--halpha-percentiles must have two values")
    if len(figsize) != 2:
        raise SystemExit("--figsize must have two values")

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
        add_scale=not args.no_scale,
        scale_kpc=args.scale_kpc,
        scale_color=args.scale_color,
        scale_fontsize=args.scale_fontsize,
        scale_x=args.scale_x,
        scale_y=args.scale_y,
        scale_dytext=args.scale_dytext,
        title=args.title,
        figsize=tuple(figsize),
        dpi=args.dpi,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
