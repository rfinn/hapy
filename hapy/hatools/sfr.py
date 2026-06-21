from pathlib import Path
import json
import numpy as np

from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo


KE12_LOG_SFR_FACTOR = 41.27
AB_ZERO_CGS = 48.60

# Cam morgan + 2024 uses Calzetti+2010 with modification to chabrier IMF
# SFR [M yr−1] = 5.01 × 10−42 L(Hα) [erg s−1].
# this gives 41.3 vs 41.27 from Kennicutt & Evans, 2012

def distance_from_velocity(v_kms):
    return cosmo.luminosity_distance(float(v_kms) / 3.0e5)


def find_sfr_input_image(cutout_dir, prefer="csgr"):
    cutout_dir = Path(cutout_dir)
    tag = cutout_dir.name

    choices = []

    if prefer.lower() == "csgr":
        choices = [
            cutout_dir / f"{tag}-CS-gr.fits",
            cutout_dir / f"{tag}-CS-ZP.fits",
        ]
    else:
        choices = [
            cutout_dir / f"{tag}-CS-ZP.fits",
            cutout_dir / f"{tag}-CS-gr.fits",
        ]

    for f in choices:
        if f.exists():
            return f

    raise FileNotFoundError(f"No CS-gr or CS-ZP image found in {cutout_dir}")


def get_header_float(header, keys, default=np.nan):
    for key in keys:
        if key in header:
            try:
                return float(header[key])
            except Exception:
                pass
    return default


def make_sfr_map(
    cutout_dir,
    velocity_kms=None,
    distance_mpc=None,
    prefer="csgr",
    output_suffix="",
    overwrite=True,
):
    """
    Make linear SFR map from HAPY continuum-subtracted Halpha image.

    Output units are Msun/yr per pixel.
    """
    cutout_dir = Path(cutout_dir)
    tag = cutout_dir.name

    metadata_file = cutout_dir / "metadata.json"
    params = json.loads(metadata_file.read_text()) if metadata_file.exists() else {}

    if distance_mpc is None:
        if velocity_kms is None:
            velocity_kms = params.get("vcosmic", params.get("vr", None))
        if velocity_kms is None:
            raise ValueError("Need velocity_kms, distance_mpc, or vr/vcosmic in metadata.json")
        dist = distance_from_velocity(velocity_kms)
        dist_cm = dist.cgs.value
        dist_mpc_value = dist.to("Mpc").value
    else:
        dist_cm = float(distance_mpc) * 3.0856775814913673e24
        dist_mpc_value = float(distance_mpc)

    ha_file = find_sfr_input_image(cutout_dir, prefer=prefer)

    with fits.open(ha_file) as hdul:
        data = hdul[0].data.astype(float)
        hdr = hdul[0].header.copy()

    hZP = get_header_float(hdr, ["PHOTZP", "MAGZP"])
    if not np.isfinite(hZP):
        raise ValueError(f"Missing/non-finite PHOTZP in {ha_file}")

    filter_width = get_header_float(
        hdr,
        ["HFILTERW", "HFILT_W", "FILTERW", "FILTER_WIDTH", "FILTWDTH"],
        default=params.get("hafilter_width_A", np.nan),
    )

    filter_center = get_header_float(
        hdr,
        ["HFILTERC", "HFILT_C", "FILTERC", "FILTER_CENTER", "FILT_CEN"],
        default=params.get("hafilter_center_A", np.nan),
    )

    if not np.isfinite(filter_width):
        filter_width = float(params.get("hafilter_width_A", np.nan))
    if not np.isfinite(filter_center):
        filter_center = float(params.get("hafilter_center_A", np.nan))

    if not np.isfinite(filter_width) or not np.isfinite(filter_center):
        raise ValueError(
            f"Missing filter width/center for {tag}; "
            "need hafilter_width_A and hafilter_center_A in metadata or header."
        )

    # AB mag/count conversion to fnu, then convert to integrated line flux
    # using dnu = c dlambda / lambda^2, then luminosity, then SFR.
    fnu_scale = 10.0 ** (-0.4 * (hZP + AB_ZERO_CGS))
    dnu = 3.0e18 * filter_width / filter_center**2
    lum_scale = 4.0 * np.pi * dist_cm**2
    sfr_scale = fnu_scale * dnu * lum_scale / (10.0 ** KE12_LOG_SFR_FACTOR)

    sfr = data * sfr_scale

    hdr["SFRMAP"] = (True, "SFR image")
    hdr["SFRUNIT"] = ("Msun/yr/pix", "Image units")
    hdr["SFRSRC"] = (ha_file.name, "Input continuum-subtracted image")
    hdr["SFRPREF"] = (prefer, "Preferred input image type")
    hdr["SFRZP"] = (float(hZP), "Input image PHOTZP")
    hdr["SFRFWID"] = (float(filter_width), "Halpha filter width, Angstrom")
    hdr["SFRFCEN"] = (float(filter_center), "Halpha filter center, Angstrom")
    hdr["SFRSCAL"] = (float(sfr_scale), "ADU to Msun/yr scale")
    hdr["SFRCAL"] = ("KennicuttEvans2012", "Halpha SFR calibration")
    hdr["DISTMPC"] = (float(dist_mpc_value), "Distance used for SFR map, Mpc")

    if velocity_kms is not None:
        hdr["VELDIST"] = (float(velocity_kms), "Velocity used for distance, km/s")

    suffix = f"-{output_suffix}" if output_suffix else ""
    outname = cutout_dir / f"{tag}-sfr{suffix}.fits"

    fits.PrimaryHDU(data=sfr, header=hdr).writeto(outname, overwrite=overwrite)

    return outname
