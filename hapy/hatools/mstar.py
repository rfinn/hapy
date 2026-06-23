from pathlib import Path
import json
import numpy as np

from astropy.io import fits
from astropy import convolution
from astropy.cosmology import WMAP9 as cosmo
from astropy.wcs import WCS

from hapy.hatools.utils import get_kpc_per_pixel

# Roediger & Courteau 2015, r-band, Marigo TP-AGB
AR_RC15 = -0.647
BR_RC15 = 1.497
LEGACY_R_ZP = 22.5
MSUN_R_AB = 4.61


def find_one(pattern, root):
    matches = sorted(Path(root).glob(pattern))
    if len(matches) == 0:
        raise FileNotFoundError(f"No file found matching {pattern} in {root}")
    return matches[0]


def distance_from_velocity(v_kms):
    return cosmo.luminosity_distance(float(v_kms) / 3.0e5)


def make_mstar_map(
    cutout_dir,
    velocity_kms=None,
    distance_mpc=None,
    smoothing=15,
    output_suffix="vcosmic",
    overwrite=True,
):
    """
    Make stellar-mass map from Legacy r-band luminosity and Legacy g-r color.

    Output image is linear stellar mass in units of 1e7 Msun per pixel.
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
    else:
        dist = distance_mpc

    legacy_dir = cutout_dir / "legacy"

    #legacy_g = find_one("*-legacy-*-g.fits", legacy_dir)

    # use the reprojected legacy images b/c we want to compare
    # with SFR map, which will be at the native telescope scale
    legacy_r = find_one("*-legacy-*-r-ha.fits", legacy_dir)



    rhdu = fits.open(legacy_r)
    rdata = rhdu[0].data.astype(float)


    rsmooth = convolution.convolve_fft(
        rdata,
        convolution.Box2DKernel(smoothing),
        allow_huge=True,
        nan_treatment="interpolate",
    )


    #mag_g = 22.5 - 2.5*np.log10(gdata)
    #mag_r = 22.5 - 2.5*np.log10(rdata)
    #mag_gr = mag_g - mag_r

    legacy_gr = find_one("*gr-ha-smooth.fits", legacy_dir)
    grhdu = fits.open(legacy_gr)
    mag_gr = grhdu[0].data.astype(float)

    good = np.isfinite(rsmooth) & (rsmooth > 0) & np.isfinite(mag_gr)

    Mr = np.full_like(rsmooth, np.nan, dtype=float)
    Mr[good] = (
        LEGACY_R_ZP
        - 2.5 * np.log10(rsmooth[good])
        - 5.0 * np.log10(dist.to("pc").value)
        + 5.0
    )

    logmstar = np.full_like(rsmooth, np.nan, dtype=float)
    logmstar[good] = AR_RC15 + BR_RC15 * mag_gr[good] + (Mr[good] - MSUN_R_AB) / (-2.5)

    # linear stellar mass image in units of Msun/pix
    mstar = 10.0 ** (logmstar)

    hdr = rhdu[0].header.copy()
    hdr["MSTAR"] = (True, "Stellar mass image")
    hdr["MSTARUN"] = ("Msun/pix", "Image units")
    hdr["MSTARZP"] = (LEGACY_R_ZP, "Legacy r-band zeropoint")
    hdr["MSUNRAB"] = (MSUN_R_AB, "Solar absolute magnitude in r, AB")
    hdr["MLREL"] = ("RoedigerCourteau2015", "Color-M/L relation")
    hdr["ML_AR"] = (AR_RC15, "log M/L intercept")
    hdr["ML_BR"] = (BR_RC15, "log M/L color coefficient")
    hdr["RIMAGE"] = (legacy_r.name, "Legacy r image used for luminosity")
    hdr["GRIMAGE"] = (legacy_gr.name, "Legacy g-r image used for color")
    hdr["SMOOTH"] = (int(smoothing), "Boxcar smoothing size in pixels")

    hdr["LUMSRC"] = ("legacy-r-ha", "Luminosity image used")
    hdr["COLSRC"] = ("legacy-gr-ha", "Color image used")


    dist_mpc = dist.to("Mpc").value if hasattr(dist, "to") else float(dist)
    mstar_hdr["DISTMPC"] = (float(dist_mpc), "Distance used, Mpc")


    if velocity_kms is not None:
        hdr["VELDIST"] = (float(velocity_kms), "Velocity used for distance, km/s")

    outname = cutout_dir / f"{tag}-mstar-{output_suffix}.fits"
    fits.PrimaryHDU(data=mstar, header=hdr).writeto(outname, overwrite=overwrite)


    # stellar mass surface density, still in units of 1e7 Msun/kpc^2
    imwcs = WCS(rhdu[0].header)
    kpc_per_pixel = get_kpc_per_pixel(imwcs, distance_mpc=dist_mpc)
    pixel_area_kpc2 = kpc_per_pixel**2

    sigma_mstar = mstar / pixel_area_kpc2

    sigma_hdr = mstar_hdr.copy()
    sigma_hdr["SIGMSTAR"] = (True, "Stellar mass surface density image")
    sigma_hdr["BUNIT"] = ("Msun/pc2", "Image units")
    sigma_hdr["KPCPIX"] = (float(kpc_per_pixel), "kpc per pixel")
    sigma_hdr["PIXARKPC"] = (float(pixel_area_kpc2), "Pixel area, kpc2")

    outname = cutout_dir / f"{tag}-sigma-mstar-{output_suffix}.fits"
    fits.PrimaryHDU(data=sigma_mstar, header=sigma_hdr).writeto(outname, overwrite=overwrite)


    rhdu.close()
    grhdu.close()



    

    return outname
