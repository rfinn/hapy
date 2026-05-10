"""
Filter color transformations for CS-gr continuum subtraction.

The fits return:

    Halpha_mag - R_mag = f(g-r)

where f is a polynomial in g-r.

Coefficients are stored in np.polyval order:

    [a2, a1, a0]

so:

    delta_mag = a2*(g-r)^2 + a1*(g-r) + a0

Use in CS-gr as:

    mag_r_to_Ha = mag_r + delta_mag
    delta_flux = 10**(-0.4 * delta_mag)
"""

from pathlib import Path

import numpy as np

from hapy.hatools.filter_transmission import (
    get_halpha_filtername,
    get_rband_filtername,
)


def _filter_stem(filtername):
    """
    Convert canonical filter filename to lookup key.
    """
    stem = Path(str(filtername)).stem

    # Remove plus signs used in some Halpha names
    stem = stem.replace("+", "")

    # BOK Halpha filter files may be named 90prime-Ha4nm,
    # while coefficient table uses BOK90prime-Ha4nm.
    if stem.startswith("90prime-Ha"):
        return "BOK" + stem

    # Historical key for INT/WFC r-band transformation table
    if stem == "WFC-SDSSr-214":
        return "WFC-SDSSr214"

    # Historical keys for INT/WFC Halpha transformation table
    if stem == "WFC-Ha-197":
        return "WFC-Ha197"

    if stem == "WFC-Ha-227":
        return "WFC-Ha227"

    return stem


# def _filter_stem(filtername):
#     """
#     Convert canonical filter filename to lookup key.

#     Examples
#     --------
#     BOK90prime-BASSr.fits -> BOK90prime-BASSr
#     WFC-SDSSr-214.fits   -> WFC-SDSSr214
#     90prime-Ha+4nm.fits  -> BOK90prime-Ha4nm
#     """
#     stem = Path(str(filtername)).stem

#     # Remove plus signs used in metadata/filter filenames
#     stem = stem.replace("+", "")

#     # Historical key for INT/WFC transformation table
#     if stem == "WFC-SDSSr-214":
#         return "WFC-SDSSr214"

#     # BOK Halpha filter files may be named 90prime-Ha4nm,
#     # while coefficient table uses BOK90prime-Ha4nm.
#     if stem.startswith("90prime-Ha"):
#         return "BOK" + stem

#     return stem


HALPHA_MINUS_R_COEFFS = {
    "BOK90prime-BASSr": {
        "BOK90prime-Ha4nm": [0.0287, -0.1685, 0.0329],
    },

    "MOS-SDSSr": {
        "MOS-Ha4nm":  [0.0353, -0.2364, 0.0430],
        "MOS-Ha8nm":  [0.0160, -0.2139, 0.0236],
        "MOS-Ha12nm": [0.0154, -0.2360, 0.0435],
        "MOS-Ha16nm": [0.0173, -0.2615, 0.0683],
    },

    "MOS-HarrisR": {
        "MOS-Ha4nm":  [0.0183, -0.1035, 0.0407],
        "MOS-Ha8nm":  [-0.0010, -0.0925, 0.0416],
        "MOS-Ha12nm": [-0.0021, -0.1056, 0.0478],
        "MOS-Ha16nm": [-0.0005, -0.1218, 0.0576],
    },

    "HDI-SDSSr": {
        "HDI-Ha4nm":  [0.0274, -0.2470, 0.0407],
        "HDI-Ha8nm":  [0.0161, -0.2355, 0.0223],
        "HDI-Ha12nm": [0.0168, -0.2664, 0.0515],
        "HDI-Ha16nm": [0.0197, -0.2979, 0.0831],
    },

    "HDI-HarrisR": {
        "HDI-Ha4nm":  [0.0107, -0.1039, 0.0308],
        "HDI-Ha8nm":  [0.0003, -0.1108, 0.0414],
        "HDI-Ha12nm": [0.0004, -0.1306, 0.0539],
        "HDI-Ha16nm": [0.0030, -0.1519, 0.0688],
    },

    "WFC-SDSSr214": {
        "WFC-Ha197": [0.1127, -0.3697, 0.1077],
        "WFC-Ha227": [0.0171, -0.2349, 0.0184],
    },
}


def get_halpha_minus_r_coeffs(instrument, rfilter, hfilter):
    """
    Return polynomial coefficients for Halpha_mag - R_mag as a function of g-r.

    Parameters
    ----------
    instrument : str
        Instrument name, e.g. BOK, MOS, HDI, INT.
    rfilter : str
        r-band FITS FILTER keyword or canonical r-band filter filename.
    hfilter : str
        Halpha filter name from metadata, e.g. Ha+4nm, Ha4nm, Ha197.

    Returns
    -------
    coeffs : list[float]
        Polynomial coefficients in np.polyval order.
    """

    # Allow already-canonical filenames from metadata.json
    if str(rfilter).endswith(".fits"):
        rband_filtername = str(rfilter)
    else:
        rband_filtername = get_rband_filtername(instrument, rfilter)

    if str(hfilter).endswith(".fits"):
        halpha_filtername = str(hfilter)
    else:
        halpha_filtername = get_halpha_filtername(instrument, hfilter)

    rkey = _filter_stem(rband_filtername)
    hkey = _filter_stem(halpha_filtername)

    try:
        return HALPHA_MINUS_R_COEFFS[rkey][hkey]
    except KeyError as exc:
        available = {
            rk: list(hdict.keys())
            for rk, hdict in HALPHA_MINUS_R_COEFFS.items()
        }

        raise KeyError(
            "No Halpha-R color transformation found for:\n"
            f"  instrument = {instrument}\n"
            f"  rfilter    = {rfilter}\n"
            f"  hfilter    = {hfilter}\n"
            f"  canonical r filter  = {rband_filtername} -> {rkey}\n"
            f"  canonical Ha filter = {halpha_filtername} -> {hkey}\n"
            f"Available transformations:\n{available}"
        ) from exc


def halpha_minus_r_color(instrument, rfilter, hfilter, gr):
    """
    Evaluate Halpha_mag - R_mag as a function of g-r color.

    Parameters
    ----------
    instrument : str
        Instrument name, e.g. BOK, MOS, HDI, INT.
    rfilter : str
        r-band FITS FILTER keyword or canonical r-band filter filename.
    hfilter : str
        Halpha filter name from metadata or canonical Halpha filter filename.
    gr : array-like
        g-r color image or scalar.

    Returns
    -------
    delta_mag : array-like
        Halpha_mag - R_mag correction in magnitudes.
    """

    coeffs = get_halpha_minus_r_coeffs(instrument, rfilter, hfilter)
    return np.polyval(coeffs, gr)


def halpha_minus_r_color_from_metadata(meta, gr):
    """
    Convenience wrapper using metadata.json contents.

    Expected metadata keys:
        telescope
        rfilter_name
        hafilter or hafilter_name
    """

    instrument = meta["telescope"]
    rfilter = meta.get("rfilter_name")

    if rfilter is None:
        raise KeyError("metadata is missing required key 'rfilter_name'")

    hfilter = meta.get("hafilter_name", meta.get("hafilter"))

    if hfilter is None:
        raise KeyError("metadata is missing required key 'hafilter' or 'hafilter_name'")

    return halpha_minus_r_color(instrument, rfilter, hfilter, gr)
