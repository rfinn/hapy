# --- Filter central wavelength and width --- #
#
# We calculate these in github/filter_transformations/filtertrans-dev.ipynb
# Becky cross checked with plots that she has


# this dictionary adds the integral of the transmission, which we need for halpha filter correction
# leaving the previous dictionar for backwards compatibility
# contains:
#   filter_center_A, filter_width_A, Rlambda_A
filter_center_width_sumtransmission = {
    '90prime-BASSr.fits':(6410.8, 1398.8, 104218.1),
    '90prime-Ha+4nm.fits':(6620.8, 83.3, 7651.4),
    'BOK90prime-BASSr.fits':(6410.8, 1398.8, 104218.1),
    'BOK90prime-Ha+4nm.fits':(6620.8, 83.3, 7651.4),
    'MOS-Ha+12nm.fits':(6698.8, 86.1, 7441.1),
    'MOS-Ha+16nm.fits':(6730.8, 83.9, 7640.8),
    'MOS-Ha+4nm.fits':(6620.8, 83.3, 7651.4),
    'MOS-Ha+8nm.fits':(6654.4, 84.1, 7614.4),
    'MOS-HarrisR.fits':(6653.9, 1551.1, 134154.0),
    'MOS-SDSSr.fits':(6287.6, 1382.5, 128347.2),
    'HDI-Ha+12nm.fits':(6701.7, 61.5, 5666.6),
    'HDI-Ha+16nm.fits':(6742.1, 59.3, 5269.4),
    'HDI-Ha+4nm.fits':(6618.5, 60.4, 5337.0),
    'HDI-Ha+8nm.fits':(6660.0, 59.8, 5564.6),
    'HDI-Ha.fits':(6580.0, 58.8, 5245.1),
    'HDI-HarrisR.fits':(6605.3, 1576.1, 133317.6),
    'HDI-SDSSr.fits':(6242.3, 1425.5, 134039.6),
    'WFC-Ha-197.fits':(6568.0, 92.9, 8421.1),
    'WFC-Ha-227.fits':(6666.2, 80.9, 7196.3),
    'WFC-SDSSr-214.fits':(6230.1, 1219.8, 109707.3),
    'panstarrs-g.fits':(4866.5, 1166.4, 58551.0),
    'panstarrs-r.fits':(6214.6, 1318.0, 90544.0),
    }

from pathlib import Path


def _norm_filter_name(fname):
    """
    Normalize filter filename keys.
    """
    fname = str(fname).strip()

    if not fname.endswith(".fits"):
        fname = fname + ".fits"

    return Path(fname).name


def get_filter_center_width_Rlambda(filtername):
    """
    Return filter center, width, and integrated transmission Rlambda.
    """
    filtername = _norm_filter_name(filtername)

    try:
        return filter_center_width_sumtransmission[filtername]
    except KeyError as exc:
        raise KeyError(
            f"Filter {filtername} not found in filter_center_width_sumtransmission"
        ) from exc


def get_continuum_oversubtraction(rfilter_name, hafilter_name):
    """
    Continuum oversubtraction correction.

    Following the Gavazzi-style correction:

        1 + Rlambda_Ha / Rlambda_r

    where Rlambda is the integral of the filter transmission.
    """

    _, _, r_Rlambda = get_filter_center_width_Rlambda(rfilter_name)
    _, _, ha_Rlambda = get_filter_center_width_Rlambda(hafilter_name)

    return 1.0 + ha_Rlambda / r_Rlambda


def get_continuum_oversubtraction_from_metadata(meta):
    """
    Compute continuum oversubtraction correction from metadata.json.

    Required metadata keys:
        rfilter_name
        hafilter_name
    """

    rfilter_name = meta["rfilter_name"]
    hafilter_name = meta["hafilter_name"]

    return get_continuum_oversubtraction(rfilter_name, hafilter_name)
