import numpy as np
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

import os
import warnings
from astropy.table import Table

import numpy as np
from hapy.imagetools.imutils import get_image_size_deg, get_image_center_deg, get_image_footprint_box_deg
from hapy.masktools.maskops import circle_pixels

def mask_radius_for_mag(mag):
    """ 
    function is from legacy pipeline  
    https://github.com/legacysurvey/legacypipe/blob/6d1a92f8462f4db9360fb1a68ef7d6c252781027/py/legacypipe/reference.py#L314-L319
    """
    # Returns a masking radius in degrees for a star of the given magnitude.
    # Used for Tycho-2 and Gaia stars.

    # This is in degrees, and is from Rongpu in the thread [decam-chatter 12099].
    return 1630./3600. * 1.396**(-mag)



def gaia_stars_in_rectangle(ra, dec, height, width, minmag=None, maxmag=None, pmsnrcut=5):
    """
    Get Gaia stars within a circular aperture.
    started by ChatGPT, rewritten by Rose Finn :)

    Parameters:
        ra (float): Right ascension in degrees.
        dec (float): Declination in degrees.
        height (float): rectangular dimension in dec direction
        width (float): rectangular dimension in ra direction
   
    Optional Parameters:
        minmag : min mag of stars to return
        maxmag : max mag of stars to return
        pmsrncut : min snr cut to use in selecting stars

    Returns:
        table: table of Gaia stars within the specified region that have pm snr > 5.


    NOTES:
    
        explanation of columns can be found here:
        https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html

    TODO: check what to use for maxmag
    """
    # Define the target coordinates
    target_coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')

    # Use DR3, which is default...
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

    # return an unlimited number of stars
    Gaia.ROW_LIMIT = -1
    
    # Query Gaia for stars within the specified rectangle
    width = u.Quantity(width, u.deg)
    height = u.Quantity(height, u.deg)

    result = Gaia.query_object_async(target_coord,height=height,width=width)
    
    #query = f"SELECT TOP 10000 * \
    #    FROM gaiadr3.gaia_source \
    #    WHERE CONTAINS(POINT(ra, dec), CIRCLE({ra}, {dec}, {radius})) = 1"

    ## Perform the Gaia query
    #job = Gaia.launch_job(query)
    #result = job.get_results()


    # Extract relevant columns (you can customize this)
    #print(result.colnames)
    selected_columns = ['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error']

    try:
        stars = result[selected_columns]
    except KeyError:
        selected_columns = ['SOURCE_ID', 'ra', 'dec', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'parallax', 'parallax_error', 'ruwe']
        stars = result[selected_columns]

    # Require SNR > 5 in proper motion
    #keepflag = np.sqrt((stars['pmra']/stars['pmra_error'])**2 + (stars['pmdec']/stars['pmdec_error'])**2) > 5
    keepflag = gaia_foreground_filter(stars)
    #keepflag = np.ones(len(stars),'bool')

    # Cut by min/max magnitude, if they are provided
    if maxmag is not None:
        keepflag = keepflag & (stars['phot_g_mean_mag'] < maxmag)
    if minmag is not None:
        keepflag = keepflag & (stars['phot_g_mean_mag'] > minmag)
    
    return stars[keepflag]



def gaia_foreground_filter(tab, pm_snr_min=5.0, plx_snr_min=5.0, ruwe_max=1.4, phot_g_max=20):
    """
    Return boolean mask selecting likely MW foreground stars.
    Keeps sources with significant proper motion OR significant parallax,
    and (optionally) reasonable RUWE.
    """
    pmra = np.array(tab["pmra"])
    pmdec = np.array(tab["pmdec"])
    pmra_err = np.array(tab["pmra_error"])
    pmdec_err = np.array(tab["pmdec_error"])

    pmag = np.array(tab["phot_g_mean_mag"])

    
    # total PM and approximate error
    pm = np.sqrt(pmra**2 + pmdec**2)
    pm_err = np.sqrt(pmra_err**2 + pmdec_err**2)
    pm_snr = pm / np.where(pm_err > 0, pm_err, np.nan)

    motion_ok = (pm_snr >= pm_snr_min) & (pmag < phot_g_max)
    motion_ok &= np.isfinite(pm_snr)
    if "parallax" in tab.colnames:
        plx = np.array(tab["parallax"])
        plx_err = np.array(tab["parallax_error"])
        plx_snr = np.abs(plx) / np.where(plx_err > 0, plx_err, np.nan)

        # allow either significant pm or significant parallax
        motion_ok = (pm_snr >= pm_snr_min) | (plx_snr >= plx_snr_min)
        motion_ok &= np.isfinite(pm_snr) | np.isfinite(plx_snr)

        
    # if "ruwe" in tab.colnames:
    #     ruwe = np.array(tab["ruwe"])
    #     motion_ok &= (ruwe < ruwe_max)

    # drop NaNs safely


    return motion_ok
def get_gaia_stars(image_name, gaiapath=None, use_cache=True,):
    """
    Retrieve Gaia bright stars within the image FOV.  
    Add radius to the returned table using the Legacy Survey magnitude-radius relation.

    Parameters
    ----------
    image_name : str
        Name of image file (used for cache filename).
    gaiapath : str, optional
        Path to local Gaia catalog fallback.
    use_cache : bool
        If True, read/write cached CSV.

    Returns
    -------
    brightstar : astropy Table or None
    x_pixels : array or None
    y_pixels : array or None
    """

    outfile = image_name.replace(".fits", "_gaia_stars.csv")

    # -------------------------------------------------
    # 1 Read cached file if available
    # -------------------------------------------------
    if use_cache and os.path.exists(outfile):
        print(f"Reading Gaia stars from {outfile}")
        brightstar = Table.read(outfile)

        # filter out foreground stars
        keepflag =  gaia_foreground_filter(brightstar)
        brightstar = brightstar[keepflag]
        return brightstar, brightstar["xpixel"], brightstar["ypixel"]


    # -------------------------------------------------
    # Get image properties
    #
    # racenter,deccenter,
    # dxdeg,dydeg,
    # image_wcs,
    # xmax, ymax,
    # -------------------------------------------------    
    image, imheader = fits.getdata(image_name,header = True)
    ymax,xmax = image.shape
    xc = xmax/2.
    yc = ymax/2.
    #image_wcs = WCS(imheader)
    # get image dimensions in deg,deg
    #dxdeg,dydeg = get_image_size_deg(image_name)


    
    # Get coord of image center.  will use when getting gaia stars
    #racenter,deccenter = get_image_center_deg(image_name) 

    racenter, deccenter, dxdeg, dydeg = get_image_footprint_box_deg(
        image_name,
        buffer_deg=0.03,
        )

    
    # -------------------------------------------------
    # 2 Try astroquery rectangle search
    # -------------------------------------------------

    print("Querying Gaia via astroquery...")

    brightstar = gaia_stars_in_rectangle(
        racenter,
        deccenter,
        dydeg,
        dxdeg,
        )    
    brightstar = gaia_stars_in_rectangle(
        racenter,
        deccenter,
        dydeg + 0.01,
        dxdeg + 0.01,
    )

    if len(brightstar) == 0:
        print("No Gaia stars found in FOV.")
        return None, None, None

    print("Found Gaia stars in FOV.")

    # Compute mask radii from legacy survey magnitude-radius relation
    mask_radius = mask_radius_for_mag(
        brightstar["phot_g_mean_mag"]
    )



    
    brightstar["radius"] = mask_radius

    # Convert to pixel coords
    starcoord = SkyCoord(
        brightstar["ra"],
        brightstar["dec"],
        frame="icrs",
        unit="deg",
    )

    xpix, ypix = image_wcs.world_to_pixel(starcoord)

    brightstar["xpixel"] = xpix
    brightstar["ypixel"] = ypix

    # convert star radius to pixels
    
    


 

    # -------------------------------------------------
    # 4 Cache result
    # -------------------------------------------------
    if use_cache:
        brightstar.write(outfile, format="csv", overwrite=True)

    return brightstar, brightstar["xpixel"], brightstar["ypixel"]




def make_gaia_mask(
    mask_array,
    x_pixels,
    y_pixels,
    pixel_scale_deg,
    gaia_table=None,
    radius_scale_factor=1):
    """
    Create a Gaia bright-star mask using a magnitude-radius relation.
    
    Parameters
    ----------
    mask_array : 2D numpy array
        Existing mask array.
    gaia_table : astropy Table
        Gaia stars with 'phot_g_mean_mag' and 'radius' (deg).
    x_pixels, y_pixels : array-like
        Star positions in pixel coordinates.
    pixel_scale_deg : float
        Pixel scale in degrees per pixel.
    
    Returns
    -------
    updated_mask : 2D numpy array
        Mask with Gaia stars added.
    gaia_mask : 2D numpy array
        Gaia-only mask.
    """

    if gaia_table is None or len(gaia_table) == 0:
        print("No bright stars on image - woo hoo!")
        return mask_array, np.zeros_like(mask_array)

    ny, nx = mask_array.shape
    gaia_mask = np.zeros_like(mask_array)

    mag = gaia_table["phot_g_mean_mag"]
    rad_deg = gaia_table["radius"]
    rad_pixels = rad_deg / pixel_scale_deg

    mask_value = np.max(mask_array) + 100
    print("mask value =", mask_value)

    for i in range(len(mag)):
        # print(
        #     f"star {i}: "
        #     f"{x_pixels[i]:.1f},"
        #     f"{y_pixels[i]:.1f},"
        #     f"{rad_pixels[i]:.1f}"
        # )

        pixel_mask = circle_pixels(
            float(x_pixels[i]),
            float(y_pixels[i]),
            float(rad_pixels[i]*radius_scale_factor),
            nx,
            ny,
        )

        gaia_mask[pixel_mask] = mask_value

    updated_mask = mask_array + gaia_mask

    return updated_mask, gaia_mask


def galaxy_overlaps_bright_star(
    ra_deg,
    dec_deg,
    gaia_table,
    mag_limit=10,
    radius_col=None,
    min_radius_arcsec=None,
):
    """
    Check whether galaxy center overlaps a bright-star mask.

    Input
    -----
    min_radius_arcsec: float
       minimum size of gaia mask.  could be e.g. 4x the image FWHM
    Returns
    -------
    overlap_flag : bool
    min_distance_arcsec : float
    nearest_radius_arcsec : float
    nearest_mag : float
    """
    if len(gaia_table) == 0:
        return False, np.nan, np.nan, np.nan

    # magnitude column
    if "phot_g_mean_mag" in gaia_table.colnames:
        mag_col = "phot_g_mean_mag"
    elif "Gmag" in gaia_table.colnames:
        mag_col = "Gmag"
    else:
        raise KeyError("Could not find Gaia magnitude column")

    # coordinates
    if "ra" in gaia_table.colnames and "dec" in gaia_table.colnames:
        ra_col, dec_col = "ra", "dec"
    elif "RA" in gaia_table.colnames and "DEC" in gaia_table.colnames:
        ra_col, dec_col = "RA", "DEC"
    else:
        raise KeyError("Could not find Gaia RA/DEC columns")

    # cut table to include bright stars
    bright = gaia_table[np.isfinite(gaia_table[mag_col]) & (gaia_table[mag_col] <= mag_limit)]

    if len(bright) == 0:
        return False, np.nan, np.nan, np.nan

    # radii in degrees
    if radius_col is not None and radius_col in bright.colnames:
        radii_deg = np.array(bright[radius_col], dtype=float)
    else:
        radii_deg = np.array(mask_radius_for_mag(bright[mag_col]), dtype=float)

    # use min radius provided, relevant for fainter stars
    if min_radius_arcsec is not None:
        min_radius_deg = min_radius_arcsec / 3600.0
        radii_deg = np.maximum(radii_deg, min_radius_deg)

    galcoord = SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs")
    #starcoord = SkyCoord(bright[ra_col] * u.deg, bright[dec_col] * u.deg, frame="icrs")
    starcoord = SkyCoord(bright[ra_col], bright[dec_col], frame="icrs")

    sep = galcoord.separation(starcoord)
    sep_arcsec = sep.arcsec

    imin = np.argmin(sep_arcsec)
    min_distance_arcsec = float(sep_arcsec[imin])
    nearest_radius_arcsec = float(radii_deg[imin] * 3600.0)
    nearest_mag = float(bright[mag_col][imin])

    overlap_flag = np.any(sep.degree < radii_deg)

    return overlap_flag, min_distance_arcsec, nearest_radius_arcsec, nearest_mag




##########################
## TESTING
##########################

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Test Gaia rectangle query")
    parser.add_argument("--ra", type=float, required=True)
    parser.add_argument("--dec", type=float, required=True)
    parser.add_argument("--dx", type=float, required=True)
    parser.add_argument("--dy", type=float, required=True)

    args = parser.parse_args()

    stars = gaia_stars_in_rectangle(args.ra, args.dec, args.dx, args.dy)
    print(stars)
    # Example usage:
    #ra = 120.0  # Example right ascension in degrees
    #dec = 45.0   # Example declination in degrees
    #radius = 0.05  # Example radius in degrees
    #result = gaia_stars_in_rectangle(ra, dec, radius,radius)
    #print(result)
