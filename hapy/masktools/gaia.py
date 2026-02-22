import numpy as np
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy import units as u

import os
import warnings
from astropy.table import Table

import numpy as np
from .geometry import circle_pixels  # wherever it lives

def mask_radius_for_mag(mag):
    """ 
    function is from legacy pipeline  
    https://github.com/legacysurvey/legacypipe/blob/6d1a92f8462f4db9360fb1a68ef7d6c252781027/py/legacypipe/reference.py#L314-L319
    """
    # Returns a masking radius in degrees for a star of the given magnitude.
    # Used for Tycho-2 and Gaia stars.

    # This is in degrees, and is from Rongpu in the thread [decam-chatter 12099].
    return 1630./3600. * 1.396**(-mag)



def gaia_stars_in_rectangle(ra, dec, height, width, minmag=None, maxmag=18, pmsnrcut=5):
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
        selected_columns = ['SOURCE_ID', 'ra', 'dec', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error']
        stars = result[selected_columns]

    # Require SNR > 5 in proper motion
    keepflag = np.sqrt((stars['pmra']/stars['pmra_error'])**2 + (stars['pmdec']/stars['pmdec_error'])**2) > 5

    # Cut by min/max magnitude, if they are provided
    if maxmag is not None:
        keepflag = keepflag & (stars['phot_g_mean_mag'] < maxmag)
    if minmag is not None:
        keepflag = keepflag & (stars['phot_g_mean_mag'] > minmag)
    
    return stars[keepflag]


def get_gaia_stars(image_name, gaiapath=None, use_cache=True,):
    """
    Retrieve Gaia bright stars within the image FOV.  
    Add radius to the returned table using the Legacy Survey magnitude-radius relation.

    Parameters
    ----------
    image_name : str
        Name of image file (used for cache filename).
    racenter, deccenter : float
        Field center in degrees.
    dxdeg, dydeg : float
        Half-width of search box in degrees.
    image_wcs : astropy.wcs.WCS
        WCS object for pixel conversion.
    xmax, ymax : int
        Image dimensions.
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
    image_wcs = WCS(imheader)
    # get image dimensions in deg,deg
    dxdeg,dydeg = imutils.get_image_size_deg(image_name)
    # Get coord of image center.  will use when getting gaia stars
    racenter,deccenter = imutils.get_image_center_deg(image_name)                


    
    # -------------------------------------------------
    # 2 Try astroquery rectangle search
    # -------------------------------------------------
    try:
        print("Querying Gaia via astroquery...")
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

    # -------------------------------------------------
    # 3 Fallback to local Gaia catalog
    # -------------------------------------------------
    except:
        print("WARNING: astroquery Gaia query failed.")
        print(e)

        if gaiapath is None:
            warnings.warn(
                "No Gaia catalog available — running without bright star masks."
            )
            return None, None, None

        try:
            brightstar = Table.read(gaiapath)

            starcoord = SkyCoord(
                brightstar["ra"],
                brightstar["dec"],
                frame="icrs",
                unit="deg",
            )

            xpix, ypix = image_wcs.world_to_pixel(starcoord)

            buffer = 0.1 * xmax
            in_bounds = (
                (xpix > -buffer)
                & (xpix < xmax + buffer)
                & (ypix > -buffer)
                & (ypix < ymax + buffer)
            )

            # Proper motion SNR cut
            pmflag = np.sqrt(
                brightstar["pmra"] ** 2
                * brightstar["pmra_ivar"]
                + brightstar["pmdec"] ** 2
                * brightstar["pmdec_ivar"]
            ) > 5

            flag = in_bounds & pmflag

            if np.sum(flag) == 0:
                return None, None, None

            brightstar = brightstar[flag]
            xpix = xpix[flag]
            ypix = ypix[flag]

            brightstar["xpixel"] = xpix
            brightstar["ypixel"] = ypix

        except FileNotFoundError:
            warnings.warn(
                f"Cannot find Gaia catalog at {gaiapath}."
            )
            return None, None, None

    # -------------------------------------------------
    # 4 Cache result
    # -------------------------------------------------
    if use_cache:
        brightstar.write(outfile, format="csv", overwrite=True)

    return brightstar, brightstar["xpixel"], brightstar["ypixel"]




def make_gaia_mask(
    mask_array,
    gaia_table = None,
    x_pixels,
    y_pixels,
    pixel_scale_deg, ):
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
        print(
            f"star {i}: "
            f"{x_pixels[i]:.1f},"
            f"{y_pixels[i]:.1f},"
            f"{rad_pixels[i]:.1f}"
        )

        pixel_mask = circle_pixels(
            float(x_pixels[i]),
            float(y_pixels[i]),
            float(rad_pixels[i]),
            nx,
            ny,
        )

        gaia_mask[pixel_mask] = mask_value

    updated_mask = mask_array + gaia_mask

    return updated_mask, gaia_mask



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
