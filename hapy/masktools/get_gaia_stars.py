import numpy as np
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy import units as u

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


if __name__ == '__main__':
    # Example usage:
    ra = 120.0  # Example right ascension in degrees
    dec = 45.0   # Example declination in degrees
    radius = 0.05  # Example radius in degrees
    result = gaia_stars_in_rectangle(ra, dec, radius,radius)
    print(result)
