#!/usr/bin/env python

from astropy.io.fits import Header
import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from astropy.stats import SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint



def estimate_sky_stats_photutils(
    data,
    weightimage=None,
    grow_radius=10,
    npixels=10,
    nsigma=2.0,
    clip_sigma=3.0,
    clip_maxiters=10,
):
    """
    Estimate (mean, median, std) of background using an object mask + sigma clipping.
    Returns floats. Does not modify data.
    """
    arr = np.asarray(data)

    # trivial cases
    if np.all(arr == 0):
        return 0.0, 0.0, 0.0
    if np.all(np.isnan(arr)):
        return np.nan, np.nan, 0.0

    user_mask = None
    if weightimage is not None:
        user_mask = (weightimage == 0)

    mask = get_object_mask_photutils(
        arr,
        grow_radius=grow_radius,
        npixels=npixels,
        nsigma=nsigma,
        user_mask=user_mask,
        sigma=clip_sigma,
        maxiters=clip_maxiters,
    )

    mean, median, std = sigma_clipped_stats(
        arr,
        sigma=clip_sigma,
        maxiters=clip_maxiters,
        mask=mask,
    )
    return float(mean), float(median), float(std)

def get_object_mask_photutils(
    data,
    grow_radius=10,
    npixels=10,
    weightimage=None,
    nsigma=2.0,
    sigma=3.0,
    maxiters=5,
):
    """
    Return boolean mask where True means 'exclude' (objects + bad pixels).
    weightimage: pixels with weight==0 are excluded.
    """
    data = np.asarray(data)

    base_mask = ~np.isfinite(data)
    if weightimage is not None:
        base_mask |= (np.asarray(weightimage) == 0)

    sigma_clip = SigmaClip(sigma=sigma, maxiters=maxiters)
    threshold = detect_threshold(data, nsigma=nsigma, sigma_clip=sigma_clip, mask=base_mask)

    segment_img = detect_sources(data, threshold, npixels=npixels, mask=base_mask)
    if segment_img is None:
        return base_mask

    footprint = circular_footprint(radius=grow_radius)
    obj_mask = segment_img.make_source_mask(footprint=footprint)

    return base_mask | obj_mask


def calculate_background_photutils(
    data,
    grow_radius=10,
    npixels=10,
    weightimage=None,
    nsigma=2.0,
    clip_sigma=3.0,
    clip_maxiters=5,
):
    """
    Estimate background mean/median/std using object mask + sigma clipping.
    """
    mask = get_object_mask_photutils(
        data,
        grow_radius=grow_radius,
        npixels=npixels,
        weightimage=weightimage,
        nsigma=nsigma,
        sigma=3.0,
        maxiters=clip_maxiters,
    )
    if weightimage is not None:
        goodmask = weightimage > 0
        badmask = ~goodmask
    else:
        badmask = np.zeros_like(data, 'bool')
    data = np.ma.array(data,mask=badmask)
    mean, median, std = sigma_clipped_stats(
        data,
        sigma=clip_sigma,
        maxiters=clip_maxiters,
        mask=mask,
    )
    return float(mean), float(median), float(std)


# def calculate_background_photutils(data,grow_radius=10, npixels=10):
#     """ from https://photutils.readthedocs.io/en/latest/user_guide/background.html """


#     mask = get_object_mask_photutils(data, grow_radius=grow_radius, npixels=npixels)
#     # calculate mean, median and std in unmasked pixels
#     mean, median, std = sigma_clipped_stats(data, sigma=5.0, mask=mask)
    
#     return mean, median, std


def estimate_and_subtract_sky(data, weightimage=None, subtract=True, **skycfg):
    """
    Estimate sky background using calculate_background_photutils and optionally subtract it.

    Returns
    -------
    data_out : ndarray
        Sky-subtracted data (or copy of original if subtract=False).
    sky_median : float
        Estimated sky median (ADU).
    sky_std : float
        Estimated sky std (ADU/pix).
    """
    arr = np.asarray(data)

    # trivial cases
    if np.all(arr == 0):
        return arr.copy(), 0.0, 0.0
    if np.all(~np.isfinite(arr)):
        return arr.copy(), np.nan, 0.0

    mean, med, std = calculate_background_photutils(arr, weightimage=weightimage, **skycfg)

    # if med is not finite, do not subtract
    if not np.isfinite(med):
        return arr.copy(), float(med), float(std)
    if weightimage is not None:
        goodmask = weightimage > 1
    else:
        goodmask = np.ones_like(arr,'bool')
    out = np.zeros_like(arr)
    if subtract:
        out[goodmask] = (arr[goodmask] - med)
    else:
        out = arr.copy()
    return out, float(med), float(std)




def subtract_median_sky(data, getstd=False, getmedian=True, subtract=True, weightimage=None,
                        **skykw):
    """
    Backwards-compatible wrapper.
    **skykw passed to estimate_sky_stats_photutils (e.g., grow_radius, npixels, nsigma).
    """
    mean, median, std = estimate_sky_stats_photutils(data, weightimage=weightimage, **skykw)

    if subtract and np.isfinite(median):
        data = data - median  # do NOT modify input in-place unless you really want that

    if getstd:
        return data, median, std
    elif getmedian:
        return data, median
    else:
        return data

# def subtract_median_sky(data,getstd=False,getmedian=True,subtract=True,weightimage=None):
#     ''' 
#     subtract median sky from image data 

#     data: 2d array to estimate median for
#     weightimage = 2d array with zero values indicating pixels to ignore


#     '''
#     # check to see if data is all zeros
#     if np.all(data == 0):
#         median = 0
#         std = 0
#         return data,median,std
#     elif np.all(np.isnan(data)): # check to see if data is all nans
#         # return NAN values for median        
#         median = np.nan
#         std = 0
#         return data,median,std

#     if weightimage is not None:
#         image_mask = weightimage == 0 # like for mosaic
#         data = np.ma.array(data, mask=image_mask)

#     try:
#         from photutils import make_source_mask
#         mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
#     except ImportError:
#         # check to see if all values in the data is zeros
#         threshold = detect_threshold(data, nsigma=3)
#         segmentation = detect_sources(data, threshold, npixels=5)        
#         #mask = segmentation.make_source_mask(data)
#         mask = segmentation.make_source_mask(size=5) # adds a dilation factor

#     #mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
#     masked_data = np.ma.array(data,mask=mask)
#     #clipped_array = sigma_clip(masked_data,cenfunc=np.ma.mean)

#     # filled masked values with nans
#     nan_filled_data = masked_data.filled(np.nan)
    
#     mean,median,std = sigma_clipped_stats(nan_filled_data,sigma=3.0)# removing this ,cenfunc=np.ma.mean)
#     if subtract:
#         data -= median
#     if getstd:
#         return data,median,std
    
#     elif getmedian:
#         return data,median
    
#     else:
#         return data
    


# def subtract_median_sky(data,getstd=False,getmedian=True,subtract=True):
#     ''' subtract median sky from image data '''
#     try:
#         from photutils import make_source_mask
#         mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
#         masked_data = np.ma.array(data,mask=mask)
#         #clipped_array = sigma_clip(masked_data,cenfunc=np.ma.mean)

#     except ImportError: # using a more recent version of photutils

#         from photutils.segmentation import SegmentationImage
#         from photutils.segmentation import detect_sources
#         from photutils.background import Background2D, MedianBackground

#         bkg_estimator = MedianBackground()

#         bkg = Background2D(data,(50, 50),filter_size=(3, 3), bkg_estimator=bkg_estimator)
#         threshold = 3 * bkg.background_rms
#         segmentation_image = detect_sources(data, threshold, npixels=10)
#         mask = segmentation_image.data > 0
#         mask = segmentation.make_source_mask(size=5) # adds a dilation factor        
#         masked_data = np.ma.array(data,mask=mask)
#     mean,median,std = sigma_clipped_stats(masked_data,sigma=3.0,cenfunc=np.ma.mean)
#     if subtract:
#         data -= median


#     if getstd:
#         return data,median,std
    
#     elif getmedian:
#         return data,median
    
#     else:
#         return data





def get_pixel_scale_chatgpt(wcs):
    """
    Compute pixel scale from WCS header.
    Assumes square pixels.
    """
    # cdelt in degrees, convert to arcsec
    if wcs.wcs.has_cd():
        cd = wcs.wcs.cd
        scale = np.sqrt(np.abs(cd[0,0]*cd[1,1] - cd[0,1]*cd[1,0])) * 3600.0
    else:
        scale = np.abs(wcs.wcs.cdelt[0]) * 3600.0
    return scale


#def get_pixel_scale(header):
#    hwcs = WCS(header)
#    pixelscale = None
#    try:
#        pixelscale = get_pixel_scale_chatgpt(hwcs)
#    except:
#        try:
#            pixelscale = np.abs(float(header['CD1_1']))*3600. # convert deg/pix to arcsec/pixel
#        except KeyError:
#            pixelscale = np.abs(float(header['PC1_1']))*3600. # Siena pipeline from astronometry.net
#    return pixelscale


def get_pixel_scale_from_filename(filename):
    ''' takes in image header and returns the pixel scale in arcsec  '''

    header = fits.getheader(filename)
    pscale = get_pixel_scale(header)

    return pscale

def get_pixel_scale(imheader):
    ''' takes in image header and returns the pixel scale in arcsec  '''
    from astropy.wcs import WCS
    import astropy.units as u
    
    # get pixel scale from image header
    # convert from degrees/pix to arcsec/pix
    
    ## making more general - not all images have CD1_1 keyword
    #self.pscale = abs(float(self.image_header['CD1_1'])*3600)

    # found a better way to get the pixel scale
    
    image_wcs = WCS(imheader)        
    pscalex,pscaley = image_wcs.proj_plane_pixel_scales()
    
    pscale = pscalex.to(u.arcsec).value # convert deg/pix to arcsec/pix

    return pscale

def get_image_size_deg(imagename):
    ''' takes in image and header and returns (sizex,sizey) dimensions in deg  '''
    from astropy.wcs import WCS
    from astropy.io import fits

    image, imheader = fits.getdata(imagename,header = True)
    image_wcs = WCS(imheader)

    # found a better way to get the pixel scale
    
    #image_wcs = WCS(imheader)        
    pscalex,pscaley = image_wcs.proj_plane_pixel_scales()
    


    nrow,ncol = image.shape

    sizex = ncol*pscalex
    sizey = nrow*pscaley
    #print()
    #print("in get_image_size_deg, sizex,sizey = ",sizex,sizey)
    #print()
    # return the size without the unit
    return sizex.value,sizey.value

def get_image_footprint_box_deg(imagename, buffer_deg=0.03):
    """
    Return center RA/Dec and rectangular footprint size in degrees
    based on the actual WCS footprint, padded by buffer_deg.
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    import numpy as np

    image, imheader = fits.getdata(imagename, header=True)
    wcs = WCS(imheader)

    fp = wcs.calc_footprint()   # shape (4, 2), columns are RA, Dec
    ra = fp[:, 0]
    dec = fp[:, 1]

    ra_min = np.min(ra)
    ra_max = np.max(ra)
    dec_min = np.min(dec)
    dec_max = np.max(dec)

    racenter = 0.5 * (ra_min + ra_max)
    deccenter = 0.5 * (dec_min + dec_max)

    width = (ra_max - ra_min) + 2.0 * buffer_deg
    height = (dec_max - dec_min) + 2.0 * buffer_deg

    return racenter, deccenter, width, height

def get_image_center_deg(imagename):
    ''' takes in image and header and returns ra and dec of center in deg  '''
    from astropy.wcs import WCS
    from astropy.io import fits
    
    image, imheader = fits.getdata(imagename,header = True)
    image_wcs = WCS(imheader)

    nrow,ncol = image.shape


    nrow_center = nrow/2
    ncol_center = ncol/2
    pixel_coord = np.array([[nrow_center,ncol_center]])

    radec = image_wcs.wcs_pix2world(pixel_coord,1,ra_dec_order=True)
    ra,dec = radec[0]
    #print("in get_image_center_deg, ra,dec = ",ra,dec)
    return ra,dec



