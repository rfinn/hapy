#!/usr/bin/env python

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
#import ccdproc

from astropy.io.fits import Header
import numpy as np


try:
    from photutils import detect_threshold, detect_sources#, make_source_mask
except ImportError:
    from photutils.segmentation import detect_threshold, detect_sources#, make_source_mask

from astropy.io.fits import Header
import numpy as np

def subtract_median_sky(data,getstd=False,getmedian=True,subtract=True,weightimage=None):
    ''' 
    subtract median sky from image data 

    data: 2d array to estimate median for
    weightimage = 2d array with zero values indicating pixels to ignore


    '''
    # check to see if data is all zeros
    if np.all(data == 0):
        median = 0
        std = 0
        return data,median,std
    elif np.all(np.isnan(data)): # check to see if data is all nans
        # return NAN values for median        
        median = np.nan
        std = 0
        return data,median,std

    if weightimage is not None:
        image_mask = weightimage == 0 # like for mosaic
        data = np.ma.array(data, mask=image_mask)

    try:
        from photutils import make_source_mask
        mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
    except ImportError:
        # check to see if all values in the data is zeros
        threshold = detect_threshold(data, nsigma=3)
        segmentation = detect_sources(data, threshold, npixels=5)        
        #mask = segmentation.make_source_mask(data)
        mask = segmentation.make_source_mask(size=5) # adds a dilation factor

    #mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
    masked_data = np.ma.array(data,mask=mask)
    #clipped_array = sigma_clip(masked_data,cenfunc=np.ma.mean)

    # filled masked values with nans
    nan_filled_data = masked_data.filled(np.nan)
    
    mean,median,std = sigma_clipped_stats(nan_filled_data,sigma=3.0)# removing this ,cenfunc=np.ma.mean)
    if subtract:
        data -= median
    if getstd:
        return data,median,std
    
    elif getmedian:
        return data,median
    
    else:
        return data
    

"""
def subtract_median_sky(data,getstd=False,getmedian=True,subtract=True):
    ''' subtract median sky from image data '''
    try:
        from photutils import make_source_mask
        mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
        masked_data = np.ma.array(data,mask=mask)
        #clipped_array = sigma_clip(masked_data,cenfunc=np.ma.mean)

    except ImportError: # using a more recent version of photutils

        from photutils.segmentation import SegmentationImage
        from photutils.segmentation import detect_sources
        from photutils.background import Background2D, MedianBackground

        bkg_estimator = MedianBackground()

        bkg = Background2D(data,(50, 50),filter_size=(3, 3), bkg_estimator=bkg_estimator)
        threshold = 3 * bkg.background_rms
        segmentation_image = detect_sources(data, threshold, npixels=10)
        mask = segmentation_image.data > 0
        mask = segmentation.make_source_mask(size=5) # adds a dilation factor        
        masked_data = np.ma.array(data,mask=mask)
    mean,median,std = sigma_clipped_stats(masked_data,sigma=3.0,cenfunc=np.ma.mean)
    if subtract:
        data -= median


    if getstd:
        return data,median,std
    
    elif getmedian:
        return data,median
    
    else:
        return data

"""
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
    
    image_wcs = WCS(imheader)        
    pscalex,pscaley = image_wcs.proj_plane_pixel_scales()
    


    nrow,ncol = image.shape

    sizex = ncol*pscalex
    sizey = nrow*pscaley
    #print()
    #print("in get_image_size_deg, sizex,sizey = ",sizex,sizey)
    #print()
    # return the size without the unit
    return sizex.value,sizey.value


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

def circle_pixels(xc,yc,r,ximage,yimage):
    '''
    GOAL:
    - return pixel values that lie within a circular aperture within radius r of position (x,y)
    
    INPUT:
    - enter the center xc,yc and radius of circle in pixels
    - also enter x and y dimensions of parent image

    OUTPUT:
    - 2D boolean array with dimension the same as the input image
    - pixel values are true for pixels within circular aperture, false otherwise
    '''

    # add some checks to make sure numbers make sense
    # actually, it works even if center is outside image boundaries
    #if (xc < 0) | (xc > ximage) | (yc < 0) | (yc > yimage):
    #    print('invalid central coordinates in circle_pixels')
        
    rows,cols = np.mgrid[0:yimage,0:ximage]
    distance = np.sqrt((rows-yc)**2+(cols-xc)**2)
    pixel_flag = distance < r
    return pixel_flag
