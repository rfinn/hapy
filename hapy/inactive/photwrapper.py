#!/usr/bin/env python

'''
This is going to be the wrapper to do photometry on the detected galaxies.

Would like to build this totally on photutils.

Useful references:


https://photutils.readthedocs.io/en/stable/segmentation.html
- detecting sources
- creating a segmentation image
- getting source properties (including total flux, Gini coefficient!!!)
- defining elliptical apertures for sources

NOTES:

'''
try:
    from photutils import detect_threshold, detect_sources#, make_source_mask
except ImportError:
    from photutils.segmentation import detect_threshold, detect_sources#, make_source_mask

    
# changing to remove deprecated function source_properties
#from photutils import source_properties
from photutils.segmentation import SourceCatalog

try:
    from photutils import Background2D, MedianBackground
except ImportError:
    from photutils.background import Background2D, MedianBackground

try:
    from photutils import EllipticalAperture
    from photutils import aperture_photometry    
except ImportError:
    from photutils.aperture import EllipticalAperture
    from photutils.aperture import aperture_photometry
    
from photutils.utils import calc_total_error
from photutils.isophote import EllipseGeometry, Ellipse

from photutils.morphology import gini

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from astropy.stats import sigma_clip, SigmaClip,sigma_clipped_stats
from astropy.visualization import simple_norm
from astropy.utils import lazyproperty

import scipy.ndimage as ndi


from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

import numpy as np
import sys

import time
start_time = time.time()

import matplotlib
#matplotlib.use('Qt5Agg')


# modules in halphagui
from hapy.hatools import utils

## filter information
## from https://www.noao.edu/kpno/mosaic/filters/
central_wavelength = {'4':6620.52,'8':6654.19,'12':6698.53,'16':6730.72,'R':6513.5,'r':6292.28,'inthalpha':6568.,'intha6657':6657,'intr':6240} # angstrom
dwavelength = {'4':80.48,'8':81.33,'12':82.95,'16':81.1,'R':1511.3,'r':1475.17,'inthalpha':95.,'intha6657':80,'intr':1347} # angstrom

# define colors - need this for plotting line and fill_between in the same color
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def display_image(image,percent=99.9,lowrange=False,mask=None,sigclip=True):
    lowrange=False
    if sigclip:
        clipped_data = sigma_clip(image,sigma_lower=5,sigma_upper=5)#,grow=10)
    else:
        clipped_data = image
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')
    #v1,v2=scoreatpercentile(image,[.5,99.5])            
    #plt.imshow(image, cmap='gray_r',vmin=v1,vmax=v2,origin='lower')    


def get_M20(catalog,objectIndex):
    """ 
    calculate M20 according to Lotz+2004 for central object only
    https://iopscience.iop.org/article/10.1086/421849/fulltext/  

    PARAMS:
    * catalog - this is a photutils.SourceCatalog
    * objectIndex - the object number in the catalog

    RETURNS:
    * M20

    NOTES:
    cat.data_ma = A 2D MaskedArray cutout from the data using the minimal bounding box of the source.
    cat.segment = A 2D ndarray cutout of the segmentation image using the minimal bounding box of the source.
    cat.cutout_centroid = The (x, y) coordinate, relative to the cutout data, of the centroid within the isophotal source segment.

    """
    # total second-order moment Mtot is the flux in each pixel fi multiplied by distance^2 of pixel to center,
    # summed over all pixels assigned to the segmentation map

    objNumber = catalog.label[objectIndex]

    # data_ma returns a 2D array with the unmasked values
    # using the min bounding box that fits the galaxy    
    dat = catalog.data_ma[objectIndex]

    # create flag for pixels associated with object in segmentation map
    segflag = catalog.segment[objectIndex] == objNumber

    # get the center coordinates of the object
    xc,yc = catalog.cutout_centroid[objectIndex]

    # can't make sense of moments that photutils includes in the catalog, so recalculating here
    ymax,xmax = catalog.data_ma[objectIndex].shape

    # create a meshgrid to represent pixels in segmentation
    X,Y = np.meshgrid(np.arange(xmax),np.arange(ymax))

    # calculate distance of each point from center of galaxy
    # this is for 2nd order moment
    distsq = (X-xc)**2 + (Y-yc)**2

    # second Moment total
    # segflag ensures that we are just counting the pixels assoc with object
    second_moment_tot = np.sum(dat[segflag]*distsq[segflag])

    ##
    # getting the second moment of 20% highest pixels
    ##

    # using https://github.com/vrodgom/statmorph/blob/master/statmorph/statmorph.py#L1430    
    # Calculate threshold pixel value
    # TODONE - update to fix M20
    sorted_pixelvals = np.sort(dat.flatten())
    flux_fraction = np.cumsum(sorted_pixelvals) / np.sum(sorted_pixelvals)
    sorted_pixelvals_20 = sorted_pixelvals[flux_fraction >= 0.8]

    threshold_brightest20 = sorted_pixelvals_20[0]

    #threshold_brightest20 = scoreatpercentile(dat[segflag].flatten(),80)

    # define flag for pixels that contain top 20% of total flux
    brightest20 = dat > threshold_brightest20

    # sum the second moment of brightest 20
    second_moment_20 = np.sum(dat[segflag & brightest20]*distsq[segflag & brightest20])

    # now calculate M20 as
    # M20 = log10(Sum_Mi/Mtot)

    M20 = np.log10(second_moment_20/second_moment_tot)

    return M20

def get_fraction_masked_pixels(catalog,objectIndex):
    """ 
    get area in the segmentation image, 
    masked area, and fraction of pixels masked
    """

    objNumber = catalog.label[objectIndex]
    dat = catalog.data[objectIndex]    
    masked_dat = catalog.data_ma[objectIndex]

    # create flag for pixels associated with object in segmentation map
    goodflag = catalog.segment[objectIndex] == objNumber

    # get number of pixels in the original segmentation image
    number_total = np.sum(goodflag)

    number_masked = number_total - np.sum(goodflag & masked_dat.mask)

    return number_total, number_masked, number_masked/number_total

# read in image and mask

# identify source for photometry

# run detect to detect source

# estimate ellipse parameters from source properties

# run 

class myStatmorph(statmorph.SourceMorphology):
#class myStatmorph(statmorph.source_morphology):

    """
    add on to statmorph 

    * changed to use the same segmentation map for the gini coefficient calculation 
      as it does for the other morph calculations

    """
    
    @lazyproperty
    def _segmap_gini(self):
        '''overwriting function so that it uses the reg segmap'''
        #self._image[self._slice_stamp]        
        segmap = np.array(self._segmap.data== 1,'i')
        return segmap[self._slice_stamp]
    
    def print(self):
        ''' adding a print method to print out the instance variables '''
        for k in self.__dict__.keys():
            if k.startswith('_'):
                continue
            print(f"{k}: {self.__dict__[k]}")
            
        
        
if __name__ == '__main__':
    image = 'MKW8-18216-R.fits'
    mask = 'MKW8-18216-R-mask.fits'
    image2 = 'MKW8-18216-CS.fits'
    nsaid='18045'
    prefix = 'MKW8-'
    nsaid='110430'
    nsaid='157146'
    prefix = 'NRGs27-'
    image = prefix+nsaid+'-R.fits'
    mask = prefix+nsaid+'-R-mask.fits'
    image2 = prefix+nsaid+'-CS.fits'
    # testing on 2017 pointing 1
    # second galaxy has clear halpha but profile is not fit
    # want to make sure we record some size
    image = 'v17p01-N119230-A742747-R.fits'
    rphot_table = 'v17p01-N119230-A742747-R_phot.fits'
    image2 = 'v17p01-N119230-A742747-CS.fits'
    haphot_table = 'v17p01-N119230-A742747-CS_phot.fits'
    mask = 'v17p01-N119230-A742747-R-mask.fits'
    image = 'v17p01-N118647-A8219-R.fits'
    rphot_table = 'v17p01-N118647-A8219-R_phot.fits'
    image2 = 'v17p01-N118647-A8219-CS.fits'
    haphot_table = 'v17p01-N118647-A8219-CS_phot.fits'
    mask = 'v17p01-N118647-A8219-R-mask.fits'
    myfilter = '4'
    myratio = .0406
    
    # testing on 2019 pointing 1
    # second galaxy has clear halpha but profile is not fit
    # want to make sure we record some size
    image = 'VFID3623-CGCG118-019-v19p001-R.fits'
    
    rphot_table = 'VFID3623-CGCG118-019-v19p001-R-phot.fits'
    image2 = 'VFID3623-CGCG118-019-v19p001-CS.fits'
    haphot_table = 'VFID3623-CGCG118-019-v19p001-CS-phot.fits'
    mask = 'VFID3623-CGCG118-019-v19p001-R-mask.fits'
    
    myfilter = 'inthalpha'
    myratio = .0356
    #image = 'MKW8-18037-R.fits'
    #mask = 'MKW8-18037-R-mask.fits'
    #image = 'r-18045-R.fits'
    #mask = 'r-18045-R-mask.fits'

    prefix = 'VFID0501-UGC09556-BOK-20210315-VFID0501'
    prefix = 'VFID2772-NGC2964-HDI-20180313-p019'    
    image = prefix+'-R.fits'
    rphot_table = prefix+'-R-phot.fits'
    image2 = prefix+'-CS.fits'
    haphot_table = prefix+'-CS-phot.fits'
    mask = prefix+'-R-mask.fits'
    myfilter = '4'
    myratio = .0497427

    try:
        e = ellipse(image,mask=mask, image2=image2, use_mpl=True,image2_filter=myfilter, filter_ratio=myratio)
    except FileNotFoundError:
        print("so sorry, but no images were loaded")
        print("try e = ellipse(imagename) to start")
    ## print('detect objects')
    ## e.detect_objects()
    ## print('find central')
    ## e.find_central_object()
    ## print('get guess')
    ## e.get_ellipse_guess()
    ## print('draw guess')
    ## e.draw_guess_ellipse_mpl()
    ## print('fit ellipse')
    ## e.fit_ellipse()
    ## print('plot results')
    ## e.draw_fit_results_mpl()

    print("--- %s seconds ---" % (time.time() - start_time))

