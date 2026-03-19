
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
import scipy.ndimage as ndi
import numpy as np
import os
import warnings
warnings.simplefilter("always", RuntimeWarning)
from pathlib import Path
try:
    from photutils import detect_threshold, detect_sources#, make_source_mask
except ImportError:
    from photutils.segmentation import detect_threshold, detect_sources#, make_source_mask

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from astropy.stats import sigma_clip, SigmaClip,sigma_clipped_stats
from astropy.visualization import simple_norm
from astropy.utils import lazyproperty
    
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

from hapy.imagetools import imutils

from hapy.hatools import morphology as morph
from hapy.imagetools.plotting import display_image
from hapy.geometry.adapters import pa_ccw_north_deg_to_photutils_theta_rad
from hapy.io.schemas import PHOT_TABLE_SCHEMA 
# This overwrites the photutils task
#from hapy.masktools.types import EllipseParams 
#import .adapters

## filter information
## from https://www.noao.edu/kpno/mosaic/filters/

translate_filter = {'ha4':'4'}

# TODONE:
# calculate central wavelength and width of these filters from new filter traces

# HDI are 4, 8,
#'ha4 s2k':(6615.45, 60.80)

# TODO - photometry should take the center wavelength and dwavelength from image header
# or pass into EllipsePhotometry
central_wavelength = {'4':6620.52,\
                          '8':6656.68,\
                          '12':6697.77,\
                          '16':6736.52,\
                          'ha4':6620.52,\
                          'ha8':6656.68,\
                          'ha12':6697.77,\
                          'ha16':6736.52,\
                          'ha4 H-alpha+4nm k1010':6620.52,\
                          'ha8 H-alpha+8nm k1011':6634.47,\
                          'ha12 H-alpha+12nm k1012':6671.64,\
                          'ha16 H-alpha+16nm k1013':6705.37,\
                          'r SDSS k1018':6292.28,\
                          'r Harris k1004':6525.5,                          
                          'HDI R':6404.63,\
                          'R':6513.5,\
                          'r':6292.28,\
                          'inthalpha':6568.,\
                          'intha6657':6657,\
                          'intr':6240,
                          } # angstrom
dwavelength = {'4':60.44,\
                   '8':59.88,\
                   '12':62.95,
                   '16':61.80,\
                   'ha4':60.44,\
                   'ha8':59.88,\
                   'ha12':62.95,\
                   'ha16':61.80,\
                   'ha4 H-alpha+4nm k1010':80.48,\
                   'Ha+4nm':80.48,\
                   'ha8 H-alpha+8nm k1011':81.33,\
                   'ha12 H-alpha+12nm k1012':82.53,\
                   'ha16 H-alpha+16nm k1013':80.77,\
                   'r SDSS k1018':1475.17,\
                   'R Harris k1004':1474.04,                   
                   'HDI R':1508.19,\
                   'R':1511.3,\
                   'r':1475.17,
                   'inthalpha':95.,\
                   'intha6657':80,\
                   'intr':1347,\
                   } # angstrom

# define colors - need this for plotting line and fill_between in the same color
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def _fraction_unmasked_pixels(cat, idx):
    """
    Return (n_total, n_unmasked, frac_unmasked) for the object's segment pixels.
    """
    label = cat.label[idx]
    seg = cat.segment[idx] == label
    #print(np.sum(seg), seg)
    npixels = np.sum(seg)
    
    # data_ma is a MaskedArray cutout; its mask marks excluded pixels
    masked = cat.data_ma[idx].mask
    
    #print(masked)
    testmask = np.ones_like(seg,'bool')
    testmask = np.array(masked, 'bool')
    nmasked = np.sum(seg & testmask)
    #print(f"DEBUG: npixels={npixels}, nmasked={nmasked}") 
        

    # plt.figure(figsize=(10,4))
    # plt.subplot(1,2,1)
    # display_image(cat.segment[idx],cmap='viridis')
    # plt.colorbar()
    # plt.subplot(1,2,2)
    # display_image(cat.data_ma[idx],cmap='viridis')
    # plt.savefig('debug_fraction_masked.png')
    n_total = int(np.sum(seg))
    n_unmasked = int(np.sum(seg & ~masked))
    frac_unmasked = (n_unmasked / n_total) if n_total > 0 else np.nan
    print("DEBUG: ",n_total, n_unmasked, frac_unmasked)
    return n_total, n_unmasked, frac_unmasked




def compute_sky_stats(data, mask=None):
    """
    Compute robust sky mean and RMS using sigma clipping.
    Ignores masked pixels if mask provided.
    """
    if mask is not None:
        good = mask == 0
        sample = data[good]
    else:
        sample = data

    mean, median, std = sigma_clipped_stats(
        sample,
        sigma=3.0,
        maxiters=5
    )

    return float(median), float(std)
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

class EllipsePhotometry():
    '''
    class to run photometry routines on image

    INPUT
    * image         - primary image (this is usually rband)
    * image2        - image2 is designed to be the Halpha image, 
                      but it can be any second image whereby you define 
                      the ellipse geometry using image 1, and
                      measure the photometry on image 1 and image 2
    * mask          - mask to apply when measuring photometry.  
                      this is usually created from the r-band image
    * image_frame   - image frame for plotting inside a gui; like if this is called from halphamain
    * use_mpl       - use mpl for testing purposes, before integrating with pyqt gui
    * napertures    - number of apertures for measuring photmetry in. default is 20.
    * apertures     - list of apetures to use, instead of generating them automatically
    * image2_filter - this is used to calculate the flux levels in the second filter.  
                      this should be one of standard halpha filters in our dictionary 
                      ('4','inthalpha','intha6657').
    * filter_ratio  - ratio of flux in image2/image1
    * psf           - used by statmorph (not sure if this is the image or the image name)
    * psf_ha        - used by statmorph



    '''
    def __init__(self, image, image2 = None, mask = None, image_frame=None, use_mpl=False, napertures=20,apertures=None, image2_filter=None, filter_ratio=None,psf=None,psf_ha=None,objra=None,objdec=None,fixcenter=False):
        '''  inputs described above '''

        self.tag = image.replace('.fits','')
        self.image, self.header = fits.getdata(image, header=True)
        self.image_name = image

        try:
            self.magzp = float(self.header['PHOTZP'])
        except KeyError:
            warnings.warn(f"No PHOTZP keyword in image {image} header. \nAssuming ZP=22.5")
            self.magzp = 22.5



        # get image dimensions - will use this to determine the max sma to measure
        self.yimage_max, self.ximage_max = self.image.shape

        self.objra = objra
        self.objdec = objdec
        self.fixcenter = fixcenter
        
        self.pixel_scale = imutils.get_pixel_scale(self.header)        
        # check to see if obj position is passed in - need to do this for off-center objects
        if (objra is not None): # unmask central elliptical region around object
            # get wcs from mask image
            wcs = WCS(self.header)
            
            # get x and y coord of galaxy from (RA,DEC) using mask wcs
            #print(f"\nobject RA={self.objra:.4f}, DEC={self.objdec:.4f}\n")
            self.xcenter,self.ycenter = wcs.wcs_world2pix(self.objra,self.objdec,0)
            self.xcenter_ra = self.xcenter
            self.ycenter_dec = self.ycenter            
            # convert sma to pixels using pixel scale from mask wcs
            #self.pixel_scale = wcs.pixel_scale_matrix[1][1]
            #self.objsma_pixels = self.objsma/(self.pixel_scale*3600)
        
        # image 2 is designed to be the Halpha image, but it can be any second
        # image whereby you define the ellipse geometry using image 1, and
        # measure the photometry on image 1 and image 2
        #
        # self.image2_flag is True is image2 is provided
        if image2 is not None:
            print("running in two image mode")
            self.image2_name = image2
            self.image2,self.header2 = fits.getdata(image2, header=True)
            self.image2_flag = True

            try:
                self.magzp2 = float(self.header2['PHOTZP'])
            except KeyError:
                warnings.warn(f"No PHOTZP keyword in image {image2} header. \nAssuming ZP=22.5")
                self.magzp2 = 22.5

        else:
            self.image2_flag = False
            self.image2 = None
            self.image2_name = None
            self.header2 = None
        self.image2_filter = image2_filter
        self.filter_ratio = filter_ratio

        # will use the gain to calculate the noise in the image
        try:
            self.gain = self.header['GAIN']
        except KeyError:
            print("WARNING: no GAIN keyword in header. Setting gain=1")
            self.gain = 1.
        self.psf = psf
        self.psf_ha = psf_ha

        if self.psf is not None:
            self.psf_data = fits.getdata(self.psf)
        if self.psf_ha is not None:
            self.hpsf_data = fits.getdata(self.psf)            
        # the mask should identify all pixels in the cutout image that are not
        # associated with the target galaxy
        # these will be ignored when defining the shape of the ellipse and when measuring the photometry
        #
        # self.mask_flag is True if a mask is provided
        if mask is not None:
            self.mask_image, self.mask_header = fits.getdata(mask,header=True)
            self.mask_flag = True
            # convert to boolean array with bad pixels = True
            self.boolmask = np.array(self.mask_image,'bool')
            self.masked_image = np.ma.array(self.image, mask = self.boolmask)
            if self.image2_flag:
                self.masked_image2 = np.ma.array(self.image2, mask = self.boolmask)
        else:
            print('not using a mask')
            self.mask_flag = False
            self.mask_image = None
            self.mask_header = None
            self.masked_image = self.image
            if self.image2_flag:
                self.masked_image2 = self.image2
        # image frame for plotting inside a gui
        # like if this is called from halphamain.py
        self.image_frame = image_frame

        # alternatively, for plotting with matplotlib
        # use this if running this code as the main program
        self.use_mpl = use_mpl
        self.napertures = napertures
        # assuming a typical fwhm
        try:
            self.fwhm = self.header['SEFWHM']
        except KeyError:
            try:
                self.fwhm = self.header['FWHM'] * self.pixel_scale
            except KeyError:
                self.fwhm = np.nan
        
    def get_noise_in_aper(self, flux, area):
        ''' calculate the noise in an area '''
        if self.sky_noise is not None:
            # if flux/pixel > n * sky_noise

            
            noise_source_e = np.sqrt(np.abs(flux)*self.gain)
            
            #noise_sky_e = self.sky_noise* np.sqrt(area) *self.gain
            # variance in sky = SUM_1^npix (sky_noise_per_pixel * gain)**2 = area * (skynoise * gain)**2
            # noise in sky = sqrt(area) * skynoise * gain
            noise_sky_e = self.sky_noise * np.sqrt(area) * self.gain
            
            noise_total_e = np.sqrt(noise_source_e**2 + noise_sky_e**2)
            noise_adu = noise_total_e/self.gain
        else:
            noise_adu = np.nan
        return noise_adu

    def run_two_image_phot(self,write1=False):
        ''' 
        batch all of the functions that we run for the gui, including:

        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()
        self.measure_phot()
        self.calc_sb()
        self.convert_units()
        self.get_image2_gini()
        self.get_asymmetry()
        self.write_phot_tables()
        self.write_phot_fits_tables()
        self.get_sky_noise()
        '''

        self.measure_sky()
        # we subtract the sky in each cutout already in get_cutouts.py
        #self.subtract_sky()
        
        #print("detect objects")
        self.detect_objects()
        #print("find central")        
        self.find_central_object() 
        #print("find ellipse guess")               
        self.get_ellipse_guess()
        #print("get frac masked pixels")                
        self.set_guess_ellipse_mask_metrics()
        #print("measure phot")                
        self.measure_phot()
        #print("get M20")                
        #self.get_all_M20()
        
        #self.get_all_frac_masked_pixels()
        #print("calc sb")         
        self.calc_sb()
        #print("convert units")                
        #self.convert_units()
        #print("get asym")        
        ##self.get_image2_gini()
        #try:
        #    self.get_asymmetry()
        #except:
        #    print("WARNING: problem measuring asymmetry")
        #    self.asym = -99
        #    self.asym_err = -99
        #    self.asym2 = -99
        #    self.asym2_err = -99
        
        #print("running statmorph - please be patient...")
        #print()
        ##self.run_statmorph()
        ##self.statmorph_flag = True
        #try:    
        #    self.run_statmorph()
        #    self.statmorph_flag = True
        #except:
        #    self.statmorph_flag = False            
        #    print("WARNING: problem running statmorph")
        #print("writing phot fits tables")
        #self.write_phot_tables()
        self.write_phot_fits_table2_simple()
        #self.get_sky_noise()

        #print()
        #print("finished with photutils")
        #print()
        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()
    
    def run_for_gui(self,runStatmorphFlag=True):
        ''' 
        batch all of the functions that we run for the gui, including:

        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()
        self.measure_phot()
        self.calc_sb()
        self.convert_units()
        self.get_image2_gini()
        self.get_asymmetry()
        self.write_phot_tables()
        self.write_phot_fits_tables()
        self.get_sky_noise()
        '''

        print("measure sky")
        self.measure_sky()

        # we subtract the sky in get_cutouts already - don't redo
        #print("subtract sky")
        #self.subtract_sky()
        
        print("detect objects")
        self.detect_objects()
        print("find central")        
        self.find_central_object() 
        print("find ellipse guess")               
        self.get_ellipse_guess()
        print("get frac masked pixels")
        self.set_guess_ellipse_mask_metrics()        
        print("measure phot")                
        self.measure_phot()
        print("get M20")                
        self.get_all_M20()

        #self.get_all_frac_masked_pixels()
        print("calc sb")         
        self.calc_sb()
        print("convert units")                
        self.convert_units()
        print("get asym")        
        self.get_image2_detected_gini()
        try:
            self.get_asymmetry()
        except:
            warnings.warn(f"{self.tag}: problem measuring asymmetry")
            self.asym = -99
            self.asym_err = -99
            self.asym2 = -99
            self.asym2_err = -99
        #self.run_statmorph()

        self.get_sky_noise()
        print("writing tables")
        #self.write_phot_tables()
        self.write_phot_fits_tables()


        print()
        print("finished with photutils")
        print()
        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()

    def run_statmorph_supervisor(self):
        self.run_statmorph(save_figs=True)
        try:
            self.run_statmorph(save_figs=True)
            print("running statmorph - please be patient...")
            print()
            self.statmorph_flag = True
        except:
            self.statmorph_flag = False            
            print("WARNING: problem running statmorph")

    def run_with_galfit_ellipse(self, xc,yc,BA=1,THETA=0):
        '''
        replicating run_for_gui(), but taking input ellipse geometry from galfit

        '''
        self.measure_sky()

        # we subtract the sky in get_cutouts
        #self.subtract_sky()
        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()

        ### RESET ELLIPSE GEOMETRY USING GALFIT VALUES
        self.xcenter = float(xc)
        self.ycenter = float(yc)
        self.position = (self.xcenter, self.ycenter)
        # leave sma as it is defined in get_ellipse_guess
        # so that we measure photometry over the same approximate region
        # in practice this could be different areas if the ellipticity is very different
        # self.sma = r

        # need to reset b to be consistent with galfit ellipticity
        BA = float(BA)
        THETA = float(THETA)
        #print('THETA inside phot wrapper',THETA, BA)
        self.b = BA*self.sma
        self.eps = 1 - BA

        # THETA is galfit version - convert to photutils
        self.theta = pa_deg_ccw_from_north_to_photutils_theta_rad(THETA)

        # replace sith 
        # EllipticalAperture gives rotation angle in radians from +x axis, CCW
        self.aperture = EllipticalAperture(self.position, self.sma, self.b, theta=self.theta)
        # EllipseGeometry using angle in radians, CCW from +x axis
        self.guess = EllipseGeometry(x0=self.xcenter,y0=self.ycenter,sma=self.sma,eps = self.eps, pa = self.theta)

        ### AND NOW BACK TO OUR REGULAR PROGRAMMING
        #print('measuring phot')
        self.measure_phot()
        #print('measuring phot')
        self.calc_sb()
        #print('measuring converting units')
        self.convert_units()
        #print('writing table')
        #self.get_image2_gini()
        try:
            self.get_asymmetry()
        except:
            print("\nWARNING: problem calculating asymmetry, probably b/c image is rectangular...")
            self.asym = -99
            self.asym_err = -99
            self.asym2 = -99
            self.asym2_err = -99
            print()
        self.write_phot_fits_tables(prefix='GAL')
        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()

    def measure_sky(self):
        skymean, skymedian, skystd = imutils.calculate_background_photutils(self.image)

        self.sky_mean = skymean
        self.sky = skymedian
        self.sky_noise = skystd

        if self.image2 is not None:
            skymean, skymedian, skystd = imutils.calculate_background_photutils(self.image2)

            self.sky_mean2 = skymean
            self.sky2 = skymedian
            self.sky_noise2 = skystd

    def subtract_sky(self):
        print(f"DEBUG: subtracting {self.sky:.2e} from image1")
        self.image -= self.sky
        if self.image2 is not None:
            print(f"DEBUG: subtracting {self.sky2:.2e} from image2")
            self.image2 -= self.sky2

        #TODO add this into image header
    def detect_objects(self, snrcut=1.5,npixels=10):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''               
        if self.mask_flag:
            if np.isfinite(self.sky_noise):
                self.threshold = self.sky_noise
            else:
                self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
            #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
            self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        else:
            if np.isfinite(self.sky_noise):
                self.threshold = self.sky_noise
            else:
            
                self.threshold = detect_threshold(self.image, nsigma=snrcut)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels)
            #self.cat = source_properties(self.image, self.segmentation)
            self.cat = SourceCatalog(self.image, self.segmentation)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation)
            

    def detect_objectsv2(self, snrcut=1.5,npixels=10):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''
        # this is not right, because the mask does not include the galaxy
        # updating based on photutils documentation
        # https://photutils.readthedocs.io/en/stable/background.html
        
        # get a rough background estimate
        sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
        if self.mask_flag:
            threshold = detect_threshold(self.image, nsigma=snrcut,sigclip_sigma=3.0, mask=self.boolmask)
        else:
            threshold = detect_threshold(self.image, nsigma=snrcut,sigclip_sigma=3.0)

        segmentation = detect_sources(self.image, threshold, npixels=npixels)


        # make an object mask, expanding the area using a circular footprint
        #mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
        #mask = make_source_mask(self.image,1.5,10,dilate_size=11)        
        try:
            from photutils import make_source_mask
            mask = make_source_mask(self.image,1.5,10,dilate_size=11)
        except ImportError:
            mask = segmentation.make_source_mask(self.image,1.5,10,dilate_size=11)
        #plt.figure()
        #plt.imshow(mask,origin="lower")
        mean, median, std = sigma_clipped_stats(self.image, sigma=3.0, mask=mask)
        self.sky = mean
        self.sky_noise = std
        self.image -= self.sky

        if self.image2 is not None:
            # measure halpha properties using same segmentation image
            mean, median, std = sigma_clipped_stats(self.image2, sigma=3.0, mask=mask)            
            #self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
            self.sky2 = mean
            self.sky_noise2 = std

            # subtract sky
            self.image2 -= self.sky2
        
        
        # now make a new segmentation image based on the new noise estimate
        # default is 1.5*std
        #self.segmentation = detect_sources(self.image, snrcut*std, npixels=npixels)


        # subtract the sky, again...
        #print("\nsky value = ",self.sky)
        #self.image -= self.sky

        ##
        # old code
        ##


        if self.mask_flag:
            self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
            #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
            self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        else:
            self.threshold = detect_threshold(self.image, nsigma=snrcut)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels)
            #self.cat = source_properties(self.image, self.segmentation)
            self.cat = SourceCatalog(self.image, self.segmentation)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        
        #self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
        #self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
        #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
        #self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
        
    def detect_objects_old(self, snrcut=1.5,npixels=11):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''

        if self.mask_flag:
            self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
            #self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
            self.cat = SourceCatalog(self.image, self.segmentation, mask=self.boolmask)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
        else:
            self.threshold = detect_threshold(self.image, nsigma=snrcut)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels)
            #self.cat = source_properties(self.image, self.segmentation)
            self.cat = SourceCatalog(self.image, self.segmentation)
            if self.image2 is not None:
                # measure halpha properties using same segmentation image
                self.cat2 = SourceCatalog(self.image2, self.segmentation, mask=self.boolmask)
            
        # get average sky noise per pixel
        # threshold is the sky noise at the snrcut level, so need to divide by this
        self.sky_noise = np.mean(self.threshold)/snrcut


    def get_all_M20(self):
        # as a kludge, I am going to set all objects' M20 equal to this value
        # in the end, I will only keep the value for the central object...
        myM20 = get_M20(self.cat,self.objectIndex)
        M20 = morph.m20_from_sourcecatalog(self.cat, self.objectIndex)
        print(f"comparing myM20={myM20:.3f} vs chatgptM20={M20:.3f}")
        allM20 = M20*np.ones(len(self.cat))
        self.cat.add_extra_property('M20',allM20)
        self.M20_1 = M20

        # repeat for image2 if it's included
        if self.image2 is not None:
            #M20 = get_M20(self.cat2,self.objectIndex)
            M20 = morph.m20_from_sourcecatalog(self.cat2, self.objectIndex)            
            allM20 = M20*np.ones(len(self.cat))
            self.cat2.add_extra_property('M20',allM20)
            self.M20_2 = M20
            
    def get_all_frac_masked_pixels(self):
        # as a kludge, I am going to set all objects' masked fraction equal to this value
        # in the end, I will only keep the value for the central object...
        ntotal,nunmasked,frac_unmasked = _fraction_unmasked_pixels(self.cat,self.objectIndex)
        frac_masked = 1.0 - frac_unmasked
        allfmasked = frac_masked*np.ones(len(self.cat))
        self.cat.add_extra_property('MASKEDFRAC',allfmasked)

        
        self.masked_fraction = frac_masked
        self.pixel_area = ntotal
        self.masked_pixel_area = ntotal - nunmasked
        print("DEBUG: masked_fraction = ",self.masked_fraction)
    def get_sky_noise(self):
        '''
        * get the noise in image1 and image2 
        * noise is stored as SKYERR in image header
          - units of sky noise are erg/s/cm^2/arcsec^2
        '''

        # I already calculate this in detect_objects
        # get sky noise for image 1
        #if self.mask_flag:
        #    threshold = detect_threshold(self.image, nsigma=1,mask=self.boolmask)
        #else:
        #    threshold = detect_threshold(self.image, nsigma=snrcut)

        # add sky noise to image 1 header

        if not hasattr(self, "sky_noise"):
            self.sky_noise = np.nan
        if not hasattr(self, "sky"):
            self.sky = np.nan
        if not hasattr(self, "sky_noise2"):
            self.sky_noise2 = np.nan
        if not hasattr(self, "sky2"):
            self.sky2 = np.nan

        #print("DEBUG uconv1,uconv2,pixscale:", self.uconversion1, self.uconversion2, self.pixel_scale)
        #print("DEBUG sky_noise ADU:", self.sky_noise, self.sky_noise2)
        
        sky_noise_erg = self.sky_noise*self.uconversion1/self.pixel_scale**2

        print('r sky noise = ',sky_noise_erg)
        try:
            self.header.set('PHOT_SKY','{:.3e}'.format(self.sky),'sky in ADU')
        except AttributeError:
            print("Warning, self.sky not found, setting to zero")
            self.sky = 0
        self.header.set('SKYSTD','{:.3e}'.format(self.sky_noise),'sky noise in ADU')        
        self.header.set('SKYERR','{:.3e}'.format(sky_noise_erg),'sky noise in erg/s/cm^2/arcsec^2')
        # save files
        fits.writeto(self.image_name,self.image,header=self.header,overwrite=True)
        self.im1_skynoise = sky_noise_erg
        # get sky noise for image 2
        if self.image2 is not None and np.isfinite(self.sky_noise2):
            try:
                sky_noise_erg2 = self.sky_noise2*self.uconversion2/self.pixel_scale**2
                self.header2.set('PHOT_SKY','{:.3e}'.format(self.sky2),'sky in ADU')
                self.header2.set('SKYSTD','{:.3e}'.format(self.sky_noise2),'sky noise in ADU')        
                self.header2.set('SKYERR','{:.3e}'.format(sky_noise_erg2),'sky noise in erg/s/cm^2/arcsec^2')
                fits.writeto(self.image2_name,self.image2,header=self.header2,overwrite=True)
            
            except AttributeError:
                print("Warning, self.sky not found, setting to zero")
                self.sky2 = 0
                self.sky_noise2 = 0
                sky_noise_erg2 = 0
            

            self.im2_skynoise = sky_noise_erg2
            print('ha sky noise = ',sky_noise_erg2)

    def find_central_object(self):
        ''' 
        find the central object in the image and get its objid in segmentation image.
        object is stored as self.objectIndex
        '''

        # TODONE - need to be able to handle objects that are not at the center - should have option to pass in RA/DEC and then do like in maskwrapper
        if self.objra is not None:
            #print()
            print("getting object position from RA and DEC")
            #print()
            xc = self.xcenter_ra
            yc = self.ycenter_dec
        else:
            ydim,xdim = self.image.shape
            xc = xdim/2
            yc = ydim/2            
        distance = np.sqrt((np.ma.array(self.cat.xcentroid) - xc)**2 + (np.ma.array(self.cat.ycentroid) - yc)**2)        
        #distance = np.sqrt((self.cat.xcentroid.value - xdim/2.)**2 + (self.cat.ycentroid.value - ydim/2.)**2)
        # save object ID as the row in table with source that is closest to center

        # check to see if len(distance) is > 1

        if len(distance) > 1:
            try:
                self.objectIndex = np.arange(len(distance))[(distance == min(distance))][0]
            except IndexError:
                print("another $#@$# version change???",np.arange(len(distance))[(distance == min(distance))],len(distance))
                print('x vars: ',self.cat.xcentroid, xc)
                print('y vars: ', self.cat.ycentroid, yc)                
                print(self.cat)
                sys.exit()
        else:
            self.objectIndex = 0
            print("WARNING: only one object in the SourceCatalog!",distance)
        #print(self.objectIndex)
        if self.image2 is not None:
            # the object index in cat 2 is not necessarily the same
            # not sure if this is something I changed or if this has always been the case...
            
            distance = np.sqrt((np.ma.array(self.cat2.xcentroid) - xc)**2 + (np.ma.array(self.cat2.ycentroid) - yc)**2)        
            # save object ID as the row in table with source that is closest to center
            self.objectIndex2 = np.arange(len(distance))[(distance == min(distance))][0]
            
        if self.objra is not None:
            # check that distance of this object is not far from the original position
            xcat = self.cat.xcentroid[self.objectIndex]
            ycat = self.cat.ycentroid[self.objectIndex]

            offset = np.sqrt((xcat-self.xcenter_ra)**2 + (ycat-self.ycenter_dec)**2)
            if offset > 100:
                print()
                print("Hold the horses - something is not right!!!")
            

            #print("")            
            #print(f"comparing xcenter {xcat:.1f} and from ra {self.xcenter_ra:.1f}")
            #print(f"comparing ycenter {ycat:.1f} and from dec {self.ycenter_dec:.1f}")
            #print()



    def run_statmorph(self, save_figs = True):
        try:
            from hapy.hatools.morphology import run_statmorph_for_photometry
        except ImportError as e:
            raise RuntimeError("statmorph not installed; install extras to run morphology") from e

        object_label = int(self.cat.label[self.objectIndex])
        mask = (self.mask_image > 0) if self.mask_image is not None else None

        stem1 = str(Path(self.image_name).with_suffix(""))
        stem2 = str(Path(self.image2_name).with_suffix("")) if self.image2_name else None

        if self.r_gini_mask is None:
            self.r_gini_mask, self.r_gini_seg, self.r_gini_threshold = self.build_rband_gini_mask(snrcut=2.5, npixels=10)
        seg_for_statmorph = self.r_gini_mask.astype(int)
        res = run_statmorph_for_photometry(
            image=self.image,
            segmentation_data=seg_fot_statmorph,
            object_label=1,
            gain=float(self.gain),
            mask=mask,
            psf=getattr(self, "psf_data", None),
            image2=self.image2 if self.image2_flag else None,
            psf2=getattr(self, "hpsf_data", None),
            make_fig=save_figs,
            make_diag=True,
            diag_outfile=f"{stem1}-statmorph-diag.pdf",
        )

        self.morph = res.morph_r
        self.morph2 = res.morph_img2

        #print("res.fig_r is None?", res.fig_r is None)
        #print("res.fig_img2 is None?", res.fig_img2 is None)

        if save_figs and res.fig_r is not None:
            out1 = f"{stem1}-statmorph-r.pdf"
            print("Saving:", out1)
            res.fig_r.savefig(out1)

        if save_figs and res.fig_img2 is not None and stem2 is not None:
            out2 = f"{stem2}-statmorph-ha.pdf"
            print("Saving:", out2)
            res.fig_img2.savefig(out2)


        
 
    def get_image2_detected_gini(self, snrcut=1.5):
        ''' 
        calculate gini coefficient for image2 using pixels that are associated with r-band object ID

        this also calculates the sum and mag of the pixels associated with the central galaxy 
        (not sure why this is done together...)
        
        '''
        if self.mask_flag:
            self.threshold2 = detect_threshold(self.image2, nsigma=snrcut, mask=self.boolmask)
            self.segmentation2 = detect_sources(self.image2, self.threshold2, npixels=10,mask=self.boolmask)
            #self.cat2 = source_properties(self.image2, self.segmentation2, mask=self.boolmask)
            cat2 = SourceCatalog(self.image2, self.segmentation2, mask=self.boolmask)            
        else:
            self.threshold2 = detect_threshold(self.image2, nsigma=snrcut)
            self.segmentation2 = detect_sources(self.image2, self.threshold2, npixels=10)
            #self.cat2 = source_properties(self.image2, self.segmentation2)
            cat2 = SourceCatalog(self.image2, self.segmentation2)            

        '''
        select pixels associated with rband image in the segmentation
        AND
        pixels that are above the SNR cut in the Halpha image (image2)
        '''
        self.gini_pixels = (self.segmentation.data == self.cat.label[self.objectIndex]) & (self.segmentation2.data > 0.)

        #self.tbl = self.cat.to_table()
        self.gini2 = gini(self.image2[self.gini_pixels])
        #self.source_sum2 = np.sum(self.image2[self.gini_pixels])
        #self.source_sum2_erg = self.uconversion1*self.source_sum2
        #self.source_sum2_mag = self.magzp2 - 2.5*np.log10(self.source_sum2)

    def build_rband_gini_mask(self, snrcut=2.5, npixels=10):
        obj_label = int(self.cat.label[self.objectIndex])
        base_mask = self.segmentation.data == obj_label

        if self.mask_flag:
            threshold = detect_threshold(self.image, nsigma=snrcut, mask=self.boolmask)
            seg = detect_sources(self.image, threshold, npixels=npixels, mask=self.boolmask)
        else:
            threshold = detect_threshold(self.image, nsigma=snrcut)
            seg = detect_sources(self.image, threshold, npixels=npixels)

        if seg is None:
            return base_mask, None, threshold

        labels = seg.data[base_mask]
        labels = labels[labels > 0]

        if labels.size == 0:
            # fallback to original mask if stricter segmap misses the target
            return base_mask, seg, threshold

        keep_label = np.bincount(labels).argmax()
        gini_mask = seg.data == keep_label

        return gini_mask, seg, threshold

    def run_hapy_gini(self, nsigma=3.0, npixels=10):
        """
        Compute custom HAPY Gini metrics.

        R_HAPY_GINI:
            Gini of r-band image over the r-band segmentation region.

        H_HAPY_GINI:
            Gini of Halpha image over the r-band segmentation region,
            with Halpha pixels below nsigma * sky_sigma set to zero.
        """
        from hapy.hatools.morphology import compute_gini, plot_hapy_gini_diagnostic
        

        self.r_gini_mask, self.r_gini_seg, self.r_gini_threshold = self.build_rband_gini_mask()
        
        #rmask = (self.segmentation.data == self.cat.label[self.objectIndex])
        rmask = self.r_gini_mask
        # R-band custom Gini over the r-band segmentation region
        rvals = self.image[rmask]
        self.R_HAPY_GINI = compute_gini(rvals, allow_negative=False)
        self.R_HAPY_NPIX = int(np.sum(rmask))
        
        # Default Halpha outputs
        self.H_HAPY_GINI = np.nan
        self.H_HAPY_FILLFRAC = np.nan
        self.H_HAPY_NPIX = int(np.sum(rmask))

        if not getattr(self, "image2_flag", False):
            return

        if self.image2 is None:
            return

        # Set thresholding base on sky noise in halpha image
        sigma_sky = getattr(self, "sky_noise2", None)
        if sigma_sky is None:
            sigma_sky = getattr(self, "sky2", None)

        if sigma_sky is None:
            raise ValueError("Could not find Halpha sky RMS attribute for H_HAPY_GINI.")

        threshold = nsigma * sigma_sky
        self.ha_gini_threshold = threshold
        hvals = np.array(self.image2[rmask], dtype=float)

        # 1D detection mask inside rmask
        det_mask_1d = hvals >= threshold

        # Save detection mask for QC
        self.hapy_ha_detect = np.zeros_like(self.image2, dtype=bool)
        self.hapy_ha_detect[rmask] = det_mask_1d

        # Counts / fill fraction
        self.H_HAPY_NPIX = int(np.sum(det_mask_1d))
        self.H_HAPY_FILLFRAC = (
            self.H_HAPY_NPIX / self.R_HAPY_NPIX if self.R_HAPY_NPIX > 0 else np.nan
            )

        # Set sub-threshold pixels to zero
        hvals[~det_mask_1d] = 0.0

        # Build full-size 2D image used for HAPY Gini plotting/QC
        self.ha_gini_image = np.full_like(self.image2, np.nan, dtype=float)
        self.ha_gini_image[rmask] = hvals

        self.H_HAPY_GINI = compute_gini(hvals, allow_negative=False)

        print("H_HAPY_GINI npix:", hvals.size)
        print("nonzero frac:", np.sum(hvals > 0) / hvals.size)

        plot_hapy_gini_diagnostic(
            r_image=self.image,
            ha_image=self.image2,
            r_gini_mask=self.r_gini_mask,
            ha_detect_mask=self.hapy_ha_detect,
            ha_gini_image=self.ha_gini_image,
            r_hapy_gini=self.R_HAPY_GINI,
            ha_hapy_gini=self.H_HAPY_GINI,
            r_hapy_npix=self.R_HAPY_NPIX,
            ha_hapy_npix=self.H_HAPY_NPIX,
            ha_hapy_fillfrac=self.H_HAPY_FILLFRAC,
            ha_threshold=self.ha_gini_threshold,
            outfile=f"{self.image_name.split('.fits')[0]}-hapy-gini-diag.pdf",
            )


    
    
    def get_asymmetry(self):
        '''
        * goal is to measure the assymetry of the galaxy about its center
        * going to measure asymmetry from pixels in the segmentation image only, so

        '''
        # TODO - need to be able to handle images that are not square
        
        # for pixels in segmentation image of central object
        # (can't figure out a way to do this without looping
        # calculate delta_x and delta_y from centroid

        self.object_pixels = self.segmentation.data == self.cat.label[self.objectIndex]

        #xc = self.cat.xcentroid[self.objectIndex].value
        #yc = self.cat.ycentroid[self.objectIndex].value
        xc = self.cat.xcentroid[self.objectIndex]
        yc = self.cat.ycentroid[self.objectIndex]
        row,col = np.where(self.object_pixels)

        grid_size = 3
        sum_diff = np.zeros((grid_size,grid_size),'f')
        source_sum = np.zeros((grid_size,grid_size),'f')
        for dxc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
            for dyc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                drow = np.array((row-(yc+dyc)),'i')
                dcol = np.array((col-(xc+dxc)),'i')
                row2 = np.array(((yc+dyc) -1*drow),'i')
                col2 = np.array(((xc+dxc) -1*dcol),'i')
                sum_diff[dyc,dxc] = np.sum(np.abs(self.masked_image[row,col] - self.masked_image[row2,col2]))
                # divide by the sum of the original pixel values for object
                source_sum[dyc,dxc] = np.sum(self.image[self.object_pixels])
        asym = sum_diff/source_sum
        #print(asym)
        self.asym = np.min(asym)
        self.asym_err = np.std(asym)
        r,c = np.where(asym == np.min(asym))
        self.asym_center = np.array([r+yc,c+xc])

        print('asymmetry = {:.3f}+/-{:.3f}'.format(self.asym,self.asym_err))
        
        if self.image2_flag:
            print("getting asym for image2")
            # using the same segmentation image at for r-band
            # is this the correct thing to do?  does segmentation2.data need to be > 0?
            self.object_pixels2 = (self.segmentation.data == self.cat.label[self.objectIndex]) #& (self.segmentation2.data > 0.)

            #xc = self.cat.xcentroid[self.objectIndex].value
            #yc = self.cat.ycentroid[self.objectIndex].value

            # using the r-band centroid
            xc = self.cat.xcentroid[self.objectIndex]
            yc = self.cat.ycentroid[self.objectIndex]
            row,col = np.where(self.object_pixels2)
            sum_diff = np.zeros((grid_size,grid_size),'f')
            source_sum = np.zeros((grid_size,grid_size),'f')


            # looks like I am measuring asymmetry myself here?
            # doesn't statmorph do this? - looks like I am not using it
            for dxc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                for dyc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                    drow = np.array((row-(yc+dyc)),'i')
                    dcol = np.array((col-(xc+dxc)),'i')
                    row2 = np.array(((yc+dyc) -1*drow),'i')
                    col2 = np.array(((xc+dxc) -1*dcol),'i')
                    sum_diff[dyc,dxc] = np.sum(np.abs(self.masked_image2[row,col] - self.masked_image2[row2,col2]))
                    # divide by the sum of the original pixel values for object
                    source_sum[dyc,dxc] = np.sum(self.image2[self.object_pixels2])
            asym2 = sum_diff/source_sum
            
            #print('asym2 = ',asym2,r,c)
            # measure halpha asymmetry at pixel where R-band asymmetry is minimum
            try:
                self.asym2 = asym2[r,c][0]
                #print(self.asym2)
            except IndexError:
                try:
                    r,c = np.where(asym == np.min(asym2))
                
                    self.asym2 = asym2[r,c][0]
                except IndexError:
                    self.asym2 = np.nan
                    self.asym2_err = np.nan
                    self.asym2_center = np.nan
                    return
            self.asym2_err = np.std(asym2)
            r,c = np.where(asym == np.min(asym2))
            self.asym2_center = np.array([r+yc,c+xc])
            #print('asymmetry2 = ',self.asym2)
            print('asymmetry = {:.3f}+/-{:.3f}'.format(self.asym2,self.asym2_err))
            '''
            # use all the same images as for r-band measurement
            self.object_pixels2 = (self.segmentation.data == self.cat.id[self.objectIndex])# & (self.segmentation2.data > 0.)

            xc = self.cat.xcentroid[self.objectIndex].value
            yc = self.cat.ycentroid[self.objectIndex].value
            row,col = np.where(self.object_pixels2)

            drow = np.array((row-yc),'i')
            dcol = np.array((col-xc),'i')
            row2 = np.array((yc -1*drow),'i')
            col2 = np.array((xc -1*dcol),'i')
            sum_diff = np.sum(np.abs(self.masked_image2[row,col] - self.masked_image2[row2,col2]))
            # divide by the sum of the original pixel values for object
            source_sum = np.sum(self.image2[self.object_pixels2])
        
        
            self.asym2b = sum_diff/source_sum
            print('asymmetry2 = ',self.asym2b)
            '''
        
    def get_ellipse_guess(self, r=2.5):
        '''
        this gets the guess for the ellipse geometry from the detection catalog 
        '''
        obj = self.cat[self.objectIndex]
        #self.xcenter = obj.xcentroid.value
        #self.ycenter = obj.ycentroid.value

        if not self.fixcenter:
            self.xcenter = float(obj.xcentroid)
            self.ycenter = float(obj.ycentroid)


        #if self.objra is not None:
        #    print("")            
        #    print(f"comparing xcenter {self.xcenter:.1f} and from ra {self.xcenter_ra:.1f}")
        #    print(f"comparing ycenter {self.ycenter:.1f} and from dec {self.ycenter_dec:.1f}")
        #    print()
        self.position = (self.xcenter, self.ycenter)
        #print(self.position,self.xcenter,obj.xcentroid,self.ycenter,obj.ycentroid)
        self.sma = obj.semimajor_sigma.value * r
        self.start_size = self.sma
        self.b = obj.semiminor_sigma.value * r
        self.eps = 1 - self.b/self.sma
        self.gini = obj.gini
        self.source_sum = self.cat[self.objectIndex].segment_flux
        self.sky_centroid = obj.sky_centroid
        # orientation is angle in radians, CCW relative to +x axis
        t = obj.orientation.value
        # orientation: radians CCW from +x axis (photutils-style)
        theta = float(obj.orientation.to(u.rad).value)
        self.theta = theta % np.pi
        try:
            self.aperture = EllipticalAperture(self.position, self.sma, self.b, theta=self.theta)
        except ValueError:
            print("\nTrouble in paradise...")
            print(self.position,self.sma,self.b,self.theta)
            sys.exit()
        # EllipseGeometry using angle in radians, CCW from +x axis
        self.guess = EllipseGeometry(x0=self.xcenter,y0=self.ycenter,sma=self.sma,eps = self.eps, pa = self.theta)

        
        self.photutils_segment_flux = np.nan
        self.photutils_segment_mag = np.nan

        if np.isfinite(self.source_sum) and (self.source_sum > 0):
            self.photutils_segment_flux = float(self.source_sum)
            self.photutils_segment_mag = self.magzp - 2.5 * np.log10(self.source_sum)

    
    def measure_mask_fraction_in_guess_ellipse(self):
        """
        Measure the fraction of masked pixels inside the photutils-derived
        guess ellipse.


        Returns
        -------
        total_pix : float
            Total geometric area of the guess ellipse in pixels.
        unmasked_pix : float
            Unmasked area inside the guess ellipse in pixels.
        masked_fraction : float
            Fraction of the ellipse area that is masked.
        """
        if self.mask_image is not None:
            mask = self.mask_image > 0
        else:
            warnings.warn(f"{self.tag}: No mask provided and self.mask is not set.", RuntimeWarning)
            return np.nan, np.nan, np.nan

        if not hasattr(self, "aperture") or self.aperture is None:
            warnings.warn(f"{self.tag}: Guess ellipse aperture is not defined. Run get_ellipse_guess() first.", RuntimeWarning)
            return np.nan, np.nan, np.nan
        
        # Fractional footprint of the ellipse on the image grid
        #print("DEBUG: self.aperture = ",self.aperture)
        apermask = self.aperture.to_mask(method="exact")
        footprint = apermask.to_image(mask.shape)
        #print("DEBUG: footprint = ",footprint)
        if footprint is None:
            raise ValueError("Ellipse footprint could not be projected onto image.")

        # plt.figure()
        # plt.imshow(footprint)
        # plt.colorbar()
        # plt.savefig("debug_masked_pixels.png")
        
        # Pixels contributing to the aperture
        total_pix = np.sum(footprint)


        # Mask convention: True = masked
        masked_weight = np.sum(footprint * mask.astype(float))
        unmasked_pix = total_pix - masked_weight
        #print(f"DEBUG: total_pix = {total_pix:.1f},masked_weight = {masked_weight:.1f}")
        if total_pix > 0:
            masked_fraction = masked_weight / total_pix
        else:
            masked_fraction = np.nan
        #print(f"DEBUG: unmasked_pix = {unmasked_pix:.1f},masked_fraction = {masked_fraction:.1f}")
        return total_pix, unmasked_pix, masked_fraction

    def set_guess_ellipse_mask_metrics(self, mask=None):
        """
        Compute and store mask metrics for the photutils-derived guess ellipse.
        """
        total_pix, unmasked_pix, masked_fraction = self.measure_mask_fraction_in_guess_ellipse()

        self.area_guess_ellipse_pix = total_pix
        self.area_guess_ellipse_unmasked_pix = unmasked_pix
        self.maskfrac_guess_ellipse = masked_fraction

    def draw_guess_ellipse(self):
        ''' DRAW INITIAL ELLIPSE ON R-BAND CUTOUT '''
        #
        markcolor='magenta'
        markwidth=1
        obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,self.sma, self.sma*(1-self.eps), rot_deg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
        self.markhltag = self.image_frame.canvas.add(obj)
        self.image_frame.fitsimage.redraw()

    def draw_guess_ellipse_mpl(self):
        ''' DRAW INITIAL ELLIPSE ON R-BAND CUTOUT '''
        #
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        #plt.imshow(self.masked_image, cmap='Greys', norm=norm , origin='lower')
        display_image(self.masked_image)#, cmap='Greys', norm=norm , origin='lower')        
        plt.colorbar()
        self.aperture.plot(color='k', lw=1.)
        #plt.show()

    def fit_ellipse(self):
        ''' FIT ELLIPSE '''
        #
        # create instance of photutils.Ellipse
        # https://photutils.readthedocs.io/en/stable/isophote.html
        self.ellipse = Ellipse(self.masked_image, self.guess)
        self.isolist = self.ellipse.fit_image()#sfix_pa = True, step=.5)#, fix_eps=True, fix_center=True)
        self.table = self.isolist.to_table()
        
    def draw_fit_results(self):
        ''' DRAW RESULTING FIT ON R-BAND CUTOUT '''
        markcolor='cyan'
        if len(self.isolist) > 5:
            smas = np.linspace(np.min(self.isolist.sma), np.max(self.isolist.sma), 3)
            objlist = []
            for sma in smas:
                iso = self.isolist.get_closest(sma)
                obj = self.image_frame.dc.Ellipse(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps), rot_deg = np.degrees(iso.pa), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            self.markhltag = self.image_frame.canvas.add(self.coadd.dc.CompoundObject(*objlist))
            self.image_frame.fitsimage.redraw()
        else:
            print('problem fitting ellipse')
    def draw_fit_results_mpl(self):
        ''' draw fit results in matplotlib figure '''
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm, origin='lower')
        apertures = []
        if len(self.isolist) > 5:
            smas = np.linspace(np.min(self.isolist.sma)+2, np.max(self.isolist.sma), 12)
            objlist = []
            for sma in smas:
                iso = self.isolist.get_closest(sma)
                #print(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps),  np.degrees(iso.pa))
                apertures.append(EllipticalAperture((iso.x0,iso.y0),iso.sma, iso.sma*(1-iso.eps), theta = np.degrees(iso.pa)))
            for aperture in apertures:
                aperture.plot(color='white',lw=1.5)
        plt.title(os.path.basename(self.image_name).replace('.fits',''))
        #plt.show()
        
        
        plt.close()
    def show_seg_aperture(self,plotname=None):
        ''' matplotlib plotting to show apertures   '''
        tbl1 = self.cat.to_table()
        cat = self.cat
        r=3.
        apertures = []
        for obj in cat:
            position = np.transpose((obj.xcentroid, obj.ycentroid))
            try:
                a = obj.semimajor_axis_sigma.value * r
                b = obj.semiminor_axis_sigma.value * r
            except AttributeError:
                a = obj.semimajor_sigma.value * r
                b = obj.semiminor_sigma.value * r
                
            theta = obj.orientation.to(u.rad).value
            #print(theta)
            apertures.append(EllipticalAperture(position, a, b, theta=theta))
    
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
        ax1.imshow(self.masked_image, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        #cmap = segm_deblend.make_cmap(random_state=12345)
        ax2.imshow(self.segmentation.data, origin='lower')
        ax2.set_title('Segmentation Image')
        for aperture in apertures:
            try:
                aperture.plot(axes=ax1, color='white', lw=1.5)
                aperture.plot(axes=ax2, color='white', lw=1.5)
            except ValueError: # photutils version 2.2.0
                aperture.plot(ax=ax1, color='white', lw=1.5)
                aperture.plot(ax=ax2, color='white', lw=1.5)
                
        #plt.show()
        if plotname is not None:
            plt.savefig(plotname)
            
            plt.close()



    def _build_aperture_grid(self):
        """
        Build elliptical aperture sequence up to image edge.
        Returns:
        apertures_a, apertures_b, area_total, allellipses
        """
        ct = max(1e-6, abs(np.cos(self.theta)))
        st = max(1e-6, abs(np.sin(self.theta)))

        rmax = np.min([
            (self.ximage_max - self.xcenter) / ct,
            (self.yimage_max - self.ycenter) / st,
            self.xcenter / ct,
            self.ycenter / st,
            ])

        index = np.arange(80)
        # Becky's list of apertures
        apertures = (index + 1) * 0.5 * self.fwhm * (1 + (index + 1) * 0.1)
        apertures_a = apertures[apertures < rmax]
        apertures_b = (1.0 - self.eps) * apertures_a
        area_total = np.pi * apertures_a * apertures_b

        allellipses = [
            EllipticalAperture(
                (self.xcenter, self.ycenter),
                apertures_a[i],
                apertures_b[i],
                self.theta,
                )
            for i in range(len(apertures_a))
            ]

        return apertures_a, apertures_b, area_total, allellipses

    def _measure_single_aperture(self, image, ap, combined_mask=None):
        """
        Measure flux and masked/unmasked area for one aperture.

        Returns:
        flux, area_total, area_unmasked
        """
        area_total = ap.area

        if combined_mask is not None:
            phot_table = aperture_photometry(image, ap, mask=combined_mask)
            area_unmasked = ap.area_overlap(image, mask=combined_mask, method="exact")
        else:
            phot_table = aperture_photometry(image, ap, method="subpixel", subpixels=5)
            area_unmasked = area_total

        flux = phot_table["aperture_sum"][0]
        return float(flux), float(area_total), float(area_unmasked)

    def _compute_annulus_snr(self, flux, prev_flux, area_unmasked, prev_area_unmasked, sky_noise):
        """
        Compute annulus-based SNR metrics.

        Returns:
            dF, dA, snr_total, snr_per_pixel, snr_image_units
        """
        if prev_flux is None:
            dF = flux
            dA = area_unmasked
        else:
            dF = flux - prev_flux
            dA = area_unmasked - prev_area_unmasked

        if dA <= 0 or not np.isfinite(dA):
            return dF, dA, -np.inf, -np.inf, -np.inf

        ave_sb_adu = dF / dA
        noise_per_pixel_adu = np.sqrt((sky_noise * self.gain) ** 2 + self.gain * np.abs(ave_sb_adu)) / self.gain
        snr_per_pixel = ave_sb_adu / noise_per_pixel_adu if noise_per_pixel_adu > 0 else -np.inf

        sigma_ann = self.get_noise_in_aper(dF, dA)
        snr_total = dF / sigma_ann if (sigma_ann is not None and np.isfinite(sigma_ann) and sigma_ann > 0) else -np.inf

        snr_image_units = ave_sb_adu / sky_noise if np.isfinite(sky_noise) and sky_noise > 0 else -np.inf

        return dF, dA, snr_total, snr_per_pixel, snr_image_units

    def _truncate_phot_arrays(self, n):
        """
        Truncate all photometry arrays to length n.
        """
        self.apertures_a = self.apertures_a[:n]
        self.apertures_b = self.apertures_b[:n]
        self.area = self.area[:n]
        self.area_total = self.area_total[:n]
        self.area_unmasked = self.area_unmasked[:n]
        self.masked_fraction = self.masked_fraction[:n]
        self.flux1 = self.flux1[:n]
        self.flux1_err = self.flux1_err[:n]
        self.allellipses = self.allellipses[:n]

        if self.image2_flag:
            self.flux2 = self.flux2[:n]
            self.flux2_err = self.flux2_err[:n]


    def measure_phot(self):
        """
        Measure elliptical aperture photometry and SNR profiles.
        """
        snr_stop = 5
        snr_consecutive = 2
        low_snr_count = 0

        prev_flux1 = None
        prev_flux2 = None
        prev_area = None

        self.apertures_a, self.apertures_b, self.area, self.allellipses = self._build_aperture_grid()

        naper = len(self.apertures_a)
        print(f"\nNumber of apertures = {naper}\n")

        self.area_total = np.zeros(naper, dtype=float)
        self.area_unmasked = np.zeros(naper, dtype=float)
        self.masked_fraction = np.zeros(naper, dtype=float)

        self.flux1 = np.zeros(naper, dtype=float)
        self.flux1_err = np.zeros(naper, dtype=float)

        if self.image2_flag:
            self.flux2 = np.zeros(naper, dtype=float)
            self.flux2_err = np.zeros(naper, dtype=float)

        self.snr_total = []
        self.snr_per_pixel = []
        self.snr_image_units = []

        self.snr_total2 = []
        self.snr_per_pixel2 = []
        self.snr_image_units2 = []

        combined_mask = None
        if self.mask_flag:
            nan_mask = ~np.isfinite(self.image)
            combined_mask = self.boolmask | nan_mask

        for i, ap in enumerate(self.allellipses):
            # --- image1 ---
            flux1, area_total, area_unmasked = self._measure_single_aperture(
                self.image,
                ap,
                combined_mask=combined_mask,
            )

            self.flux1[i] = flux1
            self.area_total[i] = area_total
            self.area_unmasked[i] = area_unmasked
            self.masked_fraction[i] = 1.0 - (area_unmasked / area_total if area_total > 0 else np.nan)

            self.flux1_err[i] = self.get_noise_in_aper(self.flux1[i], self.area_unmasked[i])

            dF, dA, snr_total, snr_per_pixel, snr_image_units = self._compute_annulus_snr(
                self.flux1[i],
                prev_flux1,
                self.area_unmasked[i],
                prev_area,
                self.sky_noise,
            )

            self.snr_total.append(snr_total)
            self.snr_per_pixel.append(snr_per_pixel)
            self.snr_image_units.append(snr_image_units)

            # --- optional stopping criterion ---
            # if snr_total < snr_stop:
            #     low_snr_count += 1
            #     if low_snr_count >= snr_consecutive:
            #         self._truncate_phot_arrays(i + 1)
            #         break
            # else:
            #     low_snr_count = 0

            # --- image2 ---
            if self.image2_flag:
                flux2, _, _ = self._measure_single_aperture(
                    self.image2,
                    ap,
                    combined_mask=combined_mask,
                )
                self.flux2[i] = flux2
                self.flux2_err[i] = self.get_noise_in_aper(self.flux2[i], self.area_unmasked[i])

                dF2, dA2, snr_total2, snr_per_pixel2, snr_image_units2 = self._compute_annulus_snr(
                    self.flux2[i],
                    prev_flux2,
                    self.area_unmasked[i],
                    prev_area,
                    self.sky_noise2 if np.isfinite(self.sky_noise2) else self.sky_noise,
                )

                self.snr_total2.append(snr_total2)
                self.snr_per_pixel2.append(snr_per_pixel2)
                self.snr_image_units2.append(snr_image_units2)

                prev_flux2 = self.flux2[i]

            prev_flux1 = self.flux1[i]
            prev_area = self.area_unmasked[i]

        self.snr_total = np.array(self.snr_total)
        self.snr_per_pixel = np.array(self.snr_per_pixel)
        self.snr_image_units = np.array(self.snr_image_units)

        self.snr_total2 = np.array(self.snr_total2)
        self.snr_per_pixel2 = np.array(self.snr_per_pixel2)
        self.snr_image_units2 = np.array(self.snr_image_units2)

        if self.image2_flag:
            print(f"DEBUG: len(flux1)={len(self.flux1)}, len(flux2)={len(self.flux2)}")
        else:
            print(f"DEBUG: len(flux1)={len(self.flux1)}")

        
    def measure_phot_old(self):
        '''
        # alternative is to use ellipse from detect
        # then create apertures and measure flux

        # rmax is max radius to measure ellipse
        # could cut this off based on SNR
        # or could cut this off based on enclosed flux?
        # or could cut off based on image dimension, and do the cutting afterward
        
        #rmax = 2.5*self.sma
        '''
        
        '''
        this is how becky set the apertures
        a = [0]
        for i in range(1,500)
        a.append(a[i-1] + hwhm + (hwhm*i*.1))
        
        '''

        # initialize conditions to implement a stopping radius
        snr_stop = 5
        snr_consecutive = 2
        low_snr_count = 0

        prev_flux1 = None
        prev_area = None
        prev_flux2 = None
        
        # rmax is set according to the image dimensions
        # look for where the semi-major axis hits the edge of the image
        # could be on side (limited by x range) or on top/bottom (limited by y range)
        # 
        #print('xcenter, ycenter, theta = ',self.xcenter, self.ycenter,self.theta)
        rmax = np.min([(self.ximage_max - self.xcenter)/abs(np.cos(self.theta)),\
                       (self.yimage_max - self.ycenter)/abs(np.sin(self.theta))])
        ct = max(1e-6, abs(np.cos(self.theta)))
        st = max(1e-6, abs(np.sin(self.theta)))
        rmax = np.min([(self.ximage_max - self.xcenter) / ct,
                           (self.yimage_max - self.ycenter) / st,
                           self.xcenter / ct,
                           self.ycenter / st,
                           ])
        #print('print rmax, ximage_max, image.shape = ',rmax,self.ximage_max,self.image.shape)
        '''
        this is how becky set the apertures
        a = [0]
        for i in range(1,500):
        a.append(a[i-1] + hwhm + (hwhm*i*.1))
        
        '''

        # TODO - update apertures to make use of input apertures
        index = np.arange(80) # why do we need 80 apertures???
        apertures = (index+1)*.5*self.fwhm*(1+(index+1)*.1)
        #apertures = (index+1)*self.fwhm*(1+(index+1)*.1)
        # cut off apertures at edge of image
        self.apertures_a = apertures[apertures < rmax]
        print(f"\nNumber of apertures = {len(index)}\n")
        #print('number of apertures = ',len(self.apertures_a))
        #self.apertures_a = np.linspace(3,rmax,40)
        self.apertures_b = (1.-self.eps)*self.apertures_a
        self.area = np.pi*self.apertures_a*self.apertures_b # area of each ellipse


        self.flux1 = np.zeros(len(self.apertures_a),'f')
        self.flux1_err = np.zeros(len(self.apertures_a),'f')
        if self.image2_flag:
            self.flux2 = np.zeros(len(self.apertures_a),'f')
            self.flux2_err = np.zeros(len(self.apertures_a),'f')
        self.allellipses = []

        self.snr_total = [] # total snr in aperture
        self.snr_per_pixel = [] # includes poisson noise from galaxy
        self.snr_image_units = [] # ave/sb to sky noise
        self.snr_total2 = [] # total snr in aperture
        self.snr_per_pixel2 = [] # includes poisson noise from galaxy
        self.snr_image_units2 = [] # ave/sb to sky noise
        
        for i in range(len(self.apertures_a)):
 
            # EllipticalAperture takes rotation angle in radians, CCW from +x axis
            ap = EllipticalAperture((self.xcenter, self.ycenter),self.apertures_a[i],self.apertures_b[i],self.theta)#,ai,bi,theta) for ai,bi in zip(a,b)]
            self.allellipses.append(ap)

            if self.mask_flag:
                # check for nans, and add them to the mask
                nan_mask = ~np.isfinite(self.image)
                combined_mask =  self.boolmask | nan_mask
                self.phot_table1 = aperture_photometry(self.image, ap, mask=combined_mask)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, mask=combined_mask)
            else:
                # subpixel is the method used by Source Extractor
                self.phot_table1 = aperture_photometry(self.image, ap, method = 'subpixel', subpixels=5)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, method = 'subpixel', subpixels=5)

            
            self.flux1[i] = self.phot_table1['aperture_sum'][0]
            # calculate noise
            self.flux1_err[i] = self.get_noise_in_aper(self.flux1[i], self.area[i])

            # --- SNR-based truncation using annulus SNR (R band) ---
            if i == 0:
                dF = self.flux1[i]
                dA = self.area[i] 
            else:
                dF = self.flux1[i] - prev_flux1
                dA = self.area[i] - prev_area
            

            # what Becky did
            # median sb in annulus as signal
            # compared with noise_per_pixel = sqrt(sky_uncertainty**2 + median_sb)
            # snr_per_pixel = median_sb/noise

            ave_sb_adu = dF/dA
            noise_per_pixel_adu = np.sqrt((self.sky_noise*self.gain)**2 + self.gain*np.abs(ave_sb_adu))/self.gain # does not account properly for gain...
            snr_per_pixel = ave_sb_adu/noise_per_pixel_adu

            sigma_ann = self.get_noise_in_aper(dF, dA)
            snr_ann =  dF / sigma_ann if (sigma_ann is not None and np.isfinite(sigma_ann) and sigma_ann > 0) else -np.inf

            # testing measurements in observed units only
            # compare ave sb in image units to sky_noise (measured from image, also in image units)

            snr_image_units = ave_sb_adu/self.sky_noise
            #print(f"DEBUG: a(pix)={self.apertures_a[i]:5.1f},sma(arc)={(self.pixel_scale * self.apertures_a[i]):5.1f}, snr_image_units={snr_image_units:.1f},snr_per_pixel={snr_per_pixel:.1f},snr_ann={snr_ann:.1f}, dF={dF:.1f},sigma_ann={sigma_ann:.1e}, sky_noise={self.sky_noise:.1e}, dA={dA:.1f}")
            self.snr_total.append(snr_ann)
            self.snr_per_pixel.append(snr_per_pixel)
            self.snr_image_units.append(snr_image_units)
            # if snr_ann < snr_stop:
            # #if snr_per_pixel < snr_stop:                
            #     low_snr_count += 1
            #     if low_snr_count >= snr_consecutive:
            #         # truncate arrays to i (inclusive) and stop
            #         n = i + 1
            #         self.apertures_a = self.apertures_a[:n]
            #         self.apertures_b = self.apertures_b[:n]
            #         self.area = self.area[:n]
            #         self.flux1 = self.flux1[:n]
            #         self.flux1_err = self.flux1_err[:n]
            #         if self.image2_flag:
            #             self.flux2 = self.flux2[:n]
            #             self.flux2_err = self.flux2_err[:n]
            #         self.allellipses = self.allellipses[:n]
            #         break
            # else:
            #     low_snr_count = 0


            
            if self.image2_flag:
                self.flux2[i] = self.phot_table2['aperture_sum'][0]
                self.flux2_err[i] = self.get_noise_in_aper(self.flux2[i], self.area[i])
                if i == 0:
                    dF = self.flux2[i] 
                    dA = self.area[i] 
                else:
                    dF = self.flux2[i] - prev_flux1
                    dA = self.area[i] - prev_area

                ave_sb_adu = dF/dA
                noise_per_pixel_adu = np.sqrt((self.sky_noise*self.gain)**2 + self.gain*np.abs(ave_sb_adu))/self.gain # does not account properly for gain...
                snr_per_pixel = ave_sb_adu/noise_per_pixel_adu

                sigma_ann = self.get_noise_in_aper(dF, dA)
                snr_ann =  dF / sigma_ann if (sigma_ann is not None and np.isfinite(sigma_ann) and sigma_ann > 0) else -np.inf

                snr_image_units = ave_sb_adu/self.sky_noise
                self.snr_total2.append(snr_ann)
                self.snr_per_pixel2.append(snr_per_pixel)
                self.snr_image_units2.append(snr_image_units)
                prev_flux2 = self.flux2[i]

            prev_flux1 = self.flux1[i]
            prev_area = self.area[i]

        self.snr_total = np.array(self.snr_total)
        self.snr_per_pixel = np.array(self.snr_per_pixel)
        self.snr_image_units = np.array(self.snr_image_units)
        self.snr_total2 = np.array(self.snr_total2)
        self.snr_per_pixel2 = np.array(self.snr_per_pixel2)
        self.snr_image_units2 = np.array(self.snr_image_units2)
            
        print(f"DEBUG: len(flux1)={len(self.flux1)}, len(flux2)={len(self.flux2)}")
    def draw_phot_apertures(self,plotname=None):
        ''' matplotlib plotting to show apertures; provide a plot name to save the output figure   '''
        tbl1 = self.cat.to_table()
        cat = self.cat
        r=3.
        apertures = []
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6))
        clipped_data = sigma_clip(self.image,sigma_lower=5,sigma_upper=5)#,grow=10)
        norm = simple_norm(clipped_data, stretch='asinh',percent=99)
        
        #display_image(self.image)
        ax1.imshow(self.image, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        #cmap = segm_deblend.make_cmap(random_state=12345)
        ax2.imshow(self.segmentation.data, origin='lower')
        ax2.set_title('Segmentation Image')
        # plot a subset of apertures
        nap = len(self.allellipses)
        plotaps = np.array([0,nap//2,-1],'i')
        for i in plotaps:
            aperture = self.allellipses[i]
            try:
                aperture.plot(axes=ax1, color='c', lw=1.5)
                aperture.plot(axes=ax2, color='white', lw=1.5)
            except ValueError: # photutils version 2.2.0
                aperture.plot(ax=ax1, color='c', lw=1.5)
                aperture.plot(ax=ax2, color='white', lw=1.5)
            
        if plotname is not None:
            plt.savefig(plotname)
            plt.close()
    def calc_sb(self):
        # calculate surface brightness in each aperture

        # first aperture is calculated differently
        self.sb1 = np.zeros(len(self.apertures_a),'f')
        self.sb1_err = np.zeros(len(self.apertures_a),'f')
        print(f"DEBUG: len(apertures_a)={len(self.apertures_a)}, len(flux1)={len(self.flux1)}") 
        self.sb1[0] = self.flux1[0]/self.area[0]
        self.sb1_err[0] = self.get_noise_in_aper(self.flux1[0], self.area[0])/self.area[0]
        # outer apertures need flux from inner aperture subtracted
        for i in range(1,len(self.area)):
            self.sb1[i] = (self.flux1[i] - self.flux1[i-1])/(self.area[i]-self.area[i-1])
            self.sb1_err[i] = self.get_noise_in_aper((self.flux1[i] - self.flux1[i-1]),(self.area[i]-self.area[i-1]))/(self.area[i]-self.area[i-1])

        # calculate SNR to follow Becky's method of cutting off analysis where SNR = 2
        self.sb1_snr = np.abs(self.sb1/self.sb1_err)
        # repeat for image 2 if it is provided
        if self.image2_flag:
            self.sb2 = np.zeros(len(self.apertures_a),'f')
            self.sb2_err = np.zeros(len(self.apertures_a),'f')
            self.sb2[0] = self.flux2[0]/self.area[0]
            self.sb2_err[0] = self.get_noise_in_aper(self.flux2[0], self.area[0])/self.area[0]
            for i in range(1,len(self.area)):
                self.sb2[i] = (self.flux2[i] - self.flux2[i-1])/(self.area[i]-self.area[i-1])
                self.sb2_err[i] = self.get_noise_in_aper((self.flux2[i] - self.flux2[i-1]),(self.area[i]-self.area[i-1]))/(self.area[i]-self.area[i-1])
            self.sb2_snr = np.abs(self.sb2/self.sb2_err)

    def get_filter_properties(self):
        try:
            self.filter_cwavelength_A = float(self.header['FILTWCEN'])
            self.filter_width_A = float(self.header['FILTWID'])
            print(f"DEBUG: got cwave = {self.filter_cwavelength_A:.2f} in header")
        except KeyError:
            print("WARNING: no center wavelength in header FILTWCEN")
            self.filter_cwavelength_A = None
            self.filter_width_A = None

        if self.header2 is not None:
            try:
                self.filter2_cwavelength_A = self.header2['FILTWCEN']
                self.filter2_width_A = self.header2['FILTWID']
                print(f"DEBUG: got cwave2 = {self.filter2_cwavelength_A:.2f} in header")
            except KeyError:
                print("WARNING: no center wavelength in header FILTWCEN")
                self.filter2_cwavelength_A = None
                self.filter2_width_A = None
            
    def convert_units(self):
        '''
        ###########################################################
        ### SET UP INITIAL PARAMETERS TO CALCULATE CONVERSION
        ### FROM ADU/S TO PHYSICAL UNITS
        ###########################################################
        '''

        self.pixel_scale = imutils.get_pixel_scale(self.header)

        # -- make use of newly calculated filter centers and widths
        # -- should be stored in the image header and read in init
        self.get_filter_properties()

        if self.filter_cwavelength_A is not None:
            cwave = self.filter_cwavelength_A * 1.e-10 # convert A to m
            dwave = self.filter_width_A * 1.e-10 # convert A to m
        else:
            # fall back on filter dictionaries, but use with caution!
            print("WARNING: no filter information - using outdated dictionaries!")
            try:
                cwave = central_wavelength[self.header["FILTER"]] * 1.e-10
                dwave = dwavelength[self.header["FILTER"]] * 1.e-10
            except KeyError:
                cwave = 6500. * 1.e-10
                dwave = 1500. * 1.e-10
                                            
        # multiply by bandwidth of filter to convert from Jy to erg/s/cm^2
        bandwidth1 = 3.e8*dwave/cwave**2
        self.uconversion1 = 3631.*10**(self.magzp/-2.5)*1.e-23*bandwidth1

        
        if self.image2_filter:
            if self.filter2_cwavelength_A is not None:

                cwave = self.filter2_cwavelength_A * 1.e-10 # convert A to m
                dwave = self.filter2_width_A * 1.e-10 # convert A to m
            else:
                # fall back on filter dictionaries, but use with caution!
                print("WARNING: no filter information - using outdated dictionaries!")
                try:
                    cwave = central_wavelength[self.header2["FILTER"]] * 1.e-10
                    dwave = dwavelength[self.header2["FILTER"]] * 1.e-10
                except KeyError:
                    cwave = 6600. * 1.e-10
                    dwave = 80. * 1.e-10
                
            bandwidth2 = 3.e8*dwave/(cwave)**2
            try:
                self.magzp2 = float(self.header2['PHOTZP'])
                self.uconversion2 = 3631.*10**(self.magzp2/-2.5)*1.e-23*bandwidth2
            except:
                # use 25 as default ZP if none is provided in header
                self.uconversion2 = 3631.*10**(25/-2.5)*1.e-23*bandwidth2
                print("WARNING: no PHOTZP keyword in image2 header. \nAssuming ZP=22.5")                
                self.magzp2 = 22.5
        else:
            self.uconversion2 = None
            
 
        if self.image2_flag:
            fr = self.filter_ratio
            self.uconversion2b = fr * self.uconversion1 #if fr is not None else np.nan
        #if self.filter_ratio is not None:
        #    if self.image2_flag:
        #        self.uconversion2b = self.filter_ratio*self.uconversion1
        #else:
        #    self.uconversion2b = None
            
        ###########################################################
        ### CONVERT UNITS TO
        ### FLUX -> ERG/S/CM^2
        ### FLUX -> MAG
        ### SURFACE BRIGHTNESS -> ERG/S/CM^2/ARCSEC^2
        ### SURFACE BRIGHTNESS -> MAG/ARCSEC^2
        ###########################################################
        self.sky_noise_erg = self.sky_noise*self.uconversion1/self.pixel_scale**2
        self.flux1_erg = self.uconversion1*self.flux1
        self.flux1_err_erg = self.uconversion1*self.flux1_err
        self.source_sum_erg = self.uconversion1*self.source_sum
        self.source_sum_mag = self.magzp - 2.5*np.log10(self.source_sum)
        self.sb1_erg_sqarcsec = self.uconversion1*self.sb1/self.pixel_scale**2
        self.sb1_erg_sqarcsec_err = self.uconversion1*self.sb1_err/self.pixel_scale**2

        # limit mag calculations to values with positive values of flux/sb
        self.mag1 = np.full_like(self.flux1, np.nan)
        self.mag1_err = np.full_like(self.flux1, np.nan)
        good = self.flux1 > 0
        self.mag1[good] = self.magzp - 2.5*np.log10(self.flux1[good])
        self.mag1_err[good] = self.mag1[good] - (self.magzp - 2.5*np.log10(self.flux1[good] + self.flux1_err[good]))
        

        self.sb1_mag_sqarcsec = np.full_like(self.sb1, np.nan)
        self.sb1_mag_sqarcsec_err = np.full_like(self.sb1, np.nan)
        sb_arcsec = self.sb1/self.pixel_scale**2
        good = self.sb1 > 0
        self.sb1_mag_sqarcsec[good] = self.magzp - 2.5*np.log10(sb_arcsec[good])
        self.sb1_mag_sqarcsec_err[good] = self.sb1_mag_sqarcsec[good] - (self.magzp - 2.5*np.log10((self.sb1[good] + self.sb1_err[good])/self.pixel_scale**2))
        
        if self.image2_flag and (self.uconversion2 is not None):
            self.flux2_erg = self.uconversion2*self.flux2
            self.flux2_err_erg = self.uconversion2*self.flux2_err
            self.source_sum2 = self.cat2.segment_flux[self.objectIndex]
            self.source_sum2_erg = self.uconversion2*self.cat2.segment_flux[self.objectIndex]
            self.source_sum2_mag = self.magzp2 - 2.5*np.log10(self.source_sum2)
            
            # limit mag calculations to values with positive values of flux/sb
            self.mag2 = np.full_like(self.flux2, np.nan)
            self.mag2_err = np.full_like(self.flux2, np.nan)
            good = self.flux2 > 0
            self.mag2[good] = self.magzp2 - 2.5*np.log10(self.flux2[good])
            self.mag2_err[good] = self.mag2[good] - (self.magzp2 - 2.5*np.log10(self.flux2[good] + self.flux2_err[good]))
        

            self.sb2_mag_sqarcsec = np.full_like(self.sb2, np.nan)
            self.sb2_mag_sqarcsec_err = np.full_like(self.sb2, np.nan)
            sb2_arcsec = self.sb2/self.pixel_scale**2
            good = self.sb2 > 0
            self.sb2_mag_sqarcsec[good] = self.magzp2 - 2.5*np.log10(sb2_arcsec[good])
            self.sb2_mag_sqarcsec_err[good] = self.sb2_mag_sqarcsec[good] - (self.magzp2 - 2.5*np.log10((self.sb2[good] + self.sb2_err[good])/self.pixel_scale**2))

            self.sb2_erg_sqarcsec = self.uconversion2*self.sb2/self.pixel_scale**2
            self.sb2_erg_sqarcsec_err = self.uconversion2*self.sb2_err/self.pixel_scale**2
            
            
            # this next set uses the filter ratio and the filter 1 flux conversion to
            # convert narrow-band flux (filter 2) to physical units.
            #if self.uconversion2b:
            #    conversion = self.uconversion2b
            #    self.flux2_erg = conversion*self.flux2
            #    self.flux2_err_erg = conversion*self.flux2_err
            #    self.sb2_erg_sqarcsec = self.uconversion2*self.sb2/self.pixel_scale**2
            #    self.sb2_erg_sqarcsec_err = self.uconversion2*self.sb2_err/self.pixel_scale**2

                
            #self.sb2_mag_sqarcsec = self.magzp2 - 2.5*np.log10(conversion*self.sb2/self.pixel_scale**2)
                #self.sb2_mag_sqarcsec_err = self.sb2_mag_sqarcsec - (self.magzp2 - 2.5*np.log10(conversion*(self.sb2+self.sb2_err)/self.pixel_scale**2))
                
    def write_phot_tables(self):
        '''
        write out photometry for image and image2 in ascii format
        '''
        
        # radius enclosed flux
        outfile = open(self.image_name.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles

        #outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
        #outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
        outfile.write('# radius flux flux_err sb sb_err sb_snr flux_erg flux_erg_err mag mag_err sb_ergsqarc sb_err_ergsqarc sb_magsqarc sb_err_magsqarc \n')
        for i in range(len(self.apertures_a)):
            s='%.2f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e '% \
                          (self.apertures_a[i],self.flux1[i],self.flux1_err[i],\
                           self.sb1[i], self.sb1_err[i], \
                          self.sb1_snr[i], \
                          self.flux1_erg[i], self.flux1_err_erg[i],\
                          self.mag1[i], self.mag1_err[i], \
                          self.sb1_erg_sqarcsec[i],self.sb1_erg_sqarcsec[i], \
                          self.sb1_mag_sqarcsec[i],self.sb1_mag_sqarcsec[i])
            s=s+'\n'
            outfile.write(s)
        outfile.close()

        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            outfile = open(self.image2_name.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles
    
            #outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
            #outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
            s = '# radius flux flux_err sb sb_err sb_snr flux_erg flux_erg_err mag mag_err sb_ergsqarc sb_err_ergsqarc sb_magsqarc sb_err_magsqarc'
            if self.uconversion2b:
                s = s +' fluxb_erg fluxb_erg_err sbb_ergsqarc sbb_err_ergsqarc sbb_magsqarc sbb_err_magsqarc'
            s = s + '\n'
            outfile.write(s)
            for i in range(len(self.apertures_a)):
                s = '%.2f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e '% \
                              (self.apertures_a[i],self.flux2[i],self.flux2_err[i],\
                               self.sb2[i], self.sb2_err[i],self.sb2_snr[i],\
                              self.flux2_erg[i], self.flux2_err_erg[i],\
                              self.mag2[i], self.mag2_err[i], \
                              self.sb2_erg_sqarcsec[i],self.sb2_erg_sqarcsec[i], \
                              self.sb2_mag_sqarcsec[i],self.sb2_mag_sqarcsec[i])
                if self.uconversion2b:
                    s=s+' %.3e %.3e %.3e %.3e %.3e %.3e'% \
                      (self.flux2b_erg[i], self.flux2b_err_erg[i],\
                      self.sb2b_erg_sqarcsec[i],self.sb2b_erg_sqarcsec[i], \
                      self.sb2b_mag_sqarcsec[i],self.sb2b_mag_sqarcsec[i])
                s = s+'\n'
                outfile.write(s)

            outfile.close()


    def _build_phot_table(
        self,
        flux,
        flux_err,
        sb,
        sb_err,
        sb_snr,
        flux_erg,
        flux_err_erg,
        mag,
        mag_err,
        sb_erg_sqarcsec,
        sb_erg_sqarcsec_err,
        sb_mag_sqarcsec,
        sb_mag_sqarcsec_err,
        snr_total,
        snr_per_pixel,
        snr_image_units,
    ):
        """
        Build one photometry table (R or image2).

        Each row corresponds to one elliptical aperture ordered by increasing
        semi-major axis. Flux quantities are cumulative within the aperture.
        Surface-brightness quantities are averaged over the unmasked aperture area.
        """

        values = {
            "ap_index": np.arange(len(self.apertures_a)),
            "sma_arcsec": self.apertures_a * self.pixel_scale,
            "sma_pix": self.apertures_a,
            "area_total_pix": self.area_total,
            "area_unmasked_pix": self.area_unmasked,
            "masked_fraction": self.masked_fraction,
            "flux_cum": flux,
            "flux_cum_err": flux_err,
            "sb_avg": sb,
            "sb_avg_err": sb_err,
            "sb_avg_snr": sb_snr,
            "flux_cgs": flux_erg,
            "flux_cgs_err": flux_err_erg,
            "mag_cum": mag,
            "mag_cum_err": mag_err,
            "sb_cgs_arcsec2": sb_erg_sqarcsec,
            "sb_cgs_arcsec2_err": sb_erg_sqarcsec_err,
            "sb_mag_arcsec2": sb_mag_sqarcsec,
            "sb_mag_arcsec2_err": sb_mag_sqarcsec_err,
            "snr_total": snr_total,
            "snr_per_pixel": snr_per_pixel,
            "snr_image_units": snr_image_units,
        }

        columns = [
            Column(values[name], name=name, unit=unit, description=desc)
            for name, unit, desc in PHOT_TABLE_SCHEMA
        ]

        return Table(columns)
            


    def write_phot_fits_tables(self, prefix=None):
        """Write out photometry tables for image1 and image2 in FITS format."""

        def _outfile(image_name):
            stem = image_name.split(".fits")[0]
            if prefix is None:
                return stem + "_phot.fits"
            return stem + f"-{prefix}_phot.fits"

        # --- image1 table ---
        t1 = self._build_phot_table(
            flux=self.flux1,
            flux_err=self.flux1_err,
            sb=self.sb1,
            sb_err=self.sb1_err,
            sb_snr=self.sb1_snr,
            flux_erg=self.flux1_erg,
            flux_err_erg=self.flux1_err_erg,
            mag=self.mag1,
            mag_err=self.mag1_err,
            sb_erg_sqarcsec=self.sb1_erg_sqarcsec,
            sb_erg_sqarcsec_err=self.sb1_erg_sqarcsec_err,
            sb_mag_sqarcsec=self.sb1_mag_sqarcsec,
            sb_mag_sqarcsec_err=self.sb1_mag_sqarcsec_err,
            snr_total=self.snr_total,
            snr_per_pixel=self.snr_per_pixel,
            snr_image_units=self.snr_image_units,
        )

        t1.meta["IMAGE"] = self.image_name
        t1.meta["BAND"] = "image1"
        t1.meta["PIXELSCL"] = self.pixel_scale

        t1.write(_outfile(self.image_name), format="fits", overwrite=True)
        self.photfile = _outfile(self.image_name)
        # --- image2 table ---
        if self.image2_flag:
            t2 = self._build_phot_table(
                flux=self.flux2,
                flux_err=self.flux2_err,
                sb=self.sb2,
                sb_err=self.sb2_err,
                sb_snr=self.sb2_snr,
                flux_erg=self.flux2_erg,
                flux_err_erg=self.flux2_err_erg,
                mag=self.mag2,
                mag_err=self.mag2_err,
                sb_erg_sqarcsec=self.sb2_erg_sqarcsec,
                sb_erg_sqarcsec_err=self.sb2_erg_sqarcsec_err,
                sb_mag_sqarcsec=self.sb2_mag_sqarcsec,
                sb_mag_sqarcsec_err=self.sb2_mag_sqarcsec_err,
                snr_total=self.snr_total2,
                snr_per_pixel=self.snr_per_pixel2,
                snr_image_units=self.snr_image_units2,
            )

            t2.meta["IMAGE"] = self.image2_name
            t2.meta["BAND"] = "image2"
            t2.meta["PIXELSCL"] = self.pixel_scale

            t2.write(_outfile(self.image2_name), format="fits", overwrite=True)
            self.photfile2 = _outfile(self.image2_name)


    def write_phot_fits_table1_simple(self, prefix=None):
        ''' write out photometry for image and image2 in fits format '''

        if prefix is None:
             outfile = self.image_name.split('.fits')[0]+'_phot.fits'
        else:
             outfile = self.image_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
        print('photometry outfile = ',outfile)

        data = [self.apertures_a*self.pixel_scale,self.apertures_a, \
             self.flux1,self.flux1_err,\
             self.sb1, self.sb1_err, \
             self.sb1_snr]
             #self.flux1_erg, self.flux1_err_erg,\
             #self.mag1, self.mag1_err, \
             #self.sb1_erg_sqarcsec,self.sb1_erg_sqarcsec_err, \
             #self.sb1_mag_sqarcsec,self.sb1_mag_sqarcsec_err]

        names = ['sma_arcsec','sma_pix','flux','flux_err',\
                 'sb', 'sb_err', \
                 'sb_snr', ]

        units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '']


        #self.sky_noise,self.sky_noise_erg]
        #'sky_noise_ADU_sqpix','sky_noise_erg_sqarcsec']
        #u.adu/u.s/u.pixel**2,u.erg/u.s/u.cm**2/u.arcsec**2]        
        columns = []
        for i in range(len(data)):
            columns.append(Column(data[i],name=names[i],unit=units[i]))
        
        t = Table(columns)
        t.write(outfile, format='fits', overwrite=True)
    def write_phot_fits_table2_simple(self, prefix=None):
        """ write out phot for second image only - don't want to overwrite R phot """
        if prefix is None:
             outfile = self.image_name.split('.fits')[0]+'_phot.fits'
        else:
             outfile = self.image_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
        
        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            if prefix is None:
                outfile = self.image2_name.split('.fits')[0]+'_phot.fits'
            else:
                outfile = self.image2_name.split('.fits')[0]+'-'+prefix+'_phot.fits'
    
            data = [self.apertures_a*self.pixel_scale*3600,self.apertures_a, \
                self.flux2,self.flux2_err,\
                self.sb2, self.sb2_err, \
                self.sb2_snr]
            names = ['sma_arcsec','sma_pix','flux','flux_err',\
                'sb', 'sb_err', \
                'sb_snr']
            units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '']
            columns = []
            for i in range(len(data)):
                columns.append(Column(data[i],name=names[i],unit=units[i]))
            self.tab2_simple = Table(columns)
            self.tab2_simple.write(outfile, format='fits', overwrite=True)


    def draw_phot_results(self):
        ''' DRAW RESULTING FIT ON R-BAND CUTOUT, for gui '''
        markcolor='cyan'
        objlist=[]
        markwidth=1.5
        for sma in self.apertures_a:
            obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,sma, sma*(1-self.eps), rot_deg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
            objlist.append(obj)
            #print(self.xcenter,self.ycenter,sma, sma*(1-self.eps), self.theta, np.degrees(self.theta))
        self.markhltag = self.image_frame.canvas.add(self.image_frame.dc.CompoundObject(*objlist))
        self.image_frame.fitsimage.redraw()
        
    def draw_phot_results_mpl(self):
        ''' draw results in matplotlib figure '''
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm, origin='lower')

        apertures = []
        for sma in self.apertures_a:
            apertures.append(EllipticalAperture((self.xcenter,self.ycenter),sma, sma*(1-self.eps), theta = self.theta))
            
        for aperture in apertures:
            aperture.plot(color='white',lw=1.5)
        plt.title(os.path.basename(self.image_name).split('.fits')[0])
        plt.savefig(self.image_name.split('.fits')[0]+'_phot_apertures.png')
        plt.close()
    def plot_profiles(self):
        ''' enclosed flux and surface brightness profiles, save figure '''
        #plt.close("all")        
        plt.figure(figsize=(10,4))
        plt.subplots_adjust(wspace=.3)
        plt.subplot(2,2,1)
        #plt.plot(self.apertures_a,self.flux1,'bo')
        plt.errorbar(self.apertures_a,self.flux1,self.flux1_err,fmt='b.')
        plt.title('R-band')
        #plt.xlabel('semi-major axis (pixels)')
        plt.ylabel('Enclosed flux')
        plt.gca().set_yscale('log')
        if self.image2_flag:
            plt.subplot(2,2,2)
            plt.errorbar(self.apertures_a,self.flux2,self.flux2_err,fmt='b.')
            #plt.xlabel('semi-major axis (pixels)')
            plt.ylabel('Enclosed flux')
            plt.title('H-alpha')
            plt.gca().set_yscale('log')
        # plot surface brightness vs radius
        plt.subplot(2,2,3)
        #plt.plot(self.apertures_a,self.flux1,'bo')
        plt.errorbar(self.apertures_a,self.sb1,self.sb1_err,fmt='b.')
        plt.xlabel('semi-major axis (pixels)')
        plt.ylabel('Surface Brightess')
        plt.gca().set_yscale('log')
        if self.image2_flag:
            plt.subplot(2,2,4)
            plt.errorbar(self.apertures_a,self.sb2,self.sb2_err,fmt='b.')
            plt.xlabel('semi-major axis (pixels)')
            plt.ylabel('Surface Brightness')
            plt.gca().set_yscale('log')
        #plt.show()
        plt.savefig(self.image_name.split('.fits')[0]+'-enclosed-flux.png')
        plt.close()
    def plot_fancy_profiles(self, logx=False):
        # plot enclosed flux        
        fig = plt.figure(figsize=(10,4))
        plt.subplots_adjust(left=.15,bottom=.15,right=.95,top=.9,wspace=.3)

        labels = ['R','Halphax100']
        alphas = [1,.4,.6,.4]
        x = self.apertures_a*self.pixel_scale
        if self.image2_flag:
            fluxes = [self.flux1_erg,self.flux2_erg]
            flux_errs = [self.flux1_err_erg,self.flux2_err_erg]
            sbfluxes = [self.sb1_erg_sqarcsec,self.sb2_erg_sqarcsec]
            sbflux_errs = [self.sb1_erg_sqarcsec_err,self.sb2_erg_sqarcsec_err]
            
        else:
            fluxes = [self.flux1_erg]#,self.flux2_erg]
            flux_errs = [self.flux1_err_erg]#,self.flux2_err_erg]
            sbfluxes = [self.sb1_erg_sqarcsec]
            sbflux_errs = [self.sb1_erg_sqarcsec_err]

        plt.subplot(1,2,1)
        plotflag = self.snr_per_pixel > 2
        x = x[plotflag]
        for i,t in enumerate(fluxes):
            y0 = fluxes[i][plotflag]
            y1 = y0+flux_errs[i][plotflag]
            y2 = y0-flux_errs[i][plotflag]

            if (i == 1) + (i == 3):
                y0=y0*100
                y1 = y1*100
                y2 = y2*100
            #print(f"DEBUG: fancy plots flux, len(x)={len(x)}, len(y1)={len(y1)}, len(y2)={len(y2)}")
            plt.fill_between(x,y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(x,y0,'-',lw=2,color=mycolors[i])
            # and points
            plt.plot(x,y0, marker='o',color=mycolors[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        plt.title(os.path.basename(self.image_name).replace('.fits',''))
        plt.gca().set_yscale('log')
        if logx:
            plt.gca().set_xscale('log')
        plt.legend(loc='lower right')

        plt.subplot(1,2,2)

        for i,t in enumerate(sbfluxes):
            y0 = sbfluxes[i][plotflag]          
            y1 = y0+sbflux_errs[i][plotflag]
            y2 = y0-sbflux_errs[i][plotflag]

            if (i == 1) + (i == 3):
                y0=y0*100
                y1 = y1*100
                y2 = y2*100
            #print(f"DEBUG: fancy plots sb, len(x)={len(x)}, len(y1)={len(y1)}, len(y2)={len(y2)}")                
            plt.fill_between(x,y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(x,y0,'-',lw=2,color=mycolors[i])
            plt.plot(x,y0,marker='o',color=mycolors[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Surface Brightness',fontsize=16)

        plt.gca().set_yscale('log')
        if logx:
            plt.gca().set_xscale('log')        
        plt.legend(loc='upper right')
            
        plt.savefig(self.image_name.split('.fits')[0]+'_enclosed_flux_fancy.png')        
        plt.close(fig)


def run_ellipse_photometry(
    r_fits: str,
    cs_fits: str | None = None,
    mask_fits: str | None = None,
    image2_filter: str | None = None,
    filter_ratio: float | None = None,
    psf: str | None = None,
    psf_ha: str | None = None,
    objra: float | None = None,
    objdec: float | None = None,
    fixcenter: bool = False,
    run_statmorph: bool = False,
    write_prefix: str | None = None,
):
    """
    Headless elliptical photometry runner.
    Returns the ellipse instance so caller can inspect results.
    Writes phot tables via existing methods.
    """
    e = EllipsePhotometry(
        r_fits,
        image2=cs_fits,
        mask=mask_fits,
        use_mpl=False,
        image2_filter=image2_filter,
        filter_ratio=filter_ratio,
        psf=psf,
        psf_ha=psf_ha,
        objra=objra,
        objdec=objdec,
        fixcenter=fixcenter,
        
    )

    e.run_for_gui(runStatmorphFlag=run_statmorph)

    # keep it simple: use your existing batch call
    #if run_statmorph:
    #    e.run_for_gui(runStatmorphFlag=True)   # (name is legacy but works)
    #else:
    #    e.run_two_image_phot()

    # If you want: allow prefixing output files without changing internal logic
    # (optional; you already have write_phot_fits_tables(prefix=...))
    if write_prefix is not None:
        e.write_phot_fits_tables(prefix=write_prefix)

    return e
