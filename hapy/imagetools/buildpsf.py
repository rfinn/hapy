'''
GOAL:
- create a psf image from an R-band mosaic


PROCEDURE
- use sextractor to identify objects
- use CLASS_STAR and flags to identify stars that are not saturated
- create a table of x,y coordinates (astropy.table.Table)
- extract stars - photutils.psf.extract_stars
- build psf - photutils.EPSFBuilder
- save psf image (e.g. for input into galfit)

sounds easy enough, right?


INPUT:
- image
- saturation level

OUTPUT:
- psf image

REFERENCES:
- following instructions found here:
https://photutils.readthedocs.io/en/stable/epsf.html#build-epsf


'''
import os
import sys
import numpy as np

try:
    from photutils import EPSFBuilder
except ImportError:
    from photutils.psf import EPSFBuilder
from photutils.psf import extract_stars
from photutils.psf import IntegratedGaussianPRF

try:
    from photutils import centroid_com, centroid_1dg, centroid_2dg
except ImportError:
    from photutils.centroids import centroid_com, centroid_1dg, centroid_2dg
from astropy.nddata import NDData
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.modeling import models, fitting
from astropy.stats import gaussian_sigma_to_fwhm

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

class psf_parent_image():
    def __init__(self, image=None, max_good=None, size=25, se_config = 'default.sex.HDI', sepath=None,nstars=100, pixelscale=0.43,oversampling=None,saturate=None):
        self.image_name = image
        self.basename = os.path.basename(self.image_name).split('.fits')[0]
        self.psf_image_name = self.basename+'-psf.fits'        
        self.data, self.header = fits.getdata(self.image_name, header=True)
        self.sat_level = max_good
        self.size = size
        self.config = se_config
        self.sextractor_files =['default.sex.HDI','default.param','default.conv','default.nnw','default.sex.INT','default.sex.BOK']
        self.saturate = saturate
        # path to source extractor files
        # default is probably fine if you have github in ~/github
        if sepath == None:
            self.sepath=os.getenv('HOME')+'/github/halphagui/astromatic/'
        else:
            self.sepath = sepath
        # number of stars to use to determine the psf
        self.nstars = nstars
        # default pixelscale is set for HDI camera

        # get pixelscale from image header
        try:
            self.pixelscale = np.abs(float(self.header['PIXSCAL1'])) # for mosaic data
        except KeyError:
            try:
                self.pixelscale = np.abs(float(self.header['CD1_1']))*3600
            except KeyError:
                self.pixelscale = pixelscale
        if oversampling == None:
            self.oversampling = 2
        else:
            self.oversampling = oversampling

        # make a directory for saving plots
        if not os.path.exists('plots'):
            os.mkdir('plots')
    def link_files(self):
        # these are the sextractor files that we need
        # set up symbolic links from sextractor directory to the current working directory        
        for file in self.sextractor_files:
            os.system('ln -s '+self.sepath+'/'+file+' .')
            
    def clean_links(self):
        # clean up symbolic links to sextractor files
        # sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
        for file in self.sextractor_files:
            os.system('unlink '+file)
    def run_all(self):
        self.runse()
        self.read_se_table()
        print("running find_stars...")
        self.find_stars() 
        print("running extract_stars...")       
        self.extract_stars() 
        print("running show_stars...")
        self.show_stars()
        print("running build_psf...")       
        self.build_psf()
        self.show_psf()
        self.measure_fwhm()
        self.save_psf_image()

    def runse(self):
        self.link_files()
        if self.saturate is not None:
            s = 'sex {} -c {}  -SATUR_LEVEL {}'.format(self.image_name,self.config,self.saturate)
        else:
            s = 'sex %s -c %s  -SATUR_LEVEL 40000.0'%(self.image_name,self.config)
        #print(s)
        os.system(s)
        ###################################
        # Read in Source Extractor catalog
        ###################################

        secat = fits.getdata('test_cat.fits',2)
        ###################################
        # get median fwhm of image
        ###################################
        fwhm = np.median(secat['FWHM_IMAGE'])*self.pixelscale
        self.se_fwhm_arcsec = fwhm
        #############################################################
        # rerun Source Extractor catalog with updated SEEING_FWHM
        #############################################################
        if self.saturate is not None:
            s = 'sex {} -c {}  -SATUR_LEVEL {} -SEEING_FWHM {:.1f}'.format(self.image_name,self.config,self.saturate,fwhm)
            print('running: ',s)
        else:
            s = 'sex %s -c %s  -SATUR_LEVEL 40000.0 -SEEING_FWHM %s'%(self.image_name,self.config,str(fwhm))

        os.system(s)

        # should update the size to be a multiple of the fwhm

    def read_se_table(self):
        self.secat = fits.getdata('test_cat.fits',2)
        pass
    def find_stars(self):
        # select objects with CLASS_STAR > 0.98
        # and FLUX_MAX 
        x = self.secat['X_IMAGE']
        y = self.secat['Y_IMAGE']

        # select based on star/class separator and no flags
        flag1 = (self.secat['CLASS_STAR'] > 0.95) & (self.secat['FLAGS'] == 0)

        # remove stars that are near the edge
        # being conservative by using the 10x full image size as the buffer rather than half image size
        # in part to compensate from zeros that sometimes are seen around perimeter after swarp
        flag2 = ((x > 5.*self.size) & (x < (self.data.shape[1] -1 - 5.*self.size)) &
                (y > 5.*self.size) & (y < (self.data.shape[0] -1 - 5.*self.size)))

        # remove stars with close neighbors
        c = SkyCoord(self.secat['ALPHA_J2000']*u.deg,self.secat['DELTA_J2000']*u.deg, frame='icrs')
        # returns index of closest match (not itself), separation in deg,
        # and 3d sep (not sure what this means if I don't provide redshift
        idx, sep2d, dist3d = c.match_to_catalog_sky(c,nthneighbor=2)
        # make sure stars don't have a neighbor within 15" (picked that fairly randomly...)
        flag3 = sep2d > 20./3600.*u.deg
        star_flag = flag1 & flag2 & flag3
        # select 25 objects with FLUX_MAX closest to median value
        # in other words, select 25 most central objects
        fm = self.secat['FLUX_MAX'][star_flag]
        x = x[star_flag]
        y = y[star_flag]
        sorted_indices = fm.argsort()
        print('number of potential stars', sum(flag1), sum(flag2), np.sum(star_flag))
        # select stars within +/- nstar/2 from a threshold
        # where threshold is the percentile ranking according to peak flux
        threshold = .30
        lower = int(threshold*len(sorted_indices)) - int((self.nstars)/2)
        upper = int(threshold*len(sorted_indices)) + int((self.nstars)/2)
        #lower = int(.25*len(sorted_indices)) 
        #upper = int(.85*len(sorted_indices))
        if lower < 0: # try increasing threshold
            threshold = .50
            lower = int(threshold*len(sorted_indices)) - int((self.nstars)/2)
            upper = int(threshold*len(sorted_indices)) + int((self.nstars)/2)
            if lower < 0:
                lower = 0
        print('number of psf stars = ',upper-lower+1)
        self.xstar = x[lower:upper]
        self.ystar = y[lower:upper]
        #print(self.xstar.shape)

        self.stars_tbl = Table()
        self.stars_tbl['x'] = self.xstar
        self.stars_tbl['y'] = self.ystar

    def extract_stars(self):
        # just in case sky was not properly subtracted
        mean_val, median_val, std_val = sigma_clipped_stats(self.data, sigma=2.)

        # working on #118
        # checking to see if more stars are retained if I don't subtract median
        # RESULT: subtracting median didn't impact result
        
        self.data -= median_val
        nddata = NDData(data=self.data)  
        self.stars = extract_stars(nddata, self.stars_tbl, size=self.size)  
        # check to make sure stars don't have zeros
        # Virgo 2017 pointing-4_R.coadd.fits is giving me trouble b/c a lot of stars have zeros
        keepflag = np.ones(len(self.stars),'bool')
        
        for i,s in enumerate(self.stars):
            # some pixels have values == 0 but they don't appear to be bad?  or are they masked?
            # increasing cut from 1 to 10 to compensate for the case where maybe 1-2 pixels are bad
            # still loosing a lot of stars so bumping up cut to 100
            if len(np.where(s.data == 0)[0]) > 10:
                keepflag[i] = False

        print(f"number of stars to keep after ==0 cut = {np.sum(keepflag)}/{len(keepflag)}")
        self.keepflag2 = keepflag
        self.stars = extract_stars(nddata, self.stars_tbl[keepflag], size=self.size)  
        #self.stars = self.stars[keepflag]
        if len(self.stars) < self.nstars:
            self.nstars = len(self.stars)
    def show_stars(self):
        # round up to integer
        nrows = int(np.ceil(np.sqrt(len(self.stars))))
        ncols = int(np.ceil(np.sqrt(len(self.stars))))
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10), squeeze=True)
        print(f"\nMaking figure with {nrows} rows and {ncols} cols with {len(self.stars)} stars...\n")
        ax = ax.ravel()
        for i in range(self.nstars):
            norm = simple_norm(self.stars[i], 'log', percent=99.)
            ax[i].imshow(self.stars[i], norm=norm, origin='lower', cmap='viridis')
        #plt.show()
        plt.savefig('plots/'+self.basename+'-allstars.png')
        plt.close()
    def build_psf(self):
        mysigclip = SigmaClip(sigma=3, cenfunc='median')
        self.oversampling=2
        if self.oversampling == None:
            epsf_builder = EPSFBuilder(maxiters=12, progress_bar=False, smoothing_kernel='quadratic', recentering_func = centroid_com)#,flux_residual_sigclip=SigmaClip)  
        else:
            epsf_builder = EPSFBuilder(oversampling=self.oversampling, maxiters=20, progress_bar=False,  recentering_func = centroid_com, smoothing_kernel='quadratic')#,sigma_clip=mysigclip)  
        self.epsf, self.fitted_stars = epsf_builder(self.stars)
    def show_psf(self):
        norm = simple_norm(self.epsf.data, 'log', percent=99.)
        plt.figure()
        plt.imshow(self.epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        #plt.show()
        plt.savefig('plots/'+self.basename+'-psf.png')
        plt.close()
    def measure_fwhm(self):
        nx,ny = self.epsf.data.shape
        x,y = np.mgrid[:nx,:ny]
        # gaussian model with initial guess for center and std
        # assume the PSF is in the center of the image
        # assume a std of 4
        p_init = models.Gaussian2D(x_mean=nx/2.,y_mean=ny/2., x_stddev=3.,y_stddev=3.)
        fit_p = fitting.LevMarLSQFitter()
        t = fit_p(p_init,x,y,self.epsf.data)
        self.x_stddev = t.x_stddev.value
        self.y_stddev = t.y_stddev.value
        # compute total radial std
        # correct for oversampling of PSF
        self.std = np.sqrt(self.x_stddev**2 + self.y_stddev**2)/self.oversampling
        # the above estimate was coming out very large.  trying this instead
        # this next estimate gives an answer that is much closer to what you get from imexam moffat
        self.std = np.mean([self.x_stddev,self.y_stddev])/self.oversampling
        # convert gaussian std to fwhm
        self.fwhm = self.std*gaussian_sigma_to_fwhm
        self.fwhm_arcsec = self.fwhm*self.pixelscale
        print('image fwhm = %.2f pix (%.2f arcsec)'%(self.fwhm, self.fwhm*self.pixelscale))
    def save_psf_image(self):
        # save the psf file
        # outfile = append '-psf' to the input image name
        # just use the basename (filename), and drop any path information
        # this will have the effect of saving the psf image in the current directory,
        # rather than the directory where the image is


        print('PSF image name = ',self.psf_image_name)
        fits.writeto(self.psf_image_name,self.epsf.data, overwrite=True)

        # update image header
        # I know this is a kludge, but not sure how to access the default header that it writes
        # BEFORE it gets written

        data, header = fits.getdata(self.psf_image_name, header=True)

        header.append(card=('FWHM', float('{:.2f}'.format(self.fwhm)), 'PSF fwhm in pixels'))
        header.append(card=('SEFWHM', float('{:.2f}'.format(self.se_fwhm_arcsec)), 'PSF fwhm in arcsec from SE'))        
        header.append(card=('STD', float('{:.2f}'.format(self.std)), 'PSF STD in pixels'))
        header.append(card=('OVERSAMP', self.oversampling, 'PSF oversampling'))
        fits.writeto(self.psf_image_name, data, header=header, overwrite=True)
    def update_image_header(self):
        '''  add fwhm and psfimage name to header of parent image '''
        self.header.append(card=('FWHM', float('{:.2f}'.format(self.fwhm)), 'PSF fwhm in pixels'))
        self.header.append(card=('SEFWHM', float('{:.2f}'.format(self.se_fwhm_arcsec)), 'PSF fwhm in arcsec from SE'))                
        self.header.append(card=('STD', float('{:.2f}'.format(self.std)), 'PSF STD in pixels'))
        self.header.append(card=('PSFSTD', float('{:.2f}'.format(self.std)), 'PSF STD in pixels'))        
        self.header.append(card=('OVERSAMP', self.oversampling, 'PSF oversampling'))
        fits.writeto(self.image_name, self.data, header=self.header, overwrite=True)
        
if __name__ == '__main__':
    import argparse

    #parser = argparse.ArgumentParser(description ='create psf image from image that contains stars')
    parser = argparse.ArgumentParser(
        description='GOAL:\n- create a psf image from an R-band mosaic \n\nPROCEDURE \n- use sextractor to identify objects (run to estimate FWHM, then input FWHM on second run) \n- use CLASS_STAR and flags to identify stars that are not saturated \n- create a table of x,y coordinates (astropy.table.Table) \n- extract stars - photutils.psf.extract_stars \n- build psf - photutils.EPSFBuilder \n- save psf image (e.g. for input into galfit)\n\nsounds easy enough, right? \n\nINPUT: \n- image \n- saturation level \n\nOUTPUT: \n- psf image \n\nREFERENCES: \n- following instructions found here: \nhttps://photutils.readthedocs.io/en/stable/epsf.html#build-epsf',
        formatter_class=argparse.RawTextHelpFormatter)
   
    parser.add_argument('--image',dest = 'image', help='input image')
    parser.add_argument('--saturate',default=None,dest = 'saturate', help='saturation limit')
    parser.add_argument('--int',default=False,dest='int',action = 'store_true', help='set this for INT data')
    parser.add_argument('--bok',default=False,dest='bok',action = 'store_true', help='set this for BOK data')        
     
    args = parser.parse_args()

    
    #image = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/MKW8_R.coadd.fits'
    #image = '/Users/rfinn/research/VirgoFilaments/Halpha/virgo-coadds-2017/pointing-4_R.coadd.fits'

    basename = os.path.basename(args.image).split('.fits')[0]
    psf_image_name = basename+'-psf.fits'
    print()
    print("PSF image name = ",psf_image_name)
    print()
    if os.path.exists(psf_image_name):
        print("found psf image for ",args.image)
        sys.exit()
    
    if args.int:
        p = psf_parent_image(image=args.image, size=39, nstars=100, oversampling=2,saturate=args.saturate,se_config='default.sex.INT')
    elif args.bok:
        p = psf_parent_image(image=args.image, size=39, nstars=100, oversampling=2,saturate=args.saturate,se_config='default.sex.BOK')
    # bok config file is causing trouble for some images - but why???
    else:
        p = psf_parent_image(image=args.image, size=25, nstars=100, oversampling=2,saturate=args.saturate)
    p.runse()
    p.read_se_table()
    p.find_stars()
    p.extract_stars()
    p.show_stars()
    p.build_psf()
    p.show_psf()
    p.measure_fwhm()
    p.save_psf_image()
    p.update_image_header()
    '''
    getting some weird noise at the bottom of the psf image

    need to read more about how it handles the edge

    not sure that this will impact results significantly, so moving on...

    '''
    
