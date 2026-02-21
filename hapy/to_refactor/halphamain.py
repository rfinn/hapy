#!/usr/bin/env python 
"""
USAGE:

Virgo 2019 INT data, running on laptop

cd github/halphagui

source venv/bin/activate

cd /data-pool/Halpha/halphagui-output-20230626

%run ~/github/halphagui/testing/halphamain.py --virgo --rimage /home/rfinn/data/reduced/virgo-coadds-feb2019-int/VF-118.1817+20.9822-INT-20190205-p001-r-shifted.fits --haimage /home/rfinn/data/reduced/virgo-coadds-feb2019-int/VF-118.1817+20.9822-INT-20190205-p001-Halpha.fits --filter inthalpha --psfdir /home/rfinn/data/reduced/psf-images/ --tabledir /home/rfinn/research/Virgo/tables-north/v1/ --auto


Testing after MVC semi-implementation.  I created a testing directory on the linux laptop, and this is command I used to run the gui:

(venv) (base) rfinn@virgof:~/research/Virgo-dev/halphagui-test$ 

python ~/github/halphagui/halphamain.py --virgo --tabledir ~/research/Virgo/tables-north/v2/ --rimage VF-145.781+31.887-HDI-20180313-p019-R.fits --haimage VF-145.781+31.887-HDI-20180313-p019-ha4.fits --csimage VF-145.781+31.887-HDI-20180313-p019-ha4-CS-ZP.fits --psfdir ~/research/Virgo-dev/halphagui-test/ --filter ha4 --prefix VF-145.781+31.887-HDI-20180313-p019

2025-December

(virgo) rfinn@s64247 HDI % 

python ~/github/halphagui/halphamain.py --rimage UAT-177.865+21.004-HDI-20150418-NRGb161-h01-R.fits --haimage UAT-177.865+21.004-HDI-20150419-NRGb161-h01-ha12.fits --filter 12 --uat --prefix NRGb161-h01

"""

# TODONE - cutout directory has a trailing dash when running in uat mode
# TODO - keep testing maskwrapper - behaving oddly when user adds objects

import sys, os
sys.path.append(os.getcwd())
#sys.path.append(os.getenv('HOME')+'github/HalphaImaging/')
sys.path.append(os.getenv('HOME')+'/github/HalphaImaging/python3/')

import numpy as np
import platform

#from PyQt5 import  QtWidgets
#from PyQt5 import QtCore
#from PyQt5.Qtcore import  Qt
from PyQt5 import QtCore,QtWidgets, QtGui
#from ginga.qtw.QtHelp import QtGui #, QtCore
from halphav5 import Ui_MainWindow
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
#from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
from ginga import colors
from ginga.canvas.CanvasObject import get_canvas_types
from ginga.misc import log
from ginga.util.loader import load_data

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, FK5
import astropy.units as u
from astropy import nddata
from astropy.table import Table, Column
from astropy.visualization import simple_norm
from astropy.cosmology import WMAP9 as cosmo
from astropy.time import Time

# packages for ellipse fitting routine
# https://photutils.readthedocs.io/en/stable/isophote.html
#from photutils.isophote import EllipseGeometry
#from photutils.isophote import Ellipse


#from photutils import EllipticalAperture



import matplotlib
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from datetime import date
import time
# routines for measuring elliptical photometry
from photwrapper import ellipse

from maskwrapper import maskwindow

from galfitwrapper import galfitwindow

from buildpsf import psf_parent_image
from halphaCommon import cutout_image

from fit_profile import profile, dualprofile, rprofile, haprofile, ratio_error
# code from HalphaImaging repository
import filterratio as runse

# code for calculating redshift cutoffs of filter
# and for calculating the transmission correction for each galaxy
# based on where it falls ini the filter bandpass
from filter_transmission import filter_trace

from join_catalogs import join_cats, make_new_cats

import imutils

#from uat_mask import mask_image
# filter information

# INT filters from https://astro.ing.iac.es/filter/list.php?instrument=WFC
# lmin is center - 0.5 FWHM
# lmin is center + 0.5 FWHM
#
# Filter # 197 is INT197 = (6568, 95)
# Filter # 227 is INT227 = (6657, 80)
#
# BOK using the NOAO filter, so Halpha 4
# 
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.,'INT197':6520.5, 'INT227':6617}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.,'INT197':6615.5, 'INT227':6697}

# Force a specific toolkit on mac
macos_ver = platform.mac_ver()[0]
try:
    matplotlib.use('Qt5Agg')
except ImportError:
    print("WARNING! could not load Qt5Agg")
    
import matplotlib.pyplot as plt

### code to measure galaxy sizes from photutils segmentation image
from get_galaxy_size import getobjectsize

# default size for cutouts, multiple of NSA PETROTH90
cutout_scale = 14

# now in terms of R25 for
cutout_scale = 2.5

mask_scalefactor = 1 # number to multiple R24 by
######################################################
## FUNCTIONS
######################################################

def get_params_from_name(image_name):
    t = os.path.basename(image_name).split('-')
    #print(t)
    if len(t) == 5:
        telescope = t[2]
        dateobs = t[3]
        pointing = t[4]
    elif len(t) == 6: # meant to catch negative declinations
        telescope = t[3]
        dateobs = t[4]
        pointing = t[5]
    #else:
    #    print("ruh roh - trouble with get_params_from_name for image ",image_name, len(t))
        #print(image_name)
        #print(t)
    return telescope,dateobs,pointing

def get_params_from_name_uat(image_name):
    '''
    coadd names should be as follows.
    if float(dec) < 0:
        outfile = f'UAT-{ra:07.3f}-{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'        
    else:
        outfile = f'UAT-{ra:07.3f}+{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'    

    each pointing will have an internal '-', like ABELL1367-h01
    '''
    t = os.path.basename(image_name).split('-')
    #print(t)
    if len(t) == 7: # meant to catch negative declinations
        telescope = t[2]
        dateobs = t[3]
        pointing = t[4]+'_'+t[5]
    elif len(t) == 8:
        telescope = t[3]
        dateobs = t[4]
        pointing = t[5]+'_'+t[6]
        
    else:
        print("ruh roh - trouble getting info from ",image_name, len(t))
        print(image_name)
        print(t)
        return
    return telescope,dateobs,pointing

class psfimage():
    def __init__(self):
        fwhm = 5.6
        fwhm_arcsec = 2.0




        
    
class uco_table():
    '''
    table for collecting positions of objects that are not in NSA or AGC catalogs
    '''
    def initialize_uco_arrays(self):
        # columns: id, x, y, ra, dec
        user = os.getenv('USER')
        today = date.today()
        str_date_today = today.strftime('%Y-%b-%d')
        self.uco_output_table = 'halpha-uco-data-'+user+'-'+str_date_today+'.fits'
        if os.path.exists(self.uco_output_table):
            self.uco_table = Table(fits.getdata(self.uco_output_table))
            self.uco_prefix = self.uco_table['PREFIX'].tolist()
            self.uco_id = self.uco_table['ID'].tolist()
            self.uco_ra = self.uco_table['RA'].tolist()
            self.uco_dec = self.uco_table['DEC'].tolist()
            self.uco_x = self.uco_table['X'].tolist()
            self.uco_y = self.uco_table['Y'].tolist()

        ## if not, create table
        else:
            self.uco_prefix = []
            self.uco_id = []
            self.uco_ra = []
            self.uco_dec = []
            self.uco_x = []
            self.uco_y = []
    def create_uco_table(self):
        c0 = Column(self.uco_prefix, name='PREFIX', description='Prefix for coadd image')
        c1 = Column(np.array(self.uco_id), name='ID',dtype=np.int32, description='ID')
        c2 = Column(np.array(self.uco_ra), name='RA',dtype='f', unit=u.deg)
        c3 = Column(np.array(self.uco_dec), name='DEC',dtype='f', unit=u.deg)
        c4 = Column(np.array(self.uco_x), name='X',dtype='f', unit=u.pixel)
        c5 = Column(np.array(self.uco_y), name='Y',dtype='f', unit=u.pixel)
        self.uco_table = Table([c0,c1,c2,c3,c4,c5])
    def write_uco_table(self):
        self.create_uco_table()

        self.uco_table.write(self.uco_output_table, format='fits',overwrite=True)

class hagui_methods():
    """ class to handle Halpha model, in Model-View-Controller design"""


    def initialize_output_arrays(self): # MVC - model
        ngal = len(self.galid)

        # galfit output
        self.gal_xc = np.zeros((ngal,2),'f')
        self.gal_xc = np.zeros((ngal,2),'f')
        self.gal_mag = np.zeros((ngal,2),'f')
        self.gal_n = np.zeros((ngal,2),'f')
        self.gal_re = np.zeros((ngal,2),'f')
        self.gal_PA = np.zeros((ngal,2),'f')
        self.gal_BA = np.zeros((ngal,2),'f')
        self.gal_sky = np.zeros((ngal,2),'f')
        
    def link_files(self): # MVC- model
        # these are the sextractor files that we need
        # set up symbolic links from sextractor directory to the current working directory
        sextractor_files=['default.sex.HDI','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            if not os.path.exists(file):
                os.system('ln -s '+self.sepath+'/'+file+' .')
    def clean_links(self): # MVC - model
        # clean up symbolic links to sextractor files
        # sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
        sextractor_files=['default.sex.HDI','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('unlink '+file)

        # remove catalog
        
        
    def build_psf(self): # MVC - model
        # check to see if R-band PSF images exist
        coadd_header = fits.getheader(self.rcoadd_fname)


        basename = os.path.basename(self.rcoadd_fname)
        psf_image_name = basename.split('.fits')[0]+'-psf.fits'
        if psf_image_name.find('-shifted') > -1:
            psf_image_name = psf_image_name.replace('-shifted','')
        basename = os.path.basename(self.hacoadd_fname)
        psf_image_name_ha = basename.split('.fits')[0]+'-psf.fits'
        psf_image_name = os.path.join(self.psfdirectory,psf_image_name)
        psf_image_name_ha = os.path.join(self.psfdirectory,psf_image_name_ha).replace('-CS','')
        if self.verbose:
            print('\nPSF NAME = ',psf_image_name,'\n')
        if os.path.exists(psf_image_name):
            # get fwhm from the image header
            # and oversampling
            print("LOADING EXISTING PSF IMAGE")
            header = fits.getheader(psf_image_name)
            self.psf_data = fits.getdata(psf_image_name)

            # add attributes
            self.psf = psfimage()
            
            self.psf.fwhm = header['FWHM'] # in pixels
            self.psf.fwhm_arcsec = self.psf.fwhm*self.pixelscale
            self.oversampling = float(header['OVERSAMP'])
            # if psf is in another directory, create a link to the current directory
            # this will avoid having a long filename b/c galfit does not handle long filenames

            ##
            # GALFIT does not handle the long image name - it craps out
            #
            # that must be why I copied the psf image to r-psf.fits
            # need to figure out another way so that I can run them in parallel
            ##
            #self.psf_image_name = psf_image_name

            ##
            # changing to use the psf without copying to current directory
            ##
            outname = os.path.basename(psf_image_name)
            command = f'ln -s {psf_image_name} {outname}'
            #print('running: ',command)
            os.system(command)
            #self.psf_image_name = 'r-psf.fits'
            self.psf_image_name = outname


            
        else:
            if self.verbose:
                print('oversampling = ',self.oversampling)
            print('PSF RESULTS FOR R-BAND COADDED IMAGE')
            self.psf = psf_parent_image(image=self.rcoadd_fname, size=21, nstars=100, oversampling=self.oversampling)
            self.psf.run_all()
            self.psf_image_name = self.psf.psf_image_name
            
        print('\nPSF NAME = ',psf_image_name_ha,'\n')
        if os.path.exists(psf_image_name_ha):
            # get fwhm from the image header
            # and oversampling
            print("LOADING EXISTING HALPHA PSF IMAGE")            
            header = fits.getheader(psf_image_name_ha)
            self.hapsf = psfimage()
            self.hapsf.fwhm = header['FWHM'] # in pixels
            self.hapsf.fwhm_arcsec = self.hapsf.fwhm*self.pixelscale

            # if psf is in another directory, create a link to the current directory
            # this will avoid having a long filename b/c galfit does not handle long filenames

            #self.psf_haimage_name = psf_image_name
            #command = 'cp {} ha-psf.fits'.format(psf_image_name)
            #print('running: ',command)
            #os.system(command)
            #self.psf_haimage_name = 'ha-psf.fits'

            outname = os.path.basename(psf_image_name_ha)
            command = f'ln -s {psf_image_name_ha} {outname}'
            #print('running: ',command)
            os.system(command)
            self.psf_haimage_name = outname
            
        else:
            print('PSF RESULTS FOR HA COADDED IMAGE')
            self.hapsf = psf_parent_image(image=self.hacoadd_fname, size=21, nstars=100, oversampling=self.oversampling)
            self.hapsf.run_all()
            self.psf_haimage_name = self.hapsf.psf_image_name

    def add_psf_to_table(self): # MVC - model
        fields = ['R_FWHM','H_FWHM']
        values = [self.psf.fwhm_arcsec,self.hapsf.fwhm_arcsec]
        for i,f in enumerate(fields):
            for j in range(len(self.table)):
                self.table[f][j]=values[i]
        pass
    def run_galfit(self, ncomp=1, asym=0, ha=0): # MVC - model?
        if self.psf_image_name is None:
            print('WARNING: psf could not be found')
            print('Please run build_psf')
            return
        if self.psf is not None:
            self.add_psf_to_table()
        self.ncomp = ncomp
        self.asym=asym
        self.galha = ha
        print('running galfit with ',ncomp,' components')
        #self.gwindow = QtWidgets.QWidget()
        '''
        if self.testing:
            self.ncomp = ncomp
            self.galfit = galfitwindow(self.gwindow, self.logger, image = 'MKW8-18037-R.fits', mask_image = 'MKW8-18037-R-mask.fits', psf='MKW8_R.coadd-psf.fits', psf_oversampling=2, ncomp=ncomp)
        else:
            self.galfit = galfitwindow(self.gwindow, self.logger, image = self.cutout_name_r, mask_image = self.mask_image_name, psf=self.psf.psf_image_name, psf_oversampling = self.oversampling, ncomp=ncomp)
        self.galfit.model_saved.connect(self.galfit_save)        
        self.galfit.setupUi(self.gwindow)

        self.gwindow.show()
        '''
        try:
            if ha:
                self.galimage = self.cutout_name_ha
                # setup psfimage
                psf = self.psf_haimage_name
                psf_oversampling = self.oversampling
            else:
                self.galimage = self.cutout_name_r
                # setup psfimage
                psf = self.psf_image_name
                psf_oversampling = self.oversampling
        except AttributeError:
            print('make sure you selected a galaxy')
            return
        try:
            if not self.auto:
                self.gwindow = QtWidgets.QWidget()
            #self.gwindow.aboutToQuit.connect(self.galfit_closed)

            
            if (ncomp == 1) & (asym == 0):
                #self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=self.nsa.rmag[self.igal], BA = self.nsa.cat.SERSIC_BA[self.igal], PA=self.nsa.cat.SERSIC_PHI[self.igal],nsersic=self.nsa.cat.SERSIC_N[self.igal], convolution_size=80)
                print('GALFIT psf image = ',psf)
                if self.auto:
                    self.galfit = galfitwindow(None, None, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, rad=self.gradius[self.igal],mag=14, BA = .8, PA=0,nsersic=2, convolution_size=80,auto=self.auto)
                else:
                    self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=14, BA = .8, PA=0,nsersic=2, convolution_size=80,auto=self.auto)                
            elif (ncomp == 1) & (asym == 1):
                #self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=self.nsa.rmag[self.igal], BA = self.nsa.cat.SERSIC_BA[self.igal], PA=self.nsa.cat.SERSIC_PHI[self.igal],nsersic=self.nsa.cat.SERSIC_N[self.igal], convolution_size=80,asym=1)
                self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=14, BA = .8, PA=5,nsersic=2, convolution_size=80,asym=1)                
            elif ncomp == 2:
                # use results from 1 component fit as input
                try:
                    mag = self.table['GAL_MAG'][self.igal]
                    re = self.table['GAL_RE'][self.igal]
                    BA = self.table['GAL_BA'][self.igal]
                    PA = self.table['GAL_BA'][self.igal]
                
                except KeyError:
                    print('WARNING!!!!')
                    print('trouble reading galfit results from data table')
                    print('make sure you run 1 component fit first')
                    return
                ########################
                # assume bulge contains 20% of light for initial guess
                ########################
                mag_disk = mag+.25
                mag_bulge = mag + 1.75
        
                ########################
                # require n=1 for disk, n=4 for bulge (allow bulge to vary)
                # also start PA=0 and BA=1 for bulge
                ########################
                nsersic_disk=1
                nsersic_bulge=4
        
                ########################
                # set re=1.5*re_initial for disk
                # set re = 0.5*re_initial for bulge
                ########################
                re1=1.2*re
                re2=.5*re
                    
                self.galfit = galfitwindow(self.gwindow, self.logger, image = galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, rad=re1, mag=mag_disk, BA=BA, PA=PA,nsersic=nsersic_disk,nsersic2=nsersic_bulge,mag2=mag_bulge, rad2=re2, fitn=False, fitn2=True, convolution_size=80)

            if not self.auto:
                self.galfit.model_saved.connect(self.galfit_save)        
                self.galfit.setupUi(self.gwindow)
                self.gwindow.show()
            else:
                self.galfit_save(None)
        except ValueError:
            print('WARNING - ERROR RUNNING GALFIT!!!')
            print('Make sure you have measured the PSF and made a mask!')
            print("error:", sys.exc_info()[0])
        #except:
        #    print("Unexpected error:", sys.exc_info()[0])
        #    raise
    def galfit_save(self,msg): # MVC - model?
        #print('galfit model saved!!!',msg)
        if self.testing:
            self.ncomp = int(msg)
            print('ncomp = ',self.ncomp)
        if self.galha:
            prefix = 'GAL_H'
        else:
            prefix = 'GAL_'
        if (self.ncomp == 1) & (self.asym == 0):
            if self.galha:
                self.galfit_hresults = self.galfit.galfit_results
            else:
                self.galfit_results = self.galfit.galfit_results
            fields = ['XC','YC','MAG','RE','N','BA','PA']
            values = np.array(self.galfit.galfit_results[:-2])[:,0].tolist()
            for i,f in enumerate(fields):
                colname = prefix+f
                if i == 3: # multiply radius by pixel scale
                    self.table[colname][self.igal]=values[i]*self.pixelscale
                else:
                    self.table[colname][self.igal]=values[i]
            fields = ['XC','YC','MAG','RE','N','BA','PA']
            print("testing, galfit_results = ", self.galfit.galfit_results[:-2])
            values = np.array(self.galfit.galfit_results[:-2])[:,1].tolist()
            for i,f in enumerate(fields):
                colname = prefix+f+'_ERR'
                if i == 3: # convert radius from pixels to arcsec
                    self.table[colname][self.igal]=values[i]*self.pixelscale
                else:
                    self.table[colname][self.igal]=values[i]
            fields = ['SKY','CHISQ']
            values = [self.galfit.galfit_results[-2],self.galfit.galfit_results[-1]]
            for i,f in enumerate(fields):
                colname = prefix+f
                self.table[colname][self.igal]=values[i]
            wcs = WCS(self.galimage)
            ra,dec = wcs.wcs_pix2world(self.galfit.galfit_results[0][0],self.galfit.galfit_results[1][0],0)
            self.table[prefix+'RA'][self.igal]=ra
            self.table[prefix+'DEC'][self.igal]=dec
            #self.update_gui_table()

            # save results to output table
            #self.table['GAL_SERSIC'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,0]
            #self.table['GAL_SERSIC_ERR'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,1]
            #self.table['GAL_SERSIC_SKY'][self.igal] = (self.galfit.galfit_results[-2])
            #self.table['GAL_SERSIC_CHISQ'][self.igal] = (self.galfit.galfit_results[-1])
            #print(self.table[self.igal])
        elif (self.ncomp == 1) & (self.asym == 1):
            if self.galha:
                self.galfit_haasym_results = self.galfit.galfit_results
            else:
                self.galfit_asym_results = self.galfit.galfit_results
            #print('writing results for galfit w/asymmetry')
            self.table[prefix+'SERSASYM'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,0]
            self.table[prefix+'SERSASYM_ERR'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,1]
            self.table[prefix+'SERSASYM_ERROR'][self.igal] = (self.galfit.galfit_results[-2])
            self.table[prefix+'SERSASYM_CHISQ'][self.igal] = (self.galfit.galfit_results[-1])
            wcs = WCS(self.galimage)
            ra,dec = wcs.wcs_pix2world(self.galfit.galfit_results[0][0],self.galfit.galfit_results[1][0],0)
            self.table[prefix+'SERSASYM_RA'][self.igal]=ra
            self.table[prefix+'SERSASYM_DEC'][self.igal]=dec
            self.update_gui_table()
        elif self.ncomp == 2:
            self.galfit_results2 = self.galfit.galfit_results
            #print(self.galfit_results2)
            self.table['GAL_2SERSIC'][self.igal] = np.array(self.galfit_results2[:-2])[:,0]
            self.table['GAL_2SERSIC_ERR'][self.igal] = np.array(self.galfit_results2[:-2])[:,1]
            self.table['GAL_2SERSIC_ERROR'][self.igal] = (self.galfit_results2[-2])
            self.table['GAL_2SERSIC_CHISQ'][self.igal] = (self.galfit_results2[-1])
            #print(self.table[self.igal])
        if not self.auto:
            self.update_gui_table()
        

    def galfit_ellip_phot(self): # MVC - model
        '''
        use galfit ellipse parameters as input for photutils elliptical photometry

        '''
        ### CLEAR R-BAND CUTOUT CANVAS
        #self.rcutout.canvas.delete_all_objects()

        ### MAKE SURE GALFIT 1 COMP MODEL WAS RUN

        try:
            xc,yc,mag,re,n,BA,pa = np.array(self.galfit_results[:-3])[:,0].tolist()
            #print('GALFIT PA = ',xc,yc,mag,re,n,BA,pa )
        except AttributeError:
            print('Warning - galfit 1 comp fit results not found!')
            print('Make sure you run galfit, then try again.')
            return
        
        ### FIT ELLIPSE
        #
        if self.auto:
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = None,image2_filter=self.hafilter, filter_ratio=self.filter_ratio,psf=self.psf_image_name,psf_ha=self.psf_haimage_name)
        else:
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = self.rcutout,image2_filter=self.hafilter, filter_ratio=self.filter_ratio,psf=self.psf_image_name,psf_ha=self.psf_haimage_name)
        #fields = ['XC','YC','MAG','RE','N','BA','PA']

        # TRANSFORM THETA
        # GALFIT DEFINES THETA RELATIVE TO Y AXIS
        # PHOT UTILS DEFINES THETA RELATIVE TO X AXIS
        THETA = pa + 90 # in degrees
        #print('THETA = ',THETA)
        self.e.run_with_galfit_ellipse(xc,yc,BA=BA,THETA=THETA)
        self.e.plot_profiles()
        #os.chdir(current_dir)

        '''
        fields = ['ASYM','ASYM_ERR','ASYM2','ASYM2_ERR']
        values = [self.e.asym, self.e.asym_err, self.e.asym2,self.e.asym2_err]
        for i,f in enumerate(fields):
            colname = 'GAL_'+f
            self.table[colname][self.igal]=values[i]
        self.update_gui_table()
        '''

        # fit profiles
        self.fit_profiles(prefix='GAL')
        # save results
        self.write_profile_fits(prefix='GAL_')

        # TODO - this is calling a view function - should not do that.  need to separate model and view, or call both from controller
        if not self.auto:
            self.draw_ellipse_results(color='cyan')

    def photutils_ellip_phot(self):
        #current_dir = os.getcwd()
        #image_dir = os.path.dirname(self.rcoadd_fname)
        #os.chdir(image_dir)

        ### CLEAR R-BAND CUTOUT CANVAS
        #self.rcutout.canvas.delete_all_objects()

        ### FIT ELLIPSE
        #
        if self.verbose:
            print()
            print("running photutils_ellip_phot")
            print()
        if self.auto:
            ra = self.objparams[0]
            dec = self.objparams[1]            
            
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = None,image2_filter='16', filter_ratio=self.filter_ratio, psf=self.psf_image_name,psf_ha=self.psf_haimage_name,objra=ra,objdec=dec)
        else:
            ra = self.defcat.cat['RA'][self.igal]
            dec = self.defcat.cat['DEC'][self.igal]            
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = self.rcutout,image2_filter='16', filter_ratio=self.filter_ratio, psf=self.psf_image_name,psf_ha=self.psf_haimage_name,objra=ra,objdec=dec )
        self.e.run_for_gui(runStatmorphFlag=True)

        if self.verbose:
            print("calling plot_profiles\n")
        self.e.plot_profiles()
        if self.verbose:
            print("finished plot_profiles\n")
        
        #os.chdir(current_dir)
        if self.verbose:
            print("saving data\n")

        ### SAVE DATA TO TABLE
        fields = ['BADGAL','XCENTROID','YCENTROID','EPS','THETA','GINI','HGINI',\
                  'M20','HM20',
                  'UNMASKED_AREA','TOTAL_AREA',\
                  'SUM','SUM_MAG','ASYM','ASYM_ERR',\
                  'HSUM','HSUM_MAG','HASYM','HASYM_ERR']#,'SUM_ERR']
        values = [self.bad_galaxy,self.e.xcenter, self.e.ycenter,self.e.eps, np.degrees(self.e.theta), \
                  self.e.cat.gini[self.e.objectIndex],self.e.cat2.gini[self.e.objectIndex],\
                  self.e.M20_1,self.e.M20_2,\
                  self.e.cat[self.e.objectIndex].area.value*self.pixelscale*self.pixelscale,\
                  #self.e.cat[self.e.objectIndex].segment_area.value*self.pixelscale*self.pixelscale,\
                  self.e.masked_pixel_area*self.pixelscale*self.pixelscale,\
                  self.e.source_sum_erg, self.e.source_sum_mag,self.e.asym, self.e.asym_err, \
                  self.e.source_sum2_erg,self.e.source_sum2_mag,self.e.asym2,self.e.asym2_err]
        for i,f in enumerate(fields):
            if i == 0:
                colname = f
            else:
                colname = 'ELLIP_'+f
            #print(colname)
            #self.table[colname][self.igal]=float('%.2e'%(values[i]))            
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
        # update sky noise
        fields = ['R_SKYNOISE','H_SKYNOISE']
        values = [self.e.im1_skynoise/1.e-17,self.e.im2_skynoise/1.e-17]
        #print("before writing skynoise: ",values)
        for i,f in enumerate(fields):
            #print(values[i])
            self.table[f][self.igal] = values[i]
        fields = ['R_SKY','H_SKY']
        values = [self.e.sky,self.e.sky2]
        #print("before writing sky values: ",values)
        for i,f in enumerate(fields):
            #print(values[i])
            self.table[f][self.igal] = values[i]

        
        # what are we doing here?
        wcs = WCS(self.cutout_name_r)
        ra,dec = wcs.wcs_pix2world(self.e.xcenter,self.e.ycenter,0)
        
        ra2,dec2 = wcs.wcs_pix2world(self.e.cat2.xcentroid[self.e.objectIndex],self.e.cat2.ycentroid[self.e.objectIndex],0)        
        self.table['ELLIP_RA'][self.igal]=ra
        self.table['ELLIP_DEC'][self.igal]=dec
        self.table['ELLIP_HRA'][self.igal]=ra2
        self.table['ELLIP_HDEC'][self.igal]=dec2
        self.write_fits_table()
        # TODONE - write out phot table
        colnames = ['area',
                    'background_mean',
                    'bbox_xmax',
                    'bbox_xmin',
                    'bbox_ymax',
                    'bbox_ymin',
                    'cxx',
                    'cxy',
                    'cyy',
                    'eccentricity',
                    'ellipticity',
                    'elongation',
                    'equivalent_radius',
                    'fwhm',
                    'gini',
                    'inertia_tensor',
                    'kron_flux',
                    'kron_fluxerr',
                    'kron_radius',
                    'local_background',
                    'moments', 
                    'moments_central',
                    'orientation',
                    'perimeter',
                    'segment_flux',
                    'segment_fluxerr',
                    'semimajor_sigma',
                    'semiminor_sigma',
                    'xcentroid',
                    'ycentroid']


        # this is hanging when trying to calculate for cat2 - not sure why
        # skipping this for now.
        # if self.verbose:
        #     print("Calculating fluxfrac_radius 30 for cat2")
            
        # r30 = self.e.cat2.fluxfrac_radius(0.3)*self.pixelscale*u.arcsec/u.pixel
        # if self.verbose:
        #     print("Calculating fluxfrac_radius 50 for cat2")
        
        # r50 = self.e.cat2.fluxfrac_radius(0.5)*self.pixelscale*u.arcsec/u.pixel
        # if self.verbose:
        #     print("Calculating fluxfrac_radius 90 for cat2")
        
        # r90 = self.e.cat2.fluxfrac_radius(0.9)*self.pixelscale*u.arcsec/u.pixel

        # if self.verbose:
        #     print()
        #     print("adding r30 to e.cat2")
        #     print()
        

        # self.e.cat2.add_extra_property('PHOT_R30',r30)

        # if self.verbose:
        #     print()
        #     print("adding r50 to e.cat2")
        #     print()
        
        # self.e.cat2.add_extra_property('PHOT_R50',r50)

        # if self.verbose:
        #     print()
        #     print("adding r90 to e.cat2")
        #     print()
        
        # self.e.cat2.add_extra_property('PHOT_R90',r90)

        if self.verbose:
            print("Calculating fluxfrac_radius")
        # calculate fractional radii, but these are circular, and in pixels
        r30 = self.e.cat.fluxfrac_radius(0.3)*self.pixelscale*u.arcsec/u.pixel
        r50 = self.e.cat.fluxfrac_radius(0.5)*self.pixelscale*u.arcsec/u.pixel
        r90 = self.e.cat.fluxfrac_radius(0.9)*self.pixelscale*u.arcsec/u.pixel

        if self.verbose:
            print()
            print("done calculating fluxfrac_radius")
        # calculate fractional radii, but these are circular, and in pixels

        if self.verbose:
            print()
            print("adding extra properties to e.cat")
            print()
        self.e.cat.add_extra_property('PHOT_R30',r30)
        self.e.cat.add_extra_property('PHOT_R50',r50)
        self.e.cat.add_extra_property('PHOT_R90',r90)



        if self.verbose:
            print("writing radii to main table")
        # write these out to the main table
        fields = ['R30','R50','R90',\
                  'HR30','HR50','HR90']
        
        #values = [self.e.cat.PHOT_R30[self.e.objectIndex].value,\
        #          self.e.cat.PHOT_R50[self.e.objectIndex].value,\
        #          self.e.cat.PHOT_R90[self.e.objectIndex].value,\
        #          self.e.cat2.PHOT_R30[self.e.objectIndex].value,\
        #          self.e.cat2.PHOT_R50[self.e.objectIndex].value,\
        #          self.e.cat2.PHOT_R90[self.e.objectIndex].value]
        values = [self.e.cat.PHOT_R30[self.e.objectIndex].value,\
                  self.e.cat.PHOT_R50[self.e.objectIndex].value,\
                  self.e.cat.PHOT_R90[self.e.objectIndex].value,\
                  0,0,0]
        for i,f in enumerate(fields):
            colname = 'ELLIP_'+f
            #print(colname,values[i])
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i].value))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except AttributeError:
                print("AttributeError: ",colname, values[i])
                print("sorry for the shit show...")
                self.table[colname][self.igal]=float('%.4e'%(values[i]))                
            #except:
            #    print("problem writing table element",colname,values[i])

        if self.verbose:
            print("writing fits table\n")
        self.write_fits_table()
        
        if self.e.statmorph_flag:
            print()
            print("running statmorph")
            self.write_statmorph()
            if self.verbose:
                print("writing fits table after running statmorph")
            self.write_fits_table()
        #c1 = Column(data=np.array(r30[self.e.objectIndex]),name='PHOTR30',unit='arcsec',description='photutils fluxfrac_radius')
        #c2 = Column(data=np.array(r50[self.e.objectIndex]),name='PHOTR50',unit='arcsec',description='photutils fluxfrac_radius')
        #c3 = Column(data=r90[self.e.objectIndex],name='PHOTR90',unit='arcsec',description='photutils fluxfrac_radius')
        #qtable.add_columns([c1,c2,c3])

        if self.verbose:
            print("setting up photutil_tab.fits")
        qtable = self.e.cat[self.e.objectIndex].to_table(colnames)
        
        phot_table_name = self.cutout_name_r.replace('.fits','-photutil_tab.fits')
        qtable = Table(qtable)
        qtable.write(phot_table_name,format='fits',overwrite=True)
        
       
        if not self.auto:
            self.update_gui_table()

        # convert theta to degrees, and subtract 90 to get angle relative to y axis
        #self.e.theta = np.degrees(self.e.theta) - 90
        # fit profiles

        if args.verbose:
            print("in photutils_ellip_phot, fitting profiles")
        self.fit_profiles()
        # save results
        self.write_profile_fits()

        
        if not self.auto:
            self.draw_ellipse_results(color='magenta')
    
    def write_statmorph(self):
        #########################################################
        ## ADD STATMORPH PARAMETERS
        #########################################################

        
        # write these out to the main table
        fields = ['XCENTROID','YCENTROID',\
                  'RPETRO_CIRC','RPETRO_ELLIP','RHALF_ELLIP',\
                  'R20','R80',\
                  'GINI','M20','F_GM20','S_GM20',\
                  'C','A','S','FLAG']
        
        values = [self.e.morph.xc_centroid,\
                  self.e.morph.yc_centroid,\
                  self.e.morph.rpetro_circ*self.pixelscale,\
                  self.e.morph.rpetro_ellip*self.pixelscale,\
                  self.e.morph.rhalf_ellip*self.pixelscale,\
                  self.e.morph.r20*self.pixelscale,\
                  self.e.morph.r80*self.pixelscale,\
                  self.e.morph.gini,\
                  self.e.morph.m20,\
                  self.e.morph.gini_m20_bulge,\
                  self.e.morph.gini_m20_merger,\
                  self.e.morph.concentration,\
                  self.e.morph.asymmetry,\
                  self.e.morph.smoothness,\
                  self.e.morph.flag]
                  
        for i,f in enumerate(fields):
            colname = 'SMORPH_'+f
            #print(colname)
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()

        if self.e.statmorph_flag2:
            ## Add Halpha values
            values = [self.e.morph2.xc_centroid,\
                      self.e.morph2.yc_centroid,\
                      self.e.morph2.rpetro_circ*self.pixelscale,\
                      self.e.morph2.rpetro_ellip*self.pixelscale,\
                      self.e.morph2.rhalf_ellip*self.pixelscale,\
                      self.e.morph2.r20*self.pixelscale,\
                      self.e.morph2.r80*self.pixelscale,\
                      self.e.morph2.gini,\
                      self.e.morph2.m20,\
                      self.e.morph2.gini_m20_bulge,\
                      self.e.morph2.gini_m20_merger,\
                      self.e.morph2.concentration,
                      self.e.morph2.asymmetry,
                      self.e.morph2.smoothness,\
                      self.e.morph2.flag]

            for i,f in enumerate(fields):
                colname = 'SMORPH_H'+f
                #print(colname)
                try:
                    self.table[colname][self.igal]=float('%.4e'%(values[i]))
                except KeyError:
                    print("KeyError: ",colname)
                    print("\ntable column names: \n",self.table.colnames)
                    sys.exit()
                except TypeError:
                    print("TypeError: ",colname, values[i])
                    print("sorry for the shit show...")
                    print("\ntable column names: \n",self.table.colnames)
                    sys.exit()

        
    def fit_profiles(self,prefix=None):
        #current_dir = os.getcwd()
        #image_dir = os.path.dirname(self.rcoadd_fname)
        #os.chdir(image_dir)
        if prefix is None:
            rphot_table = self.cutout_name_r.split('.fits')[0]+'_phot.fits'
            haphot_table = self.cutout_name_ha.split('.fits')[0]+'_phot.fits'
        else:
            rphot_table = self.cutout_name_r.split('.fits')[0]+'-'+prefix+'_phot.fits'
            haphot_table = self.cutout_name_ha.split('.fits')[0]+'-'+prefix+'_phot.fits'

        self.rfit = rprofile(self.cutout_name_r, rphot_table, label='R')
        self.rfit.becky_measurements()
        self.hafit = haprofile(self.cutout_name_ha, haphot_table, label=r"$H\alpha$")
        #if self.hafit.fit_flag:
        self.hafit.becky_measurements()
        self.hafit.get_r24_stuff(self.rfit.iso_radii[self.rfit.isophotes == 24.][0][0])

        # TODO - is this mixing model and view????  ARGGGHHHHHH
        both = dualprofile(self.rfit,self.hafit)
        try:
            both.make_3panel_plot()
        except:
            print()            
            print("problem making 3panel plot - weird values...")
            print()
    def write_profile_fits(self,prefix=None): # MVC - model
        fields = ['R24','R25','R26','R24V','R25V',\
                  'R_F25','R_F50','R_F75',\
                  'M24','M25','M26',\
                  'F_30R24','F_R24','C30',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        d = self.rfit
        values = [d.iso_radii[0],d.iso_radii[1],d.iso_radii[2],d.iso_radii[3],d.iso_radii[4],\
                  d.flux_radii[0],d.flux_radii[1],d.flux_radii[2],\
                  d.iso_mag[0],d.iso_mag[1],d.iso_mag[2],\
                  d.flux_30r24,d.flux_r24,d.c30,\
                  d.petrorad,d.petroflux_erg,d.petror50,d.petror90,d.petrocon,d.petromag
                  ]
        for i,f in enumerate(fields):
            if prefix is None:
                colname = f
            else:
                colname = prefix+f
            #print(colname, values[i])
            self.table[colname][self.igal]=float('%.4e'%(values[i][0]))
            self.table[colname+'_ERR'][self.igal]=float('%.4e'%(values[i][1]))
            
        fields = ['R16','R17','R_F25','R_F50','R_F75','M16','M17','F_30R24','F_R24','C30','R_F95R24','F_TOT',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        d = self.hafit
        values = [d.iso_radii[0],d.iso_radii[1],\
                  d.flux_radii[0],d.flux_radii[1],d.flux_radii[2],\
                  d.iso_mag[0],d.iso_mag[1],\
                  d.flux_30r24,d.flux_r24,d.c30,d.flux_95r24, d.total_flux,\
                  d.petrorad,d.petroflux_erg,d.petror50,d.petror90,d.petrocon,d.petromag
                  ]
        for i,f in enumerate(fields):
            if prefix is None:
                colname = 'H'+f
            else:
                colname = prefix+'H'+f

            #print(colname,values[i])
            self.table[colname][self.igal]=float('{:.4e}'.format(values[i][0]))
            self.table[colname+'_ERR'][self.igal]=float('{:.4e}'.format(values[i][1]))

        # SFR conversion from Kennicutt and Evans (2012)
        # log (dM/dt/Msun/yr) = log(Lx) - logCx
        logCx = 41.27
        #print(len(self.hafit.total_flux),len(self.gzdist))
        L = self.hafit.total_flux*(4.*np.pi*cosmo.luminosity_distance(self.gzdist[self.igal]).cgs.value**2)
        #print(L)
        detect_flag = L > 0
        self.sfr = np.zeros(len(L),'d')
        self.sfr[detect_flag] = np.log10(L[detect_flag]) - logCx
        if prefix is None:
            colname='LOG_SFR_HA'
        else:
            colname=prefix+'LOG_SFR_HA'
        #print('sfr = ',self.sfr)
        #print(self.sfr[0], self.sfr[1])
        self.table[colname][self.igal]=float('%.4e'%(self.sfr[0]))
        self.table[colname+'_ERR'][self.igal]=float('%.4e'%(self.sfr[1]))
        self.table[colname+'_FLAG'][self.igal]=detect_flag[0]
        # inner ssfr
        a = self.hafit.flux_30r24
        b = self.rfit.flux_30r24
        self.inner_ssfr = a[0]/b[0]
        self.inner_ssfr_err = ratio_error(a[0],b[0],a[1],b[1])
        if prefix is None:
            colname='SSFR_IN'
        else:
            colname = prefix+'SSFR_IN'
        self.table[colname][self.igal]=float('%.4e'%(self.inner_ssfr))
        self.table[colname+'_ERR'][self.igal]=float('%.4e'%(self.inner_ssfr_err))
        # outer ssfr
        c = self.hafit.flux_r24
        d = self.rfit.flux_r24
        self.outer_ssfr = (c[0] - a[0])/(d[0] - b[0])
        self.outer_ssfr_err = ratio_error(c[0] - a[0],d[0] - b[0],np.sqrt(a[1]**2 + c[1]**2),np.sqrt(b[1]**2 + d[1]**2))
        if prefix is None:
            colname='SSFR_OUT'
        else:
            colname=prefix+'SSFR_OUT'
        self.table[colname][self.igal]=float('%.4e'%(self.outer_ssfr))
        self.table[colname+'_ERR'][self.igal]=float('%.4e'%(self.outer_ssfr_err))
        self.write_fits_table()        
        if not self.auto:
            self.update_gui_table()

class hagui_interactive():
    """ 
    class to handle all gui setup and visualizations 

    I think this needs to inherit the model class, because buttons are connected to model functions

    Or the functions that are connecting the buttons should be called from the controller class

    """

    def setup_gui(self): # view
        #print(MainWindow)

        self.ui = Ui_MainWindow()        
        self.ui.setupUi(MainWindow)
        #self.ui.setGeometry(0,0,400,300)
        #self.ui.setFont(QtGui.Qfont('Arial',10))
        self.logger = logger
        self.drawcolors = colors.get_colors()
        self.dc = get_canvas_types()
        self.add_coadd_frame(self.ui.leftLayout)
        self.add_cutout_frames()
        self.connect_setup_menu()
        self.connect_ha_menu()
        self.connect_halpha_type_menu()
        self.connect_comment_menu()
        #self.connect_buttons()
        #self.add_image(self.ui.gridLayout_2)
        #self.add_image(self.ui.gridLayout_2)
        self.connect_buttons()
    def connect_buttons(self): #view
        self.ui.wmark.clicked.connect(self.find_galaxies)
        #self.ui.editMaskButton.clicked.connect(self.edit_mask)
        self.ui.makeMaskButton.clicked.connect(self.make_mask)
        #self.ui.saveCutoutsButton.clicked.connect(self.write_cutouts)
        self.ui.fitEllipseGalfitButton.clicked.connect(self.galfit_ellip_phot)
        self.ui.fitEllipsePhotutilsButton.clicked.connect(self.photutils_ellip_phot)
        self.ui.wfratio.clicked.connect(self.get_filter_ratio)
        self.ui.resetRatioButton.clicked.connect(self.reset_cutout_ratio)
        self.ui.resetSizeButton.clicked.connect(self.reset_cutout_size)
        self.ui.prefixLineEdit.textChanged.connect(self.set_prefix)
        #self.ui.fitEllipseButton.clicked.connect(self.fit_ellipse_phot)
        self.ui.galfitButton.clicked.connect(lambda: self.run_galfit(ncomp=1))
        self.ui.galfitAsymButton.clicked.connect(lambda: self.run_galfit(ncomp=1,asym=1))
        self.ui.galfitHaButton.clicked.connect(lambda: self.run_galfit(ncomp=1,ha=1))
        self.ui.galfitHaAsymButton.clicked.connect(lambda: self.run_galfit(ncomp=1,asym=1,ha=1))
        self.ui.galfit2Button.clicked.connect(lambda: self.run_galfit(ncomp=2))
        self.ui.psfButton.clicked.connect(self.build_psf)
        self.ui.saveButton.clicked.connect(self.write_fits_table)
        self.ui.clearCutoutsButton.clicked.connect(self.clear_cutouts)
        self.ui.filterRatioLineEdit.returnPressed.connect(self.set_filter_ratio)
        self.ui.cutoutSizeLineEdit.returnPressed.connect(self.set_cutout_size)
    def add_coadd_frame(self,panel_name): # view
        logger = log.get_logger("example1", log_stderr=True, level=40)
        self.coadd = image_panel(panel_name, self.ui,logger)
        self.coadd.key_pressed.connect(self.key_press_func)
        #self.coadd.add_cutouts()

    def add_cutout_framesv4(self):# with halphav4 # view
        # r-band cutout
        self.rcutout_label = QtWidgets.QLabel('r-band')
        self.ui.cutoutsLayout.addWidget(self.rcutout_label, 0, 0, 1, 1)
        a = QtWidgets.QLabel('CS Halpha')
        self.ui.cutoutsLayout.addWidget(a, 0, 1, 1, 1)
        a = QtWidgets.QLabel('Mask')
        self.ui.cutoutsLayout.addWidget(a, 0, 2, 1, 1)

        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.rcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, 8, 1,autocut_params='histogram')
        self.hacutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 1, 8, 1,autocut_params='stddev')
        self.maskcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,1, 2, 8, 1,autocut_params='stddev')
    def add_cutout_frames(self): # view
        # r-band cutout
        self.rcutout_label = QtWidgets.QLabel('r-band')
        drow = 20
        self.ui.cutoutsLayout.addWidget(self.rcutout_label, 0, 0, 1, 2)
        temp_label = QtWidgets.QLabel('')
        self.ui.cutoutsLayout.addWidget(temp_label, 0, 2, 1, 2)
        self.nsa_label = QtWidgets.QLabel('NSA ID')
        self.ui.cutoutsLayout.addWidget(self.nsa_label, 0, 4, 1, 2)
        temp_label = QtWidgets.QLabel('')
        self.ui.cutoutsLayout.addWidget(temp_label, 0, 6, 1, 2)
        self.agc_label = QtWidgets.QLabel('AGC Number')
        self.ui.cutoutsLayout.addWidget(self.agc_label, 0, 8, 1, 2)
        self.rcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, drow, 10,autocut_params='stddev')
        a = QtWidgets.QLabel('CS Halpha')
        self.ui.cutoutsLayout.addWidget(a, drow+1, 0, 1, 2)
        self.hacutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, drow+2, 0, drow, 10,autocut_params='stddev')

        a = QtWidgets.QLabel('Mask')
        self.ui.cutoutsLayout.addWidget(a, 2*drow+2, 0, 1, 2)
        self.maskcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,2*drow+3, 0, drow, 10)
        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)



    def clear_cutouts(self): # view - should this be in the controller class?
        self.rcutout.canvas.delete_all_objects()
        self.hacutout.canvas.delete_all_objects()
    def connect_setup_menu(self): # view or controller - these are calling model functions/quantities
        self.ui.actionR_coadd.triggered.connect(self.get_rcoadd_file)
        self.ui.actionHa_coadd_2.triggered.connect(self.get_hacoadd_file)
        self.ui.actionNSA_catalog_path.triggered.connect(self.getnsafile)
        self.ui.actionAGC_catalog_path.triggered.connect(self.getagcfile)
        
    def connect_ha_menu(self): # view
        #print('working on this')
        #extractAction.triggered.connect(self.close_application)

        self.ui.actionhalpha4.triggered.connect(lambda: self.set_hafilter('4'))
        self.ui.actionhalpha8.triggered.connect(lambda: self.set_hafilter('8'))
        self.ui.actionhalpha12.triggered.connect(lambda: self.set_hafilter('12'))
        self.ui.actionhalpha16.triggered.connect(lambda: self.set_hafilter('16'))
        self.ui.actioninthalpha.triggered.connect(lambda: self.set_hafilter('inthalpha'))
        self.ui.actionintha6657.triggered.connect(lambda: self.set_hafilter('intha6657'))
        self.ui.actionsienaha.triggered.connect(lambda: self.set_hafilter('sienaha'))
        
    def connect_halpha_type_menu(self): # view
        ha_types = ['Ha Emission','No Ha']
        for name in ha_types:
            self.ui.haTypeComboBox.addItem(str(name))
        self.ui.haTypeComboBox.activated.connect(self.set_halpha_type)
    def connect_comment_menu(self): # view
        comment_types = ['Cont Sub Prob','merger/tidal','scat light','asym R', 'asym Ha','fore. star', 'fore. gal','edge-on','part cov','nuc ha']
        for name in comment_types:
            self.ui.commentComboBox.addItem(str(name))
        self.ui.commentComboBox.activated.connect(self.set_comment)
    def set_prefix(self,prefix): # MVC - probably model? - no this is from the gui
        self.prefix = prefix
        #print('prefix for output files = ',self.prefix)
    def set_prefix_on_gui(self,prefix): # MVC - probably model? - no this is from the gui
        """  use this if the prefix is provided by the user - this will fill in the box with the provided prefix """
        self.ui.prefixLineEdit.setText(prefix)        
        
    def mark_galaxies(self): # MVC - view or controller, b/c this relies on model quantities
        #
        # using code in TVMark.py as a guide for adding shapes to canvas
        #
        #

        objlist = []
        markcolor='cyan'
        markwidth=1
        #size = cutout_scale*self.gradius
        size = self.cutout_sizes
        #size[size > self.global_max_cutout_size] = self.global_max_cutout_size
        #size[size < self.global_min_cutout_size] = self.global_min_cutout_size
        #print(f"DEBUGGING: in mark_galaxies: len(gximage) = {len(self.gximage)}")
        for i,x in enumerate(self.gximage):
            #print(f"{i}, {self.galid[i]}, cutout_size = {size[i]:.2f}, gximage[i]={self.gximage[i]:.1f}, agcximage={self.agcximage[i]:.1f}") 
            obj = self.coadd.dc.Box(
                x=x, y=self.gyimage[i], \
                xradius=size[i]/2,yradius=size[i]/2, \
                #xradius=100,yradius=100, \
                color=markcolor, linewidth=markwidth)
            glabel = self.coadd.dc.Text(x-50,self.gyimage[i]+60,\
                                        str(self.galid[i]), color=markcolor)
            objlist.append(obj)
            objlist.append(glabel)
        if self.agcflag:
            for i,x in enumerate(self.agcximage):
                #print(x,self.agcyimage[i],self.agc.cat.AGCNUMBER[i])
                obj = self.coadd.dc.Box(
                    x=x, y=self.agcyimage[i], xradius=75,\
                    yradius=75, color='purple', linewidth=markwidth)
                glabel = self.coadd.dc.Text(x-40,self.agcyimage[i]+40,\
                                        str(self.agc.cat['AGCnr'][i]), color='purple')
                objlist.append(obj)
                objlist.append(glabel)
            
        self.markhltag = self.coadd.canvas.add(self.coadd.dc.CompoundObject(*objlist))
        self.coadd.fitsimage.redraw()
    def setup_ratio_slider(self): # MVC - view
        self.ui.ratioSlider.setRange(0,100)
        self.ui.ratioSlider.setValue(50)
        self.ui.ratioSlider.setSingleStep(1)
        #self.ui.ratioSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        #self.ui.ratioSlider.setFocusPolicy(QtCore.StrongFocus)
        self.ui.ratioSlider.valueChanged.connect(self.ratio_slider_changed)
    def setup_cutout_slider(self): # MVC - view
        self.ui.cutoutSlider.setRange(0,100)
        self.ui.cutoutSlider.setValue(50)
        self.ui.cutoutSlider.setSingleStep(1)
        #self.ui.ratioSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        #self.ui.ratioSlider.setFocusPolicy(QtCore.StrongFocus)
        self.ui.cutoutSlider.valueChanged.connect(self.cutout_slider_changed)
            
    def clear_comment_field(self): # MVC - view?
        self.ui.commentLineEdit.clear()
        
    def display_cutouts(self): # MVC - view

        # TODONE  - split into separate model and view functions
        #################################################
        # this part needs to move to the view class
        #################################################
        # only do this when running gui (as opposed to automatically)
        self.ui.cutoutSizeLineEdit.setText(str(self.cutout_size))
        self.reset_cutout_size()
        self.reset_cutout_ratio()
        
        ############################################################
        # TODONE: this is the only part that should be in view
        #cutoutHa = Cutout2D(self.ha.data, position, self.size, wcs=self.coadd_wcs, mode = 'trim')
        ((ymin,ymax),(xmin,xmax)) = self.cutoutR.bbox_original
        bbox = '[{:d}:{:d},{:d}:{:d}]'.format(int(xmin),int(xmax),int(ymin),int(ymax))        
        #print(ymin,ymax,xmin,xmax)
        
        ############################################################                

        self.rcutout.load_image(self.r[ymin:ymax,xmin:xmax])
        self.hacutout.load_image(self.halpha_cs[ymin:ymax,xmin:xmax])
        #cutoutR.plot_on_original(color='white')
        
        # ###################################################################################
        # TODO - this should not be called in the model b/c model does not interact with view
        # ###################################################################################
        self.update_gui_table_cell(self.igal,'BBOX',str(bbox))
        self.write_fits_table()
    def draw_ellipse_results(self, color='cyan'): # MVC - view
        # mark r24
        markcolor=color#, 'yellow', 'cyan']
        markwidth=1
        #print('inside draw_ellipse_results')
        image_frames = [self.rcutout, self.hacutout]
        radii = self.rfit.iso_radii[:,0][0:2]
        objlist = []
        for i,im in enumerate(image_frames):
            for r in radii:
                #print('r = ',r)
                r = r/self.pixelscale
                #print('r = ',r)
                obj =im.dc.Ellipse(self.e.xcenter,self.e.ycenter,r, r*(1-self.e.eps), rot_deg = np.degrees(self.e.theta), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            if i == 1: # add R17 for Halpha image
                r = self.hafit.iso_radii[:,0][1]/self.pixelscale
                obj =im.dc.Ellipse(self.e.xcenter,self.e.ycenter,r, r*(1-self.e.eps), rot_deg = np.degrees(self.e.theta), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            self.markhltag = im.canvas.add(im.dc.CompoundObject(*objlist))
            im.fitsimage.redraw()
        # mark R17 in halpha image
        
    def display_mask(self, mask_image_name): # MVC - controller or view??? putting in view for now
        t = self.cutout_name_r.split('.fit')
        self.mask_image_name=t[0]+'-mask.fits'
        
        if not self.auto:
            if os.path.exists(self.mask_image_name):
                self.maskcutout.load_file(self.mask_image_name)
                self.mask_image_exists = True
                
            else:
                # clear mask frame
                self.maskcutout.fitsimage.clear()

        
        #self.mask_image = mask_image_name


        
class hacontroller():
    """ class to handle interactions with the user and send requests to model and view  """
    def get_rcoadd_file(self): # view or controller? controller b/c interact w/user
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            self.rcoadd_fname = fname[0]
            #print(self.rcoadd_fname)
            self.load_rcoadd()
            #self.le.setPixmap(QPixmap(fname))
    def get_hacoadd_file(self): # view or controller? controller
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            
            self.hacoadd_fname = fname[0]
            #print(self.hacoadd_fname)
            self.load_hacoadd()
            #self.le.setPixmap(QPixmap(fname))
    def getnsafile(self): # view or controller?
        """ get NSA file location from the user """
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            self.nsa_fname = fname[0]
            self.nsa = galaxy_catalog(self.nsa_fname,nsa=True)
        print('Got NSA catalog with {} lines!'.format(len(self.nsa.cat)))
        self.defcat = self.nsa
        #self.le.setPixmap(QPixmap(fname))
    def getagcfile(self): # controller
        """ get AGC file location from the user """        
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            self.agc_fname = fname[0]
            self.agc = galaxy_catalog(self.agc_fname,agc=True)
            self.agcflag = True
            self.defcat = self.agc
    def get_instrument(self):
        instruments = ['INT','BOK','HDI','MOS']
        for ii in instruments:
            if ii in self.rcoadd_fname:
                self.instrument = ii
                break
    def set_hafilter(self,filterid): # model or controller? I think this is a mix...
        #print('setting ha filter to ',filterid)

        self.get_instrument()
        
        # this is in controller
        self.hafilter = filterid
        print('halpha filter = ',self.hafilter)
        
        # this part should be in the model, so this should be a method in model
        # and then call model.set_filter() which would execute the following lines

        self.filter_trace = filter_trace(self.hafilter,instrument=self.instrument)
        self.zmin = self.filter_trace.minz_trans10
        self.zmax = self.filter_trace.maxz_trans10
        #self.get_zcut()
        
    def set_halpha_type(self,hatype): # controller
        self.halpha_type = hatype
        if args.verbose:
            print(f"in set_halpha_type, hatype={hatype}")
        #if int(hatype) == 0:
        #    self.haflag[self.igal]=True
        #    self.update_gui_table_cell(self.igal,'HA_FLAG',str(True))
        # why am I calling this again?
        try:
            if int(hatype) == 0:
                self.haflag[self.igal]=True
                self.update_gui_table_cell(self.igal,'HA_FLAG',str(True))
        except AttributeError:
            print('make sure you selected a galaxy')

        # add command here to write the fits table
        self.write_fits_table()
    def set_filter_ratio(self,ratio): # controller?
        try:
            self.filter_ratio = float(ratio)
            self.subtract_images()
            self.ui.filterRatioLineEdit.setText(str(ratio))
        except:
            print('did not understand.  try again')
    def set_cutout_size(self,size): # controller?
        try:
            self.cutout_size_arcsec = int(size)*u.arcsec
            self.cutout_size_arcsec = int(size)            
            self.display_cutouts()
            self.ui.cutoutSizeLineEdit.setText(str(size))            
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
            self.ui.cutoutSizeLineEdit.setText(str(self.cutout_size))            
    def key_press_func(self,key): # MVC - is this controller or view? controller
        """ define keys to that control behavior of image display """
        print(key)
        if key == 'r':
            z = self.coadd.fitsimage.settings.get_setting('zoomlevel')
            print('zoom = ',z)
            p = self.coadd.fitsimage.settings.get_setting('pan')
            print('pan = ',p)

            self.coadd.fitsimage.set_data(self.r)
            self.coadd.fitsimage.zoom_to(z.value)
            self.coadd.fitsimage.panset_xy(p.value[0],p.value[1])
            self.coadd.canvas.redraw()
        elif key == 'h':
            try:
                self.coadd.fitsimage.set_data(self.halpha_cs)
            except AttributeError:
                print('no continuum subtracted image yet - get filter ratio')
                self.coadd.fitsimage.set_data(self.ha)
        elif key == 'u': # unidentified object!
            x,y = self.coadd.fitsimage.get_last_data_xy()
            x = float(x)
            y = float(y)
            self.uco_x.append(x)
            self.uco_y.append(y)
            ra,dec = self.coadd_wcs.wcs_pix2world(x,y,0)
            self.uco_ra.append(ra)
            self.uco_dec.append(dec)
            if len(self.uco_id) == 0:
                self.uco_id.append(1)
            else:
                self.uco_id.append(np.max(self.uco_id)+1)
            self.uco_prefix.append(self.prefix)
            self.write_uco_table()
        elif key == 'down':
            '''
            up arrow will go to previous galaxy in the list
            '''
            #print(self.igal)
            if self.igal == (len(self.ra)-1):
                self.igal = 0
            else:
                self.igal += 1
                self.ui.wgalid.setCurrentIndex(self.igal)
                self.select_galaxy(self.igal)
            #print(self.igal)
        elif key == 'up':
            '''
            down arrow will go to previous galaxy in the list
            '''
            if self.igal == 0: 
                self.igal = len(self.ra)-1
            else:
                self.igal -= 1
                self.ui.wgalid.setCurrentIndex(self.igal)
                self.select_galaxy(self.igal)
    def ratio_slider_changed(self, value): # MVC - controller
        #print(self.minfilter_ratio, self.maxfilter_ratio, self.filter_ratio)
        delta = self.maxfilter_ratio - self.minfilter_ratio
        self.filter_ratio = self.minfilter_ratio + (delta)/100.*self.ui.ratioSlider.value()
        #print(value,' ratio slider changed to', round(self.filter_ratio,4))
        try:
            self.subtract_images(overwrite=True)
            self.display_cutouts()
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
    def cutout_slider_changed(self, value): # MVC - controller
        #print(self.minfilter_ratio, self.maxfilter_ratio, self.filter_ratio)
        delta = self.maxcutout_size - self.mincutout_size
        #self.cutout_size = self.mincutout_size + (delta)/100.*value
        self.cutout_size = self.cutout_sizes[self.igal]
        self.cutout_size_arcsec = self.cutout_sizes_arcsec[self.igal]        
        #print(value,' ratio slider changed to', round(self.filter_ratio,4))
        try:
            self.display_cutouts()
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
    def select_galaxy(self,id): # MVC - view or controller?
        print()
        print()
        print('selecting a galaxy')
        self.igal = self.ui.wgalid.currentIndex()
        if self.virgo:
            self.rcutout_label.setText('r-band '+str(self.defcat.cat['VFID'][self.igal]))
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]
            
            #print("compare lengths of catalogs ",len(self.defcat.cat),len(self.BA))
            print()
        elif self.uat:
            # removed this an we resolved the issue where we got the wrong galaxy in the viewer.
            #self.igal = self.igal-1 # why do we need this???
            self.rcutout_label.setText('r-band '+str(self.defcat.cat['AGCnr'][self.igal]))
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]
            print("new galaxy params = ",self.objparams)
            #print("compare lengths of catalogs ",len(self.defcat.cat),len(self.BA))
            print()
            
        else:
            self.rcutout_label.setText('r-band '+str(self.nsa2['NSAID'][self.igal]))
            self.objparams = None
        self.rcutout_label.show()
                                   
                                   
        print('active galaxy = ',self.igal)
        # when galaxy is selected from list, trigger
        # cutout imaages
        self.get_galaxy_cutout()

        # will need to also call new view functions
        self.display_cutouts()
        #self.display_mask()

        # the following two lines are the same as self.clear_cutouts
        # seems odd that this is called here - why are we clearing the cutouts after we displayed the cutouts???
        # need to test this when actually running gui, as opposed to running in auto mode
        self.clear_cutouts()
        #self.rcutout.canvas.delete_all_objects()
        #self.hacutout.canvas.delete_all_objects()            
        self.clear_comment_field()

        # TODO clear other flags, like no Halpha emission
        #################################################



    def reset_cutout_size(self): # MVC - controller
        self.cutout_size = self.reset_size
        self.update_images()
        #self.ui.cutoutSlider.setValue(50)
    def reset_cutout_ratio(self): # MVC - controller
        self.filter_ratio = self.reset_ratio
        self.update_images()
        #self.ui.ratioSlider.setValue(50)
        
    def update_images(self): # MVC - controller b/c calls model and view
        self.subtract_images(overwrite=True)
        #self.display_cutouts()
    def make_mask(self,objparams=None): # MVC - is this controller, or view?
        # TODO - break off view functions into a method within the haview class
        #current_dir = os.getcwd()
        #image_dir = os.path.dirname(self.rcoadd_fname)
        #os.chdir(image_dir)
        try:
            print()
            print("using ellipe parameters to unmask central region - woo hoo!")
            print()
            objparams = self.objparams
            #print("\t ",self.objparams)
        except AttributeError:
            print()
            print("problem getting objparams for masking routing")
            print()
            pass
    
    
        try:
            self.write_cutouts()
        except AttributeError:
            print('are you rushing to make a mask w/out selecting galaxies?')
            print('try selecting filter, then selecting galaxies')
            return
        #try:
        self.mwindow = QtWidgets.QWidget()
        print()
        print("initiating mask window")
        print("\t object params = ",objparams)
        try:
            self.mui = maskwindow(self.mwindow, self.logger, image = self.cutout_name_r, haimage=self.cutout_name_ha, sepath='~/github/halphagui/astromatic/',objparams=objparams)
        
            self.mui.mask_saved.connect(self.display_mask)
            self.mui.setupUi(self.mwindow)
            self.mwindow.show()
        except AttributeError:
            print('Hey - make sure you selected a galaxy!')
        #os.chdir(current_dir)
        
    def set_comment(self,comment): # MVC - controller or model??
        """ what is this doing? """
        col_names = ['CONTSUB_FLAG','MERGER_FLAG','SCATLIGHT_FLAG','ASYMR_FLAG','ASYMHA_FLAG','OVERSTAR_FLAG','OVERGAL_FLAG','EDGEON_FLAG','PARTIAL_FLAG','NUC_HA']
        self.table[col_names[int(comment)]][self.igal] = not(self.table[col_names[int(comment)]][self.igal])
        if not self.auto:
            if self.table[col_names[int(comment)]][self.igal]:
                pass_str = str(True)
            else:
                pass_str = str(False)
            self.update_gui_table_cell(self.igal,col_names[int(comment)],pass_str)
            self.write_fits_table()
        
        

    def closeEvent(self, event):
        # send signal that window is closed

        # write the data table
        print('writing fits table')
        self.table.write_fits_table()
        #event.accept()

class hafunctions(Ui_MainWindow, create_output_table, uco_table, hagui_methods, hagui_interactive, hacontroller):
    """ Main class for the halpha image analysis  """
    def __init__(self,MainWindow, logger, sepath=None, args=None):
        #testing=False,nebula=False,virgo=False,laptop=False,pointing=None,prefix=None,auto=False,obsyear=None,psfdir=None,rimage=None,haimage=None,csimage=None,filter=None,tabledir=None):
        super(hafunctions, self).__init__()
        self.auto = args.auto
        self.cscoadd_fname = None
        #if self.auto:
        #    matplotlib.use('TkAgg')
        #self.obsyear = args.obsyear

        # initialize psf as None
        self.psf = None
        if not(self.auto):
            self.setup_gui()
        self.prefix = args.prefix
        if (self.prefix is not None) and not(self.auto):
            self.set_prefix_on_gui(self.prefix)
        self.testing = args.testing
        self.draco = args.draco
        self.nebula = args.nebula        
        self.laptop = args.laptop
        self.virgo = args.virgo
        self.uat = args.uat
        self.verbose = args.verbose
        # this is the oversampling that I use when creating the PSF images
        self.oversampling = 2        
        if sepath == None:
            self.sepath = os.getenv('HOME')+'/github/halphagui/astromatic/'
        else:
            self.sepath = sepath
        if args.psfdir is None:
            self.psfdirectory = os.getcwd()
        else:
            self.psfdirectory = args.psfdir
        self.igal = None
        ############################################################
        ### CHECK TO SEE IF IMAGE NAMES ARE SPECIFIED
        ############################################################
        if args.rimage is not None:
            self.rcoadd_fname = args.rimage
        if args.haimage is not None:
            self.hacoadd_fname = args.haimage
        if args.csimage is not None:
            self.cscoadd_fname = args.csimage
        if args.filter is not None:
            self.filter = args.filter
            self.set_hafilter(args.filter)
        if args.tabledir is not None:
            self.tabledir = args.tabledir

        if (args.rimage is None):
            ############################################################
            ### CONFIGURATION SETUP FOR RUNNING ON DIFFERENT COMPUTERS
            ############################################################
        
            if self.draco & self.virgo:
                self.setup_draco_virgo()
                self.setup_virgo(pointing=pointing)
            elif self.laptop & self.virgo:
                self.setup_laptop_virgo()
                self.setup_virgo(pointing=pointing)
            elif self.nebula & self.virgo: 
                self.setup_nebula_virgo()
                if pointing is not None:
                    print('GOT A POINTING NUMBER FOR VIRGO FIELD')
                    self.setup_virgo(pointing=pointing)
                else:
                    self.setup_virgo()
            elif self.nebula:
                self.setup_nebula()
            elif self.testing:
                self.setup_testing()
        else:

            if not self.auto:
                # load rband image
                self.load_rcoadd()
            
                if args.haimage is not None:
                    # load Halpha
                    self.load_hacoadd()
            
            #print('running build_psf')
            #self.build_psf()

            
        ############################################################
        ### SET PARAMETERS FOR VIRGO VS HALPHA GROUPS PROJECT
        ### (WHICH USES THE NSA)
        ############################################################

        if self.virgo:
            self.setup_virgo_catalogs()
            
            self.defcat = self.vf
            self.def_label = 'VF.v1.'
            self.def_label = 'VF.v2.'            
            self.radius_label = 'radius'
            self.global_min_cutout_size = 60*u.arcsec
            self.global_max_cutout_size = 480*u.arcsec # 9 arcmin

        # this is for halpha groups analysis
        # commenting this out for now
        '''
        else:
            self.defcat = self.nsa
            self.def_label = 'NSAID'
            self.radius_label = 'PETROTH90'
            self.global_min_cutout_size = 100
            self.global_max_cutout_size = 250
        '''

        if self.uat:
            self.setup_uat_catalogs()
            self.defcat = self.agc
            self.global_min_cutout_size = 60*u.arcsec
            self.global_max_cutout_size = 360*u.arcsec # 6 arcmin
            
        self.initialize_uco_arrays()
        if self.auto:
            self.auto_run()

    def auto_run(self):
        # run the analysis without starting the gui
        # read in coadded images
        self.read_rcoadd()
        self.read_hacoadd()
        # set filter
        self.set_hafilter(self.filter)
        
        # get filter ratio
        # this will look for filter ratio in r-band image header
        # will run SE if ratio is not found in header
        self.get_filter_ratio()
        
        # subtract images
        # this will look for CS-ZP.fits image
        # if not found, it will subtract images according to the ZP ratio of r and halpha images
        self.subtract_images(overwrite=True) # checks out ok
        
        # measure psf
        # function will check if psf image already exists
        self.build_psf()

        # get galaxies
        self.find_galaxies()

        self.write_fits_table()
        # analyze galaxies
        # the default catalog has been trimmed to only include
        # galaxies in FOV

        # skipping for now
        # just making cutouts and getting galaxies in FOV

        # need to get a list of RA, DEC, SMA, BA, PA to feed into masking routine
        if self.verbose:
            print("starting processing of each galaxy ",len(self.gximage))
        for i in range(len(self.gximage)):
        #for i in [2]: # for testing
            self.igal = i
            # get cutouts
            self.auto_gal()
            if self.verbose:
                print("writing fits table after running auto_gal")
            self.write_fits_table()
            if self.verbose:
                print(f"##########################\nFinished galaxy {i+1}/{len(self.gximage)}")
        
    def auto_gal(self):
        # run the analysis on an individual galaxy
        if self.verbose:
            print("in auto_gal \n")
        self.bad_galaxy = False
        # create cutout
        if self.verbose:
            print("\ngetting galaxy cutouts\n")
        self.get_galaxy_cutout()

        if self.verbose:
            print("\nfinished cutouts\n")
        
        if self.bad_galaxy:
            print("\nzero std in sky - this is not real, so skipping galaxy ",self.cutout_name_r)
            print()
            return
        # create mask
        if args.uat:
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]
        else:
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]

        if self.verbose:
            print("initiating maskwindow\n")

        self.mui = maskwindow(None, None, image = self.cutout_name_r, haimage=self.cutout_name_ha, \
                              sepath='~/github/halphagui/astromatic/',auto=self.auto,\
                              objparams=self.objparams,unmaskellipse=True)


        # get SE parameters back from 
        # subtract the sky, using the mask image, and resave cutouts
        
        # run galfit

        # try skipping galfit for Messier 109, VFID1167...
        if self.igal == 1167:
            # skip galfit, galaxy is just too darn big and parameters are ridiculous
            print("skipping galfit b/c output is nonsense...")
        else:
            
            try:
                self.run_galfit(ncomp=1,ha=False)
            except:
                print('WARNING: problem running galfit ellip phot',self.cutout_name_r)

            # don't need this - the geometry is so similar
            self.galfit_ellip_phot()            
            try:
                if self.verbose:
                    print("running galfit_ellip_phot ")
                self.galfit_ellip_phot()
            except:
                print("##################################")
                print('WARNING: problem running galfit ellip phot',self.cutout_name_r)
                print("##################################")                

        # run galfit ellip phot
        # use try in case fit fails
        #self.galfit_ellip_phot()

        
        # run phot util ellip phot
        # use try in case fit fails
        #self.photutils_ellip_phot()

        # RF - debugging on 07/07/2024 - not getting size measurements in output file
        # taking command out of try/except
        self.photutils_ellip_phot()
        #try:
        #    self.photutils_ellip_phot()
        #except:
        #    print('\nWARNING: problem running photutils ellip phot\n',self.cutout_name_r)
        
            
      
    def setup_testing(self):
        #self.hacoadd_fname = os.getenv('HOME')+'/research/halphagui_test/MKW8_ha16.coadd.fits'
        #self.hacoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_ha16.coadd.fits'
        self.hacoadd_fname = os.getenv('HOME')+'/research/VirgoFilaments/Halpha/virgo-coadds-2017/pointing-1_ha4.coadd.fits'

        #self.ha, self.ha_header = fits.getdata(self.hacoadd_fname, header=True)
        self.haweight = self.hacoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.haweight_flag = True
        #self.rcoadd_fname = os.getenv('HOME')+'/research/halphagui_test/MKW8_R.coadd.fits'
        #self.rcoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_R.coadd.fits'
        self.rcoadd_fname = os.getenv('HOME')+'/research/VirgoFilaments/Halpha/virgo-coadds-2017/pointing-1_R.coadd.fits'
        #self.r, self.r_header = fits.getdata(self.rcoadd_fname, header=True)

        #self.rweight = self.rcoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.rweight_flag = True
        #self.pixelscale = abs(float(self.r_header['CD1_1']))*3600. # in deg per pixel
        self.nsa_fname = os.getenv('HOME')+'/research/NSA/nsa_v0_1_2.fits'
        self.nsa = galaxy_catalog(self.nsa_fname,nsa=True)
        self.agc_fname = os.getenv('HOME')+'/research/AGC/agcnorthminus1.2019Sep24.fits'
        self.agc = galaxy_catalog(self.agc_fname,agc=True)
        self.agcflag = True
        #self.coadd.load_file(self.rcoadd_fname)
        #self.filter_ratio = 0.0416
        self.filter_ratio = 0.0406
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
        self.load_hacoadd()
        self.load_rcoadd()        
        self.subtract_images()
        #self.setup_ratio_slider()
        #self.setup_cutout_slider()
    def setup_nebula(self):
        self.hacoadd_fname = '/mnt/astrophysics/reduced/20150418/MKW8_ha16.coadd.fits'
        #self.hacoadd_fname = '/mnt/qnap_home/rfinn/Halpha/reduced/virgo-coadds-2017/pointing-1_ha4.coadd.fits'
        #self.hacoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_ha16.coadd.fits'
        self.load_hacoadd()
        #self.ha, self.ha_header = fits.getdata(self.hacoadd_fname, header=True)
        #self.haweight = self.hacoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.haweight_flag = True
        self.rcoadd_fname ='/mnt/astrophysics/reduced/20150418/MKW8_R.coadd.fits'
        #self.rcoadd_fname = '/mnt/qnap_home/rfinn/Halpha/reduced/virgo-coadds-2017/pointing-1_R.coadd.fits'
        #self.rcoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_R.coadd.fits'
        #self.r, self.r_header = fits.getdata(self.rcoadd_fname, header=True)
        self.load_rcoadd()
        #self.rweight = self.rcoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.rweight_flag = True
        #self.pixelscale = abs(float(self.r_header['CD1_1']))*3600. # in deg per pixel
        #self.nsa_fname = '/mnt/qnap_home/share/catalogs/nsa_v0_1_2.fits'
        self.nsa_fname = '/mnt/astrophysics/catalogs/nsa_v0_1_2.fits'
        self.nsa = galaxy_catalog(self.nsa_fname,nsa=True)
        #self.agc_fname = '/mnt/qnap_home/share/catalogs/agcnorthminus1.2019Sep24.fits'
        self.agc_fname = '/mnt/astrophysics/catalogs/agcnorthminus1.2019Sep24.fits'
        self.agc = galaxy_catalog(self.agc_fname,agc=True)
        self.agcflag = True
        #self.coadd.load_file(self.rcoadd_fname)
        self.filter_ratio = 0.0422 #MKW8
        #self.filter_ratio = 0.0406
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
        self.subtract_images()
        #self.setup_ratio_slider()
        #self.setup_cutout_slider()
    def setup_nebula_virgo(self):
        self.imagedir =  '/mnt/astrophysics/reduced/virgo-coadds-2017/'
        self.tabledir= '/mnt/astrophysics/catalogs/Virgo/tables-north/v1/'
    def setup_laptop_virgo(self):
        if self.obsyear == '2018':
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-2018/'
        elif self.obsyear == '2020':
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-feb2020/'
        elif self.obsyear == '2019':
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-feb2019-int/'
        else:
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-2017/'
        self.tabledir= '/home/rfinn/research/Virgo/tables-north/v1/'
    def setup_draco_virgo(self):
        self.imagedir =  '/data-pool/Halpha/coadds/all-virgo-coadds/'
        self.tabledir= '/home/rfinn/research/Virgo/tables-north/v2/'
    def setup_virgo(self,pointing=None):
        ''' construct image names from input.  only works for data reduced before 2021 '''
        if pointing is None:
            self.hacoadd_fname = self.imagedir+'pointing-3_ha4.coadd.fits'
            self.rcoadd_fname = self.imagedir+'pointing-3_R.coadd.fits'
        else:
            #print('got a pointing')
            if self.obsyear == '2020':
                self.hacoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_ha4.coadd.fits'
                self.rcoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_r.coadd.fits'
            else:
                self.hacoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_ha4.coadd.fits'
                self.rcoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_R.coadd.fits'
            

        if not self.auto:
            self.load_hacoadd()
            self.load_rcoadd()
            try:
                self.subtract_images()
                #self.setup_ratio_slider()
                #self.setup_cutout_slider()
            except:
                pass
        self.setup_virgo_catalogs()

    def setup_virgo_catalogs(self):
        ## UPDATES TO USE VIRGO FILAMENT MASTER TABLE
        #self.vf_fname = self.tabledir+'vf_north_v1_main.fits'
        #self.nsa_fname = self.tabledir+'vf_north_v1_nsa_v0.fits'
        # RF 2023-03-24: updating to use the v2 catalogs
        self.vf_fname = os.path.join(self.tabledir,'vf_v2_main.fits')
        self.nsa_fname = os.path.join(self.tabledir,'vf_v2_nsa_v0.fits')
        ephot_fname = os.path.join(self.tabledir,'vf_v2_legacy_ephot.fits')        
        self.vf = galaxy_catalog(self.vf_fname,virgo=True)

        self.nsa = galaxy_catalog(self.nsa_fname,virgo=True)

        ##
        # get sizes for galaxies - will use this to unmask central region
        # need to cut this catalog based on keepflag
        ##

        ephot = Table.read(ephot_fname)

        
        #self.radius_arcsec = ephot['SMA_SB24']

        bad_sb25 = ephot['SMA_SB25'] == 0

        self.radius_arcsec = ephot['SMA_SB25']*(~bad_sb25) + 1.35*ephot['SMA_SB24']*bad_sb25
        # OK, I know what you are thinking, I can't possibly be changing this again...

        # use SMA_SB25 instead of SB24 - this should work better for both high and low SB galaxies
        # if SMA_SB25 is not available use 1.35*SMA_SB24

        # for galaxies with SMA_SB24=0, set radius to value in main table 
        noradius_flag = self.radius_arcsec == 0
        self.radius_arcsec[noradius_flag] = self.vf.cat['radius'][noradius_flag]

        # also save BA and PA from John's catalog
        # use the self.radius_arcsec for the sma
        self.BA = np.ones(len(self.radius_arcsec))
        self.PA = np.zeros(len(self.radius_arcsec))
        
        self.BA[~noradius_flag] = ephot['BA_MOMENT'][~noradius_flag]
        self.PA[~noradius_flag] = ephot['PA_MOMENT'][~noradius_flag]
        
        self.RA = self.vf.cat['RA']
        self.DEC = self.vf.cat['DEC']        
        
        self.agcflag = False
        self.nsaflag = False

        
        #self.coadd.load_file(self.rcoadd_fname)
        #self.filter_ratio = 0.0422 #MKW8
        # filter ratio for ha4
        self.filter_ratio = 0.0426
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
        

    def setup_uat_catalogs(self):
        """
        use the latest AGC as the source catalog
        """

        # on RF laptop, tabledir is /Users/rfinn/research/
        #self.agc_fname = os.path.join(self.tabledir,'AGC/agc.allsky.210720.fits')
        self.agc_fname = os.path.join(self.tabledir,'AGC/agcnorthminus1.full200617.fits')
            

        self.agc = galaxy_catalog(self.agc_fname,virgo=False,agc=True)

        self.nsa_fname = os.path.join(self.tabledir,'NSA/nsa_v1_0_1.fits')
        if not os.path.exists(self.nsa_fname):
            self.nsa_fname = os.path.join(self.tabledir,'NSA/nsa_v0_1_2.fits')            
            
        self.nsa = galaxy_catalog(self.nsa_fname,virgo=False, agc=False, nsa=True)        
        #self.agc.check_ra_colname() # should check automatically with agc=True
        ##
        # get sizes for galaxies - will use this to unmask central region
        # need to cut this catalog based on keepflag
        ##

        #agc['a'] is the blue semi-major diameter in arcmin
        self.radius_arcsec = self.agc.cat['a']*60
        
        noradius_flag = self.radius_arcsec == 0
        self.radius_arcsec[noradius_flag] = 60 # set size of galaxies with no A value to 60 arcsec

        # also save BA and PA from John's catalog
        # use the self.radius_arcsec for the sma
        self.BA = np.ones(len(self.radius_arcsec))
        self.PA = np.zeros(len(self.radius_arcsec))
        
        self.BA[~noradius_flag] = self.agc.cat['b'][~noradius_flag]/self.agc.cat['a'][~noradius_flag]

        print("setting PA to posang from AGC")
        self.PA[~noradius_flag] = self.agc.cat['posang'][~noradius_flag]
        
        self.RA = self.agc.cat['RA']
        self.DEC = self.agc.cat['DEC']        
        
        self.agcflag = True
        self.nsaflag = True

        
        # filter ratio for ha4
        self.filter_ratio = 0.0426
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
		

        
            
        
if __name__ == "__main__":
    ## RUNNING AS THE MAIN PROGRAM
    


    #####################################
    ## SETUP COMMAND-LINE PARAMETERS
    #####################################
    import argparse    
    parser = argparse.ArgumentParser(description ='Run gui for analyzing Halpha images')

    parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    
    parser.add_argument('--rimage',dest = 'rimage', default=None,help='r-band image')
    parser.add_argument('--haimage',dest = 'haimage', default=None,help='Halpha image')
    parser.add_argument('--csimage',dest = 'csimage', default=None,help='Continuum-subtracted Halpha image')    
    parser.add_argument('--filter',dest = 'filter', default=None,help='filter. options are 4, 8, 12, 16, inthalpha, or intha6657')
    parser.add_argument('--tabledir',dest = 'tabledir', default=None,help='table directory. something like /home/rfinn/research/Virgo/tables-north/v1/')
    parser.add_argument('--psfdir',dest = 'psfdir', default=None,help='set this to the directory containing PSF images')        
    parser.add_argument('--prefix',dest = 'prefix', default='v17p03',help='prefix associated with the coadded image.  Default is v17p03. required when running auto.')
    parser.add_argument('--auto',dest = 'auto', action='store_true',default=False,help='set this to process the images automatically, without the gui')
    
    parser.add_argument('--virgo',dest = 'virgo', action='store_true',default=False,help='set this if running on virgo data.  The virgo filaments catalog will be used as input.')
    parser.add_argument('--uat',dest = 'uat', action='store_true',default=False,help='set this if running on uat halpha groups.  The AGC (210720) will be used as the parent catalog.')     
    parser.add_argument('--draco',dest = 'draco', action='store_true',default=False,help='set this if running on draco.')   
    parser.add_argument('--nebula',dest = 'nebula', action='store_true',default=False,help='set this if running on open nebula virtual machine.  catalog paths will be set accordingly.')
    parser.add_argument('--laptop',dest = 'laptop', action='store_true',default=False,help="custom setting for running on Rose's laptop. catalog paths will be set accordingly.")
    
    parser.add_argument('--obsyear',dest = 'obsyear', default=None,help='year that data were taken.  this finds the right image directory if you are building the image name in pieces..  ')
    parser.add_argument('--pointing',dest = 'pointing', default=None,help='Pointing number that you want to load.  ONLY FOR VIRGO DATA, and only if you are buildling the image name in pieces.')
    
    parser.add_argument('--testing',dest = 'testing', action='store_true',default=False,help='set this if running on open nebula virtual machine')
    parser.add_argument('--onegal',dest = 'onegal', default=None, help='provide galaxy name to run halpha gui just on one galaxy')    
    parser.add_argument('--verbose',dest = 'verbose', action='store_true',default=False,help='set this for extra print statements')    
        
    args = parser.parse_args()
    
    logger = log.get_logger("example1", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)


    sepath = os.getenv('HOME')+'/github/halphagui/astromatic/'



    #################################
    ## UPDATED TO USE ARGPARSE
    #################################
    if not args.auto:
        ui = hafunctions(MainWindow, logger, sepath = sepath, args=args)        

        MainWindow.show()
        sys.exit(app.exec_())
    else:
        # run functions non-interactively
        #ui = hafunctions(MainWindow, logger, sepath = sepath, testing=args.testing,nebula=args.nebula,virgo=args.virgo,laptop=args.laptop,pointing=args.pointing,auto=args.auto,prefix=args.prefix,obsyear=args.obsyear,rimage=args.rimage,haimage=args.haimage,csimage=args.csimage,filter=args.filter,tabledir=args.tabledir,psfdir=args.psfdir)
        ui = hafunctions(MainWindow, logger, sepath = sepath, args=args)
        pass
