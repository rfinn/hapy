#!/usr/bin/env python

'''
PURPOSE:

The goal of the program is to create a mask for a galaxy image to mask
out other objects within the cutout area.

USAGE:


you just need to run this on R-band images.


PROCEDURE:


REQUIRED MODULES:
   os
   astropy
   numpy
   argsparse
   matplotlib
   scipy

USAGE:

* if running on wise images, try:

python ~/github/halphagui/maskwrapper.py --image AGC006015-unwise-1640p454-w1-img-m.fits --ngrow 1 --sesnr 2 --minarea 5 --auto


NOTES:
- rewrote using a class

# TODO: 2023-02-09: this program relies on source extractor.  I should rewrite to use photutils instead.

TESTING

/home/rfinn/research/Virgo-dev/maskwrapper-test/VFID0610-rectangle//VFID0610-NGC5985-INT-20190530-p040


objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]


python ~/github/halphagui/maskwrapper.py --image VFID0610-NGC5985-INT-20190530-p040-R.fits --haimage VFID0610-NGC5985-INT-20190530-p040-CS.fits --sepath ~/github/halphagui/astromatic/ --gaiapath /home/rfinn/research/legacy/gaia-mask-dr9.virgo.fits --objra 234.90448 --objdec 59.33198 --objsma 139.25 --objBA .496 --objPA 104.646


'''

import os
import sys
import numpy as np
import warnings

from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Tophat2DKernel, convolve
from astropy.convolution.kernels import CustomKernel
from astropy.table import Table
from astropy.coordinates import SkyCoord
import numpy as np
#import argparse
#import pyds9
from scipy.stats import scoreatpercentile


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import patches

#from maskGui import Ui_maskWindow
from hapy.maskgui.maskWidget import Ui_Form as Ui_maskWindow
from hapy.hagui.cutout_window import CutoutImage
from hapy.masktools.maskops import circle_pixels
# import gaia function to get stars within region
from hapy.masktools.gaia import gaia_stars_in_rectangle

import hapy.imagetools.imutils as imutils

try:
    from photutils.segmentation import detect_threshold, detect_sources
    #from photutils import source_properties
    from photutils.segmentation import SourceCatalog    
    from photutils.segmentation import deblend_sources
    
except ModuleNotFoundError:
    warnings.warn("Warning - photutils not found")
except ImportError:
    print("got an import error with photutils - check your version number")

from PyQt5 import QtCore, QtGui, QtWidgets
try:
    from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
    from ginga.mplw.ImageViewMpl import ImageView
    from ginga import colors
    from ginga.canvas.CanvasObject import get_canvas_types

    from ginga.misc import log
    from ginga.util.loader import load_data
    gingaflag = True
except ModuleNotFoundError:
    print("Warning - ginga was not found.  this will be a problem if running interactively")
    gingaflag = False
import timeit


#####################################
###  FUNCTIONS
#####################################



class LegacyMaskEngine():

    def update_mask(self):
        self.add_user_masks()
        print("starting add_gaia_masks...")
        self.add_gaia_masks()
        self.write_mask()

    def write_mask(self):
        """ write out mask image """

        # add ellipse params to imheader
        if self.ellipseparams is not None:
            #print("HEY!!!")
            #print()
            #print("Writing central ellipse parameters to header")
            #print(self.ellipseparams)
            #print()
            if hasattr(self.objsma,"__len__"):
                xc,yc,r,BA,PA = self.ellipseparams[0]
            else:
                xc,yc,r,BA,PA = self.ellipseparams
            self.imheader.set('ELLIP_XC',float(xc),comment='XC of mask ellipse')
            self.imheader.set('ELLIP_YC',float(yc),comment='YC of mask ellipse')
            self.imheader.set('ELLIP_A',r,comment='SMA of mask ellipse')
            self.imheader.set('ELLIP_BA',BA,comment='BA of mask ellipse')
            self.imheader.set('ELLIP_PA',np.degrees(PA),comment='PA (deg) of mask ellipse')
        else:
            print("HEY!!! writing mask, but no parameters for central ellipse!")

            
        fits.writeto(self.mask_image,self.maskdat,header = self.imheader,overwrite=True)
        invmask = self.maskdat > 0.
        invmask = np.array(~invmask,'i')
        fits.writeto(self.mask_inv_image,invmask,header = self.imheader,overwrite=True)
        if not self.auto:
            self.mask_saved.emit(self.mask_image)
            self.display_mask()
    def add_gaia_masks(self):
        # check to see if gaia stars were already masked
        if self.add_gaia_stars:
            if self.gaia_mask is None :
                try:
                    self.get_gaia_stars()
                except:
                    print("\nWARNING: problem getting gaia stars - do you have an internet connection?\nI will not add gaia stars to mask\n")
                    self.gaia_mask = None
                    return
                
                self.make_gaia_mask()
            else:
                self.maskdat += self.gaia_mask

 

    def show_mask_mpl(self):
        # plot mpl figure
        # this was for debugging purposes
        print("plotting mask and central ellipse")
        self.fig = plt.figure(1,figsize=self.figure_size)
        plt.clf()
        plt.subplots_adjust(hspace=0,wspace=0)
        plt.subplot(1,2,1)
        plt.imshow(self.image,cmap='gray_r',vmin=self.v1,vmax=self.v2,origin='lower')
        plt.title('image')
        plt.subplot(1,2,2)
        #plt.imshow(maskdat,cmap='gray_r',origin='lower')
        plt.imshow(self.maskdat,cmap=self.cmap,origin='lower',vmin=np.min(self.maskdat),vmax=np.max(self.maskdat))
        plt.title('mask')
        plt.gca().set_yticks(())
        #plt.draw()
        #plt.show(block=False)
        #print("in show_mask_mpl: objsma = ",self.objsma)        
        try:
            
            if hasattr(self.objsma, "__len__"):
                #print("working with multiple galaxies")
                # add ellipse for each galaxy if there is more than one
                for e in self.ellipseparams:
                    xc,yc,r,BA,PA = e
                    PAdeg = np.degrees(PA)
                    #print(f"BA={BA},PA={PAdeg} deg")        
                    #print("just checking - adding ellipse drawing ",self.ellipseparams)
                    ellip = patches.Ellipse((xc,yc),2*r,2*r*BA,angle=PAdeg,alpha=.2)
                    plt.gca().add_patch(ellip)
            else:
                xc,yc,r,BA,PA = self.ellipseparams
                PAdeg = np.degrees(PA)
                #print(f"BA={BA},PA={PAdeg} deg")        
                #print("just checking - adding ellipse drawing ",self.ellipseparams)
                ellip = patches.Ellipse((xc,yc),r,r*BA,angle=PAdeg,alpha=.2)
                plt.gca().add_patch(ellip)

        except:
            print("problem plotting ellipse with mask")
        # outfile
        outfile = self.mask_image.replace('.fits','.png')
        plt.savefig(outfile)
        
        #plt.show()
        




if __name__ == "__main__":
    #catalog = '/Users/rfinn/research/NSA/nsa_v0_1_2.fits'
    #gcat = galaxy_catalog(catalog)
    #from halphamain import cutout_image
    #from halphamain import cutout_image
    import argparse    
    parser = argparse.ArgumentParser(description ='Run gui for making an mask.  You can specify the RA and DEC of galaxy, which is useful if galaxy is not at the center of the cutout image.  You can also provide an elliptical region around the galaxy to unmask.  This is useful for galaxies that are shredded by source extractor.')
    parser.add_argument('--image',dest = 'image', default=None,help='r-band image')
    parser.add_argument('--haimage',dest = 'haimage', default=None,help='this is typically the continuum-subtracted Halpha image.  If no image is provided, the middle panel is left blank.')
    parser.add_argument('--sepath',dest = 'sepath', default=None,help='path to source extractor config files (e.g. ~/github/HalphaImaging/astromatic/ - this is default if no path is given.)')
    parser.add_argument('--gaiapath',dest = 'gaiapath', default=None,help='full pathname of gaia mask file from legacy dr9.')    
    parser.add_argument('--config',dest = 'config', default=None,help='source extractor config file.  default is default.sex.HDI.mask')
    parser.add_argument('--objra',dest = 'objra', default=None,help='RA of target galaxy. default is none, then object is assumed to be at center of image.')
    parser.add_argument('--objdec',dest = 'objdec', default=None,help='DEC of target galaxy')
    parser.add_argument('--objsma',dest = 'objsma', default=None,help='SMA of elliptical region to unmask around galaxy.')
    parser.add_argument('--objBA',dest = 'objBA', default=None,help='BA of elliptical region to unmask around galaxy.')
    parser.add_argument('--objPA',dest = 'objPA', default=None,help='PA of elliptical region to unmask around galaxy, measure CCW from +x axis')
    parser.add_argument('--ngrow',dest = 'ngrow', default=3,help='number of times to run grow the masked regions in auto mode.  default is 7, which is reasonable for an optical image.  try 1 if running on WISE images.')
    parser.add_argument('--auto',dest = 'auto', default=False,action='store_true',help='set this to run the masking software automatically.  the default is false, meaning that the gui window will open for interactive use.')
    parser.add_argument('--sesnr',dest = 'sesnr', default=10,help='adjust the SE SNR for detection.  Default is 10.')
    parser.add_argument('--minarea',dest = 'minarea', default=5,help='adjust the SE detection area.  Default is 10.')
    parser.add_argument('--weightim',dest = 'weightim', default=None,help='weight image to feed into source extractor.  You can use this to ignore large regions of the image by creating a mask with good regions = 1 and bad regions = 0.  then set weight_thresh to 1.')
    parser.add_argument('--weight_thresh',dest = 'weight_thresh', default=1,help='source extractor weight_thresh.  pixels in the weight images with values less than this threshold will be ignored.  Default is 1.')                        
        
    args = parser.parse_args()
    if (args.objra is not None) and (args.objBA is not None):
        objparams = [float(args.objra),float(args.objdec),float(args.objsma),float(args.objBA),float(args.objPA)]
    else:
        objparams = None
    if gingaflag:
        logger = log.get_logger("masklog", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()

    if args.objsma is not None:
        unmaskellipse = True
    else:
        unmaskellipse = False
    print("testing - unmaskellipse = ",unmaskellipse)
    if args.image is not None:
        if not args.auto:
            
            #print('got here 1')
            MainWindow = QtWidgets.QWidget()
            ui = maskwindow(MainWindow, logger,image=args.image,haimage=args.haimage,sepath=args.sepath,gaiapath=args.gaiapath,config=args.config,auto=args.auto,objparams=objparams,unmaskellipse = unmaskellipse,snr=args.sesnr,minarea=args.minarea,weightim=args.weightim,weight_threshold=args.weight_thresh)
        else:
            #print('got here 2')
            ui = maskwindow(None, None,image=args.image,haimage=args.haimage,sepath=args.sepath,gaiapath=args.gaiapath,config=args.config,auto=args.auto,objparams=objparams,unmaskellipse=unmaskellipse,snr=args.sesnr,minarea=args.minarea,ngrow=args.ngrow,weightim=args.weightim,weight_threshold=args.weight_thresh)
    else:
        #print('got here 3')
        MainWindow = QtWidgets.QWidget()
        ui = maskwindow(MainWindow, logger)
    #ui.setupUi(MainWindow)
    #ui.test()
    if not args.auto: 
        MainWindow.show()
        sys.exit(app.exec_())

    
