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
from hapy.imagetools.imutils import circle_pixels
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

defaultcat='default.sex.HDI.mask'
#####################################
###  FUNCTIONS
#####################################

def remove_central_objects(mask, sma=20, BA=1, PA=0, xc=None,yc=None):
    """ 
    find any pixels within central ellipse and set their values to zero 

    PARAMS:
    mask = 2D array containing masked pixels, like from SE segmentation image
    sma = semi-major axis in pixels
    BA = ratio of semi-minor to semi-major axes
    PA = position angle, measured in degree counter clockwise from +x axis

    OPTIONAL ARGS:
    xc = center of ellipse in pixels; assumed to be center of image if xc is not specified
    yc = center of ellipse in pixels; assumed to be center of image if yc is not specified
    
    RETURNS:
    newmask = copy of input mask, with pixels within ellipse set equal to zero

    """
    # changing the xmax and ymax - if the ellipse looks wrong, then swap back
    ymax,xmax = mask.shape
    # set center of ellipse as the center of the image
    if (xc is None) and (yc is None):
        xc,yc = xmax//2,ymax//2
    
    a = sma
    b = BA*sma
    phirad = np.radians(PA)

    X,Y = np.meshgrid(np.arange(xmax),np.arange(ymax))
    
    p1 = ((X-xc)*np.cos(phirad)+(Y-yc)*np.sin(phirad))**2/a**2
    p2 = ((X-xc)*np.sin(phirad)-(Y-yc)*np.cos(phirad))**2/b**2
    flag2 = p1+p2 < 1
    newmask = np.copy(mask)
    newmask[flag2] = 0
    # we could also get all the unique values associated with flag2, and then remove them
    ellipse_params = [xc,yc,sma,BA,phirad]
    return newmask,ellipse_params


class buildmask():
    def link_files(self):
        # TODO: replace sextractor with photutils
        # these are the sextractor files that we need
        # set up symbolic links from sextractor directory to the current working directory
        sextractor_files=['default.sex.HDI.mask','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            if os.path.exists(file):
                os.remove(file)
            os.system('ln -sf '+self.sepath+'/'+file+' .')
            #os.copy(self.sepath+'/'+file, file)
    def clean_links(self):
        # clean up symbolic links to sextractor files
        # sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
        sextractor_files=['default.sex.HDI.mask','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('unlink '+file)

    def read_se_cat(self):
        #print("reading in SE catalog")
        sexout=fits.getdata(self.catname)
        self.se_number = sexout['NUMBER']
        self.xsex = sexout['XWIN_IMAGE']
        self.ysex = sexout['YWIN_IMAGE']
        self.fwhm = sexout['FWHM_IMAGE']
        self.A_IMAGE = sexout['A_IMAGE']
        self.B_IMAGE = sexout['B_IMAGE']
        self.BA = 1./sexout['ELONGATION']
        self.THETA_IMAGE = sexout['THETA_IMAGE']        
        # also keep the ellipse params
        
        dist=np.sqrt((self.yc-self.ysex)**2+(self.xc-self.xsex)**2)
        #   find object ID

        # some objects are rturning an empty sequence - how to handle this?
        # I guess this means that the object itself wasn't detected?
        # or nothing was detected?
        if len(dist) < 1:
            # set objnumb to nan
            objnumb = np.nan
        else:
            objIndex=np.where(dist == min(dist))
            objNumber=sexout['NUMBER'][objIndex]
            objnumb = objNumber[0] # not sure why above line returns a list
        return objnumb

    def runse(self,galaxy_id = None,weightim=None,weight_threshold=1):
        # TODO update this to implement running SE with two diff thresholds
        # TODO make an alternate function that creates segmentation image from photutils
        # this is already done in ell
        print('using a deblending threshold = ',self.threshold)
        print("image = ",self.image_name)
        self.catname = self.image_name.replace('.fits','.cat')
        self.segmentation = self.image_name.replace('.fits','-segmentation.fits')

        print("segmentation image = ",self.segmentation)
        sestring = f"sex {self.image_name} -c {self.config} -CATALOG_NAME {self.catname} -CATALOG_TYPE FITS_1.0 -DEBLEND_MINCONT {self.threshold} -DETECT_THRESH {self.snr} -ANALYSIS_THRESH {self.snr_analysis} -CHECKIMAGE_NAME {self.segmentation} -DETECT_MINAREA {self.minarea}"
        if weightim is not None:
            sestring += r" -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE {weightim} -WEIGHT_THRESH {weight_threshold}"
        print(sestring)
        os.system(sestring)
        self.maskdat = fits.getdata(self.segmentation)
        # grow masked areas
        bool_array = np.array(self.maskdat.shape,'bool')
        #for i in range(len(self.xsex)):
        # check to see if the object is not centered in the cutout

        
    def get_photutils_mask(self,galaxy_id = None):
        # TODO make an alternate function that creates segmentation image from photutils
        from astropy.stats import sigma_clipped_stats
        from photutils import make_source_mask

        # create mask to cut low SNR pixels based on SNR in SFR image
        mask = make_source_mask(imdat,nsigma=self.snr,npixels=self.minarea,dilate_size=5)
        masked_data = np.ma.array(imdat,mask=mask)
        
        
        self.catname = self.image_name.replace('.fits','.cat')
        self.segmentation = self.image_name.replace('.fits','-segmentation.fits')

        
        self.maskdat = fits.getdata(self.segmentation)
        # grow masked areas
        bool_array = np.array(self.maskdat.shape,'bool')
        #for i in range(len(self.xsex)):
        # check to see if the object is not centered in the cutout
    def remove_center_object(self):
        """ this removes the object in the center of the mask, which presumably is the galaxy """
        # need to replace this with a function that will remove any objects within the specificed central ellipse
        if self.remove_center_object_flag:
            if self.off_center_flag:
                print('setting center object to objid ',self.galaxy_id)
                self.center_object = self.galaxy_id
            else:
                self.center_object = self.read_se_cat()
            if self.center_object is not np.nan:
                self.maskdat[self.maskdat == self.center_object] = 0


        if self.objsma is not None:
            if hasattr(self.objsma, "__len__"):
                # loop over objects in fov
                self.ellipseparams = []                
                for i in range(len(self.objsma)):

                    #print(f"sma={self.objsma_pixels[i]},BA={self.objBA[i]}, PA={self.objPA[i]},xc={self.xpixel[i]},yc={self.ypixel[i]}")
                    self.maskdat,eparams = remove_central_objects(self.maskdat, sma=self.objsma_pixels[i],\
                                                                             BA=self.objBA[i], PA=self.objPA[i], \
                                                                             xc=self.xpixel[i],yc=self.ypixel[i])
                    self.ellipseparams.append(eparams)
                
                pass
            else:
                # remove central objects within elliptical aperture
                #print("ellipse params in remove_central_object :",self.xpixel,self.ypixel,self.objsma_pixels,self.objBA,self.objPA)
                self.maskdat,self.ellipseparams = remove_central_objects(self.maskdat, sma=self.objsma_pixels, \
                                                                         BA=self.objBA, PA=self.objPA, \
                                                                         xc=self.xpixel,yc=self.ypixel)
        else:
            
            print("no ellipse params")
            self.ellipseparams = None
        self.update_mask()


   def run_photutil(self, snrcut=1.5,npixels=10):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''
        self.threshold = detect_threshold(self.image, nsigma=snrcut)
        segment_map = detect_sources(self.image, self.threshold, npixels=npixels)
        # deblind sources a la source extractor
        # tried this, and the deblending is REALLY slow
        # going back to source extractor
        self.segmentation = deblend_sources(self.image, segment_map,
                               npixels=10, nlevels=32, contrast=0.001)        
        self.maskdat = self.segmentation.data
        #self.cat = source_properties(self.image, self.segmentation)
        self.cat = SourceCatalog(self.image, self.segmentation)        
        # get average sky noise per pixel
        # threshold is the sky noise at the snrcut level, so need to divide by this
        self.sky_noise = np.mean(self.threshold)/snrcut
        #self.tbl = self.cat.to_table()

        if self.off_center_flag:
            print('setting center object to objid ',self.galaxy_id)
            self.center_object = self.galaxy_id
        else:
            distance = np.sqrt((self.cat.xcentroid - self.xc)**2 + (self.cat.ycentroid - self.yc)**2)
            # save object ID as the row in table with source that is closest to center
            objIndex = np.arange(len(distance))[(distance == min(distance))][0]
            # the value in shown in the segmentation image is called 'label'
            self.center_object = self.cat.label[objIndex]

        self.maskdat[self.maskdat == self.center_object] = 0
        self.update_mask()

    def update_mask(self):
        self.add_user_masks()
        print("starting add_gaia_masks...")
        self.add_gaia_masks()
        self.write_mask()
    def add_user_masks(self):
        """ this adds back in the objects that the user has masked out """
        # add back the square masked areas that the user added
        self.maskdat = self.maskdat + self.usr_mask
        # remove objects that have already been deleted by user
        if len(self.deleted_objects) > 0:
            for objID in self.deleted_objects:
                self.maskdat[self.maskdat == objID] = 0.
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

    def grow_mask(self, size=7):

        """
        Convolution: one way to grow the mask is to convolve the image with a kernel

        however, this does not preserve the pixels value of the original
        object, which come from the sextractor segmentation image.

        if the user wants to remove an object, it's much easier to do this
        by segmentation number rather than by pixels (as in the reverse of how we add objects
        to mask).

        Alternative: is to loop over just the masked pixels, and replace all pixels
        within a square region with the masked value at the central pixel.
        This will preserve the numbering from the segmentation image.

        Disadvantage: everything assumes a square shape after repeated calls.

        Alternative is currently implemented.
        
        """
        # convolve mask with top hat kernel
        # kernel = Tophat2DKernel(5)
        #mykernel = np.ones([5,5])
        #kernel = CustomKernel(mykernel)
        #self.maskdat = convolve(self.maskdat, kernel)
        #self.maskdat = np.ceil(self.maskdat)

        # we don't want to grow the size of the gaia stars, do we???
        if self.gaia_mask is not None:
            self.maskdat -= self.gaia_mask
        nx,ny = self.maskdat.shape
        masked_pixels = np.where(self.maskdat > 0.)
        for i,j in zip(masked_pixels[0], masked_pixels[1]):
            rowmin = max(0,i-int(size/2))
            rowmax = min(nx,i+int(size/2))
            colmin = max(0,j-int(size/2))
            colmax = min(ny,j+int(size/2))
            if rowmax <= rowmin:
                # something is wrong, return without editing mask
                continue
            if colmax <= colmin:
                # something is wrong, return without editing mask
                continue
            #print(i,j,rowmin, rowmax, colmin, colmax)
            self.maskdat[rowmin:rowmax,colmin:colmax] = self.maskdat[i,j]*np.ones([rowmax-rowmin,colmax-colmin])
        # add back in the gaia star masks
        if self.gaia_mask is not None:
            self.maskdat += self.gaia_mask
        if not self.auto:
            self.display_mask()
        # save convolved mask as new mask
        self.write_mask()

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
        


class my_cutout_image(QtCore.QObject):#QtCore.QObject):
    #mouse_clicked = QtCore.pyqtSignal(str)
    key_pressed = QtCore.pyqtSignal(str)
    def __init__(self,panel_name,ui, logger, row, col, drow, dcol,autocut_params='stddev',auto=False):
        #super(image_panel, self).__init__(panel_name,ui, logger, row, col, drow, dcol)
        # enable some user interaction
        #fi.get_bindings.enable_all(True)
        #super(my_cutout_image, self).__init__(panel_name,ui, logger, row, col, drow, dcol)
        #super(my_cutout_image,self).__init__(self)
        QtCore.QObject.__init__(self)
        self.logger = logger
        self.ui = ui
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params(autocut_params)
        fi.enable_autozoom('on')
        #fi.set_callback('drag-drop', self.drop_file)
        #fi.set_callback('none-move',self.cursor_cb)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        self.fitsimage = fi
        fi.show_focus_indicator(True)
        
        bd = fi.get_bindings()
        bd.enable_all(True)


        self.cutout = fi.get_widget()
        self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.fitsimage = fi
        
        self.fitsimage.set_callback('none-move',self.cursor_cb)
        #self.fitsimage.set_callback('cursor-down',self.cursor_down)
        self.readout = QtWidgets.QLabel('')
        ui.readoutGridLayout.addWidget(self.readout, 1, col, 1, 1)
        #self.ui.readoutLabel.setText('this is another test')
        self.fitsimage.set_callback('key-press',self.key_press_cb)


        # adding lines to allow drawing on canvas?
        self.dc = get_canvas_types()        
        canvas = self.dc.DrawingCanvas()
        canvas.enable_draw(True)
        canvas.enable_edit(True)
        canvas.set_drawtype('rectangle', color='lightblue')
        canvas.set_surface(fi)
        #canvas.rectangle(.5,.5,10,0)
        self.canvas = canvas
        # add canvas to view
        #fi.add(canvas)
        private_canvas = fi.get_canvas()
        private_canvas.add(canvas)
        canvas.register_for_cursor_drawing(fi)
        #canvas.add_callback('draw-event', self.draw_cb)
        canvas.set_draw_mode('draw')
        canvas.ui_set_active(True)
        self.canvas = canvas

        self.drawtypes = canvas.get_drawtypes()
        self.drawtypes.sort()

        
    def load_image(self, imagearray):
        #self.fitsimage.set_image(imagearray)
        self.fitsimage.set_data(imagearray)

    def load_file(self, filepath):
        image = load_data(filepath, logger=self.logger)

        self.fitsimage.set_image(image)
        #self.setWindowTitle(filepath)

        # adding cursor readout so we can identify the pixel values to change

    def cursor_cb(self, viewer, button, data_x, data_y):
        """This gets called when the data position relative to the cursor
        changes.
        """
        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + viewer.data_off),
                                    int(data_y + viewer.data_off))

        except Exception:
            value = None

        fits_x, fits_y = data_x + 1, data_y + 1
        
        try:
            text = "X: %.1f  Y: %.1f  Value: %.2f" % (fits_x, fits_y, float(value))
            #self.readout.setText('this is another test')
            self.readout.setText(text)
        except:
            pass
    '''    
    def cursor_down(self, fitsimage, event, data_x, data_y):
        print('clicked at ',data_x, data_y)
        #self.mouse_clicked.emit(str(data_x)+','+str(data_y))
    '''    
    def key_press_cb(self, canvas, keyname):
        print('key pressed! ',keyname)
        data_x, data_y = self.fitsimage.get_last_data_xy()
        self.key_pressed.emit(keyname+','+str(data_x)+','+str(data_y))
        #return self.imexam_cmd(self.canvas, keyname, data_x, data_y, func)

        
class maskwindow(Ui_maskWindow, QtCore.QObject,buildmask):
    mask_saved = QtCore.pyqtSignal(str)
    def __init__(self, MainWindow, logger, image=None, haimage=None, sepath=None, gaiapath=None, config=None, threshold=0.005,snr=10,cmap='gist_heat_r',objparams=None,auto=False,unmaskellipse=False,minarea=10,ngrow=3,weightim=None,weight_threshold=None):
        """

        ngrow : number of times to run grow when running in auto mode
        """
        self.auto = auto
        if MainWindow is None:
            self.auto = True
        if not self.auto:
            super(maskwindow, self).__init__()
        
            self.ui = Ui_maskWindow()
            self.ui.setupUi(MainWindow)
            MainWindow.setWindowTitle('makin a mask...')
            self.MainWindow = MainWindow


            self.logger = logger
            print("in maskwindow, I get objparams = ",objparams)
        # define the position of the target galaxy, as well as the shape and size of elliptical region to unmask around galaxy.
        #print("inside maskwrapper.init, objparams = ",objparams)
        if objparams is not None:
            self.objra = objparams[0]
            self.objdec = objparams[1]
            self.objsma = objparams[2]
            self.objBA = objparams[3]
            self.objPA = objparams[4]

        else:
            self.objra  = None  
            self.objdec = None  
            self.objsma = None  
            self.objBA  = None  
            self.objPA  = None

        if (self.objsma is not None): # unmask central elliptical region around object
            # get wcs from mask image
            wcs = WCS(fits.getheader(image))
            
            # get x and y coord of galaxy from (RA,DEC) using mask wcs
            #print(f"\nobject RA={self.objra:.4f}, DEC={self.objdec:.4f}\n")
            self.xpixel,self.ypixel = wcs.wcs_world2pix(self.objra,self.objdec,0)
            
            # convert sma to pixels using pixel scale from mask wcs
            self.pixel_scale = wcs.pixel_scale_matrix[1][1]
            self.objsma_pixels = self.objsma/(self.pixel_scale*3600)
            

        self.weightim = weightim
        self.weight_threshold = weight_threshold
        ###  The lines below are for testing purposes
        ###  and should be removed before release.
        #if image is None:
        #    image='MKW8-18216-R.fits'
        #if haimage == None:
        #    haimage='MKW8-18216-CS.fits'
        if sepath is None:
            sepath=os.getenv('HOME')+'/github/halphagui/astromatic/'
        if gaiapath is None:
            gaiapath = os.getenv("HOME")+'/research/legacy/gaia-mask-dr9.virgo.fits'

        if config is None:
            config='default.sex.HDI.mask'
        self.image_name = image
        self.haimage_name = haimage
        print(self.image_name)
        print(self.haimage_name)
        print(sepath)
        self.sepath = sepath
        self.gaiapath = gaiapath
        self.gaia_mask = None
        self.add_gaia_stars = True        
        self.config = config
        self.threshold = threshold
        self.snr = snr
        self.snr_analysis = snr
        self.minarea = minarea
        self.cmap = cmap
        self.xcursor_old = -99
        self.xcursor = -99
        self.mask_size = 20.
        # create name for output mask file
        t = self.image_name.split('.fit')
        self.mask_image=t[0]+'-mask.fits'
        self.mask_inv_image=t[0]+'-inv-mask.fits'
        #print('saving mask image as: ',self.mask_image)

        self.remove_center_object_flag = True
        
        # read in image and define center coords
        self.image, self.imheader = fits.getdata(self.image_name,header = True)
        self.ymax,self.xmax = self.image.shape
        self.xc = self.xmax/2.
        self.yc = self.ymax/2.
        self.image_wcs = WCS(self.imheader)
        self.pscalex,self.pscaley = self.image_wcs.proj_plane_pixel_scales() # appears to be degrees/pixel
        
        # get image dimensions in deg,deg
        self.dxdeg,self.dydeg = imutils.get_image_size_deg(self.image_name)
        

        # Get coord of image center.  will use when getting gaia stars
        self.racenter,self.deccenter = imutils.get_image_center_deg(self.image_name)                


        self.v1,self.v2=scoreatpercentile(self.image,[5.,99.5])
        self.adjust_mask = True
        self.figure_size = (10,5)
        self.mask_size = 20. # side of square to mask out when user clicks on a pixel

        # set up array to store the user-created object masks

        self.usr_mask = np.zeros_like(self.image)
        #print(self.image.shape, self.usr_mask.shape)
        # set off center flag as false by default
        self.off_center_flag = False

        # keep track of extra objects that the user deletes from mask

        self.deleted_objects = []

        if not self.auto:
            self.add_cutout_frames()

        # time how long it takes to run SE
        self.runse_flag = True
        runphot = False
        if self.runse_flag:
            self.link_files()
            t_0 = timeit.default_timer()        
            self.runse(weightim=self.weightim,weight_threshold=self.weight_threshold)
            self.remove_center_object()
            #self.remove_central_objects(xc=self.xpixels,yc=self.ypixels)
            t_1 = timeit.default_timer()
            #print("HELLO!!!")
            print(f"\ntime to run se: {round((t_1-t_0),3)} sec\n")
        if runphot:
            self.usephot = True
            t_1 = timeit.default_timer()
            self.run_photutil()
            t_2 = timeit.default_timer()
            print(f"\ntime to run photutils: {round((t_2-t_1),3)} sec\n")
        #self.update_mask()
        if self.auto:
            for i in range(int(ngrow)):
                # grow mask 7x when running in auto mode
                self.grow_mask()
        try:
            self.show_mask_mpl()
        except TypeError:
            print("WARNING: could not display mask")
        if not self.auto:
            self.display_cutouts()
            self.connect_buttons()
            
    def connect_buttons(self):
        #self.ui.msaveButton.clicked.connect(self.write_mask)
        self.ui.mquitButton.clicked.connect(self.quit_program)
        self.ui.mhelpButton.clicked.connect(self.print_help_menu)
        self.ui.mrunSEButton.clicked.connect(self.runse)
        #self.ui.msaveButton.clicked.connect(self.save_mask)
        #self.ui.mremoveButton.clicked.connect(self.remove_object)
        self.ui.boxSizeLineEdit.textChanged.connect(self.set_box_size)
        self.ui.seThresholdLineEdit.textChanged.connect(self.set_threshold)
        self.ui.seSNRLineEdit.textChanged.connect(self.set_sesnr)
        self.ui.seSNRAnalysisLineEdit.textChanged.connect(self.set_sesnr_analysis)
    def close_window(self):
        print('click red x to close window')
        #sys.exit()
    def add_cutout_frames(self):
        # r-band cutout
        a = QtWidgets.QLabel('r-band')
        self.ui.cutoutsLayout.addWidget(a, 0, 0, 1, 1)
        a = QtWidgets.QLabel('CS Halpha')
        self.ui.cutoutsLayout.addWidget(a, 0, 1, 1, 1)
        a = QtWidgets.QLabel('Mask')
        self.ui.cutoutsLayout.addWidget(a, 0, 2, 1, 1)

        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.rcutout = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, 4, 1)
        self.hacutout = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 1, 4, 1)
        self.maskcutout = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,1, 2, 4, 1)
        #self.maskcutout.mouse_clicked.connect(self.add_object)

        # this allows the user to press editing keys in any of the 3 image panels
        # not just in the mask panel
        self.maskcutout.key_pressed.connect(self.key_press_func)
        self.rcutout.key_pressed.connect(self.key_press_func)
        self.hacutout.key_pressed.connect(self.key_press_func)

    def display_cutouts(self):
        self.rcutout.load_file(self.image_name)
        self.rcutout.fitsimage.set_autocut_params('stddev')
        if self.haimage_name is not None:
            self.hacutout.load_file(self.haimage_name)
        self.display_mask()
    def display_mask(self):
        self.maskcutout.load_file(self.mask_image)
        self.draw_central_ellipse()
    def show_mask(self):
        if self.nods9 & (not self.auto):
            plt.close('all')
            self.fig = plt.figure(1,figsize=self.figure_size)
            plt.clf()
            plt.subplots_adjust(hspace=0,wspace=0)
            plt.subplot(1,2,1)
            plt.imshow(self.image,cmap='gray_r',vmin=self.v1,vmax=self.v2,origin='lower')
            plt.title('image')
            plt.subplot(1,2,2)
            #plt.imshow(maskdat,cmap='gray_r',origin='lower')
            plt.imshow(self.maskdat,cmap=self.cmap,origin='lower')
            plt.title('mask')
            plt.gca().set_yticks(())
            #plt.draw()
            #plt.show(block=False)
            self.draw_central_ellipse()
    def draw_central_ellipse(self, color='cyan'): # MVC - view
        # mark r24
        markcolor=color#, 'yellow', 'cyan']
        markwidth=1
        #print('inside draw_ellipse_results')
        image_frames = [self.rcutout,self.hacutout,self.maskcutout]
        if self.ellipseparams is None:
            print("")
            print("no parameters found for central ellipse")
            print()
            return
        xc,yc,r,BA,PA = self.ellipseparams
        #print("just checking - adding ellipse drawing ",self.ellipseparams)
        objlist = []
        for i,im in enumerate(image_frames):
            obj =im.dc.Ellipse(xc,yc,r,r*BA, rot_deg = np.degrees(PA), color=markcolor,linewidth=markwidth)

            objlist.append(obj)
            self.markhltag = im.canvas.add(im.dc.CompoundObject(*objlist))
            im.fitsimage.redraw()
            #print("did you see anything???")
        # mark R17 in halpha image

    def key_press_func(self,text):
        key, x, y = text.split(',')
        self.xcursor = float(x)
        self.ycursor = float(y)
        try:
            self.cursor_value = self.maskdat[int(self.ycursor),int(self.xcursor)]
        except IndexError:
            print('out of bounds, try again')
        #print('cursor value = ',self.cursor_value, key)
        if key == 'c':
            self.add_circ_object()
        elif key == 'b':
            self.add_box_object()
        elif key == 'r': 
            print('removing object')
            self.remove_object(int(self.cursor_value))
        elif key == 'o': 
            self.off_center()
        elif key == 'g': 
            self.grow_mask()
        #elif key == 't': 
        #    self.set_threshold()
        #elif key == 'n': 
        #    self.set_sesnr()
        elif key == 'h': 
            self.print_help_menu()
        elif key == 'w': 
            self.save_mask()
        elif key == 'q':
            self.quit_program()
        else:
            print('did not understand that.  \n Try again!')
        
    def print_help_menu(self):
        print('Click on mask or r/ha image, then enter:\n \t r = remove object in mask at the cursor position;'
              '\n \t c = add CIRCULAR mask at cursor position;'
              '\n \t b = add BOX mask at cursor position;'
              '\n \t g = grow the size of the current masks;'              
              '\n \t o = if target is off center (and program is removing the wrong object);'
              #'\n \t s to change the size of the mask box;'
              #'\n \t t to adjust SE threshold (0=lots, 1=no deblend );'
              #'\n \t n to adjust SE SNR; '
              '\n \t h = print this menu; '
              '\n \t w = write the mask image;'
              #'\n \t q to quit \n \n'
              '\n\n'
              'Display shortcuts (click on image to adjust):'
              '\n \t scroll  = zoom'
              '\n \t `       = zoom to fit'
              '\n \t space+s = enable contrast adjustment, click+drag, scroll wheel'
              '\n \t space   = exit contrast adjustment'              
              '\n \t a       = automatically set contrast'              
              #'\n \t ALT-right_click = adjust contrast \n \n'
              '\n\nClick Red X to close window')


    def add_box_object(self):
        '''
        this adds a square region
        '''
        print('adding pixels to the mask')
        # mask out a rectangle around click
        # size is given by mask_size
        xmin = int(self.xcursor) - int(0.5*self.mask_size)
        ymin = int(self.ycursor) - int(0.5*self.mask_size)
        xmax = int(self.xcursor) + int(0.5*self.mask_size)
        ymax = int(self.ycursor) + int(0.5*self.mask_size)
        
        # make sure cursor click is not outside of the image
        if ((self.xcursor >= self.xmax) or (self.xcursor <= 0) or (self.ycursor >= self.ymax) or (self.ycursor <= 0)):
            print('you clicked outside the image area')
            return
        
        # make sure mask dimensions are not outside of the image
        xmin = max(0,xmin)
        xmax = min(self.xmax,xmax)
        ymin = max(0,ymin)
        ymax = min(self.ymax,ymax)

        #print('xcursor, ycursor = ',self.xcursor, self.ycursor)
        mask_value = np.max(self.maskdat) + 1
        #print(xmin,xmax,ymin,ymax,self.mask_size)
        self.usr_mask[ymin:ymax,xmin:xmax] = mask_value*np.ones([ymax-ymin,xmax-xmin])
        self.maskdat = self.maskdat + self.usr_mask
        self.save_mask()
        print('added mask object '+str(mask_value))

    def add_circ_object(self):
        print('adding circular obj to the mask, with radius = ',self.mask_size)
        # mask out a rectangle around click
        # size is given by mask_size
        pixel_mask = circle_pixels(float(self.xcursor),float(self.ycursor),float(self.mask_size/2.),self.xmax,self.ymax)

        #print('xcursor, ycursor = ',self.xcursor, self.ycursor)
        mask_value = int(np.max(self.maskdat)) + 1
        #print(f"adding circular mask with value {mask_value}"
        #print(xmin,xmax,ymin,ymax,self.mask_size)
        self.usr_mask[pixel_mask] = mask_value*np.ones_like(self.usr_mask)[pixel_mask]
        self.maskdat = self.maskdat + self.usr_mask
        self.save_mask()
        print(f'added mask object {mask_value}')
        
    def remove_object(self, objID):
        '''
        this will remove masked pixels near the cursor

        '''
        #objID = int(input('enter pixel value to remove object in mask'))
        xmin = int(self.xcursor) - int(0.5*self.mask_size)
        ymin = int(self.ycursor) - int(0.5*self.mask_size)
        xmax = int(self.xcursor) + int(0.5*self.mask_size)
        ymax = int(self.ycursor) + int(0.5*self.mask_size)

        if objID == 0:
            return
        else:
            self.maskdat[self.maskdat == objID] = 0.
            self.deleted_objects.append(objID)

            # remove object from user mask
            self.usr_mask[self.usr_mask == objID] = 0
        self.save_mask()
        self.display_mask()
    def set_threshold(self,t):
        '''
        adjust threshold used in SE deblending
         (0=lots, 1=no deblend )
        '''
        print('Adjust threshold for SE deblending')
        print('0=lots, 1=no deblend')
        #t = raw_input('enter new threshold')
        try:
            self.threshold = float(t)
            if self.runse_flag:
                self.runse(weightim=self.weightim,weight_threshold=self.weight_threshold)
            else:
                self.run_photutil()

        except ValueError:
            pass
        
    def set_sesnr(self,t):
        #t = raw_input('enter new SNR')
        try:
            self.snr = float(t)
        except ValueError:
            pass
        #self.runse()
                
    def set_sesnr_analysis(self,t):
        #t = raw_input('enter new SNR')
        try:
            self.snr_analysis = float(t)
        except ValueError:
            pass

                
    def set_box_size(self,t):
        # change box size used for adding pixels to mask
        #print('current box size = '+str(self.mask_size))
        #t = input('enter new size for square area to be masked (in pixels)\n')
        try:
            self.mask_size = float(t)
        except:
            print('error reading input')

    def off_center(self):
        t = input('enter object number for target galaxy\n')
        self.off_center_flag = True
        self.galaxy_id = int(t)
        if self.runse_flag:
            self.runse(weightim=self.weightim,weight_threshold=self.weight_threshold)
        else:
            self.run_photil()
    def quit_program(self):
        self.clean_links()
        self.close_window()

    def save_mask(self):
        #super(maskwindow,self).mask_saved(event)
        print('saving mask: ',self.mask_image)
        fits.writeto(self.mask_image, self.maskdat, header = self.imheader, overwrite=True)
        if not self.auto:
            self.mask_saved.emit(self.mask_image)
            self.display_mask()
        
            #print(self.mask_image)
            self.mask_saved.emit(self.mask_image)


    def edit_mask(self):
        if self.runse_flag:
            self.runse(weightim=self.weightim,weight_threshold=self.weight_threshold)
        else:
            self.run_photutil()
        while self.adjust_mask:    
            self.show_mask()
            self.print_menu()
            fits.writeto(self.mask_image,self.maskdat,header = self.imheader,overwrite=True)
            self.mask_saved.emit(self.mask_image)
            
        
# run sextractor on input image
# return segmentation image with central object removed



    


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

    
