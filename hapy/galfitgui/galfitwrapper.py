#!/usr/bin/env python

'''
GOAL:
- open gui window to run galfit
- show image, model, residual

PROCEDURE

'''


import os
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Tophat2DKernel
import numpy as np
#import argparse
#import pyds9
from scipy.stats import scoreatpercentile
from PyQt5 import QtCore, QtGui, QtWidgets
#from PyQt5 import  QtWidgets
#from ginga.qtw.QtHelp import QtGui
#from PyQt5 import QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
#from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
from ginga.mplw.ImageViewMpl import ImageView
from ginga import colors
from ginga.canvas.CanvasObject import get_canvas_types

from ginga.misc import log
from ginga.util.loader import load_data

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


#from galfitWidget import Ui_Form as Ui_galfitWindow
from galfitWidgetBD import Ui_Form as Ui_galfitWindow
from maskwrapper import my_cutout_image

import rungalfit as rg

# import imutils to get pixel scale
import imutils

# add a frame on top to show
# cutout, mask, model, residual

class galfitwindow(Ui_galfitWindow, QtCore.QObject):
    model_saved = QtCore.pyqtSignal(str)
    def __init__(self, MainWindow, logger, image=None, sigma_image=None, mask_image=None, psf=None,psf_oversampling=None, xmaxfit=None, ymaxfit=None, xminfit=1, yminfit=1, ncomp=1, convflag = True, convolution_size=None, fitallflag=False,xc=None, yc=None,mag=None,rad=None,nsersic=None, BA=None,PA=None, mag2=None, nsersic2=None, rad2=None, BA2=None, PA2=None, xc2=None, yc2=None, fitn=True, fitn2=True,asym=False,auto=False):
        super(galfitwindow, self).__init__()
        self.auto = auto        
        if MainWindow is None:
            self.auto = True

        if not self.auto:
            # boiler plate gui stuff
            # ... at least I think it is...
            self.ui = Ui_galfitWindow()
            self.ui.setupUi(MainWindow)
            MainWindow.setWindowTitle('GALFIT!')
            self.MainWindow = MainWindow

        # name of input file
        self.image_name = image
        # store image data and header
        self.image_data, self.image_header = fits.getdata(self.image_name, header=True)

        # define root name for galfit output file
        self.image_rootname = self.image_name.split('.fits')[0]
        
        # define fit area as the complete size of the input image
        # this assumes that the user is providing a cutout rather than a large mosaic image
        if xmaxfit == None:
            self.xmaxfit = self.image_data.shape[1]
        else:
            self.xmaxfit = xmaxfit

        if ymaxfit == None:
            self.ymaxfit = self.image_data.shape[0]
        else:
            self.ymaxfit = ymaxfit
        self.xminfit = xminfit
        self.yminfit = yminfit
        
        # try to extract the magnitude ZP from image header
        # if not available, set to 25.
        try:
            self.magzp = self.image_header['PHOTZP']
        except KeyError:
            self.magzp = 25.

        # get pixel scale from image header
        # convert from degrees/pix to arcsec/pix


        self.pscale = imutils.get_pixel_scale(self.image_header)

        
        self.sigma_image = sigma_image
        self.psf_image = psf
        self.psf_oversampling = psf_oversampling
        
        self.mask_image = mask_image
        if self.mask_image != None:
            self.mask_data = fits.getdata(self.mask_image)
            # create a masked array to display to the user
            self.image_data = np.ma.masked_array(self.image_data, mask=np.array(self.mask_data,'bool'))

        # set to false for fitting central object only
        self.fitallflag = fitallflag
        # number of components to fit
        # set to 1 for single component sersic fit
        self.ncomp = ncomp

        # asymmetry flag
        self.asym = asym
        ###########################################################3
        # enable psf convolution in fit
        ###########################################################3
        self.convflag = convflag
        if convolution_size == None:
            self.convolution_size = min(self.xmaxfit, self.ymaxfit)
            print('convolution size = ',self.convolution_size)
        else:
            self.convolution_size = convolution_size
            
        ###########################################################3
        # set up initial guesses for sersic parameters
        ###########################################################3
        if xc == None:
            self.xc=self.image_data.shape[0]/2
        else:
            self.xc=xc
        if yc == None:
            self.yc=self.image_data.shape[1]/2
        else:
            self.yc=yc
        if mag == None:
            self.mag=15.
        else:
            self.mag=mag
        if rad == None:
            self.re=10
        else:
            self.re=rad
        if nsersic == None:
            self.nsersic=2
        else:
            self.nsersic=nsersic
        if BA == None:
            self.BA= 0.6
        else:
            self.BA= BA
        if PA == None:
            self.PA=0
        else:
            self.PA=PA

        self.fitn = fitn
        # parameters for second component B/D fit
        if mag2 == None:
            self.mag2=16
        else:
            self.mag2=mag2

        if nsersic2 == None:
            self.nsersic2=4
        else:
            self.nsersic2=nsersic2
            
        if rad2 == None:
            self.re2=0.5*self.re
        else:
            self.re2=rad2
            
        if PA2 == None:
            self.PA2=0
        else:
            self.PA2=PA2

        if BA2 == None:
            self.BA2= 1
        else:
            self.BA2=BA2

        if xc2 == None:
            self.xc2= self.xc
        else:
            self.xc2=xc2

        if yc2 == None:
            self.yc2= self.yc
        else:
            self.yc2=yc2

        self.fitn2 = fitn2
        
        self.fitBA = 1
        self.fitPA = 1

        # define galfit output image the same way that rungalfit does

        if self.asym:
            self.output_image=self.image_rootname+'-'+ str(self.ncomp) +'Comp-galfit-out-asym.fits'
        else:
            self.output_image=self.image_rootname+'-'+ str(self.ncomp) +'Comp-galfit-out.fits'
        # continue with other functions
        self.logger = logger


        ############################################################
        # set up gui
        ############################################################
        if not self.auto:
            self.add_cutout_frames()
            self.connect_buttons()
            self.display_image()
            self.display_initial_params()

        ############################################################
        # create an instance of rungalfit.galfit
        ############################################################        
        self.initialize_galfit()        

        ############################################################
        # the remaining commands come from user input through the gui
        ############################################################        

        if self.auto:
            self.run_galfit(fitBA=self.fitBA, fitPA=self.fitPA)
    def add_cutout_frames(self):
        # gui stuff
        # set up text labels for image, model, and residual
        a = QtWidgets.QLabel('Image')
        self.ui.cutoutsLayout.addWidget(a, 0, 0, 1, 1)
        a = QtWidgets.QLabel('Model')
        self.ui.cutoutsLayout.addWidget(a, 0, 1, 1, 1)
        a = QtWidgets.QLabel('Residual')
        self.ui.cutoutsLayout.addWidget(a, 0, 2, 1, 1)

        # gui stuff
        # set up image frames for image, model, and residual
        self.cutout_frame = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, 4, 1,autocut_params='sttdev')
        self.model_frame = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 1, 4, 1)
        self.residual_frame = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,1, 2, 4, 1)

        # connect key press function to each of the image frames
        # so that the user can call them with cursor over any image, including mainframe
        self.cutout_frame.key_pressed.connect(self.key_press_func)
        self.model_frame.key_pressed.connect(self.key_press_func)
        self.residual_frame.key_pressed.connect(self.key_press_func)
        
    def display_image(self):
        '''
        # I set this up initially to show the masked image
        # but this is making the displayed not so nice
        if self.mask_image == None:
            self.cutout_frame.load_file(self.image)
        else:
            # display image with masked values
            # if a mask is available
            self.cutout_frame.load_image(self.image_data)
        '''
        self.cutout_frame.load_file(self.image_name)
    def display_galfit_results(self):
        pass
    def connect_buttons(self):
        self.ui.quitButton.clicked.connect(self.quit_program)
        self.ui.helpButton.clicked.connect(self.print_help_menu)
        # this is the work horse function
        self.ui.runGalfitButton.clicked.connect(lambda: self.run_galfit(fitBA=self.fitBA, fitPA=self.fitPA))
        self.ui.xcLineEdit.textChanged.connect(self.set_xc)
        self.ui.ycLineEdit.textChanged.connect(self.set_yc)
        self.ui.magLineEdit.textChanged.connect(self.set_mag)
        self.ui.ReLineEdit.textChanged.connect(self.set_Re)
        self.ui.nLineEdit.textChanged.connect(self.set_nsersic)
        self.ui.PALineEdit.textChanged.connect(self.set_PA)
        self.ui.BALineEdit.textChanged.connect(self.set_BA)
 
    def display_initial_params(self):

        self.ui.xcLineEdit.setText(str(round(self.xc,2)))
        self.ui.ycLineEdit.setText(str(round(self.yc,2)))
        self.ui.magLineEdit.setText(str(round(self.mag,2)))
        self.ui.ReLineEdit.setText(str(round(self.re,2)))
        self.ui.nLineEdit.setText(str(round(self.nsersic,2)))
        self.ui.PALineEdit.setText(str(round(self.PA,2)))
        self.ui.BALineEdit.setText(str(round(self.BA,2)))
        if self.ncomp == 2:
            self.ui.xc2LineEdit.setText(str(round(self.xc2,2)))
            self.ui.yc2LineEdit.setText(str(round(self.yc2,2)))
            self.ui.mag2LineEdit.setText(str(round(self.mag2,2)))
            self.ui.Re2LineEdit.setText(str(round(self.re2,2)))
            self.ui.n2LineEdit.setText(str(round(self.nsersic2,2)))
            self.ui.PA2LineEdit.setText(str(round(self.PA2,2)))
            self.ui.BA2LineEdit.setText(str(round(self.BA2,2)))
            
    def display_fitted_params(self):
        self.ui.xcFitLineEdit.setText(str(round(self.xc,2)))
        self.ui.ycFitLineEdit.setText(str(round(self.yc,2)))
        self.ui.magFitLineEdit.setText(str(round(self.mag,2)))
        self.ui.ReFitLineEdit.setText(str(round(self.re,2)))
        self.ui.nFitLineEdit.setText(str(round(self.nsersic,2)))
        self.ui.PAFitLineEdit.setText(str(round(self.PA,2)))
        self.ui.BAFitLineEdit.setText(str(round(self.BA,2)))
        self.ui.skyFitLineEdit.setText(str(round(self.sky,2)))
        self.ui.errorFitLineEdit.setText(str(round(self.error,2)))
        try:
            self.ui.chiFitLineEdit.setText(str(round(self.chi2nu,2)))
        except TypeError:
            print("WARNING: could not update display for chinu ",self.chi2nu)
        if self.asym:
            self.ui.asymLineEdit.setText(str(round(self.asymmetry,2)))
            self.ui.asymPALineEdit.setText(str(round(self.asymmetry_PA,2)))
        if self.ncomp == 2:
            self.ui.xc2FitLineEdit.setText(str(round(self.xc2,2)))
            self.ui.yc2FitLineEdit.setText(str(round(self.yc2,2)))
            self.ui.mag2FitLineEdit.setText(str(round(self.mag2,2)))
            self.ui.Re2FitLineEdit.setText(str(round(self.re2,2)))
            self.ui.n2FitLineEdit.setText(str(round(self.nsersic2,2)))
            self.ui.PA2FitLineEdit.setText(str(round(self.PA2,2)))
            self.ui.BA2FitLineEdit.setText(str(round(self.BA2,2)))
            self.ui.chiFitLineEdit.setText(str(round(self.chi2nu,2)))
       
    def set_xc(self,dat):
        try:
            self.xc = float(dat)
        except ValueError:
            pass
    def set_yc(self,dat):
        try:
            self.yc = float(dat)
        except ValueError:
            pass
    def set_mag(self,dat):
        try:
            self.mag = float(dat)
        except ValueError:
            pass
    def set_Re(self,dat):
        try:
            self.re = float(dat)
        except ValueError:
            pass
    def set_nsersic(self,dat):
        try:
            self.nsersic = float(dat)
        except ValueError:
            pass
    def set_PA(self,dat):
        try:
            self.PA = float(dat)
        except ValueError:
            pass
    def set_BA(self,dat):
        try:
            self.BA = float(dat)
        except ValueError:
            pass


    def key_press_func(self,text):
        key, x, y = text.split(',')
        self.xcursor = float(x)
        self.ycursor = float(y)
        '''
        if key == 'n':
            self.galfit.set_n()
        elif key == 'r':
            self.galfit.set_r()
        '''
        if key == 'o':
            self.galfit.toggle_fitall()
        elif key == 'c':
            if ((self.xcursor > self.image_data.shape[0]) | (self.ycursor > self.image_data.shape[0]) | (self.xcursor < 1.) | (self.ycursor < 1.)):
                print('cursor value out of range.  put mouse over desired center and type c')
                
            self.galfit.set_center(self.xcursor, self.ycursor)
        elif key == 'f':
            self.galfit.print_fix_menu()
        elif key == 'a':
            self.galfit.toggle_asymmetry()
        elif key == 'R':
            self.galfit.reset_sersic_params()
        elif key == 'g':
            self.run_galfit()
        elif key == 'h': 
            self.print_help_menu()
        elif key == 'x':
            self.quit_program()
        else:
            print('did not understand that.  \n Try again!')

    def print_help_menu(self):
        print('What is wrong?\n \
        \n o = nearby object (toggle fitall)\
        \n b = B/A \n p = PA \n m = mag\
        \n c = recenter\
        \n f = hold values fixed\
        \n a = toggle asymmetry parameter\
        \n R = reset to original values\
        \n g = go (run galfit)\
        \n x=quit \n \
        \nDisplay shortcuts (click on image to adjust):'+
        '\n \t scroll = zoom'+
        '\n \t ` = zoom to fit'+
        '\n \t ALT-right_click = adjust contrast \n \n'+
        'Click Red X to close window')
    def quit_program(self):
        #self.close_window()
        print('click red x in top left corner to quit')
        print('lame...I know... but at least the program did not crash!')
        
#class galfit_galaxy():


    def initialize_galfit(self,convflag=True):
        '''
        GOAL: Preparing file to be run in galfit. Initialize galfit image parameters 

        INPUT: nsaid 

        OUTPUT: A definition of everything from galname to convflag, necessary for running galfit

        '''
        print('self.psfimage = ',self.psf_image)

        ## establish instance of rungalfit.galfit
        self.galfit = rg.galfit(galname=self.image_rootname,image=self.image_name,
                                mask_image = self.mask_image,
                                sigma_image=self.sigma_image,psf_image=self.psf_image,
                                psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,
                                yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,
                                convolution_size=self.convolution_size,magzp=self.magzp,
                                pscale=self.pscale,ncomp=self.ncomp,convflag=convflag,
                                fitallflag = self.fitallflag,asym=self.asym)
        
    def run_galfit(self,fitBA=1,fitPA=1):
        if self.ncomp == 1:
            self.galfit.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
        #print('in galfitwrapper, run_galfit: fitBA = ',fitBA)
            self.galfit.set_sky(0)
            self.galfit.run_galfit()
        elif self.ncomp == 2:
            self.galfit.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
            self.galfit.set_sersic_params_comp2(xobj=self.xc2,yobj=self.yc2,mag=self.mag2,rad=self.re2,nsersic=self.nsersic2,BA=self.BA2,PA=self.PA2,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
            self.galfit.set_sky(0)
            self.galfit.run_galfit()
            
        self.get_galfit_results(printflag=True)            
        if not self.auto:
            self.display_results()
            self.display_fitted_params()
            #self.galfit.close_input_file()
            #self.galfit.print_params()
            #self.galfit.print_galfit_results()
            
            # send message that a new model is available
            self.model_saved.emit(str(self.ncomp))
        
    def get_galfit_results(self,printflag = False):
        '''
        GOAL: Grab results from galfit (xc, yc, mag, re, nsersic, BA, PA, sky, error, chi2nu) and parse them into self.filename

        INPUT: nsaid

        OUTPUT: 1Comp file of output from galfit, header for the output file

        '''


        t = rg.parse_galfit_results(self.output_image, ncomp = self.ncomp, asymflag=self.asym)
        if printflag:
            self.galfit.print_galfit_results(self.output_image)
        self.galfit_results = t
        print(t)
        # for 1 comp fit
        # header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        # for 2 component fit
        # header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_SKY','CHI2NU']
        self.xc, self.xc_err = t[0]
        self.yc, self.yc_err = t[1]
        self.mag, self.mag_err = t[2]
        self.re, self.re_err = t[3]
        self.nsersic, self.nsersic_err = t[4]
        self.BA, self.BA_err = t[5]
        self.PA, self.PA_err = t[6]
        if self.asym:
            self.sky, self.sky_err = t[7]
            self.asymmetry, self.asymmetry_err = t[8]
            self.asymmetry_PA, self.asymmetry_PA_err = t[9]
            self.error = t[10]
            self.chi2nu = t[11]
            i_next = 12
        else:
            i_next = 7
            self.sky, self.sky_err = t[i_next]
            self.error = t[i_next+1]
            self.chi2nu = t[i_next+2]

        if self.ncomp == 2:
            self.xc2, self.xc2_err = t[7]
            self.yc2, self.yc2_err = t[8]
            self.mag2, self.mag2_err = t[9]
            self.re2, self.re2_err = t[10]
            self.nsersic2, self.nsersic2_err = t[11]
            self.BA2, self.BA2_err = t[12]
            self.PA2, self.PA2_err = t[13]
            self.sky, self.sky_err = t[14]
            self.error = t[15]
            self.chi2nu = t[16]


    def display_results(self):
        self.model_data = fits.getdata(self.output_image,2)
        self.residual_data = fits.getdata(self.output_image,3)
        self.model_frame.load_image(self.model_data)
        self.residual_frame.load_image(self.residual_data)
        pass
        ### THIS SHOULD BE MOVED TO PARENT PROGRAM
        ### BECAUSE HOW U DISPLAY RESULTS WILL VARY WITH USAGE
    def closeEvent(self, event):
        # send signal that window is closed


        event.accept()
        
if __name__ == "__main__":
    #catalog = '/Users/rfinn/research/NSA/nsa_v0_1_2.fits'
    #gcat = galaxy_catalog(catalog)
    #from halphamain import cutout_image
    #from halphamain import cutout_image
    logger = log.get_logger("galfitlog", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    MainWindow = QtWidgets.QWidget()
    ui = galfitwindow(MainWindow, logger, image='MKW8-18216-R.fits', mask_image = 'MKW8-18216-R-mask.fits',
                      psf='MKW8_R.coadd-psf.fits',psf_oversampling=2.)
    #ui.setupUi(MainWindow)
    #ui.test()

    MainWindow.show()
    sys.exit(app.exec_())

    
