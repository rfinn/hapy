from PyQt5 import QtCore, QtWidgets
from astropy.io import fits
from astropy.wcs import WCS

from .maskWidget import Ui_Form as Ui_maskWindow
from .cutout_view import my_cutout_image

#from hapy.masktools.engine import MaskEngine
from hapy.masktools.api import MaskEngine

class MaskWindow(Ui_maskWindow, QtCore.QObject):
    mask_saved = QtCore.pyqtSignal(str)
    def __init__(self, MainWindow, logger, image=None, haimage=None, sepath=None, gaiapath=None, config=None, threshold=0.005,snr=10,cmap='gist_heat_r',objparams=None,auto=False,unmaskellipse=False,minarea=10,ngrow=3,weightim=None,weight_threshold=None):
        """

        ngrow : number of times to run grow when running in auto mode
        """

        self.engine = MaskEngine(
        image_fits=image,
        ha_image_fits=haimage,
        sepath=sepath,
        gaiapath=gaiapath,
        config=config,
        threshold=threshold,
        snr=snr,
        minarea=minarea,
        add_gaia_stars=True,)

        def _progress_cb(stage: str, fraction: float, message: str = None):
            print(f"[mask] {stage} {fraction:0.2f} {message or ''}".strip())
            
        self.engine.build_initial_mask(weightim=weightim, weight_threshold=weight_threshold)
        self.maskdat = self.engine.maskdat

        
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
            
