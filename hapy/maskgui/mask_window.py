from PyQt5 import QtCore, QtWidgets
from astropy.io import fits
from astropy.wcs import WCS
import os
from scipy.stats import scoreatpercentile
import numpy as np
from pathlib import Path
from .maskWidget import Ui_Form as Ui_maskWindow
from .cutout_view import CutoutPanel

#from hapy.masktools.engine import MaskEngine
from hapy.masktools.api import MaskEngine, EllipseParams
from hapy.imagetools import imutils

from matplotlib import pyplot as plt

class NullLogger:
    def debug(self, *args, **kwargs): pass
    def info(self, *args, **kwargs): pass
    def warning(self, *args, **kwargs): pass
    def error(self, *args, **kwargs): pass

class MaskWindow(Ui_maskWindow, QtCore.QObject):
    mask_saved = QtCore.pyqtSignal(str)
    def __init__(self, MainWindow, logger,
             image=None,
             haimage=None,
             sepath="sex",
             config=None,
             threshold=0.005,
             snr=10,
             minarea=5,
             objparams=None,
             auto=False,
             ngrow=3,
             weightim=None,
             weight_threshold=None,
             addgaia=True,
             cmap="gist_heat_r"):

        """

        ngrow : number of times to run grow when running in auto mode
        """
        if logger is None:
            logger = NullLogger()
        self.logger = logger
        
        self.auto = auto
        if MainWindow is None:
            self.auto = True

        self.sepath = sepath

        self._syncing_views = False
        self._last_active_panel = None
        #########################################
        # WINDOW MAGIC
        #########################################       

        if not self.auto:
            self.MainWindow = MainWindow
            super(MaskWindow, self).__init__()
            self.ui = Ui_maskWindow()
            self.ui.setupUi(MainWindow)
            QtCore.QTimer.singleShot(0, self._debug_sizes)
            print("cutouts frameShape:", self.ui.cutouts.frameShape())
            print("dummyWidget exists?", hasattr(self.ui, "dummyWidget"))
            if hasattr(self.ui, "dummyWidget"):
                print("dummyWidget geom:", self.ui.dummyWidget.geometry())
            print("cutouts geom:", self.ui.cutouts.geometry())
            # Give the cutouts area (rows 0-7, cols 0-2) the space it needs
            self.ui.gridLayout_2.setRowStretch(0, 1)
            self.ui.gridLayout_2.setRowStretch(1, 1)
            self.ui.gridLayout_2.setRowStretch(2, 1)
            self.ui.gridLayout_2.setRowStretch(3, 1)
            self.ui.gridLayout_2.setRowStretch(4, 1)
            self.ui.gridLayout_2.setRowStretch(5, 1)
            self.ui.gridLayout_2.setRowStretch(6, 1)
            self.ui.gridLayout_2.setRowStretch(7, 1)

            #self.ui.gridLayout_2.setColumnStretch(0, 3)
            #self.ui.gridLayout_2.setColumnStretch(1, 3)
            #self.ui.gridLayout_2.setColumnStretch(2, 3)
            #self.ui.gridLayout_2.setColumnStretch(3, 0)            

            # Optional but very effective: keep control area from stealing vertical space
            self.ui.gridLayout_2.setRowStretch(11, 0)

            # also remove any minimum width those spacer widgets might impose
            for name in ["widget_6","widget_14","widget_15","widget_8","widget_10","widget_7","widget_9","widget_2","widget_5"]:
                if hasattr(self.ui, name):
                    w = getattr(self.ui, name)
                    w.setMinimumWidth(0)
                    w.setMaximumWidth(0)
                    w.hide()
            # Ensure the cutouts frame itself is allowed to expand
            self.ui.cutouts.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            self.ui.cutouts.setMinimumHeight(400)

            MainWindow.setWindowTitle('makin a mask...')
            self.MainWindow = MainWindow
            
        #########################################
        # INITIALIZE VARIABLES 
        #########################################
        self.r_overlay_mode = "off"
        self.ha_overlay_mode = "off"
        self.image_name = image



        image_path = Path(self.image_name)
        stem = image_path.stem

        if stem.endswith("-R") or stem.endswith("-r"):
            stem = stem[:-2]
            
        self.mask_image = str(image_path.with_name(f"{stem}-mask-manual.fits"))
        self.mask_inv_image = str(image_path.with_name(f"{stem}-inv-mask-manual.fits"))

        #self.mask_image = self.image_name.replace('-R.fits','-mask-manual.fits')
        #self.mask_inv_image=self.image_name.replace('-R.fits','-inv-mask-manual.fits')
        
        self.haimage_name = haimage
        self.weightim = weightim
        self.weight_threshold = weight_threshold
        
        self.gaia_mask = None
        self.add_gaia_stars = addgaia

        self.threshold = threshold
        self.snr = snr
        self.snr_analysis = snr
        self.minarea = minarea
        
        self.cmap = cmap
        self.xcursor_old = -99
        self.xcursor = -99
        self.mask_size = 20.

        if config is None:
            config = 'default.sex.HDI.mask'

        # create name for output mask file

        self.remove_center_object_flag = True

        # don't think I need this anymore - holdover from prior structure
        ## read in image and define center coords

        #self.ymax,self.xmax = self.image.shape
        #self.xc = self.xmax/2.
        #self.yc = self.ymax/2.
        #self.image_wcs = WCS()
        #self.pscalex,self.pscaley = self.image_wcs.proj_plane_pixel_scales() # appears to be degrees/pixel
        
        # get image dimensions in deg,deg
        #self.dxdeg,self.dydeg = imutils.get_image_size_deg(self.image_name)
        

        # Get coord of image center.  will use when getting gaia stars
        #self.racenter,self.deccenter = imutils.get_image_center_deg(self.image_name)                

        ###################################################
        self.image = fits.getdata(image)
        self.ymax, self.xmax = self.image.shape

        if self.haimage_name is not None:
            self.haimage = fits.getdata(self.haimage_name)

        #########################################
        # INITIALIZE VARIABLES 
        #########################################       
        if objparams is not None:
            self.objra = objparams[0]
            self.objdec = objparams[1]
            self.objsma = objparams[2]
            self.objBA = objparams[3]
            self.objPA = objparams[4]

            # unmask central elliptical region around object
            # get wcs from mask image
            wcs = WCS(fits.getheader(image))
            
            # get x and y coord of galaxy from (RA,DEC) using mask wcs
            #print(f"\nobject RA={self.objra:.4f}, DEC={self.objdec:.4f}\n")
            self.xpixel,self.ypixel = wcs.wcs_world2pix(self.objra,self.objdec,0)
            
            # convert sma to pixels using pixel scale from mask wcs
            self.pixel_scale = wcs.pixel_scale_matrix[1][1]
            self.objsma_pixels = self.objsma/(self.pixel_scale*3600)
            print("Found ellipse parameters !\n")
            self.ellipseparams = EllipseParams(
                xc = self.xpixel, 
                yc = self.ypixel,
                sma_pix = self.objsma_pixels,
                ba = self.objBA,
                theta_deg= self.objPA)
        else:
            print("DID NOT FIND ellipse parameters !\n")
            self.ellipseparams = None
            self.objra  = None  
            self.objdec = None  
            self.objsma = None  
            self.objBA  = None  
            self.objPA  = None
            

        self.mask_size = 20. # side of square to mask out when user clicks on a pixel

        # set up array to store the user-created object masks
        # don't do this here - do it in engine!
        #self.usr_mask = np.zeros_like(self.image)
        #print(self.image.shape, self.usr_mask.shape)
        # set off center flag as false by default
        self.off_center_flag = False

        # keep track of extra objects that the user deletes from mask

        self.deleted_objects = []


        #########################################        
        # START YOUR ENGINE!
        # initialize and build first-pass mask
        #########################################

        
        self.config = config
        self.add_gaia_stars = addgaia
        self.ngrow = ngrow


        self.engine = MaskEngine(
            image_fits=image,
            ha_image_fits=haimage,
            sepath=self.sepath,
            config=self.config,
            threshold=self.threshold,
            snr=self.snr,
            snr_analysis=self.snr_analysis,
            minarea=self.minarea,
            add_gaia_stars=self.add_gaia_stars,
            logger=self.logger,
        )
        

        def _progress_cb(stage, fraction, message=None):
            print(f"[mask] {stage} {fraction:0.2f} {message or ''}".strip())
        

        # IMPORTANT: pass progress_callback, and only build once
        self.maskdat = self.engine.build_initial_mask(
            weightim=weightim,
            weight_threshold=weight_threshold,
            progress_callback=_progress_cb,
            galaxy_ellipse = self.ellipseparams,
        )
        

        self.engine.write_mask(self.mask_image)
        
        #self.show_mask_mpl()

        #########################################        
        # ADD CUTOUT PANELS TO GUI
        #########################################        
        

        self.add_cutout_frames()
        
        #print("rcutout exists?", hasattr(self, "rcutout"))
        #print("maskcutout exists?", hasattr(self, "maskcutout"))
        #if hasattr(self, "rcutout"):
        #    print("rcutout parent:", self.rcutout.parent())

        #print("cutouts frame size:", self.ui.cutoutsLayout.parentWidget().size())
        #print("rcutout widget size:", self.rcutout.widget.size())
        #print("maskcutout widget size:", self.maskcutout.widget.size())
        QtCore.QTimer.singleShot(0, self.display_cutouts)
        #self.display_cutouts()
        self.connect_buttons()
        frame = self.ui.cutoutsLayout.parentWidget()   # the QFrame that contains the grid
        frame.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        frame.setMinimumHeight(200)  # pick something reasonable

        # make the grid stretch properly
        self.ui.cutoutsLayout.setColumnStretch(0, 1)
        self.ui.cutoutsLayout.setColumnStretch(1, 1)
        self.ui.cutoutsLayout.setColumnStretch(2, 1)
        self.ui.cutoutsLayout.setRowStretch(1, 1)
    def runse(self):
        """
        Legacy button name: 'Run SE'.
        Rebuild mask using current GUI settings and refresh displays.
        """
        def _progress_cb(stage, fraction, message=None):
            print(f"[mask] {stage} {fraction:0.2f} {message or ''}".strip())

        # keep engine config in sync with GUI state
        self.engine.threshold = self.threshold
        self.engine.snr = self.snr
        self.engine.snr_analysis = getattr(self, "snr_analysis", self.snr)
        self.engine.minarea = self.minarea
        self.engine.config = self.config
        self.engine.add_gaia_stars = self.add_gaia_stars

        self.maskdat = self.engine.build_initial_mask(
            weightim=self.weightim,
            weight_threshold=self.weight_threshold if self.weight_threshold is not None else 1.0,
            progress_callback=_progress_cb,
            galaxy_ellipse=self.ellipseparams,
            grow_size=int(self.mask_size),
            grow_iterations=int(self.ngrow),
        )

        self.engine.write_mask(self.mask_image)

        if not self.auto:
            self.display_mask()
        

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
    def _make_panel_label(self, text,fontsize):
        label = QtWidgets.QLabel(text)

        font = label.font()
        font.setPointSize(fontsize)
        font.setBold(True)
        label.setFont(font)

        label.setAlignment(QtCore.Qt.AlignCenter)

        # Let it expand horizontally but stay compact vertically
        label.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Fixed
        )

        return label        
    def add_cutout_frames(self):

        self.ui.cutoutsLayout.addWidget(self._make_panel_label("r-band",16), 0, 0)
        self.ui.cutoutsLayout.addWidget(self._make_panel_label("CS Halpha",16), 0, 1)
        self.ui.cutoutsLayout.addWidget(self._make_panel_label("Mask",16), 0, 2)

        self.rcutout = CutoutPanel(self.ui.cutoutsLayout,self.ui, self.logger, grid_pos=[1, 0, 4, 1])
        self.hacutout = CutoutPanel(self.ui.cutoutsLayout,self.ui, self.logger, grid_pos=[1, 1, 4, 1])
        self.maskcutout = CutoutPanel(self.ui.cutoutsLayout,self.ui, self.logger, grid_pos=[1, 2, 4, 1])
        
        #self.maskcutout.key_pressed.connect(self.key_press_func)
        #self.rcutout.key_pressed.connect(self.key_press_func)
        #self.hacutout.key_pressed.connect(self.key_press_func)
        self.rcutout.key_pressed.connect(lambda text: self.key_press_func(text, "r"))
        self.hacutout.key_pressed.connect(lambda text: self.key_press_func(text, "ha"))
        self.maskcutout.key_pressed.connect(lambda text: self.key_press_func(text, "mask"))

 

    def get_viewers(self):
        viewers = {
            "r": self.rcutout,
            "ha": self.hacutout,
            "mask": self.maskcutout,
        }
        return {k: v for k, v in viewers.items() if v is not None}


    def sync_pan_to(self, x: float, y: float, source_panel: str | None = None) -> None:
        if self._syncing_views:
            return

        self._syncing_views = True
        try:
            for name, viewer in self.get_viewers().items():
                if source_panel is not None and name == source_panel:
                    continue
                try:
                    viewer.set_pan_xy(x, y)
                except Exception as e:
                    print(f"Could not sync pan for {name}: {e}")
        finally:
            self._syncing_views = False


    def sync_zoom_from(self, source_panel: str) -> None:
        if self._syncing_views:
            return

        viewers = self.get_viewers()
        source = viewers.get(source_panel)
        if source is None:
            return

        self._syncing_views = True
        try:
            try:
                zoom = source.get_zoom_level()
            except Exception as e:
                print(f"Could not get zoom from {source_panel}: {e}")
                return

            for name, viewer in viewers.items():
                if name == source_panel:
                    continue
                try:
                    viewer.set_zoom_level(zoom)
                except Exception as e:
                    print(f"Could not sync zoom for {name}: {e}")
        finally:
            self._syncing_views = False


    def sync_view_from(self, source_panel: str) -> None:
        """
        Sync both pan and zoom from one panel to the others.
        """
        if self._syncing_views:
            return

        viewers = self.get_viewers()
        source = viewers.get(source_panel)
        if source is None:
            return

        self._syncing_views = True
        try:
            pan_x = pan_y = None
            zoom = None

            try:
                pan_x, pan_y = source.get_pan_xy()
            except Exception as e:
                print(f"Could not get pan from {source_panel}: {e}")

            try:
                zoom = source.get_zoom_level()
            except Exception as e:
                print(f"Could not get zoom from {source_panel}: {e}")

            for name, viewer in viewers.items():
                if name == source_panel:
                    continue

                if pan_x is not None and pan_y is not None:
                    try:
                        viewer.set_pan_xy(pan_x, pan_y)
                    except Exception as e:
                        print(f"Could not sync pan for {name}: {e}")

                if zoom is not None:
                    try:
                        viewer.set_zoom_level(zoom)
                    except Exception as e:
                        print(f"Could not sync zoom for {name}: {e}")
        finally:
            self._syncing_views = False
        
    def display_cutouts(self):
        self.rcutout.load_file(self.image_name)
        self.rcutout.fitsimage.set_autocut_params('stddev')

        
        if self.haimage_name is not None:
            self.hacutout.load_file(self.haimage_name)
        self.display_mask()
        print("loading r:", self.image_name)
        print("loading ha:", self.haimage_name)
        print("loading mask:", self.mask_image)

        try:
            self.rcutout.fitsimage.zoom_fit()
            self.rcutout.fitsimage.redraw(whence=0)

            self.maskcutout.fitsimage.zoom_fit()
            self.maskcutout.fitsimage.redraw(whence=0)

            if self.haimage_name:
                self.hacutout.fitsimage.zoom_fit()
                self.hacutout.fitsimage.redraw(whence=0)
        except Exception as e:
            print("redraw failed:", e)

        print("cutouts frame size:", self.ui.cutouts.size())

    
    def display_mask(self):
        self.maskcutout.load_file(self.mask_image)
        self.draw_central_ellipse()

        # NEW: refresh overlays
        self.update_overlay("r")
        self.update_overlay("ha")

    def _refresh(self):
        self.maskdat = self.engine.maskdat
        self.display_mask()
    
    def draw_central_ellipse(self, color='cyan'): # MVC - view
        # mark r24
        markcolor=color#, 'yellow', 'cyan']
        markwidth=1

        image_frames = [self.rcutout,self.hacutout,self.maskcutout]
        if self.ellipseparams is None:
            print("")
            print("no parameters found for central ellipse")
            print()
            return
        xc,yc = self.ellipseparams.xc, self.ellipseparams.yc
        r = self.ellipseparams.sma_pix
        BA = self.ellipseparams.ba
        PA = self.ellipseparams.theta_deg

        objlist = []
        for i,im in enumerate(image_frames):
            obj =im.dc.Ellipse(xc,yc,r,r*BA, rot_deg = PA, color=markcolor,linewidth=markwidth)

            objlist.append(obj)
            self.markhltag = im.canvas.add(im.dc.CompoundObject(*objlist))
            im.fitsimage.redraw()

    def print_cursor_readout(self, panel: str) -> None:

        ix = int(round(self.xcursor))
        iy = int(round(self.ycursor))


        if not (0 <= ix < self.xmax and 0 <= iy < self.ymax):
            print(f"[{panel}] x={ix} y={iy} is out of bounds")
            return


        rval = self.image[iy, ix] if self.image is not None else np.nan

        haval = np.nan
        if getattr(self, "haimage", None) is not None:
            try:
                haval = self.haimage[iy, ix]
            except Exception:
                haval = np.nan

        maskval = np.nan
        mask_area = None
        if self.engine is not None and self.engine.maskdat is not None:
            maskval = self.engine.maskdat[iy, ix]
            if np.isfinite(maskval) and maskval > 0:
                mask_area = int(np.sum(self.engine.maskdat == maskval))

        msg = (
            f"[{panel:4s}] x={ix:4d} y={iy:4d}  "
            f"r={rval:8.4g}  "
            f"ha={haval:8.4g}  "
            f"mask={maskval}"
        )
        if mask_area is not None:
            msg += f"  npix={mask_area}"
        print(msg)

    

    
    def key_press_func(self,text,panel):
        key, x, y = text.split(',')
        self.xcursor = float(x)
        self.ycursor = float(y)
        self.sync_pan_to(self.xcursor, self.ycursor, source_panel=panel)
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
        elif key == 'v':
            self.print_cursor_readout(panel)
        elif key == 'g': 
            self.grow()
        #elif key == 't': 
        #    self.set_threshold()
        #elif key == 'n': 
        #    self.set_sesnr()
        elif key == 'h': 
            self.print_help_menu()
        elif key == "m" and panel in ("r", "ha"):
            self.cycle_overlay_mode(panel)
        elif key == "z":
            self.sync_view_from(panel)
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
              '\n \t v = print pixel values at cursor position;'
              '\n \t m = toggle mask view on current image (outline, filled, none);'
              "\t z = sync pan/zoom from active panel to the other panels\n"
              #'\n \t s to change the size of the mask box;'
              #'\n \t t to adjust SE threshold (0=lots, 1=no deblend );'
              #'\n \t n to adjust SE SNR; '
              '\n \t h = print this menu; '
              '\n \t w = write the mask image;'
              '\n \t q = quit;'
              '\n\n'
              'Display shortcuts (click on image to adjust):'
              '\n \t scroll  = zoom'
              '\n \t `       = zoom to fit'
              '\n \t space+s = enable contrast adjustment, click+drag, scroll wheel'
              '\n \t space   = exit contrast adjustment'              
              '\n \t a       = automatically set contrast'              
              #'\n \t ALT-right_click = adjust contrast \n \n'
              #'\n\nClick Red X to close window')
              )
    def get_boolean_mask(self):
        """
        Return the current mask as a boolean array.
        """
        if self.engine is not None and self.engine.maskdat is not None:
            return self.engine.maskdat > 0
        if getattr(self, "maskdat", None) is not None:
            return self.maskdat > 0
        return None

    def add_box_object(self):
        '''
        this adds a square region
        '''
        print('adding pixels to the mask')
        # mask out a rectangle around click
        # size is given by mask_size

        mask_value = self.engine.add_box_mask(self.xcursor, self.ycursor, self.mask_size)



        if mask_value == 0:
            print("you clicked outside the image area")
            return

        self.maskdat = self.engine.maskdat
        self.save_mask()         # should call engine.write_mask()
        print(f"added mask object {mask_value}")


    def add_circ_object(self):
        print('adding circular obj to the mask, with radius = ',self.mask_size)
        # mask out a rectangle around click
        # size is given by mask_size

        
        #pixel_mask = circle_pixels(float(self.xcursor),float(self.ycursor),float(self.mask_size/2.),self.xmax,self.ymax)

        mask_value = self.engine.add_circular_mask(float(self.xcursor),float(self.ycursor),float(self.mask_size/2.))

        self.save_mask()
        print(f'added mask object {mask_value}')
        
    def remove_object(self, objID):
        '''
        this will remove masked pixels near the cursor

        '''
        objID = int(objID)
 
        self.engine.remove_object(objID)

        # keep GUI view in sync (optional if you always read from engine)
        self.maskdat = self.engine.maskdat

        self.save_mask()# should call engine.write_mask()
        self._refresh()
        #self.display_mask()

    def set_threshold(self, t):
        try:
            self.threshold = float(t)
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
        print("off_center() needs a refresh for the new engine flow; update the target ellipse/object and press Run SE again.")
    # def off_center(self):
    #     t = input('enter object number for target galaxy\n')
    #     self.off_center_flag = True
    #     self.galaxy_id = int(t)
    #     if self.runse_flag:
    #         self.runse(weightim=self.weightim,weight_threshold=self.weight_threshold)
    #     else:
    #         self.run_photil()
    # def quit_program(self):
    #     #self.clean_links()
    #     #self.close_window()
    #     window = self.window()
    #     if window is not None:
    #         window.close()
    #     else:
    #         self.close()

    def quit_program(self):
        if hasattr(self, "MainWindow") and self.MainWindow is not None:
            try:
                top = self.MainWindow.parent()
                if top is not None:
                    top.close()
                    return
            except Exception:
                pass

            try:
                self.MainWindow.close()
                return
            except Exception:
                pass

        try:
            self.close()
        except Exception:
            pass
    
    def save_mask(self):
        #super(maskwindow,self).mask_saved(event)
        print('saving mask: ',self.mask_image)

        # engin handles writing mask
        self.engine.write_mask(self.mask_image)

        if not self.auto:
            self.mask_saved.emit(self.mask_image)
            self._refresh()
            #self.display_mask()

    def edit_mask(self):
        """
        Legacy method from the old non-event-driven workflow.

        In the Qt GUI, we don't run a blocking edit loop.
        Use the interactive key commands (c/b/r/g/w/q) instead.
        """
        print("edit_mask() is deprecated in the Qt GUI. Use the interactive window + 'w' to save.")
            
    def _debug_sizes(self):
        print()
        #print("AFTER SHOW - Form size:", self.MainWindow.size() if hasattr(self, "MainWindow") else "n/a")
        #print("AFTER SHOW - cutouts geom:", self.ui.cutouts.geometry(), "size:", self.ui.cutouts.size())
        #print("AFTER SHOW - cutoutsLayout count:", self.ui.cutoutsLayout.count())

    def grow(self):
        """
        Grow the current mask using the engine's grow method.
        """
        print("Growing mask from GUI...")

        self.engine.grow(
            size=int(self.mask_size),
            ngrow=1,
            preserve_gaia=True,
        )

        self.maskdat = self.engine.maskdat
        self.engine.write_mask(self.mask_image)

        if not self.auto:
            self._refresh()
            #self.display_mask()
        

    # def grow_mask(self):
    #     """
    #     Grow the current mask using the engine's grow method.
    #     """
    #     print("Growing mask from GUI...")

    #     # 1. Call engine grow
    #     self.engine.grow_mask()

    #     # 2. Update maskdat from engine
    #     self.maskdat = self.engine.maskdat

    #     # 3. Write updated mask to disk
    #     self.engine.write_mask(self.mask_image)

    #     # 4. Refresh GUI display
    #     if not self.auto:
    #         self.display_mask()

    def cycle_overlay_mode(self, which):
        attr = f"{which}_overlay_mode"
        mode = getattr(self, attr)

        if mode == "off":
            mode = "contour"
        elif mode == "contour":
            mode = "filled"
        else:
            mode = "off"

        setattr(self, attr, mode)
        self.update_overlay(which)

    def update_overlay(self, which):
        if self.engine is None or self.engine.maskdat is None:
            return

        viewer = self.rcutout if which == "r" else self.hacutout
        mode = self.r_overlay_mode if which == "r" else self.ha_overlay_mode

        self.clear_overlay(which)

        if mode == "off":
            viewer.fitsimage.redraw()
            return
        elif mode == "contour":
            self.draw_mask_contours(viewer, which)
        elif mode == "filled":
            self.draw_mask_filled(viewer, which)

        self.draw_central_ellipse(viewer)
        viewer.fitsimage.redraw()

    def clear_overlay(self, which):
        viewer = self.rcutout if which == "r" else self.hacutout
        try:
            viewer.canvas.delete_object_by_tag(f"mask-overlay-{which}")
        except Exception:
            pass
    
    # def draw_mask_contours(self, viewer, which):
    #     mask = self.engine.maskdat
    #     ids = np.unique(mask[mask > 0]).astype(int)

    #     compound = []

    #     for objid in ids:
    #         yy, xx = np.where(mask == objid)
    #         if len(xx) == 0:
    #             continue

    #         xmin, xmax = xx.min(), xx.max()
    #         ymin, ymax = yy.min(), yy.max()

    #         rect = viewer.dc.Rectangle(
    #             xmin, ymin, xmax, ymax,
    #             color="red",
    #             linewidth=1,
    #         )
    #         compound.append(rect)

    #     if compound:
    #         obj = viewer.dc.CompoundObject(*compound)
    #         viewer.canvas.add(obj, tag=f"mask-overlay-{which}")

    def draw_mask_filled(self, viewer, which):
        """
        Draw the actual mask footprint as a semi-transparent overlay.
        Uses contiguous horizontal runs so the overlay matches the mask pixels,
        rather than drawing coarse object bounding boxes.
        """
        mask_bool = self.get_boolean_mask()
        if mask_bool is None:
            return

        compound = []

        for y, x0, x1 in self.iter_mask_runs(mask_bool):
            # rectangle spanning exactly the masked run in this row
            rect = viewer.dc.Rectangle(
                x0, y, x1 + 1, y + 1,
                color="red",
                fill=True,
                fillalpha=0.25,
                alpha=0.25,
                linewidth=0,
            )
            compound.append(rect)

        if compound:
            obj = viewer.dc.CompoundObject(*compound)
            viewer.canvas.add(obj, tag=f"mask-overlay-{which}")

        
    # def draw_mask_filled(self, viewer, which):
    #     mask = self.engine.maskdat
    #     ids = np.unique(mask[mask > 0]).astype(int)

    #     compound = []

    #     for objid in ids:
    #         yy, xx = np.where(mask == objid)
    #         if len(xx) == 0:
    #             continue

    #         xmin, xmax = xx.min(), xx.max()
    #         ymin, ymax = yy.min(), yy.max()

    #         rect = viewer.dc.Rectangle(
    #             xmin, ymin, xmax, ymax,
    #             color="red",
    #             fill=True,
    #             fillalpha=0.25,
    #             alpha=0.25,
    #             linewidth=1,
    #         )
    #         compound.append(rect)

    #     if compound:
    #         obj = viewer.dc.CompoundObject(*compound)
    #         viewer.canvas.add(obj, tag=f"mask-overlay-{which}")

    def iter_mask_runs(self, mask_bool):
        """
        Yield contiguous True runs in each row of a boolean mask.

        Yields
        ------
        y, x0, x1
            Row index and inclusive pixel range [x0, x1].
        """
        ny, nx = mask_bool.shape

        for y in range(ny):
            row = mask_bool[y]
            if not np.any(row):
                continue

            # pad with False so edges are detected cleanly
            padded = np.concatenate(([False], row, [False])).astype(int)
            diff = np.diff(padded)

            starts = np.where(diff == 1)[0]
            stops = np.where(diff == -1)[0] - 1

            for x0, x1 in zip(starts, stops):
                yield y, int(x0), int(x1)

    

    def draw_mask_contours(self, viewer, which):
        """
        Draw contours from the actual mask footprint.
        """
        from scipy.ndimage import binary_erosion
        mask_bool = self.get_boolean_mask()
        if mask_bool is None:
            return

        edge = mask_bool & ~binary_erosion(mask_bool)

        compound = []
        for y, x0, x1 in self.iter_mask_runs(edge):
            rect = viewer.dc.Rectangle(
                x0, y, x1 + 1, y + 1,
                color="red",
                fill=False,
                linewidth=1,
            )
            compound.append(rect)

        if compound:
            obj = viewer.dc.CompoundObject(*compound)
            viewer.canvas.add(obj, tag=f"mask-overlay-{which}")
    def show_mask_mpl(self):
        # plot mpl figure
        # this was for debugging purposes
        self.image, self.imheader = fits.getdata(self.image_name,header = True)        
        print("plotting mask and central ellipse")
        self.v1,self.v2=scoreatpercentile(self.image,[5.,99.5])
        self.adjust_mask = True
        self.figure_size = (10,5)
        
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
        
