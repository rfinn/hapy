import numpy as np

#from PyQt5 import  QtWidgets
from PyQt5 import QtCore,QtWidgets, QtGui
#from PyQt5.Qtcore import  Qt
#from ginga.qtw.QtHelp import QtGui, QtCore

#from halphav3 import Ui_MainWindow
try:
    from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
    #from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
    from ginga.mplw.ImageViewMpl import ImageView
    from ginga import colors
    from ginga.canvas.CanvasObject import get_canvas_types
    from ginga.misc import log
    from ginga.util.loader import load_data
    from astropy.wcs import WCS
except ModuleNotFoundError:
    print("Warning - ginga module not found.  ")
## filter information
## from https://www.noao.edu/kpno/mosaic/filters/
central_wavelength = {'4':6620.52,'8':6654.19,'12':6698.53,'16':6730.72,'R':6513.5,'r':6292.28} # angstrom
dwavelength = {'4':80.48,'8':81.33,'12':82.95,'16':81.1,'R':1511.3,'r':1475.17} # angstrom


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

class cutout_image():
    #key_pressed = QtCore.pyqtSignal(str)
    def __init__(self,panel_name,ui, logger, row, col, drow, dcol, autocut_params='stddev'):
        #super(image_panel, self).__init__()
        # enable some user interaction
        #fi.get_bindings.enable_all(True)

        self.logger = logger
        self.ui = ui
        self.dc = get_canvas_types()
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params(autocut_params)
        #fi.set_autocut_params('stddev')
        fi.enable_autozoom('on')
        #fi.set_callback('drag-drop', self.drop_file)
        #fi.set_callback('none-move',self.cursor_cb)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        fi.show_focus_indicator(True)
        self.fitsimage = fi
        
        bd = fi.get_bindings()
        bd.enable_all(True)
        self.cutout = fi.get_widget()
        self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.fitsimage = fi

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
        self.image_wcs = WCS(filepath)
        self.fitsimage.set_image(image)
        #self.setWindowTitle(filepath)


class mplWindow(QtWidgets.QDialog):
    
    def __init__(self, parent=None):
        super(mplWindow, self).__init__(parent)



        # a figure instance to plot on
        self.figure = Figure(figsize=(2,2))
        self.axes = self.figure.add_subplot(111)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
        self.setParent(parent)
        
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        # set the layout
        #self.plot()
        #parent.addWidget(self.canvas)
        #self.show()
        if parent:
            parent.addWidget(self.canvas,1,2,4,1)
        else:
            layout = QtWidgets.QGridLayout(parent)
            layout.addWidget(self.canvas,1,2,4,1)
        self.setLayout(layout)
  
    def show(self):
        self.axes.plot(np.arange(10))
        self.canvas.draw()

