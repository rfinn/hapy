import sys
from PyQt5 import QtWidgets  # or PySide6
from hapy.maskgui.mask_window import MaskWindow
from hapy.masktools.api import MaskEngine, EllipseParams
from hapy.geometry.adapters import pa_ccw_north_to_photutils_theta
testim = "/Users/rfinn/github/hapy/test_images/VFID2507-UGC08693-HDI-20170523-p007-R.fits"

def main():
    app = QtWidgets.QApplication(sys.argv)

    mainwin = QtWidgets.QMainWindow()
    form = QtWidgets.QWidget(mainwin)     # <-- QWidget parented to main window
    mainwin.setCentralWidget(form)        # <-- QWidget becomes the main content

    #testim = "/Users/rfinn/github/hapy/test_images/VFID0569-NGC5989-BOK-20220424-VFID0607-R.fits"
    pixscale = 0.45
    # trying to find one with obvious gaia stars
    testim = "/Users/rfinn/github/hapy/test_images/VFID2507-UGC08693-HDI-20170523-p007-R.fits"
    pixscale = 0.45
    sma = 84
    sma_pix = sma/pixscale
    ba = 0.38
    b = ba * sma
    pa_deg = 170. + 90
    # convert to 
    ra = 206.1183347
    dec = 35.1925572
    objparams = [ra,dec,sma, ba, pa_deg]
    theta_deg = pa_ccw_north_to_photutils_theta(pa_deg)
    #class EllipseParams:
    #"""Ellipse in pixel coords.  theta_deg is angle CCW from +x axis"""
    #xc: float
    #yc: float
    #sma_pix: float
    #ba: float
    #theta_deg: float # angle in deg from +x axis
    #theta_deg = pa_deg

    # work on this in the future, to pass in EllipseParams to MaskWindow and use the hapy datastructure instead of the objparams list
    #galEllipse = EllipseParams(xc=None, yc=None, sma_pix=sma_pix, ba=ba, theta_deg=theta_deg)

    mw = MaskWindow(
        form,              # <-- PASS THE QWidget HERE
        logger=None,
        image=testim,
        haimage=None,
        config="default.sex.HDI.mask",
        threshold=0.005,
        snr=10,
        minarea=5,
        objparams=objparams,
        auto=False,
    )

    mainwin.setWindowTitle("makin a mask...")
    mainwin.resize(900, 600)
    # or: mainwin.showMaximized()
    mainwin.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()


