import sys
from PyQt5 import QtWidgets  # or PySide6
from hapy.maskgui.mask_window import MaskWindow

def main():
    app = QtWidgets.QApplication(sys.argv)
    win = QtWidgets.QMainWindow()
    logger = None

    mw = MaskWindow(
        win, logger,
        image="path/to/image.fits",
        haimage=None,
        sepath="sex",
        gaiapath=None,
        config="default.sex.HDI.mask",
        auto=False,
    )

    win.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
