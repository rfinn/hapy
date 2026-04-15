from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from PyQt5 import QtCore, QtWidgets

from hapy.maskgui.maskWidget import Ui_Form as Ui_maskWindow
from hapy.masktools.api import MaskEngine, EllipseParams


class MaskWindow(QtWidgets.QWidget):
    """
    Thin GUI wrapper around MaskEngine.

    GUI responsibilities:
    - load/display image products
    - expose a few user-editable controls
    - call the shared masking engine
    - update status/readout widgets

    Engine responsibilities:
    - source detection / segmentation
    - gaia masking
    - mask growth / cleanup
    - writing masks
    """

    mask_saved = QtCore.pyqtSignal(str)

    def __init__(
        self,
        *,
        image_fits: str | Path,
        haimage: str | Path | None = None,
        sepath: str | None = None,
        config: str | None = None,
        gaiapath: str | None = None,
        galaxy_ellipse: EllipseParams | None = None,
        threshold: float = 0.005,
        snr: float = 10.0,
        minarea: int = 5,
        grow_size: int = 7,
        grow_iterations: int = 1,
        weightim: str | None = None,
        weight_threshold: float = 1.0,
        add_gaia_stars: bool = True,
        logger=None,
        parent=None,
    ):
        super().__init__(parent)

        self.logger = logger
        self.image_fits = Path(image_fits)
        self.haimage_fits = Path(haimage) if haimage is not None else None
        self.sepath = sepath
        self.config = config
        self.gaiapath = gaiapath
        self.galaxy_ellipse = galaxy_ellipse

        self.threshold = float(threshold)
        self.snr = float(snr)
        self.minarea = int(minarea)
        self.grow_size = int(grow_size)
        self.grow_iterations = int(grow_iterations)
        self.weightim = weightim
        self.weight_threshold = float(weight_threshold)
        self.add_gaia_stars = bool(add_gaia_stars)

        self.ui = Ui_maskWindow()
        self.ui.setupUi(self)

        self.engine: MaskEngine | None = None
        self.mask: np.ndarray | None = None
        self.mask_path: Path | None = None

        self.image_data: np.ndarray | None = None
        self.image_header = None
        self.haimage_data: np.ndarray | None = None
        self.haimage_header = None

        self.image_viewer = None
        self.mask_viewer = None
        self.readout_label = None

        self._load_inputs()
        self._populate_controls()
        self._init_engine()
        self._setup_viewer_area()
        self._setup_readout_area()
        self._connect_signals()
        self._initialize_display()

    # ------------------------------------------------------------
    # setup helpers
    # ------------------------------------------------------------

    def _log(self, msg: str) -> None:
        if self.logger is not None:
            self.logger.info(msg)
        else:
            print(msg)

    def _load_inputs(self) -> None:
        self._log(f"Loading image: {self.image_fits}")
        self.image_data, self.image_header = fits.getdata(self.image_fits, header=True)

        if self.haimage_fits is not None and self.haimage_fits.exists():
            self._log(f"Loading Halpha image: {self.haimage_fits}")
            self.haimage_data, self.haimage_header = fits.getdata(
                self.haimage_fits, header=True
            )

        stem = self.image_fits.name
        if stem.lower().endswith(".fits"):
            stem = stem[:-5]
        self.mask_path = self.image_fits.with_name(f"{stem}-mask.fits")

    def _populate_controls(self) -> None:
        """
        Populate line edits from current state.

        The generated UI exposes:
        - boxSizeLineEdit
        - seThresholdLineEdit
        - seSNRLineEdit
        - seSNRAnalysisLineEdit
        """
        self.ui.boxSizeLineEdit.setText(str(self.grow_size))
        self.ui.seThresholdLineEdit.setText(str(self.threshold))
        self.ui.seSNRLineEdit.setText(str(self.snr))
        self.ui.seSNRAnalysisLineEdit.setText(str(self.snr))

    def _read_controls(self) -> None:
        """
        Read updated values from the GUI before rebuilding.
        """
        self.grow_size = self._safe_int(self.ui.boxSizeLineEdit.text(), self.grow_size)
        self.threshold = self._safe_float(
            self.ui.seThresholdLineEdit.text(), self.threshold
        )
        self.snr = self._safe_float(self.ui.seSNRLineEdit.text(), self.snr)

    @staticmethod
    def _safe_float(value: str, default: float) -> float:
        try:
            return float(value)
        except Exception:
            return float(default)

    @staticmethod
    def _safe_int(value: str, default: int) -> int:
        try:
            return int(float(value))
        except Exception:
            return int(default)

    def _init_engine(self) -> None:
        self.engine = MaskEngine(
            image_fits=str(self.image_fits),
            sepath=self.sepath,
            config=self.config,
            threshold=self.threshold,
            snr=self.snr,
            minarea=self.minarea,
            weightim=self.weightim,
            weight_threshold=self.weight_threshold,
            add_gaia_stars=self.add_gaia_stars,
        )

    def _connect_signals(self) -> None:
        self.ui.mrunSEButton.clicked.connect(self.build_mask)
        self.ui.mhelpButton.clicked.connect(self.show_help)
        self.ui.mquitButton.clicked.connect(self.close)

    # ------------------------------------------------------------
    # viewer / readout
    # ------------------------------------------------------------

    def _setup_viewer_area(self) -> None:
        """
        maskWidget.py shows that cutoutsLayout is already empty.
        No dummy widget removal is needed.
        """
        container = QtWidgets.QWidget(self.ui.cutouts)
        self.viewer_layout = QtWidgets.QGridLayout(container)
        self.viewer_layout.setContentsMargins(0, 0, 0, 0)

        # Temporary placeholders.
        # Replace with Ginga / CutoutImage widgets once the shell is stable.
        self.image_viewer = QtWidgets.QLabel("Main image viewer")
        self.image_viewer.setAlignment(QtCore.Qt.AlignCenter)
        self.image_viewer.setMinimumSize(300, 300)
        self.image_viewer.setStyleSheet("border: 1px solid gray;")

        self.mask_viewer = QtWidgets.QLabel("Mask viewer")
        self.mask_viewer.setAlignment(QtCore.Qt.AlignCenter)
        self.mask_viewer.setMinimumSize(300, 300)
        self.mask_viewer.setStyleSheet("border: 1px solid gray;")

        self.viewer_layout.addWidget(self.image_viewer, 0, 0)
        self.viewer_layout.addWidget(self.mask_viewer, 0, 1)

        self.ui.cutoutsLayout.addWidget(container, 0, 0)

    def _setup_readout_area(self) -> None:
        self.readout_label = QtWidgets.QLabel("Ready")
        self.ui.readoutGridLayout.addWidget(self.readout_label, 0, 0)

    def _set_readout(self, text: str) -> None:
        if self.readout_label is not None:
            self.readout_label.setText(text)
        self._log(text)

    def _initialize_display(self) -> None:
        self._display_main_image()
        self._display_mask()

    def _display_main_image(self) -> None:
        if self.image_data is None:
            self.image_viewer.setText("No image loaded")
            return

        self.image_viewer.setText(
            f"Main image\n{self.image_fits.name}\nshape={self.image_data.shape}"
        )

    def _display_mask(self) -> None:
        if self.mask is None:
            self.mask_viewer.setText("No mask yet")
            return

        nmask = int(np.sum(self.mask > 0))
        self.mask_viewer.setText(f"Mask ready\nmasked pixels={nmask}")

    # ------------------------------------------------------------
    # actions
    # ------------------------------------------------------------

    def build_mask(self) -> None:
        if self.engine is None:
            raise RuntimeError("Mask engine is not initialized")

        self._read_controls()
        self._init_engine()

        self._set_readout("Building mask...")

        self.mask = self.engine.build_initial_mask(
            galaxy_ellipse=self.galaxy_ellipse,
            grow_size=self.grow_size,
            grow_iterations=self.grow_iterations,
        )

        self._display_mask()
        self.save_mask()
        self._set_readout("Mask build complete")

    def save_mask(self) -> None:
        if self.engine is None:
            raise RuntimeError("Mask engine is not initialized")
        if self.mask_path is None:
            raise RuntimeError("Mask output path is not set")
        if self.mask is None:
            self._set_readout("No mask to save")
            return

        self.engine.write_mask(self.mask_path)
        self._set_readout(f"Saved mask: {self.mask_path}")
        self.mask_saved.emit(str(self.mask_path))

    def show_help(self) -> None:
        msg = (
            "Mask GUI\n\n"
            "Run SE builds or rebuilds the automatic mask.\n"
            "Help shows this dialog.\n"
            "Quit closes the window.\n\n"
            "Current control mapping:\n"
            "- Box Size -> grow_size\n"
            "- SE threshold -> engine threshold\n"
            "- SE SNR detect -> engine snr\n"
            "- SE SNR analysis -> not yet used separately in this thin wrapper"
        )
        QtWidgets.QMessageBox.information(self, "Mask Help", msg)

    # ------------------------------------------------------------
    # optional helpers
    # ------------------------------------------------------------

    def set_galaxy_ellipse(self, ellipse: EllipseParams | None) -> None:
        self.galaxy_ellipse = ellipse
        self._set_readout("Updated galaxy ellipse")

    def set_growth(self, grow_size: int, grow_iterations: int) -> None:
        self.grow_size = int(grow_size)
        self.grow_iterations = int(grow_iterations)
        self.ui.boxSizeLineEdit.setText(str(self.grow_size))
        self._set_readout(
            f"Updated growth: size={self.grow_size}, iterations={self.grow_iterations}"
        )

    def set_detection_params(self, threshold: float, snr: float, minarea: int) -> None:
        self.threshold = float(threshold)
        self.snr = float(snr)
        self.minarea = int(minarea)
        self.ui.seThresholdLineEdit.setText(str(self.threshold))
        self.ui.seSNRLineEdit.setText(str(self.snr))
        self._set_readout(
            f"Updated detection params: threshold={self.threshold}, "
            f"snr={self.snr}, minarea={self.minarea}"
        )
