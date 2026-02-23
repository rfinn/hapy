# hapy/maskgui/cutout_view.py
"""
Ginga + PyQt cutout display panel.

This is *view-layer* code: it owns widgets, callbacks, and readouts.
Keep it in maskgui so masktools can remain importable in headless mode.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from PyQt5 import QtCore, QtWidgets

# Ginga imports
try:
    from ginga.qtw.ImageViewQt import CanvasView
    from ginga.canvas.CanvasObject import get_canvas_types
    from ginga.util.loader import load_data
except ModuleNotFoundError as e:
    raise ModuleNotFoundError(
        "ginga is required for CutoutPanel. Install it (e.g. `pip install ginga`)."
    ) from e


@dataclass
class CutoutReadout:
    x: float
    y: float
    value: Optional[float]


class CutoutPanel(QtCore.QObject):
    """
    A single image viewer panel with:
      - ginga CanvasView widget
      - cursor readout label
      - key-press signal emitted with key and data coords

    Parameters
    ----------
    parent_layout:
        Layout to which the image widget is added.
    readout_layout:
        Layout to which the readout label is added.
    logger:
        Ginga logger (e.g., ginga.misc.log.get_logger(...)).
    grid_pos:
        (row, col, row_span, col_span) for placing the image widget.
    readout_pos:
        (row, col, row_span, col_span) for placing the readout label.
    autocut_params:
        Ginga autocut preset (e.g. "stddev").
    enable_drawing:
        If True, enable rectangle drawing on a canvas overlay.
    """

    key_pressed = QtCore.pyqtSignal(str)  # emits "keyname,data_x,data_y"

    def __init__(
        self,
        parent_layout,
        readout_layout,
        logger,
        grid_pos: Tuple[int, int, int, int],
        readout_pos = None,
        autocut_params: str = "stddev",
        enable_drawing: bool = True,
    ):
        super().__init__()

        self.logger = logger
        self._autocut_params = autocut_params

        # ---- Create ginga viewer ----
        fi = CanvasView(self.logger, render="widget")
        fi.enable_autocuts("on")
        fi.set_autocut_params(autocut_params)
        fi.enable_autozoom("on")
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        fi.show_focus_indicator(True)

        bd = fi.get_bindings()
        bd.enable_all(True)

        self.fitsimage = fi
        self.widget = fi.get_widget()

        r, c, rs, cs = grid_pos
        parent_layout.addWidget(self.widget, r, c, rs, cs)
        self.widget.show()
        # ---- Readout label ----
        self.readout = QtWidgets.QLabel("")
        # commenting for now
        #rr, rc, rrs, rcs = readout_pos
        #readout_layout.addWidget(self.readout, rr, rc, rrs, rcs)

        # ---- Callbacks ----
        self.fitsimage.set_callback("none-move", self._cursor_cb)
        self.fitsimage.set_callback("key-press", self._key_press_cb)

        # ---- Optional drawing canvas overlay ----
        self.dc = get_canvas_types()
        self.canvas = None
        self.drawtypes = []

        enable_drawing = True
        if enable_drawing:
            canvas = self.dc.DrawingCanvas()
            canvas.enable_draw(True)
            canvas.enable_edit(True)
            canvas.set_drawtype("rectangle", color="lightblue")
            canvas.set_surface(fi)

            private_canvas = fi.get_canvas()
            private_canvas.add(canvas)

            canvas.register_for_cursor_drawing(fi)
            canvas.set_draw_mode("draw")
            canvas.ui_set_active(True)

            self.canvas = canvas
            self.drawtypes = sorted(canvas.get_drawtypes())

    # ----------------------
    # Public API
    # ----------------------
    def load_file(self, filepath: str) -> None:
        """Load a FITS file into the viewer."""
        image = load_data(filepath, logger=self.logger)
        self.fitsimage.set_image(image)
        # Force display update
        self.fitsimage.zoom_fit()
        self.fitsimage.redraw(whence=0)
        
    def load_data(self, image_array) -> None:
        """Load a numpy array into the viewer."""
        self.fitsimage.set_data(image_array)
        self.fitsimage.zoom_fit()
        self.fitsimage.redraw(whence=0)
        
    def set_autocut_params(self, autocut_params: str) -> None:
        self._autocut_params = autocut_params
        self.fitsimage.set_autocut_params(autocut_params)

    def get_last_data_xy(self) -> Tuple[float, float]:
        return self.fitsimage.get_last_data_xy()

    # ----------------------
    # Callbacks
    # ----------------------
    def _cursor_cb(self, viewer, button, data_x, data_y):
        """Update cursor readout label."""
        try:
            value = viewer.get_data(int(data_x + viewer.data_off), int(data_y + viewer.data_off))
        except Exception:
            value = None

        fits_x, fits_y = data_x + 1, data_y + 1

        if value is None:
            text = f"X: {fits_x:.1f}  Y: {fits_y:.1f}  Value: --"
        else:
            try:
                text = f"X: {fits_x:.1f}  Y: {fits_y:.1f}  Value: {float(value):.2f}"
            except Exception:
                text = f"X: {fits_x:.1f}  Y: {fits_y:.1f}  Value: --"

        self.readout.setText(text)

    def _key_press_cb(self, canvas, keyname):
        """Emit key press + current data coords as a single string."""
        data_x, data_y = self.fitsimage.get_last_data_xy()
        self.key_pressed.emit(f"{keyname},{data_x},{data_y}")
