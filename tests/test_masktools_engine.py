import numpy as np
from hapy.masktools.api import MaskEngine

def test_engine_imports_without_qt():
    # If this import works in an env without Qt, you're headless.
    from hapy.masktools.api import MaskEngine  # noqa

def test_progress_callback_is_called(tmp_path, monkeypatch):
    # This test assumes you can run build_initial_mask without actually running sextractor
    # If sextractor is required, see the "mocking" section below.
    calls = []

    def cb(stage, fraction, message=None):
        calls.append(stage)

    # Use a small synthetic FITS you write for the test
    from astropy.io import fits
    img = np.zeros((20, 20), dtype=float)
    fpath = tmp_path / "img.fits"
    fits.writeto(fpath, img, overwrite=True)

    engine = MaskEngine(image_fits=str(fpath), sepath=None, gaiapath=None)
    # If your engine always runs sextractor, you’ll need to monkeypatch that function.
    # See next section.
    try:
        engine.build_initial_mask(progress_callback=cb)
    except Exception:
        pass

    assert "start" in calls
