from hapy.masktools.api import MaskEngine

def progress_cb(stage, fraction, message=None):
    print(stage, fraction, message)

engine = MaskEngine(
    image_fits="path/to/image.fits",
    sepath="sex",                 # or None if you handle default
    gaiapath=None,
    config="default.sex.HDI.mask",
    threshold=0.005,
    snr=10,
    minarea=5,
)

mask = engine.build_initial_mask(progress_callback=progress_cb)

print(mask.shape, mask.dtype, mask.min(), mask.max())
