from hapy.masktools.api import MaskEngine, EllipseParams

def progress_cb(stage, fraction, message=None):
    print(stage, fraction, message)

testim = "/Users/rfinn/github/hapy/test_images/VFID0569-NGC5989-BOK-20220424-VFID0607-R.fits"
testim = "/Users/rfinn/github/hapy/test_images/VFID2463-NGC5240-INT-20190209-p106-R.fits"
testim = "/Users/rfinn/github/hapy/test_images/VFID2486-UGC08318-BOK-20210417-VFID2488-R.fits"
pixscale = 0.45
# trying to find one with obvious gaia stars
testim = "/Users/rfinn/github/hapy/test_images/VFID2507-UGC08693-HDI-20170523-p007-R.fits"
pixscale = 0.45
sma = 84
sma_pix = sma/pixscale
ba = 0.38
b = ba * sma
pa_deg = 170. + 90

galEllipse = EllipseParams(None, None, sma_pix, ba, pa_deg)

    
engine = MaskEngine(
    image_fits=testim,
    sepath="sex",                 # or None if you handle default
    gaiapath=None,
    config="default.sex.HDI.mask",
    threshold=0.005,
    snr=10,
    minarea=5,
    add_gaia_stars=True
)

mask = engine.build_initial_mask(progress_callback=progress_cb)

print(mask.shape, mask.dtype, mask.min(), mask.max())


# testing some more

engine.build_initial_mask()

engine.write_mask(testim.replace('.fits','-mask.fits'))
