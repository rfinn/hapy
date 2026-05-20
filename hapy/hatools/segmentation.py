def make_simple_photutils_segmentation(rfile, segfile, maskfile=None,
                                       nsigma=2.0, npixels=20):


    from photutils.segmentation import detect_sources
    from photutils.background import Background2D, MedianBackground
    from astropy.stats import SigmaClip
    from astropy.io import fits
    import numpy as np

    data, hdr = fits.getdata(rfile, header=True)

    bad = ~np.isfinite(data)

    if maskfile is not None and os.path.exists(maskfile):
        m = fits.getdata(maskfile)
        bad |= m > 0

    sigma_clip = SigmaClip(sigma=3.0, maxiters=5)
    bkg = Background2D(
        data,
        box_size=64,
        filter_size=3,
        sigma_clip=sigma_clip,
        bkg_estimator=MedianBackground(),
        mask=bad,
    )

    threshold = bkg.background + nsigma * bkg.background_rms

    segm = detect_sources(
        data,
        threshold,
        npixels=npixels,
        mask=bad,
    )

    if segm is None:
        segdata = np.zeros(data.shape, dtype=np.int16)
    else:
        segdata = segm.data.astype(np.int16)

    hdr["SEGSRC"] = ("photutils", "Segmentation source")
    hdr["SEGNSIG"] = (float(nsigma), "Detection threshold")
    hdr["SEGNPX"] = (int(npixels), "Minimum connected pixels")

    fits.PrimaryHDU(data=segdata, header=hdr).writeto(segfile, overwrite=True)

    return str(segfile)
