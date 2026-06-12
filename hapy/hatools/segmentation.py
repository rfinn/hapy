# def make_simple_photutils_segmentation(rfile, segfile, maskfile=None,
#                                        nsigma=2.0, npixels=20):


#     from photutils.segmentation import detect_sources
#     from photutils.background import Background2D, MedianBackground
#     from astropy.stats import SigmaClip
#     from astropy.io import fits
#     import numpy as np

#     data, hdr = fits.getdata(rfile, header=True)

#     bad = ~np.isfinite(data)

#     if maskfile is not None and os.path.exists(maskfile):
#         m = fits.getdata(maskfile)
#         bad |= m > 0

#     sigma_clip = SigmaClip(sigma=3.0, maxiters=5)
#     bkg = Background2D(
#         data,
#         box_size=64,
#         filter_size=3,
#         sigma_clip=sigma_clip,
#         bkg_estimator=MedianBackground(),
#         mask=bad,
#     )

#     threshold = bkg.background + nsigma * bkg.background_rms

#     segm = detect_sources(
#         data,
#         threshold,
#         npixels=npixels,
#         mask=bad,
#     )

#     if segm is None:
#         segdata = np.zeros(data.shape, dtype=np.int16)
#     else:
#         segdata = segm.data.astype(np.int16)

#     hdr["SEGSRC"] = ("photutils", "Segmentation source")
#     hdr["SEGNSIG"] = (float(nsigma), "Detection threshold")
#     hdr["SEGNPX"] = (int(npixels), "Minimum connected pixels")

#     fits.PrimaryHDU(data=segdata, header=hdr).writeto(segfile, overwrite=True)

#     return str(segfile)


def make_simple_photutils_segmentation(rfile, segfile, maskfile=None,
                                       nsigma=2.0, npixels=20):

    import os
    import numpy as np
    from astropy.io import fits
    from astropy.stats import SigmaClip, sigma_clipped_stats
    from photutils.background import Background2D, MedianBackground
    from photutils.segmentation import detect_sources

    data, hdr = fits.getdata(rfile, header=True)
    data = np.asarray(data, dtype=float)

    bad = ~np.isfinite(data)

    if maskfile is not None and os.path.exists(maskfile):
        m = fits.getdata(maskfile)
        bad |= m > 0

    ny, nx = data.shape
    min_dim = min(nx, ny)

    # Choose a background box size appropriate for small cutouts.
    # Keep it large enough to estimate sky, but not larger than ~1/4 image size.
    if min_dim < 128:
        box_size = max(16, min_dim // 4)
    elif min_dim < 256:
        box_size = 32
    else:
        box_size = 64

    # Ensure npixels is not too large for tiny cutouts.
    npixels_use = min(npixels, max(5, int(0.002 * data.size)))

    sigma_clip = SigmaClip(sigma=3.0, maxiters=5)

    try:
        bkg = Background2D(
            data,
            box_size=box_size,
            filter_size=3,
            sigma_clip=sigma_clip,
            bkg_estimator=MedianBackground(),
            mask=bad,
            exclude_percentile=50.0,
        )

        threshold = bkg.background + nsigma * bkg.background_rms
        bkg_method = "Background2D"

    except Exception as err:
        print(
            f"WARNING: Background2D failed for {rfile} with box_size={box_size}: {err}"
        )
        print("WARNING: falling back to sigma_clipped_stats background")

        good = np.isfinite(data) & ~bad
        if np.sum(good) > 0:
            mean, med, std = sigma_clipped_stats(data[good], sigma=3.0, maxiters=5)
        else:
            med, std = 0.0, np.nan

        if not np.isfinite(std) or std <= 0:
            std = np.nanstd(data[good]) if np.sum(good) > 0 else 1.0

        if not np.isfinite(std) or std <= 0:
            std = 1.0

        threshold = med + nsigma * std
        bkg_method = "sigma_clipped_stats"

    segm = detect_sources(
        data,
        threshold,
        npixels=npixels_use,
        mask=bad,
    )

    if segm is None:
        segdata = np.zeros(data.shape, dtype=np.int16)
    else:
        segdata = segm.data.astype(np.int16)

    hdr["SEGSRC"] = ("photutils", "Segmentation source")
    hdr["SEGMETH"] = (bkg_method, "Background method for segmentation")
    hdr["SEGNSIG"] = (float(nsigma), "Detection threshold")
    hdr["SEGNPX"] = (int(npixels_use), "Minimum connected pixels")
    hdr["SEGBOX"] = (int(box_size), "Background2D box size")

    fits.PrimaryHDU(data=segdata, header=hdr).writeto(segfile, overwrite=True)

    return str(segfile)
