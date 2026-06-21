#!/usr/bin/env python


"""
GOAL:
* subtract the continuum using color-correction from legacy images

USEAGE:

python subtract_continuum.py dirname

this assumes the file structure and naming convention used in the Virgo Filament Survey, with 

dirname-R.fits
dirname-Ha.fits
 
mask image:
dirname-R-mask.fits

subdirectory called
legacy/

that has g and r images that have been reprojected onto the halpha image pixel scale


PROCEDURE:
* create a g-r image
   2.5 log10(flux_r/flux_g)

* calculate rms in r-band outside the masked regions

* apply color correction to continuum-subtracted images
  

Filter transformations are in hapy/hatools/filter_transformations.py


"""
import sys
import os
import numpy as np
from astropy.io import fits
from astropy import stats, convolution
from astropy import wcs
from astropy.table import Table

import glob
from reproject import reproject_interp


from hapy.hatools.filter_transformations import halpha_minus_r_color_from_metadata
from hapy.hatools.filter_properties import get_continuum_oversubtraction_from_metadata
from hapy.hatools.segmentation import make_simple_photutils_segmentation
from hapy.imagetools.imutils import calculate_background_photutils
#from hapy.hatools.segmentation import make_simple_photutils_segmentation

#import warnings
#warnings.filterwarnings('ignore')


######################################################################
###  FILTER DEFINITIONS
######################################################################
# # Halpha filter width in angstrom
# filter_width_AA = {'BOK':80.48,'HDI':80.48,'INT':95,'MOS':80.48,'INT6657':80}

# # central wavelength in angstroms
# filter_lambda_c_AA = {'BOK':6620.52,'HDI':6620.52,'INT':6568,'MOS':6620.52,'INT6657':6657}


# # integral of filter transmission
# # calculated in Halpha-paper1.ipynb
# filter_Rlambda = {"KPNO_Ha+4nm": 78.58, "WFC_Ha": 84.21, "WFC_Ha6657": 71.96,\
#                   "KPNO_R" : 1341.54, "KPNO_r" : 1283.47, "BASS_r": 1042.18, "WFC_r": 1097.07}

# # from oversubtraction of continuum, a la Gavazzi+2006
# # take telescope as the key
# halpha_continuum_oversubtraction = {'BOK':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["BASS_r"]),\
#                             'HDI':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["KPNO_r"]),\
#                             'INT':(1 +filter_Rlambda["WFC_Ha"]/filter_Rlambda["WFC_r"]),\
#                             'MOS':(1 +filter_Rlambda["KPNO_Ha+4nm"]/filter_Rlambda["KPNO_R"]),\
#                             'INT6657':(1 +filter_Rlambda["WFC_Ha6657"]/filter_Rlambda["WFC_r"])}



def estimate_extra_continuum_scale(
    ha_data,
    rcont_data,
    galaxy_mask=None,
    bad_mask=None,
    min_pixels=100,
    clip_range=(0.75, 1.15),
    ratio_range=(0.7, 1.3),
    scale_percentile=30.0,
):
    if galaxy_mask is None:
        print("WARNING: galaxy_mask is None; using all valid pixels")
        galaxy_mask = np.ones_like(ha_data, dtype=bool)
    else:
        galaxy_mask = np.asarray(galaxy_mask, dtype=bool)

        if np.count_nonzero(galaxy_mask) == 0:
            print("WARNING: galaxy_mask is empty; using all valid pixels")
            galaxy_mask = np.ones_like(ha_data, dtype=bool)

    good = (
        galaxy_mask
        & np.isfinite(ha_data)
        & np.isfinite(rcont_data)
        & (rcont_data > 0)
    )

    if bad_mask is not None:
        bad_mask = np.asarray(bad_mask, dtype=bool)
        good &= ~bad_mask

    if np.count_nonzero(good) < min_pixels:
        print("WARNING: not enough valid pixels for auto continuum scaling")
        return 1.0
    
    ratio = ha_data[good] / rcont_data[good]
    ratio = ratio[np.isfinite(ratio)]

    if ratio_range is not None:
        lo, hi = ratio_range
        ratio = ratio[(ratio > lo) & (ratio < hi)]

    if len(ratio) < min_pixels:
        print("WARNING: not enough pixels after ratio clipping for auto continuum scaling")
        return 1.0

    raw_scale = np.nanpercentile(ratio, scale_percentile)
    scale = raw_scale

    if clip_range is not None:
        scale = np.clip(scale, clip_range[0], clip_range[1])

    print(f"auto continuum raw scale = {raw_scale:.4f}")
    print(f"auto continuum clipped scale = {scale:.4f}")
    print(f"auto continuum percentile = {scale_percentile:.1f}")
    print(f"pixels used for auto continuum scale = {len(ratio)}")

    return float(scale)




def estimate_scale_from_negative_tail(
    ha_data,
    rcont_data,
    galaxy_mask,
    sky_sigma,
    bad_mask=None,
    target=-1.5,
    percentile=5,
    scale_range=(0.75, 1.15),
    ngrid=81,
    min_pixels=100,
):
    good = (
        galaxy_mask
        & np.isfinite(ha_data)
        & np.isfinite(rcont_data)
        & np.isfinite(galaxy_mask)
    )

    if bad_mask is not None:
        good &= ~bad_mask

    if np.count_nonzero(good) < min_pixels:
        return 1.0

    ha = ha_data[good]
    rc = rcont_data[good]

    scales = np.linspace(scale_range[0], scale_range[1], ngrid)
    scores = []

    for s in scales:
        cs = ha - s * rc
        p = np.nanpercentile(cs, percentile)
        scores.append(abs((p / sky_sigma) - target))

    best = np.nanargmin(scores)
    return float(scales[best])

def estimate_scale_from_negative_tail_bisect(
    ha_data,
    rcont_data,
    galaxy_mask,
    sky_sigma,
    bad_mask=None,
    target=-1.5,
    percentile=5,
    scale_range=(0.75, 1.15),
    min_pixels=100,
    tol=0.02,
    max_iter=25,
):
    good = (
        galaxy_mask
        & np.isfinite(ha_data)
        & np.isfinite(rcont_data)
        & (rcont_data > 0)
    )

    if bad_mask is not None:
        good &= ~bad_mask

    if np.count_nonzero(good) < min_pixels:
        print("WARNING: not enough valid pixels for negative-tail scaling")
        return 1.0

    ha = ha_data[good]
    rc = rcont_data[good]

    if not np.isfinite(sky_sigma) or sky_sigma <= 0:
        print("WARNING: invalid sky_sigma for negative-tail scaling")
        return 1.0

    def tail_value(scale):
        cs = ha - scale * rc
        return np.nanpercentile(cs, percentile) / sky_sigma

    lo, hi = scale_range
    tail_lo = tail_value(lo)
    tail_hi = tail_value(hi)

    print(f"negtail: scale_lo={lo:.4f}, tail_lo={tail_lo:.3f}")
    print(f"negtail: scale_hi={hi:.4f}, tail_hi={tail_hi:.3f}")
    print(f"negtail: target={target:.3f}")

    # larger scale should make tail more negative
    if tail_lo < target and tail_hi < target:
        print("negtail: entire range is too negative; using minimum scale")
        return float(lo)

    if tail_lo > target and tail_hi > target:
        print("negtail: entire range is not negative enough; using maximum scale")
        return float(hi)

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        tail_mid = tail_value(mid)

        if abs(tail_mid - target) < tol:
            print(f"negtail: converged scale={mid:.4f}, tail={tail_mid:.3f}")
            return float(mid)

        # larger scale -> more negative tail
        if tail_mid < target:
            # too negative / oversubtracted
            hi = mid
        else:
            # not negative enough / undersubtracted
            lo = mid

    scale = 0.5 * (lo + hi)
    tail = tail_value(scale)
    print(f"negtail: final scale={scale:.4f}, tail={tail:.3f}")

    return float(scale)

# def estimate_scale_from_negative_tail_bisect(
#     ha_data,
#     rcont_data,
#     galaxy_mask,
#     sky_sigma,
#     bad_mask=None,
#     target=-1.5,
#     percentile=5,
#     scale_range=(0.75, 1.15),
#     min_pixels=100,
#     tol=0.002,
#     max_iter=25,
# ):
#     good = (
#         galaxy_mask
#         & np.isfinite(ha_data)
#         & np.isfinite(rcont_data)
#         & (rcont_data > 0)
#     )

#     if bad_mask is not None:
#         good &= ~bad_mask

#     if np.count_nonzero(good) < min_pixels:
#         print("WARNING: not enough valid pixels for negative-tail scaling")
#         return 1.0

#     ha = ha_data[good]
#     rc = rcont_data[good]

#     def tail_value(scale):
#         cs = ha - scale * rc
#         return np.nanpercentile(cs, percentile) / sky_sigma

#     lo, hi = scale_range
#     tail_lo = tail_value(lo)
#     tail_hi = tail_value(hi)

#     # If even the minimum scale is too negative, use minimum scale.
#     if tail_lo <= target:
#         print(f"negtail: minimum scale already too negative: tail={tail_lo:.2f}")
#         return float(lo)

#     # If even the maximum scale is not negative enough, use maximum scale.
#     if tail_hi >= target:
#         print(f"negtail: maximum scale not negative enough: tail={tail_hi:.2f}")
#         return float(hi)

#     for _ in range(max_iter):
#         mid = 0.5 * (lo + hi)
#         tail_mid = tail_value(mid)

#         if abs(tail_mid - target) < tol:
#             return float(mid)

#         # larger scale makes tail more negative
#         if tail_mid > target:
#             lo = mid
#         else:
#             hi = mid

#     return float(0.5 * (lo + hi))

def get_galaxy_region_from_segmap(segfile, mask=None, label=None):
    seg = fits.getdata(segfile)

    if label is None:
        labels = np.unique(seg)
        labels = labels[labels > 0]

        if len(labels) == 0:
            return None

        # fallback: use largest segment
        label = labels[np.argmax([(seg == lab).sum() for lab in labels])]

    galaxy_region = seg == label

    if mask is not None:
        galaxy_region &= ~mask

    return galaxy_region

def getEllipseFocii(xcent, ycent, a, ba, pa):
    pa_rad = pa * np.pi/180.
    b = a*ba
    c = np.sqrt(a**2-b**2)
    x1, y1 = c*np.cos(np.pi/2+pa_rad) + xcent, c*np.sin(np.pi/2+pa_rad) + ycent
    x2, y2 = c*np.cos(3*np.pi/2+pa_rad) + xcent, c*np.sin(3*np.pi/2+pa_rad) + ycent
    return x1, y1, x2, y2

def getEllipseCriterion(x, y, x1, y1, x2, y2, a):
    tot_dist = np.sqrt((x-x1)**2+(y-y1)**2) + np.sqrt((x-x2)**2+(y-y2)**2)
    cond = tot_dist <= 2*a
    return cond

def getEllipseAll(x, y, xcent, ycent, a, ba, pa):
    x1, y1, x2, y2 = getEllipseFocii(xcent, ycent, a, ba, pa)
    cond = getEllipseCriterion(x, y, x1, y1, x2, y2, a)
    return cond

def getPixLength(fits_image):
    hdr = fits_image[0].header
    try:
        cd1_1, cd1_2, cd2_1, cd2_2 = hdr['CD1_1'], hdr['CD1_2'], hdr['CD2_1'], hdr['CD2_2']
        pixlength_x = 3600.0 * np.sqrt(cd1_1**2 + cd2_1**2)
        pixlength_y = 3600.0 * np.sqrt(cd1_2**2 + cd2_2**2)
        return pixlength_x, pixlength_y
    except:
        return 3600.0 * abs(hdr['CDELT1']), 3600.0 * abs(hdr['CDELT2'])


def get_params_from_name(image_name):
    #print(t)
    tels = ['BOK','HDI','INT','MOS']
    for t in tels:
        if t in image_name:
            telescope = t
            break
    t = os.path.basename(image_name).split('-')
    for item in t:
        if item.startswith('20'):
            dateobs = item
            break
    pointing = t[-1]

    return telescope,dateobs,pointing
    
def filter_transformation(telescope,rfilter, gr_col):
    
    """
    use Matteo's linear fits to transform r to Halpha 

    I need to get more info on what these are - are they in mag or flux?

    QFM: is it ok that I am using these transformations on legacy g-r
    when they were derived for panstarrs g-r

    RETURN:
    * ha_r_no_nan - (Ha - r) color as a function of (g-r) color
    """
    if (telescope == 'BOK') :
        # need to get updated transformation from matteo that is using the
        # BASS r filter, which seems to havesignificantly lower transmission
        # than other r-band filters
        #Ha4_KPSr = -0.1804 * (gr_col) + 0.0158

        # FROM MATTEO 4/15/2024
        #Number of stars meeting the color and brightness cuts: 54842
        #---------------------------------------------------------------------
        #Best fit linear    Ha4 - BASSr = -0.1274 * (PS1_g-PS1_r) + 0.0151
        #Best fit quadratic Ha4 - BASSr = -0.0230*(PS1_g-PS1_r)^2 + -0.0976*(PS1_g-PS1_r) + 0.0063
        #---------------------------------------------------------------------
        #Initial bias/scatter:             -0.0608/0.0276
        #Corr linear fit              :    -0.0000/0.0125
        #Corr quadratic fit           :    0.0000/0.0125
        #---------------------------------------------------------------------
        #ha_r = -0.1804*gr_col + 0.0158
        ha_r = -0.1274*gr_col + 0.0151 # Ha4 vs BASS r
    
    elif ((telescope == 'HDI') and (rfilter == 'r')):
        #Ha4_KPSr = -0.1804 * (gr_col) + 0.0158
        ha_r = -0.1804*gr_col + 0.0158
    elif telescope == 'INT': # should be another case for the redder halpha, right?
        #Intha_INTSr = -0.2334 * (gr_col) + 0.0711
        ha_r = -0.2334 * (gr_col) + 0.0711
    elif (telescope == 'MOS') | ((telescope == 'HDI') and (rfilter == 'R')):
        #Ha4_KPHr = -0.0465 * (gr_col) + 0.0012
        ha_r = -0.0465 * (gr_col) + 0.0012
    else:
        print(f"HEY - DID NOT FIND A MATCH FOR TELESOPE {telescope} AND FILTER {rfilter}")

    # replace the nans with zeros
    ha_r_no_nan = np.nan_to_num(ha_r)
    return ha_r_no_nan

def getCorrelation(Halpha, cont):
    return filter2D(Halpha, ddepth=-1, kernel=cont)
     

def get_gr(gfile,rfile,mask=None, smooth_kernel=0):
    
    """ take g and r filenames, return g-r data and save g-r color image """
    g = fits.open(gfile)
    r = fits.open(rfile)
    data_g = g[0].data
    data_r = r[0].data
    g.close()
    r.close()
    
    ###
    # the following is from Matteo Fossati
    ###
    
    # get noise in the image    
    #stat is a tuple of mean, median, sigma
    #print('\nIn get_gr \nComputing median values for g and r images')
    stat_r = stats.sigma_clipped_stats(data_r,mask=mask)
    print('Subtracting {0:3.2e} from r-band image'.format(stat_r[1]))

    # this is not going to mask out the galaxy, so the sky values will likely be skewed
    # I am going to assume that the legacy images don't need another round of sky subtraction???
    data_r -= stat_r[1]
    stat_g = stats.sigma_clipped_stats(data_g,mask=mask)
    print('Subtracting {0:3.2e} from g-band image'.format(stat_g[1]))
    data_g -= stat_g[1]

    # create a mask, where SNR > 10
    # QUESTION : why is this 3 instead of 10?
    #usemask = (data_g>3*stat_g[2])
    usemask = (data_r>5*stat_r[2]) & (data_g>5*stat_g[2])

    # calculate the g-r color 
    gr_col = -2.5*np.log10(data_g/data_r)

    # TODONE - should add masking here - we don't want stars to be in our g-r image, right?
    if mask is not None:
        gr_col[mask] = np.nan
    print('Smoothing images for color calculation')
    # changing convolution size from 20 to 10 b/c I'm wondering if it's blurring the color
    # gradients too much - specific example is

    # Testing to remove smoothing on M109

    # we are using the legacy images that are reprojected on the halpha footprint,
    # so pixel scale is slightly larger than native
    if smooth_kernel and smooth_kernel > 0:
        print(f"Smoothing g-r image with Box2DKernel({smooth_kernel})")
        gr_col = convolution.convolve_fft(gr_col, convolution.Box2DKernel(10), allow_huge=True, nan_treatment='interpolate')

        if mask is not None:
            gr_col[mask] = np.nan
    else:
        print("Not smoothing g-r image")

    bad = np.logical_not(usemask)

    if mask is not None:
        bad = bad | mask

    gr_col[bad] = np.nan

    # set the pixel with SNR < 10 to nan - don't use these for color correction
    #gr_col[np.logical_not(usemask)] = np.nan
    
    # save gr color image
    hdu = fits.PrimaryHDU(gr_col, header=r[0].header)
    outimage = rfile.replace('r-ha.fits','gr-smooth.fits')
    #print(f"name for g-r image is {outimage}")
    print(f"writing g-r color image to {outimage}")
    hdu.writeto(outimage, overwrite=True)
    #hdu.close()
    return gr_col

def plot_image(data):
    ###########################################################
    # show the continuum subtracted image
    ###########################################################
    from matplotlib import pyplot as plt
    #from scipy.stats import scoreatpercentile
    from astropy.visualization import simple_norm
    

    plt.ion()
    fig = plt.figure()    
    norm = simple_norm(data, stretch='asinh',max_percent=99,min_percent=.5)
    plt.imshow(data, norm=norm,origin='lower',interpolation='nearest')#,vmin=v1,vmax=v2)
    #plt.show()
    #plt.draw()


# def zp_scale_r_to_ha(zp_ha, zp_r):
#     """Scale factor alpha so that CS = Ha - alpha * R."""
#     if zp_ha is None or zp_r is None:
#         return np.nan
#     zp_ha = float(zp_ha)
#     zp_r = float(zp_r)
#     if not (np.isfinite(zp_ha) and np.isfinite(zp_r)):
#         print("WARNING: could not calculate the zp scale!")
#         return np.nan
#     return float(10 ** (-0.4 * (zp_r - zp_ha)))


if __name__ == '__main__':

    import json
    from pathlib import Path
    import argparse


    parser = argparse.ArgumentParser(description="Make CS-gr continuum-subtracted image for one HAPY cutout directory.")

    parser.add_argument("cutdir", help="HAPY cutout directory, e.g. cutouts/VFID2550-UGC05020-INT-20190208-p031")
    parser.add_argument("--contscale", type=float, default=1.0, help="Manual extra continuum scale factor applied to the r-continuum image.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing CS-gr products.")


    parser.add_argument("--auto-contscale", action="store_true", help="Estimate an extra continuum scale factor using the galaxy segmentation region.")
    parser.add_argument("--auto-contscale-method", choices=["ratio", "negtail"], default="ratio", help="Method for estimating the extra continuum scale.")
    parser.add_argument("--auto-contscale-min", type=int, default=100, help="Minimum number of valid segmentation pixels required for auto continuum scale.")
    parser.add_argument("--auto-contscale-min-scale", type=float, default=0.75, help="Minimum allowed auto continuum scale.")
    parser.add_argument("--auto-contscale-max-scale", type=float, default=1.15, help="Maximum allowed auto continuum scale.")

    parser.add_argument("--auto-contscale-percentile", type=float, default=30.0, help="For method='ratio': percentile of clipped Ha/Rcont ratio used for auto continuum scale.")
    parser.add_argument("--auto-contscale-ratio-min", type=float, default=0.7, help="For method='ratio': minimum Ha/Rcont ratio included.")
    parser.add_argument("--auto-contscale-ratio-max", type=float, default=1.3, help="For method='ratio': maximum Ha/Rcont ratio included.")

    parser.add_argument("--auto-contscale-negtail-percentile", type=float, default=5.0, help="For method='negtail': lower percentile of CS pixels to compare against sky noise.")
    parser.add_argument("--auto-contscale-negtail-target", type=float, default=-1.5, help="For method='negtail': target lower-tail value in units of sky sigma.")
    parser.add_argument("--auto-contscale-negtail-ngrid", type=int, default=81, help="For method='negtail': number of trial continuum scales.")


    args = parser.parse_args()

    dirname = args.cutdir
    contscale = args.contscale
    overwrite = args.overwrite

    # get current directory
    topdir = os.getcwd()

    # define cutout directory and load metadata.json
    cutdir = Path(dirname)
    metafile = cutdir / "metadata.json"

    if not metafile.exists():
        raise FileNotFoundError(f"Missing metadata.json: {metafile}")

    with open(metafile, "r") as f:
        meta = json.load(f)

    # Use metadata.json as the source of truth
    tag = meta.get("tag", cutdir.name)
    objid = meta.get("objid", tag)
    vfid = objid.split("-")[0]

    telescope = meta["telescope"]
    dateobs = meta.get("dateobs", "")
    pointing = meta.get("pointing", "")
    hafilter = meta.get("hafilter", "")


    
    # Metadata-driven filter correction and filter properties
    halpha_filter_cor = float(meta.get("filter_correction", 1.0))
    if halpha_filter_cor == 0:
        print("resetting filter correction to 1")
        halpha_filter_cor = 1.0

    # This replaces header FLTRATIO / old lookup logic
    rscale = float(meta["filter_ratio"])

    # This replaces filter_lambda_c_AA[telescope] and filter_width_AA[telescope]
    hfilter_center_A = float(meta["hafilter_center_A"])
    hfilter_width_A = float(meta["hafilter_width_A"])

    # Not currently in metadata.json. Keep as neutral correction.
    halpha_extinction_correction = 1.0

    print(f"{tag}: telescope = {telescope}")
    print(f"{tag}: dateobs = {dateobs}")
    print(f"{tag}: pointing = {pointing}")
    print(f"{tag}: hafilter = {hafilter}")
    print(f"{tag}: halpha filter correction = {halpha_filter_cor:.4f}")
    print(f"{tag}: filter ratio = {rscale:.6f}")
    print(f"{tag}: halpha filter center = {hfilter_center_A:.2f} A")
    print(f"{tag}: halpha filter width = {hfilter_width_A:.2f} A")
    print(f"{tag}: continuum scale factor = {contscale:.3f}")
    print(f"{tag}: MW extinction correction = {halpha_extinction_correction:.3f}")


    cont_oversub = get_continuum_oversubtraction_from_metadata(meta)
    print(f"{tag}: continuum oversubtraction correction = {cont_oversub:.4f}")
    
    # move to subdirectory specified in the command line
    os.chdir(cutdir)


    ############################################################
    ## Define image names
    ############################################################    

    # define the file names
    Rfile = f"{tag}-R.fits"       # r-band image taken with same telescope as halpha
    Hfile = f"{tag}-Ha.fits"      # halpha image
    outname = f"{tag}-CS-gr.fits"


    segfile = f"{tag}-R-phot-segmentation.fits"

    
    # get legacy images that are reprojected to the halpha image
    # these are in the legacy subdirectory
    legacy_path = os.path.join("legacy", vfid + "*r-ha.fits")
    rfiles = glob.glob(legacy_path)

    if len(rfiles) < 1:
        print("problem getting r-ha.fits legacy image", len(rfiles))
        os.chdir(topdir)
        sys.exit(1)
    else:
        leg_rfile = rfiles[0]  # legacy r-band image

    # legacy g-band image, shifted to match halpha footprint and pixel scale
    gfiles = glob.glob(os.path.join("legacy", vfid + "*g-ha.fits"))

    if len(gfiles) < 1:
        print("problem getting g-ha.fits legacy image")
        os.chdir(topdir)
        sys.exit(2)
    else:
        leg_gfile = gfiles[0]  # legacy g-band image

    # define the mask file
    # HAPY convention:
    #   manual mask: <tag>-mask-manual.fits
    #   auto mask:   <tag>-mask.fits
    #   csgr sky mask fallback: <tag>-R-phot-segmentation.fits
    manual_mask = f"{tag}-mask-manual.fits"
    auto_mask = f"{tag}-mask.fits"
    segfile = f"{tag}-R-phot-segmentation.fits"

    if os.path.exists(manual_mask):
        maskfile = manual_mask
        print(f"Using manual mask: {maskfile}")

    elif os.path.exists(auto_mask):
        maskfile = auto_mask
        print(f"Using auto mask: {maskfile}")

    else:
        print("WARNING: no HAPY mask found; creating photutils segmentation mask")
        maskfile = make_simple_photutils_segmentation(
            Rfile,
            segfile,
            maskfile=None,
            nsigma=2.0,
            npixels=20,
        )
        print(f"Using photutils segmentation mask: {maskfile}")

    mask = fits.getdata(maskfile)
    mask = mask > 0

    if np.sum(mask) == 0:
        print(f"WARNING: mask file exists but has no masked pixels: {maskfile}")
        mask = None
        


    #overwrite = True

    """
    reproject infile to reffile image

    PARAMS:
    Rfile : r-band image taken with halpha, to be used for continuum
    Hfile : halpha image filename
    gfile : g-band filename to be used for calculating g-r color (legacy image)
    rfile : r-band filename to be used for calculating g-r color (legacy image)
    outname: output name

    RETURN:
    nothing, but save the CS subtracted image that uses the g-r color in the current directory
    """

    if os.path.exists(outname) & (not overwrite):
        print("continuum-subtracted image exists - not redoing it")
        os.chdir(topdir)
        sys.exit()




    ############################################################
    ## Get g-r image
    ############################################################
    outimage = leg_rfile.replace("r-ha.fits", "gr-smooth.fits")
    
    if os.path.exists(outimage) & (not overwrite):
        print("found g-r image. not remaking this")
        hdu = fits.open(outimage)
        gr_col = hdu[0].data
        hdu.close()
    else:
        gr_col = get_gr(leg_gfile, leg_rfile, mask=mask, smooth_kernel=5)

    # usemask should be all the values in the color image that are not equal to np.nan
    usemask = ~np.isnan(gr_col)  # these are the good values in the g-r color

    # this should be the text describing the galaxy
    # like : VFID0569-NGC5989-INT-20190530-p002
    fileroot = Rfile.replace("-R.fits", "")

    # read in *our* r-band and halpha images
    hhdu = fits.open(Hfile)
    rhdu = fits.open(Rfile)

    # get photometric ZP for each image
    rZP = rhdu[0].header["PHOTZP"]
    hZP = hhdu[0].header["PHOTZP"]

    # get filter names
    rfilter = rhdu[0].header["FILTER"]
    hfilter = hhdu[0].header["FILTER"]

    # get the pixel scale in the halpha image
    wcs_NB = wcs.WCS(Hfile)
    pscale_NB = wcs.utils.proj_plane_pixel_scales(wcs_NB) * 3600.0

    rscale_zp = zp_scale_r_to_ha(hZP, rZP)
    #print(f"scaling r-band continuum by {rscale:.6f}, ratio of ZP={zp_fratio}")

    rscale_meta = float(meta.get("filter_ratio", np.nan))
    rscale_zp = zp_scale_r_to_ha(hZP, rZP)

    print(f"{tag}: metadata filter_ratio = {rscale_meta:.6f}")
    print(f"{tag}: ZP r-to-Halpha scale = {rscale_zp:.6f}")
    print(f"{tag}: metadata / ZP scale = {rscale_meta / rscale_zp:.3f}")

    rscale = rscale_zp

    
    ##
    # The following is from Matteo Fossati
    ##

    print("\nGenerate NET image\n")

    # ############################################################
    # ## Subtract local sky from cutouts
    # ############################################################

    # HAPY cutouts should already be sky-subtracted when
    # metadata["cutout_sky_subtracted"] == True.
    data_r = rhdu[0].data.astype(float)
    data_NB = hhdu[0].data.astype(float)

    data_r_to_Ha = data_r * rscale
    
    # subtract local sky
    mean_sky_r, median_sky_r, std_sky_r = calculate_background_photutils(
        data_r_to_Ha,
        grow_radius=10,
        npixels=10,
        weightimage=r_weight,
        nsigma=2.0,
        clip_sigma=3.0,
    )
    data_r = data_r - median_sky_r

    mean_sky_h, median_sky_h, std_sky_h = calculate_background_photutils(
        data_NB,
        grow_radius=10,
        npixels=10,
        weightimage=r_weight,
        nsigma=2.0,
        clip_sigma=3.0,
    )
    data_NB = data_NB - median_sky_h
    

    print("Using HAPY cutout images with additional local sky subtraction")
    print(f"\tr-band sky subtracted = {median_sky_r:.2e}")
    print(f"\thalpha sky subtracted = {median_sky_h:.2e}")    

    # # TODONE - revisit this and examine the masking. - skipping sky subtraction here b/c already done when making cutouts
    # # the mask we are currently using does not mask the central galaxy
    # # also, we already subtract the sky from each continuum image when making cutouts...


    
    # # subtract sky from r-band image
    # print("Computing median values for r and halpha images")
    
    # print("subtracting these values from the image...")

    # stat_r = stats.sigma_clipped_stats(rhdu[0].data, mask=mask)
    # print("Subtracting {0:3.2e} from r-band image".format(stat_r[1]))

    #data_r = rhdu[0].data - stat_r[1]
    # data_r_to_Ha = data_r * rscale

    # # sky subtracted r-band image
    # skysub_r_name = Rfile.replace("-R.fits", "-R-sky.fits")
    # hdu = fits.PrimaryHDU(data_r, header=rhdu[0].header)
    # hdu.writeto(skysub_r_name, overwrite=True)

    # # subtract sky from Halpha image
    # stat_h = stats.sigma_clipped_stats(hhdu[0].data, mask=mask)
    # print("Subtracting {0:3.2e} from halpha image".format(stat_h[1]))
    # data_NB = hhdu[0].data - stat_h[1]


    ############################################################
    ## Transform images
    ############################################################
    
    ##
    # These comments are from Matteo's program
    ##
    # Generate the r band mag image and the r band calibrated to Halpha wave
    # This works only for positive flux pixels. Take this into account

    mag_r_to_Ha = -2.5 * np.log10(data_r_to_Ha) + rZP
    mag_NB = -2.5 * np.log10(data_NB) + hZP

    # now calc fluxes using the same ZP
    #data_r_ZP30 = 10.0 ** (-0.4 * (mag_r - 30))
    #data_NB_ZP30 = 10.0 ** (-0.4 * (mag_NB - 30))

    # Transform the mag_r image to the observed Halpha filter
    #
    # mag_r in AB mags
    # g-r color is in AB mags

    # going to stay in counts to avoid nans
    # create an image with delta needed to correct for color term
    # this is the fit to
    # delta_mag = (halpha - r) = f(g-r)
    delta_mag = halpha_minus_r_color_from_metadata(meta, gr_col)
    #delta_mag = np.zeros_like(gr_col, dtype=float)

    delta_mag_name = f"{tag}-CS-gr-delta-mag.fits"
    hdu = fits.PrimaryHDU(delta_mag, header=hhdu[0].header)
    hdu.header.set("IMTYPE", "DELTMAG", "Halpha - R color correction")
    hdu.header.set("CSTYPE", "CS-gr", "Used for CS-gr continuum subtraction")
    hdu.writeto(delta_mag_name, overwrite=True)
    print(f"Wrote {delta_mag_name}")
    
    # mag_r_to_Ha and mag_r should 
    mag_r_to_Ha = mag_r_to_Ha + delta_mag
    
    # convert to flux units
    delta_flux = 10.0 ** (-0.4 * delta_mag)

    # use the color correction for pixels with sufficient SNR
    data_r_to_Ha[usemask] = data_r_to_Ha[usemask] * delta_flux[usemask]


    # check for existing segmentation map
    if not os.path.exists(segfile):
        print(f"WARNING: missing photutils segmentation; creating {segfile}")
        segfile = make_simple_photutils_segmentation(
            Rfile,
            segfile,
            maskfile=maskfile,
        )

    # get additional scale factor to make 
    galaxy_region = get_galaxy_region_from_segmap(segfile, mask=mask)

    if args.auto_contscale_method == "ratio":
        extra_scale = estimate_extra_continuum_scale(
            ha_data=data_NB,
            rcont_data=data_r_to_Ha,
            galaxy_mask=galaxy_region,
            bad_mask=mask,
            min_pixels=args.auto_contscale_min,
            clip_range=(args.auto_contscale_min_scale, args.auto_contscale_max_scale),
            ratio_range=(args.auto_contscale_ratio_min, args.auto_contscale_ratio_max),
            scale_percentile=args.auto_contscale_percentile,
        )

    elif args.auto_contscale_method == "negtail":

        cs0 = data_NB - data_r_to_Ha

        stat_cs = stats.sigma_clipped_stats(cs0, mask=mask)
        sky_sigma = stat_cs[2]

        header_cskystd = hhdu[0].header.get("CSKYSTD", np.nan)
        print(f"Header CSKYSTD = {header_cskystd}")
        print(f"Measured CS sigma = {sky_sigma:.3f}")


        extra_scale = estimate_scale_from_negative_tail_bisect(
            ha_data=data_NB,
            rcont_data=data_r_to_Ha,
            galaxy_mask=galaxy_region,
            sky_sigma=sky_sigma,
            bad_mask=mask,
            target=args.auto_contscale_negtail_target,
            percentile=args.auto_contscale_negtail_percentile,
            scale_range=(args.auto_contscale_min_scale, args.auto_contscale_max_scale),
            min_pixels=args.auto_contscale_min,
        )
    
    # if args.auto_contscale:
    #     extra_scale = estimate_extra_continuum_scale(
    #         ha_data=data_NB,
    #         rcont_data=data_r_to_Ha,
    #         galaxy_mask=galaxy_region,
    #         bad_mask=mask,
    #         min_pixels=args.auto_contscale_min,
    #         clip_range=(args.auto_contscale_min_scale, args.auto_contscale_max_scale),
    #         ratio_range=(args.auto_contscale_ratio_min, args.auto_contscale_ratio_max),
    #         scale_percentile=args.auto_contscale_percentile,
    #         )

    # else:
    #     extra_scale = 1.0
    print(f"\nauto CS extra_scale = {extra_scale:.4f}\n")
    csgr_data = data_NB - extra_scale * data_r_to_Ha
    
    ##
    # Matteo Comment: Go to cgs units
    ##

    # TODONE - make sure filter quantities are correct
    # need:
    #   - center wavelength in A
    
    fnu_NB = 3.631E3 * data_NB * 1E-12
    flam_NB = 2.99792458E-5 * fnu_NB / (hfilter_center_A**2) * 1E18

    cnu_NB = 3.631E3 * data_r_to_Ha * 1E-12
    clam_NB = 2.99792458E-5 * cnu_NB / (hfilter_center_A**2) * 1E18

    # need to multiply by width of filter to convert from flux/A to flux
    flam_net = hfilter_width_A * (flam_NB - contscale * clam_NB)

    # correct CS flux for variations in filter transmission and oversubtraction due to halpha in r-band filter

    # TODONE - need to update continuum oversubtraction terms for full filters
    # skipping for now
    # flam_net = flam_net * halpha_continuum_oversubtraction[telescope] * halpha_filter_cor

    cont_oversub = get_continuum_oversubtraction_from_metadata(meta)
    flam_net = flam_net * cont_oversub * halpha_filter_cor
    
    # MW extinction correction is currently neutral because it is not in metadata.json
    flam_net = flam_net * halpha_extinction_correction

    # Save a version in AB/count-like units for compatibility with HAPY photometry programs

    # why are we using contscale again here when data_r_to_Ha is already scaled already
    # here, contscale is an extra factor the user can input to tweak the continuum subtraction


    # correct for filter transmission variations and for halpha emission in the continuum filter

    # TODONE - need to update oversubtraction terms using the correct filter traces
    # skipping for now
    # NB_ABmag = NB_ABmag * halpha_continuum_oversubtraction[telescope] * halpha_filter_cor
    csgr_data = csgr_data * cont_oversub * halpha_filter_cor


    # MW extinction correction is currently neutral because it is not in metadata.json
    csgr_data = csgr_data * halpha_extinction_correction

    ############################################################
    # Write CS-gr image
    ############################################################

    hhdu[0].header.set("CSTYPE", "CS-gr", "Continuum subtraction type")
    hhdu[0].header.set("CONSCALE", float(f"{contscale:.4f}"), "Continuum scale factor")
    hhdu[0].header.set("FILT_COR", float(f"{halpha_filter_cor:.4f}"), "Filter transmission correction")
    hhdu[0].header.set("FLTRATIO", float(f"{rscale:.8f}"), "r-to-Halpha continuum scale")
    hhdu[0].header.set("CONTOSUB", float(f"{cont_oversub:.4f}"), "CONT OVERSUB COR")
    hhdu[0].header.set("MWEXTCOR", float(f"{halpha_extinction_correction:.4f}"), "MW extinction correction")
    hhdu[0].header.set("HACEN_A", float(f"{hfilter_center_A:.2f}"), "Halpha filter center Angstrom")
    hhdu[0].header.set("HAWID_A", float(f"{hfilter_width_A:.2f}"), "Halpha filter width Angstrom")
    hhdu[0].header.set("SRCMETA", "metadata.json", "Source of correction metadata")

    hhdu[0].header.set("AUTOCONT", bool(args.auto_contscale), "Auto continuum scale used")
    hhdu[0].header.set("CONTSCL", float(extra_scale), "Extra continuum scale")
    #hhdu[0].header.set("CONTQ", float(args.auto_contscale_q), "Percentile used for auto continuum scale")
    hdu = fits.PrimaryHDU(csgr_data, header=hhdu[0].header)
    hdu.writeto(outname, overwrite=True)

    print(f"Wrote {outname}")

    # The rest are different versions of the CS image that Matteo saves.
    # Keeping calculations for diagnostics, but not writing by default.

    hdu = fits.PrimaryHDU(flam_NB, header=hhdu[0].header)

    hdu = fits.PrimaryHDU(clam_NB, header=hhdu[0].header)

    # Calculate clipped statistic
    # stat is a tuple of mean, median, sigma
    stat = stats.sigma_clipped_stats(flam_net, mask=mask)

    print("Unbinned SB limit 1sigma {0:3.2e} e-18".format(stat[2] / (pscale_NB[0] ** 2)))

    # This is the continuum-subtracted image in physical flux units
    hdu = fits.PrimaryHDU(flam_net, header=hhdu[0].header)

    # convert image to surface brightness units
    sblam_net = flam_net / (pscale_NB[0] ** 2)

    hdu = fits.PrimaryHDU(sblam_net, header=hhdu[0].header)

    print("Smoothing net image")
    flam_net_smooth = convolution.convolve_fft(
        flam_net,
        convolution.Box2DKernel(10),
        allow_huge=True,
        nan_treatment="interpolate",
    )

    hdu = fits.PrimaryHDU(flam_net_smooth, header=hhdu[0].header)

    stat_sm = stats.sigma_clipped_stats(flam_net_smooth, mask=mask)

    print(
        "Smoothed {1}x{1} SB limit 1sigma {0:3.2e} e-18".format(
            stat_sm[2] / (pscale_NB[0] ** 2), 15
        )
    )

    # close hdu files
    hhdu.close()
    rhdu.close()

    # move back to the top directory
    os.chdir(topdir)
