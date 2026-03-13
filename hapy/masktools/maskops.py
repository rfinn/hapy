
# hapy/masktools/maskops.py
from __future__ import annotations
from typing import Iterable, List
import numpy as np

from photutils.aperture import EllipticalAperture
from .types import EllipseParams, MaskFractionResult

defaultcat='default.sex.HDI.mask'

#from hapy.imagetools.imutils import circle_pixels


def circle_pixels(xc,yc,r,ximage,yimage):
    '''
    GOAL:
    - return pixel values that lie within a circular aperture within radius r of position (x,y)
    
    INPUT:
    - enter the center xc,yc and radius of circle in pixels
    - also enter x and y dimensions of parent image

    OUTPUT:
    - 2D boolean array with dimension the same as the input image
    - pixel values are true for pixels within circular aperture, false otherwise
    '''

    # add some checks to make sure numbers make sense
    # actually, it works even if center is outside image boundaries
    #if (xc < 0) | (xc > ximage) | (yc < 0) | (yc > yimage):
    #    print('invalid central coordinates in circle_pixels')
        
    rows,cols = np.mgrid[0:yimage,0:ximage]
    distance = np.sqrt((rows-yc)**2+(cols-xc)**2)
    pixel_flag = distance < r
    return pixel_flag


def remove_central_objects(mask, ellipse_params=None):
    """ 
    find any pixels within central ellipse and set their values to zero 

    PARAMS:
    mask = 2D array containing masked pixels, like from SE segmentation image
    ellipse_params = EllipseDataStructure 
      - sma = semi-major axis in pixels
      - BA = ratio of semi-minor to semi-major axes
      - PA = position angle, measured in degree counter clockwise from +x axis
      - xc
      - yc
    OPTIONAL ARGS:
    xc = center of ellipse in pixels; assumed to be center of image if xc is not specified
    yc = center of ellipse in pixels; assumed to be center of image if yc is not specified
    
    RETURNS:
    newmask = copy of input mask, with pixels within ellipse set equal to zero

    """

    if ellipse_params is None:
        objid = find_central_objid(mask)
        newmask = remove_object_by_id(mask, objid)
        return newmask, ellipse_params
    else:
        # changing the xmax and ymax - if the ellipse looks wrong, then swap back
        ymax,xmax = mask.shape

        if ellipse_params is not None:
            xc = ellipse_params.xc
            yc = ellipse_params.yc
            a = ellipse_params.sma_pix
            b = ellipse_params.ba * a
            phirad = np.radians(ellipse_params.theta_deg)
        else:
            # set center of ellipse as the center of the image
            a = None
            b = None
            phirad = None
            xc = None
            yc = None

        if xc is None:
            xc = xmax//2
            yc = ymax//2

        X,Y = np.meshgrid(np.arange(xmax),np.arange(ymax))
    
        p1 = ((X-xc)*np.cos(phirad)+(Y-yc)*np.sin(phirad))**2/a**2
        p2 = ((X-xc)*np.sin(phirad)-(Y-yc)*np.cos(phirad))**2/b**2
        flag2 = p1+p2 < 1
        newmask = np.copy(mask)
        newmask[flag2] = 0
        # we could also get all the unique values associated with flag2, and then remove them
        #ellipse_params = [xc,yc,sma,BA,phirad]
        return newmask,ellipse_params



def find_central_objid(mask):
    ymax,xmax = mask.shape

    xc = xmax//2
    yc = ymax//2

    return mask[yc,xc]

def grow_mask_square(maskdat: np.ndarray, size: int = 7, ngrow=1) -> np.ndarray:
    """
    Grow labeled regions by painting a size×size square around each masked pixel,
    preserving the pixel's label value.
    """
    out = np.array(maskdat, copy=True)
    nx, ny = out.shape
    half = int(size / 2)

    for k in range(ngrow):
        masked = np.where(out > 0)
        for i, j in zip(masked[0], masked[1]):
            r0 = max(0, i - half)
            r1 = min(nx, i + half)
            c0 = max(0, j - half)
            c1 = min(ny, j + half)
            if r1 <= r0 or c1 <= c0:
                continue
            out[r0:r1, c0:c1] = out[i, j]
    return out




def apply_user_masks(maskdat, usr_mask, deleted_objects=None):
    """
    Apply user-created masks and remove deleted object IDs.

    Parameters
    ----------
    maskdat : ndarray
        Base segmentation/mask image.
    usr_mask : ndarray
        User mask layer (same shape). Nonzero values are added.
    deleted_objects : iterable of int
        Object IDs to delete (set to zero).

    Returns
    -------
    newmask : ndarray
    """
    newmask = maskdat + usr_mask
    if deleted_objects is not None:
        for obj_id in deleted_objects:
            newmask[newmask == int(obj_id)] = 0
    return newmask



def remove_object_by_id(mask_array, object_id):
    """
    Remove a segmentation object from mask by ID.

    Parameters
    ----------
    mask_array : 2D numpy array
    object_id : int or None

    Returns
    -------
    new_mask : 2D numpy array
    """

    if object_id is None or np.isnan(object_id):
        return mask_array.copy()

    mask = mask_array.copy()
    mask[mask == object_id] = 0

    return mask


def remove_objects_in_ellipse(mask_array, ellipse_params):
    """
    Remove objects within one or more ellipses.

    Parameters
    ----------
    mask_array : 2D numpy array
    ellipse_params : dict or list of dicts

        Each dict must contain:
        {
            "sma": float,
            "BA": float,
            "PA": float,
            "xc": float,
            "yc": float
        }

    Returns
    -------
    new_mask : 2D numpy array
    used_ellipses : list of returned ellipse parameters
    """

    mask = mask_array.copy()

    if ellipse_params is None:
        return mask, None

    if isinstance(ellipse_params, dict):
        ellipse_params = [ellipse_params]

    used_ellipses = []

    for e in ellipse_params:
        mask, eparams = remove_central_objects(
            mask,
            sma=e["sma"],
            BA=e["BA"],
            PA=e["PA"],
            xc=e["xc"],
            yc=e["yc"],
        )
        used_ellipses.append(eparams)

    return mask, used_ellipses





def ellipse_mask_fraction(mask, ell: EllipseParams):
    """
    Compute fraction of masked pixels inside an ellipse.

    Parameters
    ----------
    mask : 2D ndarray
        Boolean or integer mask image. Nonzero / True = masked.
    ell : EllipseParams
        Ellipse definition in pixel coordinates.

    Returns
    -------
    frac_masked : float
        Fraction of ellipse pixels that are masked.
    n_total : int
        Total number of pixels in ellipse.
    n_masked : int
        Number of masked pixels in ellipse.
    n_unmasked : int
        Number of unmasked pixels in ellipse.
    """
    mask_bool = np.asarray(mask) > 0

    if not np.isfinite(ell.xc) or not np.isfinite(ell.yc):
        return np.nan, 0, 0, 0
    if not np.isfinite(ell.sma_pix) or ell.sma_pix <= 0:
        return np.nan, 0, 0, 0
    if not np.isfinite(ell.ba) or ell.ba <= 0:
        return np.nan, 0, 0, 0
    if not np.isfinite(ell.theta_deg):
        return np.nan, 0, 0, 0

    smb_pix = ell.sma_pix * ell.ba
    theta_rad = np.deg2rad(ell.theta_deg)

    aper = EllipticalAperture(
        (ell.xc, ell.yc),
        a=ell.sma_pix,
        b=smb_pix,
        theta=theta_rad,
    )

    aper_mask = aper.to_mask(method="center")
    aper_image = aper_mask.to_image(mask_bool.shape)

    if aper_image is None:
        return np.nan, 0, 0, 0

    inside = aper_image > 0
    n_total = int(np.sum(inside))

    if n_total == 0:
        return np.nan, 0, 0, 0


    n_masked = int(np.sum(mask_bool[inside]))
    n_unmasked = n_total - n_masked
    frac_masked = n_masked / n_total

    if n_total == 0:
        return MaskFractionResult(
            frac_masked=np.nan,
            n_total=0,
            n_masked=0,
            n_unmasked=0
            )
    else:
        return MaskFractionResult(
            frac_masked=frac_masked,
            n_total=n_total,
            n_masked=n_masked,
            n_unmasked=n_unmasked
            )

import numpy as np


def distance_to_nearest_mask(mask, xc, yc):
    """
    Distance from (xc, yc) to nearest masked pixel.

    Parameters
    ----------
    mask : 2D ndarray
        Boolean or integer mask image. True/nonzero = masked.

    xc, yc : float
        Pixel coordinates of galaxy center.

    Returns
    -------
    dist_pix : float
        Distance in pixels to nearest masked pixel.
        Returns np.nan if no masked pixels exist.
    """

    mask_bool = np.asarray(mask) > 0

    if not mask_bool.any():
        return np.nan

    y_mask, x_mask = np.where(mask_bool)

    dx = x_mask - xc
    dy = y_mask - yc

    dist2 = dx**2 + dy**2

    return np.sqrt(dist2.min())


from scipy.ndimage import label


def largest_mask_region(mask, ellipse_mask=None):
    """
    Size of largest contiguous masked region.

    Parameters
    ----------
    mask : 2D ndarray
        Boolean or integer mask image.

    ellipse_mask : ndarray, optional
        Boolean mask defining region of interest (e.g., ellipse).

    Returns
    -------
    max_pixels : int
        Size of largest connected masked region.
    """

    mask_bool = np.asarray(mask) > 0

    if ellipse_mask is not None:
        mask_bool = mask_bool & ellipse_mask

    if not mask_bool.any():
        return 0

    labeled, nlab = label(mask_bool)

    sizes = np.bincount(labeled.ravel())

    sizes[0] = 0  # background

    return sizes.max()
