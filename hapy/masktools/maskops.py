
# hapy/masktools/maskops.py
from __future__ import annotations
from typing import Iterable, List
import numpy as np

defaultcat='default.sex.HDI.mask'


def remove_central_objects(mask, sma=20, BA=1, PA=0, xc=None,yc=None):
    """ 
    find any pixels within central ellipse and set their values to zero 

    PARAMS:
    mask = 2D array containing masked pixels, like from SE segmentation image
    sma = semi-major axis in pixels
    BA = ratio of semi-minor to semi-major axes
    PA = position angle, measured in degree counter clockwise from +x axis

    OPTIONAL ARGS:
    xc = center of ellipse in pixels; assumed to be center of image if xc is not specified
    yc = center of ellipse in pixels; assumed to be center of image if yc is not specified
    
    RETURNS:
    newmask = copy of input mask, with pixels within ellipse set equal to zero

    """
    # changing the xmax and ymax - if the ellipse looks wrong, then swap back
    ymax,xmax = mask.shape
    # set center of ellipse as the center of the image
    if (xc is None) and (yc is None):
        xc,yc = xmax//2,ymax//2
    
    a = sma
    b = BA*sma
    phirad = np.radians(PA)

    X,Y = np.meshgrid(np.arange(xmax),np.arange(ymax))
    
    p1 = ((X-xc)*np.cos(phirad)+(Y-yc)*np.sin(phirad))**2/a**2
    p2 = ((X-xc)*np.sin(phirad)-(Y-yc)*np.cos(phirad))**2/b**2
    flag2 = p1+p2 < 1
    newmask = np.copy(mask)
    newmask[flag2] = 0
    # we could also get all the unique values associated with flag2, and then remove them
    ellipse_params = [xc,yc,sma,BA,phirad]
    return newmask,ellipse_params




def grow_mask_square(maskdat: np.ndarray, size: int = 7) -> np.ndarray:
    """
    Grow labeled regions by painting a size×size square around each masked pixel,
    preserving the pixel's label value.
    """
    out = np.array(maskdat, copy=True)
    nx, ny = out.shape
    half = int(size / 2)

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




def apply_user_masks(maskdat: np.ndarray, usr_mask: np.ndarray, deleted_objects: Iterable[int]) -> np.ndarray:
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




