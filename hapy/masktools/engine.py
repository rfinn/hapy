# hapy/masktools/engine.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Tuple, List, Union, Callable

import numpy as np
from astropy.io import fits
#from astropy.wcs import WCS
from astropy import wcs
from .maskops import (
    remove_central_objects,
    apply_user_masks,
    grow_mask_square,
    circle_pixels,
)
from .gaia import get_gaia_stars, make_gaia_mask
from .sextractor import run_sextractor

from .types import EllipseParams


class MaskEngine:
    """
    Headless mask builder.

    Owns:
      - image data + header + WCS
      - segmentation mask (maskdat)
      - user masks / deleted IDs
      - optional Gaia mask
    """

    def __init__(
        self,
        image_fits: str,
        ha_image_fits: Optional[str] = None,
        sepath: Optional[str] = None,
        gaiapath: Optional[str] = None,
        config: str = "default.sex.HDI.mask",
        threshold: float = 0.005,
        snr: float = 10.0,
        snr_analysis: Optional[float] = None,
        minarea: int = 5,
        add_gaia_stars: bool = True,
        #ellipse_params=None,  
        verbose: bool = True,
        logger: Optional[Any] = None,
            ):
        self.image_fits = str(image_fits)
        self.ha_image_fits = str(ha_image_fits) if ha_image_fits else None

        #self.sepath = sepath
        #self.gaiapath = gaiapath
        self.config = config

        self.threshold = float(threshold)
        self.snr = float(snr)
        self.snr_analysis = float(snr_analysis) if snr_analysis is not None else float(snr)
        self.minarea = int(minarea)

        self.add_gaia_stars = bool(add_gaia_stars)
        #self.ellipse_params = ellipse_params
        self.verbose = verbose
        self.logger = logger

        # loaded image state
        self.image: np.ndarray
        self.header: fits.Header
        self.wcs: WCS

        self.xmax: int
        self.ymax: int
        self.xc: float
        self.yc: float

        # mask state
        self.maskdat: Optional[np.ndarray] = None
        self.usr_mask: Optional[np.ndarray] = None
        self.deleted_objects: List[int] = []
        self.gaia_mask: Optional[np.ndarray] = None

        # ellipse bookkeeping (for GUI to draw)
        self.ellipse_params: Optional[Union[EllipseParams, List[EllipseParams]]] = None

        self._load_image()


    def _progress(
        self,
        progress_callback: Optional[Callable[..., None]],
        *,
        stage: str,
        fraction: float,
        message: Optional[str] = None,
    ) -> None:
        """
        Headless progress hook. Never imports Qt. GUI can adapt this callback.
        """
        if progress_callback is not None:
            progress_callback(stage=stage, fraction=float(fraction), message=message)

    # ---------- core image load ----------
    def _load_image(self) -> None:
        self.image, self.header = fits.getdata(self.image_fits, header=True)
        self.wcs = wcs.WCS(self.header)
        pscale = wcs.utils.proj_plane_pixel_scales(self.wcs)
        self.pixel_scale_deg = pscale[0]
        self.ymax, self.xmax = self.image.shape
        self.xc = self.xmax / 2.0
        self.yc = self.ymax / 2.0
        self.usr_mask = np.zeros_like(self.image, dtype=float)

    # ---------- “pipeline” actions ----------
    def build_initial_mask(
        self,
        weightim: Optional[str] = None,
        weight_threshold: float = 1.0,
        remove_center_object: bool = True,
        center_object_id = None,
        galaxy_ellipse = None,
        grow_size: int = 5,
        grow_iterations: int = 3,
        output_prefix: Optional[str] = None,
        progress_callback=None,
        
        
    ):
        """
        Run detection (SE for now), optionally remove galaxy, then apply Gaia + user masks.
        Returns maskdat.
        """
        self._progress(progress_callback, stage="start", fraction=0.0, message="Starting mask build")

        self._progress(progress_callback, stage="sextractor", fraction=0.1, message="Running Source Extractor")
      
        segdata, catname, segmentation = run_sextractor(
            image_name=self.image_fits,
            #sepath=self.sepath,
            config=self.config,
            threshold=self.threshold,
            snr=self.snr,
            snr_analysis=self.snr_analysis,
            minarea=self.minarea,
            weight_image=weightim,
            weight_threshold=weight_threshold,
            
            #verbose=self.verbose,
            #output_prefix=output_prefix,
        )
        self.maskdat = segdata

        #print("after sextractor", self.maskdat.shape)
 
        self._progress(progress_callback, stage="cleanup", fraction=0.45, message="Applying object removals / user masks")
        # 2) remove center galaxy object(s)

        if remove_center_object:
            print("removing central object")
            self.maskdat, ellipse_params = remove_central_objects(
                self.maskdat,
                #center_object_id=center_object_id,
                ellipse_params=galaxy_ellipse,
             )
            self.ellipse_params = ellipse_params


        if self.usr_mask is not None:
            self.maskdat = apply_user_masks(self.maskdat, self.usr_mask)


        # GROW MASK
        self._progress(progress_callback, stage="grow", fraction=0.7, message="Growing mask")
        if grow_iterations > 0:
            self._progress(progress_callback, stage="grow", fraction=0.7,
                   message=f"Growing mask ({grow_iterations} iterations)")
            self.grow(size=grow_size, ngrow=grow_iterations, preserve_gaia=True)

        # ADD GAIA STARS
        if self.add_gaia_stars:# and self.gaiapath is not None:

            self._progress(progress_callback, stage="gaia", fraction=0.85, message="Adding Gaia stars")
            brightstar, x_pixels, y_pixels = get_gaia_stars(self.image_fits)            
            self.gaia_mask, star_masks = make_gaia_mask(self.maskdat, x_pixels, y_pixels, self.pixel_scale_deg, gaia_table=brightstar)
            #print("after gaia", self.maskdat.shape, " and gaia mask shape = ",self.gaia_mask.shape)
            self.maskdat = apply_user_masks(self.maskdat, self.gaia_mask)
            #print("after gaia - add back user masks", self.maskdat.shape)
        #print("after gaia", self.maskdat.shape)
        self._progress(progress_callback, stage="done", fraction=1.0, message="Mask build complete")
        return self.maskdat

        
    def add_circular_mask(self, x: float, y: float, radius_pix: float) -> int:
        """
        Add a circular user mask centered at (x, y) with radius in pixels.
        Returns the integer mask value used.
        """
        # lazily create a user mask if you don't already have one
        if getattr(self, "usr_mask", None) is None:
            self.usr_mask = np.zeros_like(self.maskdat, dtype=float)

        ymax, xmax = self.maskdat.shape

        # compute boolean pixel mask
        pixel_mask = circle_pixels(float(x), float(y), float(radius_pix), xmax, ymax)

        # choose next id
        mask_value = int(np.max(self.maskdat)) + 1

        self.usr_mask[pixel_mask] = mask_value
        self.maskdat = self.maskdat + self.usr_mask
        
        return mask_value
    # ---------- galaxy removal ----------
    def remove_object(self, objid):
        """Remove an object from the mask by zeroing all pixels with value obj_id."""
        if objid == 0:
            return

        # main mask
        self.maskdat[self.maskdat == objid] = 0.0

        # user mask bookkeeping (if you keep one)
        if getattr(self, "usr_mask", None) is not None:
            self.usr_mask[self.usr_mask == objid] = 0.0

        # optional: keep a list in the engine too
        if not hasattr(self, "deleted_objects"):
            self.deleted_objects = []
        self.deleted_objects.append(int(objid))

    
    def remove_galaxy_from_mask(
        self,
        center_object_id: Optional[int] = None,
        galaxy_ellipses: Optional[Union[EllipseParams, Sequence[EllipseParams]]] = None,
    ) -> None:
        """
        Remove the target galaxy from segmentation mask:
          - either by deleting a label (center_object_id)
          - and/or by clearing everything inside specified ellipses
        """
        if self.maskdat is None:
            raise RuntimeError("maskdat is None; run build_initial_mask() first.")

        # A) remove a specific segmentation label
        if center_object_id is not None:
            self.maskdat[self.maskdat == center_object_id] = 0

        # B) remove anything inside ellipses (good for shredded galaxies)
        if galaxy_ellipses is not None:
            if isinstance(galaxy_ellipses, EllipseParams):
                ell_list = [galaxy_ellipses]
            else:
                ell_list = list(galaxy_ellipses)

            self.ellipse_params = ell_list if len(ell_list) > 1 else ell_list[0]

            for e in ell_list:
                self.maskdat, _ = remove_central_objects(
                    self.maskdat,
                    sma=e.sma_pix,
                    BA=e.ba,
                    PA=e.pa_deg,
                    xc=e.xc,
                    yc=e.yc,
                )
        else:
            self.ellipse_params = None

    # ---------- Gaia ----------
    def apply_gaia_masks(self) -> None:
        """
        Compute and add Gaia masks (cached if already present).
        """
        if self.maskdat is None:
            raise RuntimeError("maskdat is None; run build_initial_mask() first.")

        if self.gaia_mask is None:
            gaia_tbl = get_gaia_stars(
                image_fits=self.image_fits,
                #gaiapath=self.gaiapath,
                cache_csv=True,
                verbose=self.verbose,
            )
            self.gaia_mask = make_gaia_mask(
                mask_shape=self.maskdat.shape,
                gaia_table=gaia_tbl,
                wcs=self.wcs,
                mask_value=int(np.max(self.maskdat)) + 100,
            )

        self.maskdat = self.maskdat + self.gaia_mask

    # ---------- user edits ----------
    def apply_user_masks_and_deletions(self) -> None:
        if self.maskdat is None or self.usr_mask is None:
            raise RuntimeError("maskdat/usr_mask missing; run build_initial_mask() first.")

        self.maskdat = apply_user_masks(self.maskdat, self.usr_mask, self.deleted_objects)

    def add_box_mask(self, x: float, y: float, size: float) -> int:
        """
        Adds a square user mask. Returns the new mask label/value used.
        """
        if self.maskdat is None or self.usr_mask is None:
            raise RuntimeError("maskdat/usr_mask missing; run build_initial_mask() first.")

        x = float(x); y = float(y); size = float(size)
        xmin = max(0, int(x - 0.5 * size))
        xmax = min(self.xmax, int(x + 0.5 * size))
        ymin = max(0, int(y - 0.5 * size))
        ymax = min(self.ymax, int(y + 0.5 * size))

        new_val = int(np.max(self.maskdat)) + 1
        self.usr_mask[ymin:ymax, xmin:xmax] = new_val
        self.apply_user_masks_and_deletions()
        return new_val

    def remove_object_id(self, obj_id: int) -> None:
        if self.maskdat is None:
            raise RuntimeError("maskdat missing; run build_initial_mask() first.")

        obj_id = int(obj_id)
        if obj_id <= 0:
            return
        self.deleted_objects.append(obj_id)
        self.apply_user_masks_and_deletions()


    def grow(self,
        size: int = 7,
        ngrow: int = 1,
        preserve_gaia: bool = True,
        ) -> None:
        """
        Grow mask by applying square dilation.
        
        Parameters
        ----------
        size : int
        Width of square dilation kernel.
        ngrow : int
        Number of growth iterations.
        preserve_gaia : bool
        If True, Gaia mask is not grown.
        """
        if self.maskdat is None:
            raise RuntimeError("maskdat missing; run build_initial_mask() first.")

        base = self.maskdat

        if preserve_gaia and self.gaia_mask is not None:
            base = base - self.gaia_mask

        grown = grow_mask_square(base, size=size, ngrow=ngrow)

        if preserve_gaia and self.gaia_mask is not None:
            grown = grown + self.gaia_mask

        self.maskdat = grown
        self.apply_user_masks_and_deletions()
        


    # ---------- outputs ----------
    def write_mask(self, mask_fits: str, inv_mask_fits: Optional[str] = None, overwrite: bool = True) -> None:
        if self.maskdat is None:
            raise RuntimeError("maskdat missing; run build_initial_mask() first.")
        #print("in write_mask, shape = ",self.maskdat.shape)
        fits.writeto(mask_fits, self.maskdat, header=self.header, overwrite=overwrite)

        if inv_mask_fits is not None:
            inv = np.array(~(self.maskdat > 0), dtype=int)
            fits.writeto(inv_mask_fits, inv, header=self.header, overwrite=overwrite)
    def auto_mask(
        self,
        mask_path,
        ngrow=3,
        weightim=None,
        weight_threshold=None,
        progress_callback=None,
    ):
        """
        Fully automatic mask creation:
          - build initial segmentation mask
          - remove central object (if configured)
          - grow mask `ngrow` times
          - add Gaia stars (optional)
          - write mask to disk

        PARAMS:
        - mask_path: filename for output mask
        - ngrow: number of times to grow the masked regions from SE (before adding stars)
        - weightim: weight image

        RETURNS
        - mask_path: name of output mask file that is written to disk
        """

        # 1. Build initial mask
        self.build_initial_mask(
            weightim=weightim,
            weight_threshold=weight_threshold,
            progress_callback=progress_callback,
        )

        # 2. Remove central object (if engine supports it)
        if hasattr(self, "remove_center_object_flag") and self.remove_center_object_flag:
            self.remove_center_object()

        # 3. Grow mask
        for i in range(ngrow):
            if progress_callback:
                progress_callback("grow", i / max(1, ngrow), f"Growing mask {i+1}/{ngrow}")
            self.grow_mask()

        # 4. Add Gaia stars
        if self.add_gaia_stars:
            self.add_gaia_stars_to_mask()

        # 5. Write mask
        self.write_mask(mask_path)

        if progress_callback:
            progress_callback("done", 1.0, "Mask build complete")

        return mask_path
