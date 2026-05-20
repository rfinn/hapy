import logging
log = logging.getLogger(__name__)

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u


import os
import numpy as np
from pathlib import Path


#from . import utils
from hapy.imagetools.imutils import get_pixel_scale
#from hapy.imagetools.imutils import calculate_background_photutils
from hapy.imagetools.imutils import estimate_and_subtract_sky
from hapy.hatools.filter_transmission import (
    get_halpha_filtername,
    get_rband_filtername,
    get_filter_wavelength_info,
    filter_center_width,
)
instruments = ['BOK','INT','HDI','MOS']


# INT filters from https://astro.ing.iac.es/filter/list.php?instrument=WFC
# lmin is center - 0.5 FWHM
# lmin is center + 0.5 FWHM
#
# Filter # 197 is INT197 = (6568, 95)
# Filter # 227 is INT227 = (6657, 80)
#
# BOK using the NOAO filter, so Halpha 4
# 



def _safe_set_float_header(header, key, value, comment, ndec=6):
    if value is not None and np.isfinite(value):
        header[key] = (round(float(value), ndec), comment)
    else:
        header[key] = ("NaN", f"{comment} (not measurable)")
        
def fix_header_gain(input_header):
    header = input_header.copy()
    # if FIXGAIN is in header, return
    if "FIXGAIN" in header:
        return None

    # otherwise get EXPTIME and GAIN from header
    exptime = header.get("EXPTIME")
    gain = header.get("GAIN")

    # Missing keywords → warn and mark skipped
    if exptime is None or gain is None:
        log.warning(f"fix_gain: missing EXPTIME or GAIN (EXPTIME={exptime}, GAIN={gain}); leaving GAIN unchanged")
        header["FIXGAIN"] = (False, "GAIN not scaled (missing EXPTIME/GAIN)")
        return None

    # Coerce to floats
    try:
        exptime_f = float(exptime)
        gain_f = float(gain)
    except Exception:
        log.warning(f"fix_gain: non-numeric EXPTIME/GAIN (EXPTIME={exptime}, GAIN={gain}); leaving GAIN unchanged")
        header["FIXGAIN"] = (False, "GAIN not scaled (non-numeric EXPTIME/GAIN)")
        return

    # Invalid EXPTIME → warn and mark skipped
    if exptime_f <= 0:
        log.warning(f"fix_gain: invalid EXPTIME={exptime_f}; leaving GAIN unchanged")
        header["FIXGAIN"] = (False, "GAIN not scaled (invalid EXPTIME)")
        return

    # mv GAIN to GAINORIG
    header["GAINORIG"] = (gain, "Original detector gain (e-/ADU)")

    # set newgain to gain * exptime
    new_gain = gain * exptime
    header["GAIN"] = (new_gain, "Eff gain=GAINORIG*EXPTIME")

    header["FIXGAIN"] = (True, "GAIN was multiplied by EXPTIME")
    return header

def fix_header_exptime(input_header):
    header = input_header.copy()
    # if FIXGAIN is in header, return
    if "FIXEXPT" in header:
        return None

    # otherwise get EXPTIME and GAIN from header
    exptime = header.get("EXPTIME")
    
    # Missing keywords → warn and mark skipped
    if exptime is None:
        log.warning(f"fix_exptime: missing EXPTIME (EXPTIME={exptime}); leaving EXPTIME unchanged")
        header["FIXEXPT"] = (False, "EXPTIME not set to 1")
        return None

    # Coerce to floats
    try:
        exptime_f = float(exptime)
    except Exception:
        log.warning(f"fix_exptime: non numerical EXPTIME (EXPTIME={exptime}); leaving EXPTIME unchanged")
        header["FIXEXPT"] = (False, "EXPTIME not set to 1")
        return None

    # Invalid EXPTIME → warn and mark skipped
    if exptime_f <= 0:
        log.warning(f"fix_exptime: negative EXPTIME (EXPTIME={exptime}; leaving EXPTIME unchanged")
        header["FIXEXPT"] = (False, "EXPTIME not set to 1")
        return
        

    # mv GAIN to GAINORIG
    header["EXPTORIG"] = (exptime_f, "Original exposure time of coadd (sec)")

    # set newgain to gain * exptime
    header["EXPTIME"] = (1, "new exptime=1 to match PHOTZP")

    return header


def zp_scale_r_to_ha(zp_ha, zp_r):
    """Scale factor α so that CS = Ha - α * R."""
    if zp_ha is None or zp_r is None:
        return np.nan
    zp_ha = float(zp_ha)
    zp_r = float(zp_r)
    if not (np.isfinite(zp_ha) and np.isfinite(zp_r)):
        print("WARNING: could not calculate the zp scale!")
        return np.nan
    return float(10 ** (-0.4 * (zp_r - zp_ha)))




def _weight_ok(weight_cutout, central_frac=0.25, min_center_frac=0.5):
    """
    Decide whether the cutout is safely on valid data using the weight map.

    This checks the central region rather than requiring the exact central
    pixel to have positive weight, because isolated bad/zero weight pixels
    should not reject otherwise valid galaxy cutouts.
    """
    if weight_cutout is None:
        return True, {"center_weight": np.nan, "central_good_frac": np.nan}

    w = np.asarray(weight_cutout, dtype=float)
    good = np.isfinite(w) & (w > 0)

    ny, nx = w.shape
    yc, xc = ny // 2, nx // 2

    center_weight = w[yc, xc] if np.isfinite(w[yc, xc]) else np.nan

    halfx = max(1, int(0.5 * central_frac * nx))
    halfy = max(1, int(0.5 * central_frac * ny))

    x1 = max(0, xc - halfx)
    x2 = min(nx, xc + halfx + 1)
    y1 = max(0, yc - halfy)
    y2 = min(ny, yc + halfy + 1)

    central_good = good[y1:y2, x1:x2]
    central_good_frac = np.mean(central_good)

    #ok = (center_weight > 0) and (central_good_frac >= min_center_frac)
    # dropping some good galaxies, so trying this cut instead of one above
    ok = central_good_frac >= min_center_frac

    stats = {
        "center_weight": float(center_weight) if np.isfinite(center_weight) else np.nan,
        "central_good_frac": float(central_good_frac),
        "min_center_frac": float(min_center_frac),
    }
    return ok, stats

def ellipse_missing_fraction(weight_data, xc, yc, a, b, theta):
    """
    Fraction of pixels inside ellipse with zero weight.
    """
    from photutils.aperture import EllipticalAperture
    import numpy as np

    aper = EllipticalAperture((xc, yc), a, b, theta=theta)
    amask = aper.to_mask(method="center")
    mask_image = amask.to_image(weight_data.shape)

    if mask_image is None:
        return np.nan, 0

    inside = mask_image > 0
    npix = np.sum(inside)
    if npix == 0:
        return np.nan, 0

    missing = np.sum(weight_data[inside] <= 0)
    return missing / npix, npix

class CoaddImage:
    """
    Class to handle coadded images and create cutouts around galaxies.
    """
    def __init__(self, image_file, verbose=True):
        self.image_file = image_file
        self.pixelscale = None
        self.verbose = verbose

        self.data = None
        self.header = None
        self.wcs = None

        self.psf_image_name = None

        self.weight_image = None
        self.weight_flag = False

    def load_image(self):
        """
        Load FITS image and WCS.
        """
        self.image_file = str(self.image_file)

        self.header = fits.getheader(self.image_file)
        self.data = fits.getdata(self.image_file)
        self.wcs = WCS(self.header)
        self.pixelscale = get_pixel_scale(self.header)

        if self.verbose:
            print(
                f"Loaded image {self.image_file} with shape {self.data.shape} "
                f"and pixel scale {self.pixelscale:.3f} arcsec/pix"
            )

        weight_image = self.image_file.replace(".fits", ".weight.fits")

        if os.path.exists(weight_image):
            self.weight_image = weight_image
            self.weight_flag = True
        else:
            self.weight_image = None
            self.weight_flag = False


            
    def check_for_psf(self, psfdir=None): # MVC - model
        """
        check for 
        - psf file
        """
        # check for rband psf
        
        basename = os.path.basename(self.image_file)
        psf_image_name = basename.split('.fits')[0]+'-psf.fits'
        if psfdir is not None:
            psf_image_name = os.path.join(psfdir,psf_image_name)
        if os.path.exists(psf_image_name):
            self.psf_image_name = psf_image_name
        else:
            self.psf_image_name = None

    def get_instrument(self):
        for ii in instruments:
            if ii in self.image_file:
                self.instrument = ii
                break
    def get_filter(self):
        self.filter = self.header['FILTER']
    def get_filter_center_width(self):
        if self.filter in ["r","R","r SDSS k1018", "R Harris k1004"]:
            filter_file = get_rband_filtername(self.instrument, self.filter)
            fcenter_wave, fwidth = filter_center_width[filter_file]
        else:
            filter_file = get_halpha_filtername(self.instrument,self.filter)
            print("filter file = ",filter_file)
            fcenter_wave, fwidth = filter_center_width[filter_file]
        # add center wavelength and width to header
        self.header['FILTNAME'] = (filter_file, "Normalized filter filename")
        self.header["FILTWCEN"] = (fcenter_wave,"Filter Center Wavelength (A)")
        self.header["FILTWID"] = (fwidth, "Filter Width (A)")
        self.filter_file = filter_file
        self.filter_center = fcenter_wave
        self.filter_width = fwidth
    def get_target(self):
        self.target = self.header['OBJECT']
        
    def get_fwhm(self):
        try:
            self.fwhm_se_arcsec = float(self.header['SEFWHM'])
        except KeyError:
            self.fwhm_se_arcsec = None
        try:
            self.fwhm_psf_arcsec = float(self.header['FWHM'])*self.pixelscale
        except KeyError:
            self.fwhm_psf_arcsec = None

        # clean up header

        # replace SEFWHM with FWHM_SE
        self.header.rename_keyword("SEFWHM","FWHM_SE")#] = (self.fwhm_se_arcsec, "SE FWHM in arcsec")
        # add FWHM_PSF
        fwhm = round(float(self.fwhm_psf_arcsec),2) if self.fwhm_psf_arcsec is not None else None
        self.header["FWHM_PSF"] = (fwhm, "PSF FWHM in arcsec")
    def cutout_region_is_valid(self, ra, dec, size_arcsec, central_frac=0.25, min_center_frac=0.8):
        """
        Quick validity check using the weight image at a proposed cutout location.
        Returns (ok, stats_dict).
        """
        if not getattr(self, "weight_flag", False) or not getattr(self, "weight_image", None):
            return True, {"reason": "no_weight"}

        size_pix = int(size_arcsec / self.pixelscale)
        position = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

        try:
            wdata = fits.getdata(self.weight_image)
            wcut = Cutout2D(wdata, position=position, size=(size_pix, size_pix), wcs=self.wcs).data
        except Exception:
            return False, {"reason": "weight_cutout_failed"}

        ok, stats = _weight_ok(wcut, central_frac=central_frac, min_center_frac=min_center_frac)
        return ok, stats



    def get_ellipse_missing_fraction(self, xc, yc, a_arcsec, b_arcsec, theta_deg):
        """
        Measure the fraction of the full intended ellipse that is missing.

        Missing includes:
          1) ellipse pixels that fall off the image
          2) ellipse pixels on the image with weight == 0

        Parameters
        ----------
        xc, yc : float
            Ellipse center in image pixel coordinates.
        a_arcsec, b_arcsec : float
            Semi-major and semi-minor axes in arcseconds.
        theta_deg : float
            Ellipse position angle in degrees, in the convention expected by
            EllipticalAperture.

        Sets
        ----
        self.frac_missing : float
            Fraction of the full ellipse area that is missing.
        self.npix_total : int
            Total number of pixels in the full intended ellipse.
        self.npix_onimage : int
            Number of ellipse pixels that overlap the image.
        self.npix_good : int
            Number of ellipse pixels on the image with weight > 0.
        """
        from photutils.aperture import EllipticalAperture
        import numpy as np
        from astropy.io import fits

        self.frac_missing = np.nan
        self.npix_total = 0
        self.npix_onimage = 0
        self.npix_good = 0

        if not getattr(self, "weight_flag", False) or not getattr(self, "weight_image", None):
            return

        try:
            wdata = fits.getdata(self.weight_image)
        except Exception:
            return

        ny, nx = wdata.shape
        a = a_arcsec/self.pixelscale
        b = b_arcsec/self.pixelscale        
        if not np.all(np.isfinite([xc, yc, a, b, theta_deg])) or a <= 0 or b <= 0:
            return

        # Build a bounding box large enough to contain the full ellipse.
        # Add a small margin so edge pixels are not clipped.
        halfsize_x = int(np.ceil(a)) + 3
        halfsize_y = int(np.ceil(a)) + 3

        x0 = int(np.floor(xc)) - halfsize_x
        x1 = int(np.floor(xc)) + halfsize_x + 1
        y0 = int(np.floor(yc)) - halfsize_y
        y1 = int(np.floor(yc)) + halfsize_y + 1

        # Ellipse center in local bounding-box coordinates
        xc_local = xc - x0
        yc_local = yc - y0

        aper = EllipticalAperture((xc_local, yc_local), a, b, theta=np.deg2rad(theta_deg))
        amask = aper.to_mask(method="center")

        bbox_shape = (y1 - y0, x1 - x0)
        mask_image = amask.to_image(bbox_shape)

        if mask_image is None:
            return

        ellipse_full = mask_image > 0
        npix_total = int(np.sum(ellipse_full))
        self.npix_total = npix_total

        if npix_total == 0:
            return

        # Figure out overlap of the bounding box with the real image
        ix0 = max(0, x0)
        ix1 = min(nx, x1)
        iy0 = max(0, y0)
        iy1 = min(ny, y1)

        # No overlap with image at all
        if ix0 >= ix1 or iy0 >= iy1:
            self.npix_onimage = 0
            self.npix_good = 0
            self.frac_missing = 1.0
            return

        # Matching slices in the local ellipse-mask image
        mx0 = ix0 - x0
        mx1 = mx0 + (ix1 - ix0)
        my0 = iy0 - y0
        my1 = my0 + (iy1 - iy0)

        ellipse_onimage = ellipse_full[my0:my1, mx0:mx1]
        weight_onimage = wdata[iy0:iy1, ix0:ix1]

        npix_onimage = int(np.sum(ellipse_onimage))
        self.npix_onimage = npix_onimage

        if npix_onimage == 0:
            self.npix_good = 0
            self.frac_missing = 1.0
            return

        # Good coverage means ellipse pixel is on-image and weight > 0
        good = ellipse_onimage & np.isfinite(weight_onimage) & (weight_onimage > 0)
        npix_good = int(np.sum(good))
        self.npix_good = npix_good

        self.frac_missing = 1.0 - (npix_good / npix_total)


    def make_cutout(self, ra=None, dec=None, size_arcsec=None, output_name=None,subtract_sky=False, skycfg=None, return_cutout=True,fix_gain=True, fix_exptime=True, overwrite=False,return_slices=False, slices_original=None):
        """
        Create a science cutout either from RA/Dec/size or from an explicit
        parent-image slice.
  

        If a weight image exists, reject cutouts whose central region falls in
        a chip gap or off the valid image area.

        If a weight image exists, write a matched weight cutout using the same
        pixel slice and reject cutouts whose central weight region is invalid.

        """

        if output_name is not None:
            outpath = Path(output_name)
            if outpath.is_file() and not overwrite:
                print(f"Skipping existing cutout: {outpath}")
                if return_slices:
                    return ("skipped", None, None, None)
                return ("skipped", None, None)

        if self.data is None or self.wcs is None:
            self.load_image()

        # --------------------------------------------------
        # Define science cutout and slice
        # --------------------------------------------------
        if slices_original is None:
            if ra is None or dec is None or size_arcsec is None:
                raise ValueError("ra, dec, and size_arcsec are required when slices_original is not provided")

            size_pix = int(size_arcsec / self.pixelscale)
            position = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

            cutout = Cutout2D(
                self.data,
                position=position,
                size=(size_pix, size_pix),
                wcs=self.wcs,
            )
            cutout_data = cutout.data.copy()
            cutout_wcs = cutout.wcs
            slices_original = cutout.slices_original

        else:
            yslice, xslice = slices_original
            cutout_data = self.data[yslice, xslice].copy()
            cutout_wcs = self.wcs.slice((yslice, xslice))

        # --------------------------------------------------
        # Output name
        # --------------------------------------------------
        if output_name is None:
            base = os.path.basename(self.image_file).replace(".fits", "")
            output_name = f"{base}-cutout.fits"


        # --------------------------------------------------
        # Matched weight cutout
        # --------------------------------------------------
        wcut = None 
        weight_header = None
        weight_output_name = None

        if getattr(self, "weight_flag", False) and getattr(self, "weight_image", None):
            try:
                wdata, wheader = fits.getdata(self.weight_image, header=True)
                wcut_obj = Cutout2D(wdata, position=position, size=(size_pix, size_pix), wcs=self.wcs)
                wcut = wcut_obj.data.copy()
                weight_header = wheader.copy()
                weight_header.update(wcut_obj.wcs.to_header())
            except Exception:
                wcut = None
                weight_header = None

        # --- reject edge / chip-gap cases based on weight image ---
        weight_ok, weight_stats = _weight_ok(wcut, central_frac=0.25, min_center_frac=0.8)
        if not weight_ok:
            print(f"Rejecting cutout due to invalid central weight region: {output_name}")
            return ("invalid", None, None)


        # --------------------------------------------------
        # Optional sky subtraction
        # --------------------------------------------------
        skycfg = skycfg or {}
        med = None
        std = None
        if subtract_sky:
            cutout_data, med, std = estimate_and_subtract_sky(
                cutout_data,
                weightimage=wcut,
                subtract=subtract_sky,
                **skycfg,
            )

        if output_name is None:
            base = os.path.basename(self.image_file).replace('.fits', '')
            output_name = f"{base}-cutout.fits"

        # weight cutout filename
        if wcut is not None:
            if output_name.endswith(".fits"):
                weight_output_name = output_name.replace(".fits", ".weight.fits")
            else:
                weight_output_name = output_name + ".weight.fits"

        outheader = self.header.copy()
        outheader.update(cutout.wcs.to_header())

        if self.psf_image_name is not None:
            outheader.set('PSFIMAGE', self.psf_image_name)


        if getattr(self, "weight_image", None) is not None:
            outheader.set("WGTIMAGE", str(self.weight_image))
        

        # record weight diagnostics
        outheader["WGTCUTOK"] = (bool(weight_ok), "Central weight region valid")
        _safe_set_float_header(outheader, "WGTCNTR", weight_stats["center_weight"], "Central weight value")
        _safe_set_float_header(outheader, "WGTCFRAC", weight_stats["central_good_frac"], "Central good-weight fraction")

        if subtract_sky:
            parent_hdr = self.header
            if "SKYMED" in parent_hdr:
                outheader["PSKYMED"] = (self.header["SKYMED"], "Parent coadd SKYMED")
            if "SKYSTD" in parent_hdr:
                outheader["PSKYSTD"] = (self.header["SKYSTD"], "Parent coadd SKYSTD")

            outheader["CUTSKY"] = (True, "Sky median subtracted from cutout")
            outheader["CSKYSRC"] = ("CUTOUT", "Sky measured on cutout")
            _safe_set_float_header(outheader, "CSKYMED", med, "Median sky (ADU) subtracted")
            _safe_set_float_header(outheader, "CSKYSTD", std, "Std sky (ADU)")
            outheader["CSKYMETH"] = ("PHOTUTILS", "Background estimation method")
            outheader["CSKYOK"] = (bool(np.isfinite(med) and np.isfinite(std)), "Sky estimate finite")

        if fix_gain:
            newheader = fix_header_gain(outheader)
            if newheader is not None:
                outheader = newheader

        if fix_exptime:
            newheader = fix_header_exptime(outheader)
            if newheader is not None:
                outheader = newheader


        # print("DEBUG before write:",
        #     output_name,
        #     np.nanmin(cutout_data),
        #     np.nanmax(cutout_data),
        #     np.count_nonzero(cutout_data))
        # write science cutout
        fits.PrimaryHDU(data=cutout_data, header=outheader).writeto(output_name, overwrite=True)

        # write weight cutout if available
        if wcut is not None and weight_header is not None and weight_output_name is not None:
            weight_header["WGTCUT"] = (True, "This file is a weight-image cutout")
            weight_header["SCIIM"] = (os.path.basename(output_name), "science-image cutout")
            fits.PrimaryHDU(data=wcut, header=weight_header).writeto(weight_output_name, overwrite=True)

        if self.verbose:
            print(f"Cutout saved to {output_name}")
            if weight_output_name is not None:
                print(f"Weight cutout saved to {weight_output_name}")

        if return_cutout:
            if return_slices:
                return ("ok", cutout_data, outheader, slices_original)
            return ("ok", cutout_data, outheader)
        return None

    # #def make_cutout(self, ra, dec, size_arcsec, output_name=None):
    # def make_cutout(self, ra, dec, size_arcsec, output_name=None,
    #                 subtract_sky=False, skycfg=None, return_cutout=True, fix_gain=True, overwrite=False):

    #     """
    #     Create a cutout centered act (ra, dec) with size in arcsec.
   #     """

    #     # --- check if output_name exists, return if it does
    #     if output_name is not None:
    #         outpath = Path(output_name)
    #         if outpath.is_file() and not overwrite:
    #             print(f"Skipping existing cutout: {outpath}")
    #             return None
    
    #     if self.data is None or self.wcs is None:
    #         self.load_image()

    #     # convert size from arcsec to pixels
    #     size_pix = int(size_arcsec / self.pixelscale)
    #     #position = (ra, dec)
    #     position = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame="icrs")
    #     cutout = Cutout2D(self.data, position=position, size=(size_pix,size_pix), wcs=self.wcs)


    #     # --- optional weight cutout (for background masking) ---
    #     wcut = None
    #     if getattr(self, "weight_flag", False) and getattr(self, "weight_image", None):
    #         try:
    #             wdata = fits.getdata(self.weight_image)
    #             wcut = Cutout2D(wdata, position=position, size=(size_pix, size_pix), wcs=self.wcs).data
    #         except Exception:
    #             wcut = None

    #     # --- optional sky subtraction on the cutout ---
    #     skycfg = skycfg or {}
    #     if subtract_sky:
 
    #         cutout_data, med, std = estimate_and_subtract_sky(
    #             cutout.data,
    #             weightimage=wcut,
    #             subtract=subtract_sky
    #             )
        
    #     if output_name is None:
    #         base = os.path.basename(self.image_file).replace('.fits', '')
    #         output_name = f"{base}-cutout.fits"



    #     outheader = self.header.copy()
    #     outheader.update(cutout.wcs.to_header())
    #     #outheader = cutout.wcs.to_header()
        
    #     if self.psf_image_name is not None:
    #         outheader.set('PSFIMAGE',self.psf_image_name)

    #     # TODO - propagate header keywords to the cutout

        
    #     # record sky stats

    #     if subtract_sky:
    #         parent_hdr = self.header
    #         # check for previous values of SKYMED in header
    #         if "SKYMED" in parent_hdr:
    #             outheader["PSKYMED"] = (self.header["SKYMED"], "Parent coadd SKYMED")
    #         if "SKYSTD" in parent_hdr:
    #             outheader["PSKYSTD"] = (self.header["SKYSTD"], "Parent coadd SKYSTD")

    #         outheader["CUTSKY"] = (True, "Sky median subtracted from cutout")
    #         outheader["CSKYSRC"] = ("CUTOUT", "Sky measured on cutout")
    #         _safe_set_float_header(outheader, "CSKYMED", med, "Median sky (ADU) subtracted")
    #         _safe_set_float_header(outheader, "CSKYSTD", std, "Std sky (ADU)")
    #         #_safe_set_float_header(outheader, "CSKYMEA", mean, "Mean sky (ADU)")
    #         #outheader["CSKYMED"] = (float(med), "Median sky (ADU) subtracted")
    #         #outheader["CSKYSTD"] = (float(std), "Sigma-clipped sky std (ADU/pix)")
    #         outheader["CSKYMETH"] = ("PHOTUTILS", "Background estimation method")
    #         outheader["CSKYOK"] = (bool(np.isfinite(med) and np.isfinite(std)), "Sky estimate finite")
    #     # fix gain
    #     if fix_gain:
    #         newheader = fix_header_gain(outheader)
    #         if newheader is not None:
    #             outheader = newheader
    #     hdu = fits.PrimaryHDU(data=cutout_data, header=outheader)                    
    #     hdu.writeto(output_name, overwrite=True)

    #     if self.verbose:
    #         print(f"Cutout saved to {output_name}")

    #     if return_cutout:
    #         return cutout_data, outheader
    #     return None
        
class HalphaImageSet:
    def __init__(self, rcoadd_fname, hacoadd_fname, psfdir=None):
        self.rcoadd_fname = rcoadd_fname
        self.hacoadd_fname = hacoadd_fname
        self.r = CoaddImage(self.rcoadd_fname)
        self.h = CoaddImage(self.hacoadd_fname)

        self.cs_fname = hacoadd_fname.replace(".fits","-CS-ZP.fits")
        if os.path.exists(self.cs_fname):
            self.cs_flag = True
            self.cs = CoaddImage(self.cs_fname)
        else:
            print("WARNING: did not find CS image!")
            self.cs_flag = False
        self.psfdir = psfdir
        
    def load_coadds(self):

        self.r.load_image()
        self.h.load_image()

        if self.cs_flag:
            self.cs.load_image()
            
        self.r.check_for_psf(psfdir=self.psfdir)
        self.h.check_for_psf(psfdir=self.psfdir)

        self.r.get_instrument()
        self.h.get_instrument()
        self.r.get_filter()
        self.h.get_filter()    

        self.r.get_filter_center_width()
        self.h.get_filter_center_width()

        if self.cs_flag:
            self.cs.header["FILTNAME"] = self.h.header["FILTNAME"]
            self.cs.header["FILTWCEN"] = self.h.header["FILTWCEN"]
            self.cs.header["FILTWID"] = self.h.header["FILTWID"]
        # get fwhm if
        self.r.get_fwhm()
        self.h.get_fwhm()


        # fix gain
        #self.r.fix_gain()
        #self.h.fix_gain()
        #if self.cs_flag:
        #    self.cs.fix_gain()

    def cutout_location_is_valid(self, ra, dec, size_arcsec):
        """
        Require valid cutout region in both R and Halpha images.
        """
        ok_r, stats_r = self.r.cutout_region_is_valid(ra, dec, size_arcsec)
        if not ok_r:
            return False, "r_invalid"

        ok_h, stats_h = self.h.cutout_region_is_valid(ra, dec, size_arcsec)
        if not ok_h:
            return False, "ha_invalid"

        return True, "ok"

        
    def get_ellipse_coverage(self, xc, yc, a, b, theta_deg):
        self.r.get_ellipse_missing_fraction(xc, yc, a, b, theta_deg)
        self.h.get_ellipse_missing_fraction(xc, yc, a, b, theta_deg)

        self.frac_missing_r = self.r.frac_missing
        self.frac_missing_h = self.h.frac_missing
        self.max_frac_missing = max(self.frac_missing_r, self.frac_missing_h)

        self.ellipse_npix_total_r = self.r.npix_total
        self.ellipse_npix_total_h = self.h.npix_total
        self.ellipse_npix_onimage_r = self.r.npix_onimage
        self.ellipse_npix_onimage_h = self.h.npix_onimage
        self.ellipse_npix_good_r = self.r.npix_good
        self.ellipse_npix_good_h = self.h.npix_good
    
        
    def get_cutout_all_filters_old(self, ra, dec, size_arcsec, rootname):
        
        self.r.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-R.fits")
        self.h.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-Ha.fits")
        if self.cs_flag:
            self.cs.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-CS-ZP.fits")            



        
    def get_cutout_all_filters(self, ra, dec, size_arcsec, rootname, subtract_sky=False, overwrite=False):
        """
        Write R and Ha cutouts, optionally with local sky subtraction.

        Then write CS-ZP from the cutouts as:

            CS-ZP = Ha - scale * R

        where scale is derived from PHOTZP in the parent coadd headers.

        If subtract_sky=True, CS-ZP is built from the sky-subtracted cutouts.

        Existing outputs are skipped unless overwrite=True.
        """

        r_name = f"{rootname}-R.fits"
        h_name = f"{rootname}-Ha.fits"
        cs_name = f"{rootname}-CS-ZP.fits"

        r_path = Path(r_name)
        h_path = Path(h_name)
        cs_path = Path(cs_name)


        # --------------------------------------------------
        # Ha cutout
        # --------------------------------------------------
        h_status, h_data, h_hdr, h_slices = self.h.make_cutout(
            ra, dec, size_arcsec,
            output_name=h_name,
            subtract_sky=subtract_sky,
            overwrite=overwrite,
            return_slices=True
        )

        if h_status == "skipped":
            if h_path.is_file():
                h_data, h_hdr = fits.getdata(h_name, header=True)
            else:
                raise FileNotFoundError(f"Ha cutout was skipped but does not exist: {h_name}")

        elif h_status == "invalid":
            print()
            print(f"WARNING: Skipping galaxy because Ha cutout is invalid: {rootname}")
            print()
            return "invalid"

        elif h_status != "ok":
            raise RuntimeError(f"Unknown Ha make_cutout status: {h_status}")
        

        # --------------------------------------------------
        # R cutout
        # --------------------------------------------------
        r_status, r_data, r_hdr = self.r.make_cutout(
            ra, dec, size_arcsec,
            output_name=r_name,
            subtract_sky=subtract_sky,
            overwrite=overwrite,
            slices_original=h_slices,
        )

        if r_status == "skipped":
            if r_path.is_file():
                r_data, r_hdr = fits.getdata(r_name, header=True)
            else:
                raise FileNotFoundError(f"R cutout was skipped but does not exist: {r_name}")

        elif r_status == "invalid":
            print()
            print(f"WARNING: Skipping galaxy because R cutout is invalid: {rootname}")
            print()
            return "invalid"

        elif r_status != "ok":
            raise RuntimeError(f"Unknown R make_cutout status: {r_status}")


        # --------------------------------------------------
        # CS-ZP cutout: always build from R and Ha cutouts
        # --------------------------------------------------
        if cs_path.is_file() and not overwrite:
            print(f"Skipping existing cutout: {cs_path}")
            return "ok"

        zp_r = self.r.header.get("PHOTZP", np.nan)
        zp_h = self.h.header.get("PHOTZP", np.nan)
        scale = zp_scale_r_to_ha(zp_h, zp_r)

        cs_hdr = h_hdr.copy()
        cs_hdr["CSMAKE"] = (True, "Continuum-subtracted cutout written")
        cs_hdr["CSFROM"] = ("cutouts", "CS-ZP made from R/Ha cutouts")
        cs_hdr["CSFORM"] = ("Ha - a*R", "Continuum subtraction form")
        cs_hdr["CSSCALE"] = (float(scale) if np.isfinite(scale) else np.nan, "a: scale applied to R")
        cs_hdr["CSPHZPR"] = (float(zp_r) if np.isfinite(zp_r) else np.nan, "Parent R PHOTZP")
        cs_hdr["CSPHZPH"] = (float(zp_h) if np.isfinite(zp_h) else np.nan, "Parent Ha PHOTZP")
        cs_hdr["CSLSKY"] = (bool(subtract_sky), "Local sky-sub applied before CS")

        if np.isfinite(scale):
            cs_data = h_data.astype(float) - scale * r_data.astype(float)
        else:
            cs_data = np.full_like(h_data, np.nan, dtype=float)
            cs_hdr["CSMAKE"] = (False, "CS failed: missing/nonfinite PHOTZP")

        fits.PrimaryHDU(data=cs_data, header=cs_hdr).writeto(cs_name, overwrite=True)
        print(f"CS-ZP cutout saved to {cs_name}")

        return "ok"


if __name__ == "__main__":
    pass
    # initiate HalphaImageSet

    # initiate catalog

    # cull catalog to galaxies within FOV and filter

    # define size array
    
    # loop over galaxies

    # get_all_cutouts
