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



def _safe_set_float_header(header, key, value, comment, ndec=4):
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


def _weight_ok(weight_cutout, central_frac=0.25, min_center_frac=0.8):
    """
    Decide whether the cutout is safely on valid data.

    Parameters
    ----------
    weight_cutout : 2D array
        Weight map cutout.
    central_frac : float
        Fractional size of central box to test.
    min_center_frac : float
        Minimum fraction of pixels in the central box that must have weight > 0.

    Returns
    -------
    ok : bool
    stats : dict
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

    ok = (center_weight > 0) and (central_good_frac >= min_center_frac)

    stats = {
        "center_weight": float(center_weight) if np.isfinite(center_weight) else np.nan,
        "central_good_frac": float(central_good_frac),
    }
    return ok, stats

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

    def load_image(self):
        """
        Load FITS image and WCS.
        """
        self.header = fits.getheader(self.image_file)
        self.data = fits.getdata(self.image_file)
        self.wcs = WCS(self.header)
        self.pixelscale = get_pixel_scale(self.header)
        if self.verbose:
            print(f"Loaded image {self.image_file} with shape {self.data.shape} and pixel scale {self.pixelscale:.3f} arcsec/pix")

        weight_image = self.image_file.replace('.fits','.weight.fits')

        if os.path.exists(weight_image):
            self.weight_image = weight_image
            self.weight_flag = True
        else:
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
            self.fwhm_arcsec = float(self.header['SEFWHM'])
        except KeyError:
            self.fwhm_arcsec = None
        try:
            self.fwhm_pixels = float(self.header['FWHM'])
        except KeyError:
            self.fwhm_pixels = None

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
        
    def make_cutout(self, ra, dec, size_arcsec, output_name=None,
                    subtract_sky=False, skycfg=None, return_cutout=True,
                    fix_gain=True, fix_exptime=True, overwrite=False):

        """
        Create a cutout centered at (ra, dec) with size in arcsec.

        If a weight image exists, reject cutouts whose central region falls in
        a chip gap or off the valid image area.

        If a weight image exists and the cutout is valid, also write a weight-image
        cutout with filename output_name.replace(".fits", ".weight.fits").
        """

        if output_name is not None:
            outpath = Path(output_name)
            if outpath.is_file() and not overwrite:
                print(f"Skipping existing cutout: {outpath}")
                return ("skipped", None, None)

        if self.data is None or self.wcs is None:
            self.load_image()

        size_pix = int(size_arcsec / self.pixelscale)
        position = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

        cutout = Cutout2D(self.data, position=position, size=(size_pix, size_pix), wcs=self.wcs)
        cutout_data = cutout.data.copy()

        # --- optional weight cutout ---
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

        # --- optional sky subtraction ---
        skycfg = skycfg or {}
        med = None
        std = None
        if subtract_sky:
            cutout_data, med, std = estimate_and_subtract_sky(
                cutout.data,
                weightimage=wcut,
                subtract=subtract_sky
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

        if self.weight_image is not None:
            outheader.set('WGTIMAGE', str(self.weight_image))

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

    def get_cutout_all_filters_old(self, ra, dec, size_arcsec, rootname):
        
        self.r.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-R.fits")
        self.h.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-Ha.fits")
        if self.cs_flag:
            self.cs.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-CS-ZP.fits")            

    def get_cutout_all_filters(self, ra, dec, size_arcsec, rootname, subtract_sky=False, overwrite=False):
        """
        Write sky-subtracted (optional) R and Ha cutouts. Then always write CS-ZP cutout as:
            CS = Ha - scale * R
        where scale is derived from PHOTZP in the parent coadd headers.

        If subtract_sky=True, CS is built from the sky-subtracted cutouts.

        Existing outputs are skipped unless overwrite=True.
        """

        r_name = f"{rootname}-R.fits"
        h_name = f"{rootname}-Ha.fits"
        cs_name = f"{rootname}-CS-ZP.fits"

        r_path = Path(r_name)
        h_path = Path(h_name)
        cs_path = Path(cs_name)

        # --------------------------------------------------
        # R cutout
        # --------------------------------------------------

        status, r_data, r_hdr = self.r.make_cutout(
            ra, dec, size_arcsec,
            output_name=r_name,
            subtract_sky=subtract_sky,
            overwrite=overwrite
            )

        if status == "skipped":
            if r_path.is_file():
                r_data, r_hdr = fits.getdata(r_name, header=True)
            else:
                raise FileNotFoundError(f"R cutout was skipped but does not exist: {r_name}")

        elif status == "invalid":
            print(f"Skipping galaxy because R cutout is invalid: {rootname}")
            return

        elif status == "ok":
            pass

        else:
            raise RuntimeError(f"Unknown R make_cutout status: {r_status}")

        # --------------------------------------------------
        # Ha cutout
        # --------------------------------------------------

        h_status, h_data, h_hdr = self.h.make_cutout(
            ra, dec, size_arcsec,
            output_name=h_name,
            subtract_sky=subtract_sky,
            overwrite=overwrite
        )

        if h_status == "skipped":
            if h_path.is_file():
                h_data, h_hdr = fits.getdata(h_name, header=True)
            else:
                raise FileNotFoundError(f"Ha cutout was skipped but does not exist: {h_name}")

        elif h_status == "invalid":
            print(f"Skipping galaxy because Ha cutout is invalid: {rootname}")
            return "invalid"

        elif h_status != "ok":
            raise RuntimeError(f"Unknown Ha make_cutout status: {h_status}")


        # --------------------------------------------------
        # CS cutout
        # --------------------------------------------------
        if cs_path.is_file() and not overwrite:
            print(f"Skipping existing cutout: {cs_path}")
        else:
            zp_r = self.r.header.get("PHOTZP", np.nan)
            zp_h = self.h.header.get("PHOTZP", np.nan)
            scale = zp_scale_r_to_ha(zp_h, zp_r)

            cs_hdr = h_hdr.copy()
            cs_hdr["CSMAKE"] = (True, "Continuum-subtracted cutout written")
            cs_hdr["CSFORM"] = ("Ha - a*R", "Continuum subtraction form")
            cs_hdr["CSSCALE"] = (float(scale), "a: scale applied to R")
            cs_hdr["CSPHZPR"] = (float(zp_r) if np.isfinite(zp_r) else np.nan, "Parent R PHOTZP")
            cs_hdr["CSPHZPH"] = (float(zp_h) if np.isfinite(zp_h) else np.nan, "Parent Ha PHOTZP")
            cs_hdr["CSLSKY"] = (bool(subtract_sky), "Local sky-sub applied before CS")

            if np.isfinite(scale):
                cs_data = h_data - scale * r_data
            else:
                cs_data = np.full_like(h_data, np.nan, dtype=float)
                cs_hdr["CSMAKE"] = (False, "Continuum subtraction failed (missing PHOTZP)")

            fits.PrimaryHDU(data=cs_data, header=cs_hdr).writeto(cs_name, overwrite=True)

        # --------------------------------------------------
        # Optional legacy CS coadd cutout
        # --------------------------------------------------
        if self.cs_flag and not subtract_sky:
            self.cs.make_cutout(
                ra, dec, size_arcsec,
                output_name=cs_name,
                overwrite=overwrite
            )

        
  

if __name__ == "__main__":
    pass
    # initiate HalphaImageSet

    # initiate catalog

    # cull catalog to galaxies within FOV and filter

    # define size array
    
    # loop over galaxies

    # get_all_cutouts

