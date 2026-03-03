from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u


import os
import numpy as np

#from . import utils
from hapy.imagetools.imutils import get_pixel_scale
#from hapy.imagetools.imutils import calculate_background_photutils
from hapy.imagetools.imutils import estimate_and_subtract_sky
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



def zp_scale_r_to_ha(zp_ha, zp_r):
    """Scale factor α so that CS = Ha - α * R."""
    if zp_ha is None or zp_r is None:
        return np.nan
    zp_ha = float(zp_ha)
    zp_r = float(zp_r)
    if not (np.isfinite(zp_ha) and np.isfinite(zp_r)):
        return np.nan
    return float(10 ** (-0.4 * (zp_ha - zp_r)))


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

    #def make_cutout(self, ra, dec, size_arcsec, output_name=None):
    def make_cutout(self, ra, dec, size_arcsec, output_name=None,
                    subtract_sky=False, skycfg=None, return_cutout=True):

        """
        Create a cutout centered act (ra, dec) with size in arcsec.
        """
        if self.data is None or self.wcs is None:
            self.load_image()

        # convert size from arcsec to pixels
        size_pix = int(size_arcsec / self.pixelscale)
        #position = (ra, dec)
        position = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame="icrs")
        cutout = Cutout2D(self.data, position=position, size=(size_pix,size_pix), wcs=self.wcs)


        # --- optional weight cutout (for background masking) ---
        wcut = None
        if getattr(self, "weight_flag", False) and getattr(self, "weight_image", None):
            try:
                wdata = fits.getdata(self.weight_image)
                wcut = Cutout2D(wdata, position=position, size=(size_pix, size_pix), wcs=self.wcs).data
            except Exception:
                wcut = None

        # --- optional sky subtraction on the cutout ---
        skycfg = skycfg or {}
        if subtract_sky:
 
            cutout_data, med, std = estimate_and_subtract_sky(
                cutout.data,
                weightimage=wcut,
                subtract=subtract_sky
                )
        
        if output_name is None:
            base = os.path.basename(self.image_file).replace('.fits', '')
            output_name = f"{base}-cutout.fits"



        outheader = self.header.copy()
        outheader.update(cutout.wcs.to_header())
        #outheader = cutout.wcs.to_header()
        
        if self.psf_image_name is not None:
            outheader.set('PSFIMAGE',self.psf_image_name)

        # TODO - propagate header keywords to the cutout

        
        # record sky stats

        if subtract_sky:
            parent_hdr = self.header
            # check for previous values of SKYMED in header
            if "SKYMED" in parent_hdr:
                outheader["PSKYMED"] = (self.header["SKYMED"], "Parent coadd SKYMED")
            if "SKYSTD" in parent_hdr:
                outheader["PSKYSTD"] = (self.header["SKYSTD"], "Parent coadd SKYSTD")

            outheader["CUTSKY"] = (True, "Sky median subtracted from cutout")
            outheader["SKYSRC"] = ("CUTOUT", "Sky measured on cutout")
            outheader["SKYMED"] = (float(med), "Median sky (ADU) subtracted")
            outheader["SKYSTD"] = (float(std), "Sigma-clipped sky std (ADU/pix)")
            outheader["SKYMETH"] = ("PHOTUTILS", "Background estimation method")
            
        hdu = fits.PrimaryHDU(data=cutout_data, header=outheader)                    
        hdu.writeto(output_name, overwrite=True)

        if self.verbose:
            print(f"Cutout saved to {output_name}")

        if return_cutout:
            return cutout_data, outheader
        return None
        
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
        
        self.h.get_instrument()
        self.h.get_filter()    

        # get fwhm if
        self.r.get_fwhm()
        self.h.get_fwhm()        
    def get_cutout_all_filters_old(self, ra, dec, size_arcsec, rootname):
        
        self.r.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-R.fits")
        self.h.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-Ha.fits")
        if self.cs_flag:
            self.cs.make_cutout(ra, dec, size_arcsec, output_name=f"{rootname}-CS-ZP.fits")            
    
    def get_cutout_all_filters(self, ra, dec, size_arcsec, rootname, subtract_sky=False):
        """
        Write sky-subtracted (optional) R and Ha cutouts. Then always (re)write CS-ZP cutout as:
            CS = Ha - scale * R
        where scale is derived from PHOTZP in the parent coadd headers.

        If subtract_sky=True, CS is built from the sky-subtracted cutouts.
        """
        # 1) Make R and Ha cutouts (optionally sky-subtracted locally)
        r_data, r_hdr = self.r.make_cutout(
            ra, dec, size_arcsec,
            output_name=f"{rootname}-R.fits",
            subtract_sky=subtract_sky
            )
        h_data, h_hdr = self.h.make_cutout(
            ra, dec, size_arcsec,
            output_name=f"{rootname}-Ha.fits",
            subtract_sky=subtract_sky
            )


        # 3) Compute ZP scale from parent coadd headers
        # (these are loaded when you call load_coadds())
        zp_r = self.r.header.get("PHOTZP", np.nan)
        zp_h = self.h.header.get("PHOTZP", np.nan) 
        scale = zp_scale_r_to_ha(zp_h, zp_r)

        # 4) Build CS-ZP cutout (overwrite OK)
        cs_hdr = h_hdr.copy()  # use Ha cutout header/WCS as reference

        cs_hdr["CSMAKE"] = (True, "Continuum-subtracted cutout written")
        cs_hdr["CSFORM"] = ("Ha - a*R", "Continuum subtraction form")
        cs_hdr["CSSCALE"] = (float(scale), "a: scale applied to R")
        cs_hdr["CSPHZPR"] = (float(zp_r) if np.isfinite(zp_r) else np.nan, "Parent R PHOTZP")
        cs_hdr["CSPHZPH"] = (float(zp_h) if np.isfinite(zp_h) else np.nan, "Parent Ha PHOTZP")
        cs_hdr["CSLSKY"] = (bool(subtract_sky), "Local sky-sub applied before CS")

        if np.isfinite(scale):
            cs_data = h_data - scale * r_data
        else:
            # no ZP -> write NaNs, but still create the file for pipeline consistency
            cs_data = np.full_like(h_data, np.nan, dtype=float)
            cs_hdr["CSMAKE"] = (False, "Continuum subtraction failed (missing PHOTZP)")

        fits.PrimaryHDU(data=cs_data, header=cs_hdr).writeto(f"{rootname}-CS-ZP.fits", overwrite=True)

        # 5) Optional: if you still want to cut out from a pre-made CS coadd, keep it behind a flag
        if self.cs_flag and not subtract_sky:
             self.cs.make_cutout(ra, dec, size_arcsec, output_name=cs_name)


if __name__ == "__main__":
    pass
    # initiate HalphaImageSet

    # initiate catalog

    # cull catalog to galaxies within FOV and filter

    # define size array
    
    # loop over galaxies

    # get_all_cutouts

