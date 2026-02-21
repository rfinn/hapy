from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
import os
import numpy as np

from . import utils

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
        self.pixelscale = utils.get_pixel_scale(self.header)
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

    def make_cutout(self, ra, dec, size_arcsec, output_name=None):
        """
        Create a cutout centered act (ra, dec) with size in arcsec.
        """
        if self.data is None or self.wcs is None:
            self.load_image()

        # convert size from arcsec to pixels
        size_pix = int(size_arcsec / self.pixelscale)
        position = (ra, dec)

        cutout = Cutout2D(self.data, position=position, size=size_pix, wcs=self.wcs)

        if output_name is None:
            base = os.path.basename(self.image_file).replace('.fits', '')
            output_name = f"{base}-cutout.fits"

        outheader = cutout.wcs.to_header()
        if self.psf_image_name is not None:
            outheader.set('PSFIMAGE',self.psf_image_name)

        # TODO - propagate header keywords to the cutout
        
        hdu = fits.PrimaryHDU(data=cutout.data, header=outheader)
        hdu.writeto(output_name, overwrite=True)

        if self.verbose:
            print(f"Cutout saved to {output_name}")

        #return cutout.data, cutout.wcs

        
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


    def get_all_cutouts(self, ra, dec, size_arcsec, rootname):
        
        t = self.r.make_cutout(self, ra, dec, size_arcsec, output_name=f"{rootname}-R.fits")
        t = self.h.make_cutout(self, ra, dec, size_arcsec, output_name=f"{rootname}-Ha.fits")
        if self.cs_flag:
            t = self.cs.make_cutout(self, ra, dec, size_arcsec, output_name=f"{rootname}-CS-ZP.fits")            
    
    


if __name__ == "__main__":
    pass
    # initiate HalphaImageSet

    # initiate catalog

    # cull catalog to galaxies within FOV and filter

    # define size array
    
    # loop over galaxies

    # get_all_cutouts

