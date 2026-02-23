# hapy/hatools/hasession.py
import os

from hapy.masktools.api import MaskEngine  # adjust import if needed


class HaSession:
    """
    Headless Halpha processing session.

    Goal: provide a non-interactive pipeline that can also be called by the GUI.
    """

    def __init__(self, args, sepath=None):
        self.args = args

        # flags
        self.auto = bool(getattr(args, "auto", False))
        self.verbose = bool(getattr(args, "verbose", False))
        self.virgo = bool(getattr(args, "virgo", False))
        self.uat = bool(getattr(args, "uat", False))
        self.nebula = bool(getattr(args, "nebula", False))
        self.laptop = bool(getattr(args, "laptop", False))
        self.draco = bool(getattr(args, "draco", False))
        self.testing = bool(getattr(args, "testing", False))

        # inputs
        self.rcoadd_fname = getattr(args, "rimage", None)
        self.hacoadd_fname = getattr(args, "haimage", None)
        self.cscoadd_fname = getattr(args, "csimage", None)
        self.filter = getattr(args, "filter", None)
        self.prefix = getattr(args, "prefix", None)

        # paths
        self.sepath = sepath or (os.getenv("HOME") + "/github/halphagui/astromatic/")
        self.psfdirectory = getattr(args, "psfdir", None) or os.getcwd()
        self.tabledir = getattr(args, "tabledir", None)

        # per-galaxy / runtime state
        self.igal = None
        self.bad_galaxy = False

        # these will be filled later by your existing catalog logic
        self.gximage = []           # list of galaxies in FOV
        self.radius_arcsec = None   # per galaxy radius array
        self.BA = None              # per galaxy BA array
        self.PA = None              # per galaxy PA array
        self.defcat = None          # default catalog object with .cat['RA'], .cat['DEC']

        # filenames created during processing
        self.cutout_name_r = None
        self.cutout_name_ha = None

    # ------------------------
    # Top-level pipeline
    # ------------------------
    def auto_run(self):
        """
        Run the full pipeline headlessly.
        This is basically your existing hafunctions.auto_run(), moved here.
        """
        # These calls are placeholders—hook them to your existing implementations:
        self.read_rcoadd()
        self.read_hacoadd()

        if self.filter is not None:
            self.set_hafilter(self.filter)

        self.get_filter_ratio()
        self.subtract_images(overwrite=True)
        self.build_psf()
        self.find_galaxies()

        self.write_fits_table()

        if self.verbose:
            print(f"starting processing of each galaxy: {len(self.gximage)}")

        for i in range(len(self.gximage)):
            self.igal = i
            self.auto_gal()
            self.write_fits_table()

            if self.verbose:
                print(f"Finished galaxy {i+1}/{len(self.gximage)}")

    def auto_gal(self):
        """
        Process one galaxy headlessly.
        Replaces the old maskwindow(None, None, ...) call with MaskEngine.
        """
        self.bad_galaxy = False

        self.get_galaxy_cutout()
        if self.bad_galaxy:
            if self.verbose:
                print(f"Skipping galaxy {self.igal}: bad cutout/sky stats")
            return

        # --- Build objparams / ellipse for masking ---
        # This matches what you had:
        ra = float(self.defcat.cat["RA"][self.igal])
        dec = float(self.defcat.cat["DEC"][self.igal])

        mask_scalefactor = 1.0  # keep your existing value if different
        sma_arcsec = float(mask_scalefactor * self.radius_arcsec[self.igal])
        ba = float(self.BA[self.igal])
        pa_deg = float(self.PA[self.igal] + 90)  # your convention

        # Old code used list: [RA, DEC, SMA_arcsec, BA, PA_deg]
        galaxy_ellipse = [ra, dec, sma_arcsec, ba, pa_deg]

        if self.verbose:
            print("Building mask headlessly for", self.cutout_name_r)

        # Continue your pipeline here:
        # - run galfit
        # - galfit_ellip_phot
        # - photutils_ellip_phot
        # Keep those methods in HaSession and call them:
        self.run_galfit(ncomp=1, ha=False)
        self.galfit_ellip_phot()
        self.photutils_ellip_phot()


        
class HaSession:
    def __init__(self, args, sepath=None):

        # flags
        self.auto = bool(getattr(args, "auto", False))
        self.verbose = bool(getattr(args, "verbose", False))
        self.virgo = bool(getattr(args, "virgo", False))
        self.uat = bool(getattr(args, "uat", False))
        self.nebula = bool(getattr(args, "nebula", False))
        self.laptop = bool(getattr(args, "laptop", False))
        self.draco = bool(getattr(args, "draco", False))
        self.testing = bool(getattr(args, "testing", False))

        # inputs
        self.rcoadd_fname = getattr(args, "rimage", None)
        self.hacoadd_fname = getattr(args, "haimage", None)
        self.cscoadd_fname = getattr(args, "csimage", None)
        self.filter = getattr(args, "filter", None)
        self.prefix = getattr(args, "prefix", None)

        # paths
        #self.sepath = sepath or (os.getenv("HOME") + "/github/halphagui/astromatic/")
        self.psfdirectory = getattr(args, "psfdir", None) or os.getcwd()
        self.tabledir = getattr(args, "tabledir", None)

        # per-galaxy / runtime state
        self.igal = None
        self.bad_galaxy = False

        # per-galaxy / runtime state
        self.igal = None
        self.bad_galaxy = False

        # these will be filled later by your existing catalog logic
        self.gximage = []           # list of galaxies in FOV
        self.radius_arcsec = None   # per galaxy radius array
        self.BA = None              # per galaxy BA array
        self.PA = None              # per galaxy PA array
        self.defcat = None          # default catalog object with .cat['RA'], .cat['DEC']

        # filenames created during processing
        self.cutout_name_r = None
        self.cutout_name_ha = None
   
        self.oversampling = 2        

        if args.psfdir is None:
            self.psfdirectory = os.getcwd()
        else:
            self.psfdirectory = args.psfdir

        ############################################################
        ### CHECK TO SEE IF IMAGE NAMES ARE SPECIFIED
        ############################################################
        if self.filter is not None:
            self.set_hafilter(self.filter)
        if self.rcoadd_fname is not None:
            self.mask_image = self.rcoadd_fname.replace('.fits','-mask.fits')
        else:
            self.mask_image = None
            
    def auto_run(self):
        # run the analysis without starting the gui
        # read in coadded images
        self.read_rcoadd()
        self.read_hacoadd()
        # set filter
        self.set_hafilter(self.filter)
        
        # get filter ratio
        # this will look for filter ratio in r-band image header
        # will run SE if ratio is not found in header
        self.get_filter_ratio()
        
        # subtract images
        # this will look for CS-ZP.fits image
        # if not found, it will subtract images according to the ZP ratio of r and halpha images
        self.subtract_images(overwrite=True) # checks out ok
        
        # measure psf
        # function will check if psf image already exists
        self.build_psf()

        # get galaxies
        self.find_galaxies()

        self.write_fits_table()
        # analyze galaxies
        # the default catalog has been trimmed to only include
        # galaxies in FOV

        # skipping for now
        # just making cutouts and getting galaxies in FOV

        # need to get a list of RA, DEC, SMA, BA, PA to feed into masking routine
        if self.verbose:
            print("starting processing of each galaxy ",len(self.gximage))
        for i in range(len(self.gximage)):
        #for i in [2]: # for testing
            self.igal = i
            # get cutouts
            self.auto_gal()
            if self.verbose:
                print("writing fits table after running auto_gal")
            self.write_fits_table()
            if self.verbose:
                print(f"##########################\nFinished galaxy {i+1}/{len(self.gximage)}")
    def _make_mask(self):
        #########################################        
        # START YOUR MASK ENGINE!
        # initialize and build first-pass mask
        #########################################
        from hapy.masktools.api import MaskEngine, EllipseParams


        engine = MaskEngine(
            image_fits=self.cutout_name_r,
            ha_image_fits=self.cutout_name_ha,
            config="default.sex.HDI.mask",  # or track this from args
            threshold=0.005,
            snr=10,
            minarea=5,
            add_gaia_stars=True,
        )

        def _progress(stage, fraction, message=None):
            if self.verbose:
                print(f"[mask] {stage} {fraction:0.2f} {message or ''}".strip())
 
        self.maskdat = engine.build_initial_mask(
            galaxy_ellipse = self.objparams,
            remove_center_object=True,
            grow_size=5,
            grow_iterations=3,
            progress_callback=_progress,
        )
        self.engine.write_mask(self.mask_image)
    def auto_gal(self):
        # run the analysis on an individual galaxy

        self.bad_galaxy = False
        # create cutout
 
        self.get_galaxy_cutout()

        if self.bad_galaxy:
            if self.verbose:
                print("\nzero std in sky - this is not real, so skipping galaxy ",self.cutout_name_r)
                print()
            return

        # --- Build objparams / ellipse for masking ---
        # create mask
        if self.uat:
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]
        else:
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]

        self._make_mask()


  
        self.run_galfit(ncomp=1,ha=False)
        self.galfit_ellip_phot()
        self.photutils_ellip_phot()
        
            
    # ------------------------
    # Stubs to connect to your existing methods
    # ------------------------
    def read_rcoadd(self): raise NotImplementedError
    def read_hacoadd(self): raise NotImplementedError
    def set_hafilter(self, f): raise NotImplementedError
    def get_filter_ratio(self): raise NotImplementedError
    def subtract_images(self, overwrite=False): raise NotImplementedError
    def build_psf(self): raise NotImplementedError
    def find_galaxies(self): raise NotImplementedError
    def write_fits_table(self): raise NotImplementedError
    def get_galaxy_cutout(self): raise NotImplementedError
    def run_galfit(self, ncomp=1, ha=False): raise NotImplementedError
    def galfit_ellip_phot(self): raise NotImplementedError
    def photutils_ellip_phot(self): raise NotImplementedError
        
