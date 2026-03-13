#!/usr/bin/env python

import os
import numpy as np
from scipy import interpolate

import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.table import Table

wave_halpha = 6563. # angstrom
from . import utils


# --- Filter central wavelength and width --- #
#
# We calculate these in github/filter_transformations/filtertrans-dev.ipynb
# Becky cross checked with plots that she has

filter_center_width = {
    'BOK90prime-BASSr.fits':(6410.8, 1398.8),
    'BOK90prime-Ha+4nm.fits':(6620.8, 83.3),
    'MOS-Ha+12nm.fits':(6698.8, 86.1),
    'MOS-Ha+16nm.fits':(6730.8, 83.9),
    'MOS-Ha+4nm.fits':(6620.8, 83.3),
    'MOS-Ha+8nm.fits':(6654.4, 84.1),
    'MOS-HarrisR.fits':(6653.9, 1551.1),
    'MOS-SDSSr.fits':(6287.6, 1382.5),
    'HDI-Ha+12nm.fits':(6701.7, 61.5),
    'HDI-Ha+16nm.fits':(6742.1, 59.3),
    'HDI-Ha+4nm.fits':(6618.5, 60.4),
    'HDI-Ha+8nm.fits':(6660.0, 59.8),
    'HDI-Ha.fits':(6580.0, 58.8),
    'HDI-HarrisR.fits':(6605.3, 1576.1),
    'HDI-SDSSr.fits':(6242.3, 1425.5),
    'WFC-Ha-197.fits':(6568.0, 92.9),
    'WFC-Ha-227.fits':(6666.2, 80.9),
    'WFC-SDSSr-214.fits':(6230.1, 1219.8),
    'panstarrs-g.fits':(4866.5, 1166.4),
    'panstarrs-r.fits':(6214.6, 1318.0),
}

instrument_to_prefix = {'INT':'WFC','BOK':'90prime','HDI':'HDI','MOS':'MOS'}

hafilter_to_suffix = {4:'Ha+4nm',8:'Ha+8nm',12:'Ha+12nm',16:'Ha+16nm',\
                        '4':'Ha+4nm','8':'Ha+8nm','12':'Ha+12nm','16':'Ha+16nm',\
                        'ha4':'Ha+4nm','ha8':'Ha+8nm','ha12':'Ha+12nm','ha16':'Ha+16nm','ha':'Ha-197',\
                        'Ha+4nm':'Ha+4nm','Ha4nm':'Ha+4nm',\
                        'Halpha':'Ha-197','Ha6657':'Ha-227',\
                        'ha4 H-alpha+4nm k1010':'Ha+4nm',
                        'ha8 H-alpha+8nm k1011':'Ha+8nm',
                        'ha12 H-alpha+12nm k1012':'Ha+12nm',
                        'ha16 H-alpha+16nm k1013':'Ha+16nm',}

def get_halpha_filtername(instrument, hfilter):
    print("testing, self.hafilter = ",hfilter, instrument)
    halpha_filtername = f"{instrument_to_prefix[instrument]}-{hafilter_to_suffix[hfilter]}.fits"
    return halpha_filtername

def get_rband_filtername(instrument, rfilter):
    """
    Normalize instrument + FILTER keyword into the canonical r-band filter filename.
    """

    instrument = str(instrument).strip()
    rfilter = str(rfilter).strip()

    # exact / common aliases
    filter_map = {
        ("BOK", "r"): "BOK90prime-BASSr.fits",
        #("BOK", "BASSr"): "BOK90prime-BASSr.fits",
        #("BOK", "BASS-r"): "BOK90prime-BASSr.fits",

        ("HDI", "R"): "HDI-HarrisR.fits",
        ("HDI", "r"): "HDI-SDSSr.fits",
        #("HDI", "HarrisR"): "HDI-HarrisR.fits",
        #("HDI", "SDSSr"): "HDI-SDSSr.fits",

        #("MOS", "R"): "MOS-HarrisR.fits",
        ("MOS", "r SDSS k1018"): "MOS-SDSSr.fits",
        ("MOS", "R Harris k1004"): "MOS-HarrisR.fits",
        #("MOS", "SDSSr"): "MOS-SDSSr.fits",

        ("INT", "r"): "WFC-SDSSr-214.fits",
        #("INT", "SDSSr"): "WFC-SDSSr-214.fits",
        #("WFC", "r"): "WFC-SDSSr-214.fits",
        #("WFC", "SDSSr"): "WFC-SDSSr-214.fits",

        ("PANSTARRS", "r"): "panstarrs-r.fits",
        ("PS1", "r"): "panstarrs-r.fits",
    }

    key = (instrument, rfilter)
    if key in filter_map:
        return filter_map[key]

    # more permissive fallback logic
    rf = rfilter.lower()

    if instrument == "BOK":
        return "BOK90prime-BASSr.fits"

    if instrument == "HDI":
        if "harris" in rf or rf == "r":
            return "HDI-HarrisR.fits" if rfilter == "R" else "HDI-SDSSr.fits"

    if instrument == "MOS":
        if "harris" in rf or rfilter == "R":
            return "MOS-HarrisR.fits"
        if "sdss" in rf or rfilter == "r":
            return "MOS-SDSSr.fits"

    if instrument in ("INT", "WFC"):
        return "WFC-SDSSr-214.fits"

    if instrument in ("PANSTARRS", "PS1"):
        return "panstarrs-r.fits"

    raise KeyError(f"Could not determine r-band filter filename for instrument={instrument}, FILTER={rfilter}")


def get_filter_wavelength_info(filter_name):
    if filter_name not in filter_wavelengths:
        raise KeyError(f"No wavelength info for filter {filter_name}")
    return filter_wavelengths[filter_name]

class FilterTrace():
    def __init__(self,hafilter,filterpath=None,instrument=None,mintrans=10.):
        '''
        hafilter can be 4, 8, 12, 16, inthalpha, or intha6657

        filter path should point to filter_trace subdirectory of github/halphagui

        mintrans is the minimum transmission to use for selecting galaxies
        '''
        # read in filter trace
        # assume github directory is off main dir
        self.hafilter = hafilter
        if filterpath == None:
            self.filterpath = os.getenv('HOME')+'/github/halphagui/filter_traces/'
        else:
            self.filterpath = filterpath
        self.instrument = instrument

        self.read_filter()
        self.get_filter_properties()

        
    def get_halpha_filtername(self):
        print("testing, self.hafilter = ",self.hafilter, self.instrument)
        self.halpha_filtername = utils.get_filter_file(f"{instrument_to_prefix[self.instrument]}-{hafilter_to_suffix[self.hafilter]}.fits")

    def read_filter(self):
        """ updating to use the new filter curves """

        if self.hafilter == 'sienaha':
            # wavelength in nm, transmission in percent
            filterfile = self.filterpath+'/chroma-halpha-transmission-ascii.txt'
            wavescale=10
            
            tab = ascii.read(filterfile)
            self.wave = tab['col1']*wavescale # wavelength in angstrom
            self.trans = tab['col2'] # transmission percent
            self.trans = 100*self.trans
            

        else:
            self.get_halpha_filtername()
            ftab = Table.read(self.halpha_filtername)
            self.wave = ftab['wavelength']
            self.trans = ftab['transmission']
        
            
    def read_filter_old(self):
        wavescale = 1
        if hafilter == 'inthalpha':
            filterfile = self.filterpath+'/wfc-int-197-halpha.txt'
            # wavelength is in nm so scale by 10
            wavescale=10
            pass
        elif hafilter == 'intha6657':
            filterfile = self.filterpath+'/wfc-int-227-ha6657.txt'
            wavescale=10
            pass
        elif hafilter == 'sienaha':
            # wavelength in nm, transmission in percent
            filterfile = self.filterpath+'/chroma-halpha-transmission-ascii.txt'
            wavescale=10
            pass
        else:
            filterfile = self.filterpath+'/ha'+str(hafilter)+'-sim04.txt'
        tab = ascii.read(filterfile)
        self.wave = tab['col1']*wavescale # wavelength in angstrom
        self.trans = tab['col2'] # transmission percent
        if hafilter == 'sienaha':
            self.trans = 100*self.trans
            
    def get_filter_properties(self):
        self.maxtrans = np.max(self.trans)
        # get wavelengths where transmission crosses 10 percent level
        ids = np.where(self.trans > 10.)
        # calculate min and max redshifts that correspond to transmission cut
        self.minz_trans10 = (self.wave[ids[0][0]]/wave_halpha -1.)
        self.maxz_trans10 = (self.wave[ids[0][-1]]/wave_halpha -1.)

        # repeat for 30 and 50
        ids = np.where(self.trans > 30.)
        self.minz_trans30 = (self.wave[ids[0][0]]/wave_halpha -1.)
        self.maxz_trans30 = (self.wave[ids[0][-1]]/wave_halpha -1.)
        ids = np.where(self.trans > 50.)
        self.minz_trans50 = (self.wave[ids[0][0]]/wave_halpha -1.)
        self.maxz_trans50 = (self.wave[ids[0][-1]]/wave_halpha -1.)
        self.spline_fit()

        # calculate again to find where transmission is 10% of max transmission
        # get wavelengths where transmission crosses 10 percent level
        ids = np.where(self.trans/self.maxtrans > .1)
        # calculate min and max redshifts that correspond to transmission cut
        self.minz_trans10max = (self.wave[ids[0][0]]/wave_halpha -1.)
        self.maxz_trans10max = (self.wave[ids[0][-1]]/wave_halpha -1.)

        
    def spline_fit(self):
        
        # create spline fit
        self.spline_fit = interpolate.splrep(self.wave, self.trans)

    def get_transmission(self, wave):
        '''
        INTPUT:
        - wavelength, either individual value or array

        RETURNS:
        - transmission (spline fit to transmission curve at each wavelength
        - fitflag - True if wavelength is within the filter trace, false otherwise
        '''
        # make sure that the wavelength is in the right range
        if len(wave) > 1: #array
            #print('input is an array')
            transmission = np.zeros(len(wave),'f')
            self.fitflag = (wave > np.min(self.wave)) & (wave < np.max(self.wave))
            #print(self.fitflag)
            if sum(self.fitflag) == 0:
                print('all wavelengths out of range')
                return None, False
            elif sum(self.fitflag) < len(wave):
                print('WARNING: some galaxies are outside the filter window')
                print('wavelength out of range.  needs be between %.1f and %.1f Angstrom'%( np.min(self.wave),np.max(self.wave)))
            #print(self.fitflag)
            #print(wave)
            transmission[self.fitflag] = interpolate.splev(wave[self.fitflag],self.spline_fit)
            return transmission, self.fitflag
        elif (wave < np.min(self.wave)) | (wave > np.max(self.wave)): # check wavelength for a is single value
            print('wavelength out of range.  needs be between %.1f and %.1f Angstrom'%( np.min(self.wave),np.max(self.wave)))
            return None, False
        else:
            # function = return transmission for a given wavelength
            return interpolate.splev(wave,self.spline_fit), True
    def get_trans_correction(self, redshift,outfile=None):
        """
        calculate ratio of max filter transmission to transmission at obs wavelength of Halpha
        at the provided redhifts.

        This also generates a plot galaxies_in_filter.png of the filter transmission with histogram of redshifts overlaid.

        INPUT:
        redshift : float, can be a single value or an array

        RETURNS:
        correction : float, ratio of max filter transmission to transmission at that wavelength

        """
        wave = (redshift+1)*wave_halpha
        transmission, flag = self.get_transmission(wave)
        #self.test = transmission
        correction = np.zeros(len(transmission),'f')
        correction = self.maxtrans/transmission
        #plt.show()
        if outfile is not None:
            plt.figure()
            plt.plot(self.wave, self.trans/10,'k-')
            ##
            # set bin size to some constant range of min/max wavelength
            # this now fixes the odd behavior of the redshift histogram
            ##
            minwave = (self.minz_trans10max +1)*wave_halpha
            maxwave = (self.maxz_trans10max +1)*wave_halpha
            mybins = np.linspace(minwave,maxwave,20)


            plt.hist(wave, bins=mybins)

            # adding grid lines so we can better see how transmission varies with wavelength
            plt.grid(visible=True)
            plt.xlim((self.minz_trans10+1)*wave_halpha-50,(self.maxz_trans10+1)*wave_halpha+50)
            plt.xlabel('Wavelength (Angstrom)')
            plt.ylabel('Transmission %/10')
            titlestring = 'Halpha Filter = {}'.format(self.hafilter)
            plt.title(titlestring)
            
            plt.savefig(outfile)
            plt.close()
        return correction

    def get_response(self):
        """
        TODO : write this function
        integrate the filter to get the filter response
        
        integral of R(lambda) dlamba

        return integral
        """
        
        pass

