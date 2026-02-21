#!/usr/bin/env python

import os
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.table import Table

wave_halpha = 6563. # angstrom
from . import utils



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
        instrument_to_prefix = {'INT':'WFC','BOK':'90prime','HDI':'HDI','MOS':'MOS'}

        filter_to_suffix = {4:'Ha+4nm',8:'Ha+8nm',12:'Ha+12nm',16:'Ha+16nm',\
                                '4':'Ha+4nm','8':'Ha+8nm','12':'Ha+12nm','16':'Ha+16nm',\
                                'ha4':'Ha+4nm','ha8':'Ha+8nm','ha12':'Ha+12nm','ha16':'Ha+16nm','ha':'Ha-197',\
                                'Ha+4nm':'Ha+4nm','Ha4nm':'Ha+4nm',\
                                'Halpha':'Ha-197','Ha6657':'Ha-227'}
        print("testing, self.hafilter = ",self.hafilter, self.instrument)
        self.halpha_filtername = utils.get_filter_file(f"{instrument_to_prefix[self.instrument]}-{filter_to_suffix[self.hafilter]}.fits")

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

        return correction

    def get_response(self):
        """
        TODO : write this function
        integrate the filter to get the filter response
        
        integral of R(lambda) dlamba

        return integral
        """
        
        pass

