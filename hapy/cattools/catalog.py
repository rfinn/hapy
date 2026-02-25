#!/usr/bin/env python
"""
catalog.py

Core catalog-handling utilities for the halphagui project.

This module provides functionality for:

- Loading galaxy catalogs (NSA, AGC, Virgo, etc.)
- Standardizing coordinate column names (RA/DEC)
- Identifying galaxies within an image field of view using WCS
- Applying redshift cuts (e.g., matching Hα filter transmission windows)
- Writing culled catalogs to disk

This module is backend-only and contains no GUI dependencies.

Intended Usage
--------------
The galaxy_catalog class is designed to be used by:

- Standalone processing scripts
- Continuum subtraction pipelines
- The halphagui interface

Example
-------
>>> from halphagui.core.catalog import galaxy_catalog
>>> cat = galaxy_catalog("nsa_catalog.fits", nsa=True)
>>> keep = cat.galaxies_in_fov(wcs, nrow=2048, ncol=4096,
...                            zmin=0.015, zmax=0.025)
>>> cat.cull_catalog(keep, prefix="field1")

Author
------
Rose Finn
2026 Feb 20

"""

from astropy.table import Table
import numpy as np
import os


from astropy.io import fits
#from astropy.wcs import WCS
#from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
#from astropy.coordinates import ICRS, FK5
#import astropy.units as u


import re

def _safe(s):
    """
    Make string safe for filenames.
    Keeps: letters, numbers, ., _, -
    Replaces everything else with _
    """
    s = str(s).strip().replace(" ", "")
    return re.sub(r"[^A-Za-z0-9._+\-]+", "_", s)

class GalaxyCatalog():
    """ 
    A container for galaxy catalog operations.

    Parameters
    ----------
    catalog : str
        Path to FITS table containing galaxy catalog.
    nsa : bool, optional
        Set True if catalog is an NSA catalog.
    agc : bool, optional
        Set True if catalog is an AGC catalog.
    virgo : bool, optional
        Set True if catalog is a Virgo catalog.
    sizecat : optional
        Optional associated size catalog.
    verbose : bool, optional
        Enable diagnostic output.

    Notes
    -----
    This class assumes that the catalog contains RA and DEC columns.
    If AGC-style column names (radeg/decdeg) are present, they will be
    renamed automatically.
    """
    
    agcflag: bool
    

    def __init__(self,catalog,nsa=False,agc=False,virgo=False,sizecat=None, verbose=False):
        self.cat = Table.read(catalog)
        
        #self.cat = Table(self.cat)
        self.catalog_name = catalog
        self.agcflag = agc
        self.nsaflag = nsa
        self.virgoflag = virgo
        
        if self.agcflag:
            scheme = "agc"
        elif self.virgoflag:
            scheme = "virgo"
        else:
            scheme = "generic"
        self.ensure_objid(scheme=scheme)
        
        self.verbose = verbose
        if self.agcflag:
            self.check_ra_colname()
            self.get_shape_agc()
        if self.virgoflag:
            self.get_shape_virgo()
        if sizecat is not None:
            self.sizecat = sizecat
        else:
            self.sizecat = None


    def check_ra_colname(self):
        """
        GOAL:
        make sure the catalog has RA and DEC
        columns that are named RA and DEC
        
        this is set up to rename the AGC fields radeg/decdeg to the more standard RA/DEC

        PARAMS:
        * self

        METHOD:
        * will edit the column names of self.cat (the galaxy catalog)

        """
        try:
            t = self.cat['RA']
        except AttributeError:
            print('defining new catalog columns for RA/DEC')
            self.cat.rename_column('radeg','RA')
            self.cat.rename_column('decdeg','DEC')            

        except KeyError:
            print('defining new catalogs columns for RA/DEC')
            #print(self.cat.colnames)
            self.cat.rename_column('radeg','RA')
            self.cat.rename_column('decdeg','DEC')            

    def ensure_objid(self, scheme: str = "generic"):
        """
        Ensure self.cat has an 'objid' column.
        - virgo: VFID-NEDname (or VFID alone if NEDname missing)
        - agc: AGCnr
        - generic: existing objid if present, else row index
        """
        if "objid" in self.cat.colnames:
            return

        n = len(self.cat)

        if scheme == "virgo":
            if "VFID" in self.cat.colnames and "NEDname" in self.cat.colnames:
                self.cat["objid"] = [f"{_safe(v)}-{_safe(name)}"
                                         for v, name in zip(self.cat["VFID"], self.cat["NEDname"])]
                #self.cat["objid"] = [f"{v}-{n}" for v, n in zip(self.cat["VFID"], self.cat["NEDname"])]
                return
            if "VFID" in self.cat.colnames:
                self.cat["objid"] = [str(v) for v in self.cat["VFID"]]
                return

        if scheme == "agc":
            if "AGCnr" in self.cat.colnames:
                self.cat["objid"] = [str(x) for x in self.cat["AGCnr"]]
                return

        # generic fallback
        self.cat["objid"] = [str(i) for i in range(n)]

    
    def galaxies_in_fov(self,wcs,zmin=None,zmax=None,image_name = None,weight_image=None, agcflag=None,virgoflag=None):
        """
        GOAL: get galaxies in FOV

        PROCEDURE:
        * transforms catalog coords to image coords using wcs

        PARAMS:
        * wcs - of image
        
        OPTIONAL PARAMS:
        * nrow
        * ncol
        * zmin - apply redshift cut to galaxies, e.g. that fall within halpha filter window
        * zmax - apply redshift cut to galaxies, e.g. that fall within halpha filter window
        
        """
        #print('in galaxies in fov, nrow,ncol = ',nrow,ncol) # debug
        #print(f"self.nsa flag is {self.nsa}")



        ###########################################################################
        # use astropy.WCS.wcs.footprint_contains to get galaxies w/in FOV of image
        ###########################################################################        
        # 
        # this can replace the method below, where I transform all the coordinates
        # to pixels.  However, astropy now returns nans for objects that are far
        # from the field center, and this is causing errors downstream.
        # So footprint_contains should be more robust.
        #
        coords = SkyCoord(ra=self.cat['RA'],dec=self.cat['DEC'],unit='deg') 
        self.keepflag = wcs.footprint_contains(coords)
        print(f"number of galaxies based on keepflag  = {np.sum(self.keepflag)}")       


 

        # check number of galaxies in fov
        if self.keepflag is None:
            print("WARNING: found no galaxies in FOV")
            return
        else:
            print(f"found {np.sum(self.keepflag)} after RA/DEC cuts")
            #print()


                     

        ###########################################################################
        # check weight image to make sure galaxy is in good part of image
        ###########################################################################
            
        # should also check the weight image and remove galaxies with weight=0
        # this won't take care of images with partial exposures, but we can deal with that later...
        # TODO - how to handle images with partial exposures, meaning only part of galaxy is in FOV?
        
        
        imagename = image_name
        if imagename is not None:
            if imagename.find('shifted.fits') > -1:
                weightimage = imagename.replace('-r-shifted.fits','-r.weight-shifted.fits')
            else:
                weightimage = imagename.replace('.fits','.weight.fits')

            #if os.path.exists(weightimage):
            if 'MOS' not in imagename: # TODO: not sure why I am skipping MOS.  should just check to see if weightimage exists?
                if os.path.exists(weightimage):
                    px,py = wcs.wcs_world2pix(self.cat['RA'][self.keepflag],self.cat['DEC'][self.keepflag],0)
                    print()
                    print("cross checking object locations with weight image")
                    print()
                    whdu = fits.open(weightimage)
                    # just check center position?
                    int_px = np.array(px,'i')
                    int_py = np.array(py,'i')        
                    centerpixvals = whdu[0].data[int_py,int_px]
                    # weight image will have value > 0 if there is data there
                    weightflag = centerpixvals > 0
                    self.keepflag[self.keepflag] = self.keepflag[self.keepflag] & weightflag

        if (zmin is not None) & (zmax is not None):
            self.keepflag = self.apply_redshift_cut(zmin=zmin,zmax=zmax, agcflag=agcflag,virgoflag=virgoflag)
        
        return self.keepflag
    
    def apply_redshift_cut(self,zmin=None,zmax=None,agcflag=False,virgoflag=False):
        ###########################################################################
        # get redshift cut
        ###########################################################################
        print(f"\nApplying redshift cut: zmin={zmin:.4f}, zmax={zmax:.4f}\n")
        # initialize value of zFlag
        zFlag = np.zeros(len(self.cat), 'bool')
        #print(f"DEBUGGING: len(self.cat)={len(self.cat)}, len(keepflag)={len(self.keepflag)}")
        if agcflag:
            print("\t using the AGC velocities")
            zFlag1 = (self.cat['vopt']/3.e5 > zmin) & (self.cat['vopt']/3.e5 < zmax)
            zFlag2 = (self.cat['v21']/3.e5 > zmin) & (self.cat['v21']/3.e5 < zmax)
            zFlag = zFlag1 | zFlag2
            return (zFlag & self.keepflag)
        else:
            try:
                if self.virgoflag:
                    #print('virgo, right?')
                    print("\t using the Virgo velocities")                
                    zFlag = (self.cat['vr']/3.e5 > zmin) & (self.cat['vr']/3.e5 < zmax)
                elif self.nsaflag:
                    print("\t using the NSA velocities")                
                    zFlag = (self.cat.Z > zmin) & (self.cat.Z < zmax)
                print('number of galaxies on image, after z cut = ',np.sum(zFlag & self.keepflag))
                return (zFlag & self.keepflag)

            except AttributeError:
                print('AttributeError')
                print('make sure you selected the halpha filter')
                return self.keepflag


    def cull_catalog(self, keepflag,prefix):
        self.cat = self.cat[keepflag]
        if self.nsaflag:
            self.rmag = 22.5 - 2.5*np.log10(self.cat.NMGY[:,4])
        
        if self.nsaflag:
            outfile = prefix+'_nsa.fits'
            fits.writeto(outfile,self.cat, overwrite=True)
        elif self.agcflag:
            outfile = prefix+'_agc.fits'
            try:
                fits.writeto(outfile,self.cat, overwrite=True)
            except:
                self.cat.write(outfile, overwrite=True)
        elif self.virgoflag:
            #print('virgo, right???')
            outfile = prefix+'_virgo_cat.fits'
            #print('culled catalog = ',outfile)
            self.cat.write(outfile,format='fits',overwrite=True)
            # cull ephot
    def get_shape_virgo(self):

        ephot = Table.read(self.catalog_name.replace('main.fits','legacy_ephot.fits'))

        
        bad_sb25 = ephot['SMA_SB25'] == 0

        self.radius_arcsec = ephot['SMA_SB25']*(~bad_sb25) + 1.35*ephot['SMA_SB24']*bad_sb25
        # OK, I know what you are thinking, I can't possibly be changing this again...

        # use SMA_SB25 instead of SB24 - this should work better for both high and low SB galaxies
        # if SMA_SB25 is not available use 1.35*SMA_SB24

        # for galaxies with SMA_SB24=0, set radius to value in main table 
        noradius_flag = self.radius_arcsec == 0
        self.radius_arcsec[noradius_flag] = self.cat['radius'][noradius_flag]

        # also save BA and PA from John's catalog
        # use the self.radius_arcsec for the sma
        self.BA = np.ones(len(self.radius_arcsec))
        self.PA = np.zeros(len(self.radius_arcsec))
        
        self.BA[~noradius_flag] = ephot['BA_MOMENT'][~noradius_flag]
        self.PA[~noradius_flag] = ephot['PA_MOMENT'][~noradius_flag]
        
        self.RA = self.cat['RA']
        self.DEC = self.cat['DEC']        
        
    def get_shape_agc(self):
        """
        use the latest AGC as the source catalog
        """


        self.radius_arcsec = self.cat['a']*60
        
        noradius_flag = self.radius_arcsec == 0
        self.radius_arcsec[noradius_flag] = 60 # set size of galaxies with no A value to 60 arcsec

        # also save BA and PA from John's catalog
        # use the self.radius_arcsec for the sma
        self.BA = np.ones(len(self.radius_arcsec))
        self.PA = np.zeros(len(self.radius_arcsec))
        
        self.BA[~noradius_flag] = self.agc.cat['b'][~noradius_flag]/self.agc.cat['a'][~noradius_flag]

        self.PA[~noradius_flag] = self.agc.cat['posang'][~noradius_flag]
        
        self.RA = self.cat['RA']
        self.DEC = self.cat['DEC']        
        
