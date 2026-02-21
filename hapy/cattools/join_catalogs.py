#!/usr/bin/env python

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import numpy as np
from astropy.table import Table
def join_cats(ra1,dec1,ra2,dec2,maxoffset=30.):
    '''
    INPUT:
    - enter two sets of coordinates: ra1, dec1 and ra2, dec2
       - these are assumed to be in deg
    - maxoffset = max separation to consider a match, in arcsec

    OUTPUT:
    - 4 arrays, with length equal to the number of unique galaxies

      - cat1_index = index of galaxy in catalog 1
      - cat1_flag = boolean array, true if galaxy is in catalog 1
      - cat2_index = index of galaxy in catalog 2
      - cat2_flag = boolean array, true if galaxy is in catalog 2

    
    '''
    maxoffset = maxoffset*u.arcsec
    c1 = SkyCoord(ra1*u.deg,dec1*u.deg, frame='icrs')
    c2 = SkyCoord(ra2*u.deg,dec2*u.deg, frame='icrs')

    # match catalog2 to catalog1
    id2, d2d, d3d = c1.match_to_catalog_sky(c2)
    # keep matches that have offsets < maxoffset
    # (match_to_catalog_sky finds closest object,
    # so some are not real matches
    matchflag2to1 = d2d < maxoffset

    # match catalog1 to catalog2
    id1, d2d, d3d = c2.match_to_catalog_sky(c1)
    # keep matches that have offsets < maxoffset
    # (match_to_catalog_sky finds closest object,
    # so some are not real matches
    matchflag1to2 = d2d < maxoffset

    # joined catalog will contain all objects in cat1
    # plus any in cat2 that weren't matched to cat1
    join_index = np.arange(len(ra1) + sum(~matchflag1to2))
    #print(len(join_index))
    cat1_index = np.zeros(len(join_index),'i')
    
    # flag for objects in cat 1 that are in the joined catalog
    # should be all of the objects in cat 1
    cat1_flag = join_index < len(ra1)
    cat1_index[cat1_flag] = np.arange(len(join_index))[cat1_flag]
    cat2_index = np.zeros(len(join_index),'i')
    cat2_index[cat1_flag] = id2
    
    #print(cat2_index[cat1_flag][matchflag2to1],len(cat2_index[cat1_flag][matchflag2to1]))
    #print(id2[matchflag2to1])
    cat2_flag = np.zeros(len(join_index),'bool')
    cat2_flag[cat1_flag] = matchflag2to1
    # this zeros the index of anything that doesn't have a match to cat1
    # otherwise index would be showing closest match to cat 1,
    # even if offset is > maxoffset
    cat2_index[~cat2_flag]=np.zeros(sum(~cat2_flag),'bool')


    # fill in details for objects in cat2 that weren't matched to cat1
    #print('number in cat2 that are not in cat1 = ',sum(~matchflag1to2))
    c2row = np.arange(len(ra2))
    cat2_index[~cat1_flag] = np.arange(len(ra2))[~matchflag1to2]
    cat2_flag[~cat1_flag] = np.ones(sum(~cat1_flag),'bool')
    return cat1_index,cat1_flag, cat2_index, cat2_flag

def make_new_cats(cat1, cat2):
    # returns tables with one line for each unique galaxy
    # tables contain contents of original two catalogs, but now line matched.
    # For example, table1 may have a row of zeros for a galaxy
    # that is not in tab1 but it is in tab2.
    cat1_index,cat1_flag, cat2_index, cat2_flag = join_cats(cat1['RA'],cat1['DEC'],cat2['RA'],cat2['DEC'])
    newcat1 = np.zeros(len(cat1_index),dtype=cat1.dtype)
    newcat1[cat1_flag] = cat1[cat1_index[cat1_flag]]
    newcat2 = np.zeros(len(cat1_index),dtype=cat2.dtype)
    newcat2[cat2_flag] = cat2[cat2_index[cat2_flag]]
    return Table(newcat1), cat1_flag, Table(newcat2), cat2_flag
if __name__ == '__main__':
    nsa = fits.getdata('NRGs27_nsa.fits')
    agc = fits.getdata('NRGs27_agc.fits')
    t = join_cats(nsa.RA,nsa.DEC,agc.RA,agc.DEC)
    j = make_new_cats(nsa,agc)
    
    
