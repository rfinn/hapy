#!/usr/bin/env python
"""



USAGE:
python ~/github/halphagui/scripts/get_cutouts.py --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --virgo
"""
from astropy.io import fits
from astropy.table import Table
import numpy as np
from hapy.hatools import GalaxyCatalog, CoaddImage, HalphaImageSet, FilterTrace


def main(args=None):

    import argparse

    parser = argparse.ArgumentParser(description ='create psf image from image that contains stars')

    #parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    parser.add_argument('--rimage',dest = 'rimage', help='r-band image name')
    
    parser.add_argument('--psfdir',dest = 'psfdir', help='set to coadd directory')
    parser.add_argument('--catalog',dest = 'catalog', help='full path to galaxy catalog to use for cutouts.  ')    
    parser.add_argument('--agc',dest='agc', default=False, action='store_true', help='set this if using agc catalog')
    parser.add_argument('--virgo',dest='virgo', default=False, action='store_true', help='set this if using virgo catalog')
    parser.add_argument('--maxcorrection',dest='maxcorrection', default=3, help='maximum filter correction for galaxies in FOV.  default is 3, so galaxies whose redshift falls where filter transmission < 33 percent will be skipped.')        
    parser.add_argument('--oneimage',dest = 'oneimage',default=None, help='give full path to the r-band image name to run on just one image')
    
    args = parser.parse_args()

    
    try:
        rheader = fits.getheader(args.rimage)
        himage = rheader['HAIMAGE']
    except KeyError:
        print(f"WARNING: could not get HAIMAGE in rimage header {args.rimage}")
    image_set = HalphaImageSet(args.rimage, himage, psfdir=args.psfdir)
    image_set.load_coadds()


    ###################################################    
    # get galaxy catalog
    ###################################################    
    gcat = GalaxyCatalog(args.catalog,nsa=False,agc=args.agc,virgo=args.virgo,sizecat=None, verbose=False)

    
    gcat.galaxies_in_fov(image_set.h.wcs,zmin=None,zmax=None,image_name=himage, agcflag=args.agc,virgoflag=args.virgo)

    ###################################################
    # get redshift from filter trace module
    ###################################################    
    if args.virgo:
        redshift = gcat.cat['vr']/3.e5
    if args.agc:
        redshift = gcat.cat['vopt']/3.e5
        flag = gcat.cat['vopt'] == 0
        redshift[flag] = gcat.cat['v21'][flag]/3.e5
    redshift = redshift[gcat.keepflag] # redshifts for gals in FOV
    
    myfilter = FilterTrace(image_set.h.filter, instrument=image_set.h.instrument)
    corrections = myfilter.get_trans_correction(redshift,outfile=None)
    
    filter_keepflag = corrections < args.maxcorrection # this is a crazy big cut, but we can adjust with halphagui
    
    gcat.keepflag[gcat.keepflag] = filter_keepflag

    print(f"number of galaxies in FOV = {np.sum(gcat.keepflag)}")
    
    ###################################################
    # get galaxy cutouts
    ###################################################    
    gra = gcat.RA[gcat.keepflag]
    gdec = gcat.DEC[gcat.keepflag]
    gradius = gcat.radius_arcsec[gcat.keepflag]
    gBA = gcat.BA[gcat.keepflag]
    gPA = gcat.PA[gcat.keepflag]    
    #print("gra = ",gra)

    for i in range(len(gra)):
        print(f"ra={gra[i]:.6f}, dec={gdec[i]:.6f}, radius_arcsec={gradius[i]:6.2f}, BA={gBA[i]:.2f}, PA={gPA[i]:5.1f}")

if __name__ == '__main__':


    main()
