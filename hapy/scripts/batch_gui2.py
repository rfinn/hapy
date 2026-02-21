#!/usr/bin/env python

'''
GOAL: 

Run this from, e.g. /home/rfinn/research/Virgo/gui-output-2019
- the gui will create a cutout folder in this directory that has a subdirectory for each pointing

cd github/halphagui

source venv/bin/activate

cd /data-pool/Halpha/halphagui-output-20230626

python ~/github/halphagui/batch_gui2.py --coaddir /media/rfinn/hdata/coadds/BOK2021pipeline/ --bok --psfdir /media/rfinn/hdata/psf-images/


NOTES:

2023-06-26: major rewrite for running on draco with new naming convention

'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib

import argparse

parser = argparse.ArgumentParser(description ='run the halpha gui in automatic mode.  Run this from the directory where you want the output data stored.  For example: /home/rfinn/research/Virgo/gui-output-NGC5846/')
parser.add_argument('--coadd_dir',dest = 'coadd_dir', default ='/data-pool/Halpha/coadds/all-virgo-coadds/', help = 'directory for coadds. Default is /home/rfinn/data/reduced/virgo-coadds-feb2019-int/')
parser.add_argument('--psfdir',dest = 'psfdir', default='/data-pool/Halpha/psf-images/', help = "directory containing PSF images.  Default is /data-pool/Halpha/psf-images/, which is for draco.")
parser.add_argument('--tabledir',dest = 'tabledir', default='/home/rfinn/research/Virgo/tables-north/v2/', help = "directory containing VF tables. Default is /home/rfinn/research/Virgo/tables-north/v2/")

args = parser.parse_args()



matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")

working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True


# get list of r-band coadded images
a = glob.glob(args.coadd_dir+'VF*INT*-r-shifted.fits')
b = glob.glob(args.coadd_dir+'VF*HDI*-r.fits')
c = glob.glob(args.coadd_dir+'VF*HDI*-R.fits')
d = glob.glob(args.coadd_dir+'VF*BOK*-r.fits')         
flist1 = a + b + c + d


flist1.sort()
##
# update for draco - stealing from build_web_coadds2.py
##
vtabledir = args.tabledir 
vmain = fits.getdata(vtabledir+'vf_v2_main.fits')
#homedir = '/mnt/qnap_home/rfinn/'
VFFIL_PATH = vtabledir+'/vf_v2_environment.fits'
vffil = fits.getdata(VFFIL_PATH)
#psfdir = homedir+'/data/reduced/psf-images/'


# get list of r-band coadded images
a = glob.glob(args.coadd_dir+'VF*INT*-r-shifted.fits')
b = glob.glob(args.coadd_dir+'VF*HDI*-r.fits')
c = glob.glob(args.coadd_dir+'VF*HDI*-R.fits')
d = glob.glob(args.coadd_dir+'VF*BOK*-r.fits')         
flist1 = a + b + c + d

flist1.sort()


i = 0
for rimage in flist1: # loop through list

    print()
    print('##########################################')        
    print('WORKING ON IMAGE: ',rimage)
    print('##########################################')        
    print()
    # read in r-band images
    # find matching Halpha image

    # grab other coadds

    rootname = rimage.split('-r')[0]

    # haimage is stored in header of r-band image
    header = fits.getheader(rimage)
    try:
        haimage = header['HAIMAGE']
        haimage = os.path.join(os.path.dirname(rimage),haimage)
    except KeyError:
        print("can't find halpha image for ",rimage)
        print("moving to next image...")
    
    ##
    # get filter based on the image name
    ##
    if 'BOK' in rimage:
        hfilter = 4
    elif 'HDI' in rimage:
        hfilter = 4
    elif 'INT' in rimage:
        if 'Ha6657' in rimage:
            hfilter = 'inthalpha6657'
        else:
            hfilter = 'inthalpha'
    else:
        print("could not guess filter for ",rimage)
        print("moving to the next image")

    ##
    # get prefix from image name
    ##
    t = rimage.split('-')
    if 'INT' in rimage:
        prefix = t[-3]
    else:
        prefix = t[-2]
    print(rimage, haimage,hfilter, prefix)
    
    command_string = f'python  ~/github/halphagui/halphamain.py --virgo --rimage {rimage} --haimage {haimage} --filter {hfilter} --psfdir {args.psfdir} --tabledir {args.tabledir} --prefix {prefix} --auto'

    # check to see if shifted r-band image exists.  if 
    try:
        print('running : ',command_string)
        # testing
        os.system(command_string)
    except:
        print('##########################################')
        print('WARNING: problem running auto gui on ',rimage)
        print('##########################################')

    #just running on one directory for testing purposes
    #i += 1
    #if i > 0:
    #    break


