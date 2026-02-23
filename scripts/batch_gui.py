#!/usr/bin/env python

'''
GOAL: 

Run this from, e.g. /home/rfinn/research/Virgo/gui-output-2019
- the gui will create a cutout folder in this directory that has a subdirectory for each pointing


python ~/github/halphagui/batch_gui.py --coaddir /media/rfinn/hdata/coadds/BOK2021pipeline/ --bok --psfdir /media/rfinn/hdata/psf-images/


'''

import os
import sys
import shutil
import glob
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(description ='run the halpha gui in automatic mode.  Run this from the directory where you want the output data stored.  For example: /home/rfinn/research/Virgo/gui-output-NGC5846/')
parser.add_argument('--coadd_dir',dest = 'coadd_dir', default ='/data-pool/Halpha/coadds/all-virgo-coadds/', help = 'directory for coadds. Default is /data-pool/Halpha/coadds/all-virgo-coadds/')
#parser.add_argument('--hdi',dest = 'hdi', default =False, action='store_true', help = 'set this for HDI data.  it will grab the filenames according to the HDI naming convention for the coadded images.')
#parser.add_argument('--mosaic',dest = 'mosaic', default =False, action='store_true', help = "set this for Mosaic data.  it will grab the filenames according to Becky's naming convention for the coadded images.")
#parser.add_argument('--bok',dest = 'bok', default =False, action='store_true', help = "set this for Bok 90prime data.  it will grab the filenames according to naming convention for the 90Prime images.")
parser.add_argument('--psfdir',dest = 'psfdir', default='/data-pool/Halpha/psf-images/', help = "directory containing PSF images.  Default is for draco /data-pool/Halpha/psf-images/.  When running on virgo vms, set to /mnt/qnap_home/rfinn/Halpha/reduced/psf-images/")
#parser.add_argument('--getgalsonly',dest = 'psfdir', default='/home/rfinn/data/reduced/psf-images/', help = "directory containing PSF images.  Default is /home/rfinn/data/reduced/psf-images/, which is for laptop.  When running on virgo vms, set to /mnt/qnap_home/rfinn/Halpha/reduced/psf-images/")
parser.add_argument('--oneimage',dest = 'oneimage', default=None, help = "use this to run on one image only.  Specify the full path to the r-band image. ")
parser.add_argument('--testing',dest = 'testing', default=False,action='store_true', help = "use this to run on the first image only for testing purposes. ")
parser.add_argument('--verbose',dest = 'verbose', action='store_true',default=False,help='set this for extra print statements')    
args = parser.parse_args()



matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")
#telescope = 'INT'


# get list of current directory
imagedir = args.coadd_dir

working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True


##
# update for draco - stealing from build_web_coadds2.py
##
coadd_dir = '/data-pool/Halpha/coadds/all-virgo-coadds/'                        
zpdir = coadd_dir+'/plots/'
fratiodir = coadd_dir+'/plots-filterratio/'
vtabledir = homedir+'/research/Virgo/tables-north/v2/'
vmain = fits.getdata(vtabledir+'vf_v2_main.fits')
#homedir = '/mnt/qnap_home/rfinn/'
VFFIL_PATH = vtabledir+'/vf_v2_environment.fits'
vffil = fits.getdata(VFFIL_PATH)
#psfdir = homedir+'/data/reduced/psf-images/'

outpathbase = '/data-pool/Halpha/'        
psfdir = outpathbase+'psf-images/'
outdir = outpathbase+'/html_dev/coadds/'

if args.oneimage is not None:
    # add code to run on the one image image
    # easiest is to cut the flist1 to include the image
    
    # make sure that the image exists
    print()
    print(args.oneimage)
    #print()
    print(str(os.path.exists(args.oneimage.rstrip())))
    if not os.path.exists(args.oneimage):
        teststring = '/data-pool/Halpha/coadds/all-virgo-coadds/VF-118.182+20.982-INT-20190205-p001-r-shifted.fits'
        print('test string = ',os.path.exists(teststring))
        print(f"Could not find {args.oneimage} - please check the r-band coadd name you provided")
        
        sys.exit()
        
    # redefine the list to include the one image only
    flist1 = [args.oneimage]

else:
    # updating for new naming convention and for the setup on draco
    # get list of r-band coadded images
    a = glob.glob(args.coadd_dir+'VF*INT*-r-shifted.fits')
    b = glob.glob(args.coadd_dir+'VF*HDI*-r.fits')
    c = glob.glob(args.coadd_dir+'VF*HDI*-R.fits')
    #######################################################
    # updating to use the shifted r-band images
    #######################################################    
    d = glob.glob(args.coadd_dir+'VF*BOK*-r-shifted.fits')
    e = glob.glob(args.coadd_dir+'VF*MOS*-R.fits')         
    flist1 = a + b + c + d + e

    #######################################################
    # for now, just want to rerun on the BOK images after
    # aligning the r to the halpha
    #######################################################
    
flist1.sort()
print("number of r-band images = ",len(flist1))


i = 0
for rimage in flist1: # loop through list

    print()
    print('##########################################')        
    print('WORKING ON IMAGE: ',rimage)
    print('##########################################')        
    print()
    # read in r-band images
    # find matching Halpha image    
    hdu = fits.open(rimage)
    rfilter = hdu[0].header['FILTER']
    try:
        haimage = hdu[0].header['HAIMAGE']
        haimage = os.path.join(os.path.dirname(rimage),haimage)
    except KeyError:
        print()
        print("WARNING: no halpha image in header of ",rimage)
        print("\t moving to next image")
        hdu.close()
        continue
    hdu.close()

    ##
    # get the halpha filter name
    ##
    if '-Halpha.fits' in haimage:
        hfilter = 'inthalpha'
    elif '-Ha6657.fits' in haimage:
        hfilter = 'intha6657'
    else: # includes BOK, HDI, MOSAIC
        hfilter = '4'

    ##
    # get rootname of the image
    ##
    rootname = rimage.split('-'+rfilter+'-')[0] # should be VF-2018-03-16-HDI-p054
    print(rootname)        
    rweightimage = rimage.replace('.fits','.weight.fits')
    pointing = rootname.split('-')[-1]
    dirname = os.path.dirname(rimage)
        
    prefix = os.path.basename(rootname).replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")

    if args.verbose:
        command_string = f'python  ~/github/halphagui/halphamain.py --virgo --rimage {rimage} --haimage {haimage} --filter {hfilter} --psfdir {args.psfdir} --tabledir /home/rfinn/research/Virgo/tables-north/v2/ --prefix {prefix} --verbose --auto'
    else:
        command_string = f'python  ~/github/halphagui/halphamain.py --virgo --rimage {rimage} --haimage {haimage} --filter {hfilter} --psfdir {args.psfdir} --tabledir /home/rfinn/research/Virgo/tables-north/v2/ --prefix {prefix} --auto'

    try:
        print('running : ',command_string)
        os.system(command_string)
    except:
        print('##########################################')
        print('WARNING: problem running auto gui on ',rimage)
        print('##########################################')

    #just running on one directory for testing purposes
    plt.close("all")
    i += 1
    if (i > 0) and args.testing:
        break


