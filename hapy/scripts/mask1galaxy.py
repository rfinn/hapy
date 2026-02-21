#!/usr/bin/env python

"""
GOAL:
* to create a mask from the r-band image of a galaxy
* to reproject that mask onto the WISE pixel scale
* this was developed for Virgo project and is using that directory structure

USAGE:

~/github/halphagui/mask1galaxy.py VFID_dirname


NOTES:

when running on vms, the topdir is:
/mnt/astrophysics/muchogalfit-output

when running on grawp, the topdir is:
/mnt/astrophysics/rfinn/muchogalfit-output

this switch is handled automatically


* 2024-02-08
- not sure how this will handle galaxies within JM's group images - I don't think it does

"""
import os
import sys
#import glob
import numpy as np
import subprocess

overwrite = True


print("trying to find virgo catalogs")
vdirs = ["/mnt/astrophysics/rfinn/catalogs/Virgo/v2/",\
         "/mnt/astrophysics/catalogs/Virgo/v2/",\
         "/home/rfinn/research/Virgo/tables-north/v2/"]

virgotabledir = None
for vd in vdirs:
    if os.path.exists(vd):
        virgotabledir = vd
        print("found ",virgotabledir)
        break

gdirs = ["/home/siena.edu/rfinn/github/halphagui/astromatic/",\
         "/home/rfinn/github/halphagui/astromatic/",\
         "/home/rfinn/github/halphagui/astromatic/"]


sepath = None
for gd in gdirs:
    if os.path.exists(gd):
        sepath = gd
        print("found ",sepath)
        break

if virgotabledir is None:
    print("found the machine name but could not find table dir")
    sys.exit()
            
def funpack_image(input,output,nhdu=1):
    command = 'funpack -O {} {}'.format(output,input)
    print(command)
    os.system(command)

def get_galaxy_params(VFID):
    ##
    # get sizes for galaxies - will use this to unmask central region
    # need to cut this catalog based on keepflag
    ##
    from astropy.table import Table
    
    ephot = Table.read(virgotabledir+'/vf_v2_legacy_ephot.fits')
    vmain = Table.read(virgotabledir+'/vf_v2_main.fits')


    # get galaxy id

    galid = np.arange(len(vmain))[vmain['VFID']==VFID][0]
    #self.radius_arcsec = ephot['SMA_SB24']
    
    bad_sb25 = ephot['SMA_SB25'] == 0
    
    radius_arcsec = ephot['SMA_SB25']*(~bad_sb25) + 1.35*ephot['SMA_SB24']*bad_sb25
    # OK, I know what you are thinking, I can't possibly be changing this again...
    
    # use SMA_SB25 instead of SB24 - this should work better for both high and low SB galaxies
    # if SMA_SB25 is not available use 1.35*SMA_SB24
    
    # for galaxies with SMA_SB24=0, set radius to value in main table 
    noradius_flag = radius_arcsec == 0
    radius_arcsec[noradius_flag] = vmain['radius'][noradius_flag]
    
    # also save BA and PA from John's catalog
    # use the self.radius_arcsec for the sma
    BA = np.ones(len(radius_arcsec))
    PA = np.zeros(len(radius_arcsec))
    
    BA[~noradius_flag] = ephot['BA_MOMENT'][~noradius_flag]
    PA[~noradius_flag] = ephot['PA_MOMENT'][~noradius_flag]

    gRAD = radius_arcsec[galid]

    gBA = BA[galid]
    gPA = PA[galid]    
    gRA = vmain['RA'][galid]
    gDEC = vmain['DEC'][galid]
    return gRA,gDEC,gRAD,gBA,gPA

                        
    
if __name__ == '__main__':
    homedir = os.getenv("HOME")
    topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
    image_source_dir = '/mnt/astrophysics/virgofilaments-data/'
    try:
        os.chdir(topdir)
    except FileNotFoundError: # assuming that we are running on virgo vms
        
        topdir = '/mnt/astrophysics/muchogalfit-output/'
        image_source_dir = '/mnt/virgofilaments-data/'
        try:
            os.chdir(topdir)
        except FileNotFoundError:
            topdir = os.getcwd()+'/'
            os.chdir(topdir)
    # take as input the galaxy name
    galname = sys.argv[1]
    print("working on ",galname)
    # move to muchogalfit-output directory
    output_dir = topdir+galname+'/'
    print(output_dir)
    if not os.path.exists(output_dir):
        print('WARNING: {} does not exist\n Be sure to run setup_galfit.py first')
        os.chdir(topdir)
        sys.exit()
    
    os.chdir(output_dir)

    # read in input file to get galname, objname, ra, dec, and bandpass
    sourcefile = open(galname+'sourcelist','r')
    galaxies = sourcefile.readlines()
    if len(galaxies) > 1:
        # set the flag to have more than one galaxy in the galfit input file
        multiflag = True
    elif len(galaxies) == 1:
        multiflag = False
    else:
        print('Problem reading sourcelist for {}'.format(galname))
        print('Please check the setup directory')
        os.chdir(topdir)
        sys.exit()

    # parse information from file
    vfid, objname, ra, dec, bandpass = galaxies[0].rstrip().split()
    ra = float(ra)
    dec = float(dec)

    # set up path name for image directory
    # directory where galaxy images are
    data_dir = f"{image_source_dir}/{int(ra)}/{objname}/"





    # construct image name
    image = f"{objname}-custom-image-r.fits"
    # check if r-band image exists
    if not os.path.exists(image):
        # if not, copy from john's directories and funpack
        funpack_image(os.path.join(data_dir,image),os.path.join(output_dir,image.replace('.fz','')))
        # get information from the image, like the image size and position of gal
        image = image.replace('.fz','')
        
    # define the mask name created by maskwrapper.py
    mask = image.replace('.fits','-mask.fits')

    # check if the mask exists
    if os.path.exists(mask):
        if overwrite:
            pass
            print("remaking mask")
        else:
            print(f"mask already exists for {vfid}.")
            print("moving to the next galaxy")
            sys.exit()

    # get shape parameters for galaxy
    if virgotabledir is not None:
        gRA,gDEC,gRAD,gBA,gPA = get_galaxy_params(vfid)
        gPA = gPA + 90
        cmd = f"python {homedir}/github/halphagui/maskwrapper.py --image {image} --objra {ra} --objdec {dec} --objsma {gRAD:.1f} --objBA {gBA:.1f} --objPA {gPA:.1f} --sepath {sepath} --auto"
    else:
        cmd = f"python {homedir}/github/halphagui/maskwrapper.py --image {image} --auto"
    # call maskwrapper.py
    print(cmd)

    # plot mask and central ellipse
    
    #plt.figure()
    #data = fits.getdata(mask)
    #plt.imshow(data)

    # plot ellipse
    os.system(cmd)




    # reproject mask onto wise image
    
    # construct WISE W3 image name
    wimage = f"{objname}-custom-image-W3.fits"
    # check if r-band image exists
    if not os.path.exists(wimage):
        # if not, copy from john's directories and funpack
        funpack_image(os.path.join(data_dir,wimage),os.path.join(output_dir,wimage.replace('.fz','')))
        # get information from the image, like the image size and position of gal
        wimage = wimage.replace('.fz','')
        
    # call reproject_mask.py
    cmd = f"python {homedir}/github/halphagui/reproject_mask.py {mask} {wimage}"
    #print(cmd)
    print()
    os.system(cmd)
    print(f"done masking {vfid}")
    print()
