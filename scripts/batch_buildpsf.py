#!/usr/bin/env python

'''
GOAL: 
- create psf images for all the coadds

Run this from, e.g. /home/rfinn/data/reduced/psf-images/

Change the coadd_image_directory as needed before running


'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib
import time
import multiprocessing as mp


matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")
telescope = 'INT'
working_dir = os.getcwd()

# overwrite output files if they exist
overwrite = False

import argparse


image_results = []
def collect_results(result):

    global results
    image_results.append(result)


def runone(image,int=False,bok=False):
    #print("hello ",image)
    basename = os.path.basename(image).split('.fits')[0]
    psf_image_name = basename+'-psf.fits'
    print()
    print("PSF image name = ",psf_image_name)
    print()
    if os.path.exists(psf_image_name):
        print("found psf image for ",image)
        return

    print('##########################################')
    print('##########################################')        
    print('BUILDING PSF FOR: ',image)
    print('##########################################')
    print('##########################################')
    
    if int:
        command_string = 'python ~/github/halphagui/buildpsf.py --image {} --int'.format(image)
    elif bok:
        command_string = 'python ~/github/halphagui/buildpsf.py --image {} --bok'.format(image)
    else:
        command_string = 'python ~/github/halphagui/buildpsf.py --image {} '.format(image)
    try:
        print('running : ',command_string)
        os.system(command_string)
    except:
        print('##########################################')
        print('WARNING: problem running buildpsf.py for ',image)
        print('##########################################')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Run buildpsf.py for all images in the specified directory')
    parser.add_argument('--coaddir',dest = 'coaddir', default ='/home/rfinn/data/reduced/virgo-coadds-feb2019-int/', help = 'directory for coadds. Default is /home/rfinn/data/reduced/virgo-coadds/feb2019-int/')
    parser.add_argument('--int', dest = 'int', default = False,action='store_true', help = 'set this for INT data')
    parser.add_argument('--bok', dest = 'bok', default = False,action='store_true', help = 'set this for BOK data')
    parser.add_argument('--hdi', dest = 'hdi', default = False,action='store_true', help = 'set this for HDI data')
    parser.add_argument('--ngc', dest = 'ngc', default = False,action='store_true', help = "set this for Becky's NGC5846 data")
    parser.add_argument('--uat', dest = 'uat', default = False,action='store_true', help = "set this for UAT Halpha groups data")

    args = parser.parse_args()

    coadd_image_directory = args.coaddir


    filters = ['r','Halpha','Ha6657','ha4','R','Ha','Ha+4nm','Ha4','ha8','ha12','ha16']
    #saturate_level = [100,30,30]

    for i,f in enumerate(filters):
        # get list of current directory
        # this will grab the coadds but not the weight images
        if args.bok:
            print("running bok")
            flist1 = glob.glob(coadd_image_directory+'VF-*BOK*-'+f+'.fits') 
           
        elif args.int or args.hdi:
            print("running int/hdi")                       
            flist1 = glob.glob(coadd_image_directory+'VF-*-'+f+'.fits')
        elif args.ngc: 
            print("running ngc")                      
            flist1 = glob.glob(coadd_image_directory+'nNGC*'+f+'.fits')
        elif args.uat: 
            print("running uat")           
            flist1 = glob.glob(coadd_image_directory+'UAT*'+f+'.fits')
        else:
            print("running default search")
            searchstring = coadd_image_directory+'VF-*'+f+'.fits'
            print('looking for files: ',searchstring)
            flist1 = glob.glob(searchstring)

        print(i,f,len(flist1))
        # changing this to just do pointing 022 and 026
        #flistp20 = glob.glob(coadd_image_directory+'VF-*p017*-'+f+'.fits')
        #flistp26 = glob.glob(coadd_image_directory+'VF-*p072*-'+f+'.fits')
        #flist1 = flistp20+flistp26

        # checking commit
        flist1.sort()

        
        #print(flist1[0:3])


        #for rimage in rimages: # loop through list
        #image_pool = mp.Pool(mp.cpu_count())
        #myresults = [image_pool.apply_async(runone,args=(im,args.int,args.bok),callback=collect_results) for im in flist1]

        #image_pool.close()
        #image_pool.join()
        #image_results = [r.get() for r in myresults]

        # going back to serial because it is failing on some when using mp but not when running individually
        for fimage in flist1: # loop through list
            if 'weight.fits' in fimage:
                continue
            if 'INT' in fimage:
                runone(fimage, int=True,bok=False)
            elif 'BOK' in fimage:
                runone(fimage, int=False,bok=True)
            else:
                runone(fimage, int=False,bok=False)
            #start_time = time.perf_counter()
            # adding saturation limit for normalized images
            # I estimated this from the r-band image for p001
            # this is in counts/sec
            #if not overwrite:
                # check if psf image exists

            # clock time to run buildpsf
            #end_time = time.perf_counter()
            #print('\t total time = ',end_time - start_time)

            # just running on one directory for testing purposes
            #break
        #break

        #'''

