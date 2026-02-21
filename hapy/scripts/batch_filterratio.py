#!/usr/bin/env python

'''
GOAL: 
- get filter ratio for all images in the current directory

- need to account for the different naming conventions of the different Halpha filters


USAGE:

Run from the directory with the coadds /media/rfinn/hdata/coadds/all-virgo-coadds

python ~/github/halphagui/batch_filterratio.py

'''

import os
import sys
import glob
from astropy.io import fits
import matplotlib
import time
import filterratio as runse
import multiprocessing as mp

try:
    matplotlib.use("Qt5agg")
except ImportError:
    print("WARNING: could not use Qt5agg backend")
homedir = os.getenv("HOME")

# overwrite output files if they exist
overwrite = False

import argparse

image_results = []
def collect_results(result):

    global results
    image_results.append(result)


def subtract_images(rimage,himage,filter_ratio,zpflag=False):
    r = fits.open(rimage)
    ha = fits.open(himage)
    halpha_cs = ha[0].data - filter_ratio*r[0].data
    hacoadd_cs_fname = himage.split('.fits')[0]+'-CS.fits'
    if zpflag:
        hacoadd_cs_fname = himage.split('.fits')[0]+'-CS-ZP.fits'
    print("subtracting images: ",himage," -> ",hacoadd_cs_fname)
    # ratio is already written to r-band image header
    fits.writeto(hacoadd_cs_fname,halpha_cs,header=ha[0].header,overwrite=True)
    r.close()
    ha.close()
    pass




def getoneratio(rimage,instrument,plotdir):
    # find the corresponding Halpha image
    f = instrument
    if f == 'INT':
        # could end in Halpha or Ha6657
        t = rimage.split('-')
        #print(t)
        if rimage[10] == '+':
            dateobs = t[3]
            pointing = t[4]
        else:
            dateobs = t[4]
            pointing = t[5]
        if dateobs == '20190531': # ugh - new month too!
            newdate = '20190601'
        else:
            newdate = dateobs[:-2]
        hastring = 'VF-*INT-'+newdate+'*-'+pointing+'*-Halpha.fits'
        hfiles = glob.glob(hastring)
        if len(hfiles) < 1:
            hastring = 'VF-*INT-'+newdate+'*-'+pointing+'*-Ha6657.fits'
            hfiles = glob.glob(hastring)

        #print("found these ",hfiles)
        if len(hfiles) > 0:
            himage = hfiles[0]
        else:
            print()
            print('Warning: no halpha image found for ',rimage)
            print('moving to the next image...')
            return
    elif f == 'BOK':
        # should end in Ha4.fits
        testname = rimage.replace('-r.fits','-Ha4.fits')
        testname = rimage.replace('-r-shifted.fits','-Ha4.fits')
        if os.path.exists(testname):
            himage = testname
        # check for images that were taken on different dates
        elif rimage == 'VF-266.477+58.350-BOK-20220423-VFID0783-r-shifted.fits':
            himage = 'VF-266.477+58.350-BOK-20220428-VFID0783-Ha4.fits'
        #elif rimage == 'VF-246.890+21.523-BOK-20220426-VFID3598-r.fits':
        #    himage = 'VF-246.890+21.523-BOK-20220425-VFID3598-Ha4.fits'
        else:
            print()
            print('Warning: no halpha image found for ',rimage)
            print("looking for ",testname)
            print('moving to the next image...')
            return
    elif f == 'HDI':
        # should end in ha4.fits
        testname = rimage.replace('-r.fits','-ha4.fits').replace('-R.fits','-ha4.fits')

        print('testname = ',testname)
        if os.path.exists(testname):
            himage = testname

        else:
            # halpha image could have a different date
            t = rimage.split('-')

            if rimage[10] == '+':
                dateobs = t[3]
                pointing = t[4]
            else:
                dateobs = t[4]
                pointing = t[5]
            hastring = 'VF-*'+f+'-*'+dateobs[:-2]+'*-'+pointing+'*-ha4.fits'
            hfiles = glob.glob(hastring)

            #print("found these ",hfiles)
            if len(hfiles) > 0:
                himage = hfiles[0]
                print('halpha image = ',himage)                
            else:
                print()
                print('Warning: no halpha image found for ',rimage)
                print('moving to the next image...')
                return
    elif f == 'MOS':
        testname = rimage.replace('-R.fits','-Ha4.fits')
        if os.path.exists(testname):
            himage = testname
        else:
            print()
            print('Warning: no halpha image found for ',rimage)
            print('moving to the next image...')
            return

        

    print()
    print('##########################################')        
    print('GETTING FILTER RATIO FOR: ',rimage,himage)
    print('##########################################')

    start_time = time.perf_counter()

    # TODO - we should use the himage as a reference b/c it's lower snr
    # we don't want extra noise
    #ZP1,zp1flag,ZP2,zp2flag = runse.run_sextractor(rimage, himage)
    # this line might do the job, but should check to make sure I'm not missing something
    #ZP2,zp2flag,ZP1,zp1flag = runse.run_sextractor(himage, rimage)
    ZP1,zp1flag,ZP2,zp2flag = runse.run_sextractor(rimage, himage)
    
    if zp1flag and zp2flag:
        #print("got ZP ratio")
        zpargs = (ZP1,ZP2)
    else:
        zpargs = None
    t = runse.make_plot(rimage, himage, return_flag = True, plotdir = plotdir,zps = zpargs)
    if t == -1:
        print("no CS images make for ",rimage,himage)
        return

    if len(t) == 2:
        ave, std = t
        fzpratio = None
    elif len(t) == 3:
        ave, std, fzpratio = t
    #print(ave,std)

    subtract_images(rimage,himage,ave)

    if fzpratio is not None:
        subtract_images(rimage,himage,fzpratio,zpflag=True)

    # add ratio to r-band image headers
    r,header = fits.getdata(rimage,header=True)
    header.set('FLTRATIO',ave)
    header.set('FLTR_ERR',std)
    header.set('HAIMAGE',os.path.basename(himage))
    if fzpratio is not None:
        header.set('FRATIOZP',fzpratio)
        
    fits.writeto(rimage,r,header=header,overwrite=True)        
    # clock time to get filter ratio
    end_time = time.perf_counter()
    print('\t total time = ',end_time - start_time)

def getoneratio_vfs(rimage,plotdir):
    # find the corresponding Halpha image
    rheader = fits.getheader(rimage)
    himage = rheader['HAIMAGE']

    print()
    print('##########################################')        
    print('GETTING FILTER RATIO FOR: ',rimage,himage)
    print('##########################################')

    start_time = time.perf_counter()

    # TODO - we should use the himage as a reference b/c it's lower snr
    # we don't want extra noise
    #ZP1,zp1flag,ZP2,zp2flag = runse.run_sextractor(rimage, himage)
    # this line might do the job, but should check to make sure I'm not missing something
    #ZP2,zp2flag,ZP1,zp1flag = runse.run_sextractor(himage, rimage)
    ZP1,zp1flag,ZP2,zp2flag = runse.run_sextractor(rimage, himage)
    
    if zp1flag and zp2flag:
        #print("got ZP ratio")
        zpargs = (ZP1,ZP2)
    else:
        zpargs = None
    t = runse.make_plot(rimage, himage, return_flag = True, plotdir = plotdir,zps = zpargs)
    if t == -1:
        print("no CS images make for ",rimage,himage)
        return

    if len(t) == 2:
        ave, std = t
        fzpratio = None
    elif len(t) == 3:
        ave, std, fzpratio = t
    #print(ave,std)

    subtract_images(rimage,himage,ave)

    if fzpratio is not None:
        subtract_images(rimage,himage,fzpratio,zpflag=True)

    # add ratio to r-band image headers
    r,header = fits.getdata(rimage,header=True)
    header.set('FLTRATIO',ave)
    header.set('FLTR_ERR',std)
    #header.set('HAIMAGE',os.path.basename(himage))
    if fzpratio is not None:
        header.set('FRATIOZP',fzpratio)
        
    fits.writeto(rimage,r,header=header,overwrite=True)        
    # clock time to get filter ratio
    end_time = time.perf_counter()
    print('\t total time = ',end_time - start_time)

def getoneratio_uat(rimage,plotdir):
    # find the corresponding Halpha image

    t = rimage.split('-')
    if len(t) == 8:
        instrument = t[3]
        dateobs = t[4]
        pointing = t[6]+'-'+t[7]
    elif len(t) == 7: # positive declinations
        instrument = t[2]
        dateobs = t[3]
        pointing = t[4]+'-'+t[5]

    # get halpha image from r-band header
    rheader = fits.getheader(rimage)
    himage = rheader['HAIMAGE']
    print()
    print('##########################################')        
    print('GETTING FILTER RATIO FOR: ',rimage,himage)
    print('##########################################')

    start_time = time.perf_counter()

    # TODO - we should use the himage as a reference b/c it's lower snr
    # we don't want extra noise
    #ZP1,zp1flag,ZP2,zp2flag = runse.run_sextractor(rimage, himage)
    # this line might do the job, but should check to make sure I'm not missing something
    #ZP2,zp2flag,ZP1,zp1flag = runse.run_sextractor(himage, rimage)
    ZP1,zp1flag,ZP2,zp2flag = runse.run_sextractor(rimage, himage)
    
    if zp1flag and zp2flag:
        #print("got ZP ratio")
        zpargs = (ZP1,ZP2)
    else:
        zpargs = None
    t = runse.make_plot(rimage, himage, return_flag = True, plotdir = plotdir,zps = zpargs)
    if t == -1:
        print("no CS images make for ",rimage,himage)
        return

    if len(t) == 2:
        ave, std = t
        fzpratio = None
    elif len(t) == 3:
        ave, std, fzpratio = t
    #print(ave,std)

    subtract_images(rimage,himage,ave)

    if fzpratio is not None:
        subtract_images(rimage,himage,fzpratio,zpflag=True)

    # add ratio to r-band image headers
    r,header = fits.getdata(rimage,header=True)
    header.set('FLTRATIO',ave)
    header.set('FLTR_ERR',std)
    #header.set('HAIMAGE',os.path.basename(himage))
    if fzpratio is not None:
        header.set('FRATIOZP',fzpratio)
        
    fits.writeto(rimage,r,header=header,overwrite=True)        
    # clock time to get filter ratio
    end_time = time.perf_counter()
    print('\t total time = ',end_time - start_time)
    

    
########################################################
# loop through telescopes and get r-band images
########################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Run filterratio.py for all images in the specified directory')
    parser.add_argument('--plotdir',dest = 'plotdir', default ='plots-filterratio', help = 'directory for to store the plots in')
    parser.add_argument('--prefix', dest = 'prefix', default = 'VF-',help = "prefix for coadds.  the default is VF-")
    parser.add_argument('--uat', dest = 'uat', default = False, action='store_true',help = "set this for running with UAT data.")
    #parser.add_argument('--vfs', dest = 'uat', default = False, action='store_true',help = "set this for running with UAT data.")        
    parser.add_argument('--oneimage',dest = 'oneimage',default=None, help='give full path to the r-band image name to run on just one image')    
    args = parser.parse_args()


    plotdir = os.path.join(os.getcwd(),args.plotdir)
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)

    # INT
    # r-shifted.fits
    # Halpha or Ha6657

    # BOK
    # r-shifted.fits
    # Ha4

    # HDI
    # -r.fits -R.fits
    # ha4

    # will need to add more for mosaic images

    # telescope/instrument names

    # why no MOS here?

    #inames = ["BOK","HDI"]
    #inames = ["HDI"]
    #inames = ["BOK"]
    #inames = ["INT","HDI"]
        
    if args.oneimage is not None:
        if not os.path.exists(args.oneimage):
            print(f"Could not find {args.oneimage} - please check the r-band coadd name you provided")
            sys.exit()
        
        if args.uat:
            getoneratio_uat(args.oneimage,plotdir)
            
        else:
 
            getoneratio_vfs(args.oneimage,plotdir)

    else:
        for i,f in enumerate(inames):
            # get list of current directory
            # this will grab the coadds but not the weight images
            if f == 'INT':
                #print("match string = ",'VF-*'+f+'*r-shifted.fits')
                rimages = glob.glob('VF-*'+f+'*r-shifted.fits')
            elif f == 'BOK':
                rimages = glob.glob('VF-*'+f+'*r-shifted.fits')
            elif f == 'HDI':
                rimages1 = glob.glob('VF-*'+f+'*r.fits')
                rimages2 = glob.glob('VF-*'+f+'*R.fits')
                rimages = rimages1 + rimages2
            elif f == 'MOS':
                rimages = glob.glob('VF-*'+f+'*R.fits')

            print(f"found {len(rimages)} rband images for {f}")

            # sort the file list
            rimages.sort()

        
            #for rimage in rimages: # loop through list
            image_pool = mp.Pool(mp.cpu_count())
            myresults = [image_pool.apply_async(getoneratio,args=(im,f,plotdir),callback=collect_results) for im in rimages]
    
            image_pool.close()
            image_pool.join()
            image_results = [r.get() for r in myresults]




