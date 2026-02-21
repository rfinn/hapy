#!/usr/bin/env python

'''
GOAL:
  The goal of this program is to run SExtractor in two image mode, to detect objects on the R-band image and
  apply the same apertures to the Halpha image.

PROCEDURE:
  - copy default setup files to current directory
  - Run sextractor on each image

EXAMPLE:

EXTRA NOTES:
- updated 8/14/19 to encase code in two functions.  this makes it easier to import into other programs, like the Halpha gui.

WRITTEN BY:
Rose Finn, 04 Jan 2017

'''
import os
from astropy.io import fits
import argparse
import numpy as np
from astropy.stats.sigma_clipping import sigma_clip

from astropy.stats import sigma_clip
from astropy.modeling import models, fitting


catdir = 'SEcats_filterratio'
if not os.path.exists(catdir):
    os.mkdir(catdir)



###############################################
### FUNCTIONS
###############################################

def run_sextractor(image1,image2, default_se_dir = '/Users/rfinn/github/halphagui/astromatic'):
    # get magnitude zeropoint for image 1 and image 2
    header1 = fits.getheader(image1)
    header2 = fits.getheader(image2)
    try:
        ZP1 = header1['PHOTZP']
        print('got ZP = ',ZP1)
        zp1flag = True
    except KeyError:
        print('no PHOTZP found in image 1 header.  Too bad :(')
        print('did you run getzp.py?')
        zp1flag = False
        ZP1 = None
    try:
        ZP2 = header2['PHOTZP']
        print('got ZP = ',ZP2)
        zp2flag = True
    except KeyError:
        print('no PHOTZP found in image 2 header.  Too bad :(')
        print('did you run getzp.py?')
        zp2flag = False
        ZP2 = None

    base = os.path.basename(image1)
    froot1 = os.path.splitext(base)[0]
    base = os.path.basename(image2)
    froot2 = os.path.splitext(base)[0]
    
    # check if output catalogs exist - if they do, don't rerun SE
    secatalog1 = f"{catdir}/{froot1}.cat"
    secatalog2 = f"{catdir}/{froot2}.cat"
    #print(secatalog1,secatalog2)
    if 'BOK' in image1:
        defaultcat = "default.sex.BOK"
    elif 'INT' in image1:
        defaultcat = "default.sex.INT"
    else:
        defaultcat = "default.sex.HDI"
    if os.path.exists(secatalog1) and os.path.exists(secatalog2):
        print("FOUND SE CATALOGS - NOT RERUNNING")
    else:
        print('RUNNING SOURCE EXTRACTOR')
        #s = f"sex {image1},{image2} -c default.sex.HDI -CATALOG_NAME {os.path.join(catdir,f'{froot2}.cat}'"
        if zp1flag: # why do I need to run image 1 in two image mode???
            #s ='sex ' + image1+','+image1 + ' -c default.sex.HDI -CATALOG_NAME ' + froot1 + '.cat -MAG_ZEROPOINT '+str(ZP1)
            s = f"sex {image1},{image1} -c {defaultcat} -CATALOG_NAME {catdir}/{froot1}.cat  -MAG_ZEROPOINT {ZP1}"
            os.system(s)
        else:
            s = f"sex {image1},{image1} -c {defaultcat} -CATALOG_NAME {catdir}/{froot1}.cat "            
            os.system(s)
        if zp2flag:
            s = f"sex {image1},{image2} -c {defaultcat} -CATALOG_NAME {catdir}/{froot2}.cat -MAG_ZEROPOINT {ZP2}"
            os.system(s)
        else:
            s = f"sex {image1},{image2} -c {defaultcat} -CATALOG_NAME {catdir}/{froot2}.cat "
            os.system(s)
    #print('in run_sextractor, returning for ZP and flags: ',ZP1, zp1flag, ZP2, zp2flag)
    return ZP1, zp1flag, ZP2, zp2flag

def make_plot(image1, image2, return_flag = False, plotdir = './', zps=None):
    from matplotlib import pyplot as plt
    import matplotlib
    matplotlib.use('Agg') # or 'Qt5Agg', 'GTKAgg', etc
    from scipy.stats import scoreatpercentile
    base = os.path.basename(image1)
    froot1 = os.path.splitext(base)[0]
    #cat1 = fits.getdata(froot1+'.cat')
    cat1name = os.path.join(catdir,froot1+'.cat')
    #print(froot1+'.cat')    
    # updating for draco - not sure what else this will break
    if not os.path.exists(cat1name):
        cat1name = os.path.join(catdir,image1.replace('.fits','.cat'))
    cat1 = fits.getdata(cat1name,2)

    base = os.path.basename(image2)
    froot2 = os.path.splitext(base)[0]
    cat2name = os.path.join(catdir,froot2+'.cat')
    if not os.path.exists(cat2name):
        cat2name = os.path.join(catdir,image1.replace('.fits','.cat'))
    
    try:
        cat2 = fits.getdata(cat2name,2)
    except IndexError:
        print("problem getting SE catalog for ",image2, " and catalog ",froot2," aborting make_plot")
        return None
    plt.figure(figsize=(6,6))
    plt.subplots_adjust(bottom=.2,left=.15,right=.95,hspace=.5)
    plt.subplot(2,1,1)

    # cut out extreme outliers using scoreatpercentile
    # May 2023 - changing apertures to APERTURE magnitudes
    # note - this will be a problem
    #
    # from getzp.py - I am using aperture mags here
    # y = self.matchedarray1['MAG_APER'][:,self.naper][flag]
    # yerr = self.matchedarray1['MAGERR_APER'][:,self.naper][flag]
    
    xmax = scoreatpercentile(cat1.FLUX_AUTO,95)
    ymax = scoreatpercentile(cat2.FLUX_AUTO,95)

    keepflag = (cat1.FLUX_AUTO < xmax) & (cat1.FLUX_AUTO > 0) & (cat2.FLUX_AUTO < ymax) & (cat2.FLUX_AUTO > 0)
    
    # TODONE - need to use SE flags to avoid things that are contaminated by nearby neighbors
    # also, use bright unsaturated sources to fit the filter ratio
    # also, compare measured ratio to diff in ZP
    keepflag = keepflag & (cat1.FLAGS < 1) & (cat2.FLAGS < 1)
    
    plt.plot(cat1.FLUX_AUTO[keepflag],cat2.FLUX_AUTO[keepflag],'k.',alpha=.4)

    # TODO - implement sigma clipping in fit, like in getzp.py
    c = np.polyfit(cat1.FLUX_AUTO[keepflag],cat2.FLUX_AUTO[keepflag],1,cov=True)
    print("results from polyfit = ",c)

    
    xline = np.linspace(0,xmax,100)
    plt.plot(xline,np.polyval(c[0],xline),ls='--')
    #plt.title("med ratio from np.polyfit (no clipping)")

    if zps is not None:
        print("adding zp line to filter ratio plot!")
        # we are plotting f2/f1 - ratio of Halpha to r
        ZP1,ZP2 = zps
        # get expected flux ratio from difference in ZP
        dm = ZP2-ZP1
        fratiozp = 10**(dm/2.5) # f2/f1
        #print('fratiozp = ',fratiozp)
        plt.plot(xline,fratiozp*xline,ls=':',c='r')
        plt.text(0.05,.8,'$ZP\ fratio = %.4f \ np.polyfit$'%(fratiozp),transform=plt.gca().transAxes,fontsize=8)
    
    print()
    #plt.xlim(0,xmax)
    #plt.ylim(0,ymax)
    plt.xlabel("r-band FLUX_AUTO")
    plt.ylabel("NB FLUX_AUTO")

    filename = os.path.basename(image2)
    t = filename.replace('.fits','')    
    plt.title(t,fontsize=12)
    plt.text(0.05,.9,'$med ratio = %.4f (%.4f)$'%(c[0][0],np.sqrt(c[1][0][0])),transform=plt.gca().transAxes,fontsize=8)    


    
    plt.subplot(2,1,2)
    x = cat2.FLUX_AUTO
    y = cat2.FLUX_AUTO/cat1.FLUX_AUTO
    flag = keepflag & (cat1.FLAGS == 0)
    plt.plot(x[flag],y[flag],'k.',alpha=.4)
    #plt.axis([0,500,.02,.08])
    x1,x2 = plt.xlim()
    #plt.xlim(.2,x2)
    #  TODO this is not working well for BOK
    #  need to do sigma clipping first?
    data = y[flag]
    clipped_data = sigma_clip(data,cenfunc='median',stdfunc='mad_std',sigma=2)
    ave = np.ma.median(clipped_data)
    # use the MAD instead
    std = np.ma.std(clipped_data)
    plt.axhline(y=ave,ls='--',label='SE flux ratios',lw=2)
    plt.ylim(-0.5*ave,3*ave)
    plt.title("med ratio with sigma clipping")
    #print('%.4f (%.4f)'%(ave,std))

    ##
    # Add line for ratio of zps 
    ##
    if zps is not None:
        plt.axhline(y = fratiozp,ls=':',c='r',label='ZP ratios')
        plt.text(0.05,.8,'$ZP\ fratio = %.4f$'%(fratiozp),transform=plt.gca().transAxes,fontsize=8)
    plt.ylabel('Flux(Halpha)/Flux(R)')
    plt.xlabel('Flux(Halpha) (ADU)')
    plt.text(0.05,.9,'$med ratio = %.4f (%.4f)$'%(ave,std),transform=plt.gca().transAxes,fontsize=8)
    plt.legend()
    #plt.gca().set_xscale('log')

    #plt.show()
    #plt.axis([0,5000,-.06,.06])
    plt.savefig(plotdir+'/'+t+'-filter-ratio.png')

    if return_flag:
        if zps is not None:
            return ave, std,fratiozp
        else:
            return ave, std
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ="Run sextractor in two-image mode.  \n To run from within ipython:\n %run ~/github/HalphaImaging/uat_sextractor_2image.py --image1 pointing-1_R.coadd.fits --image2 pointing-1_ha4.coadd.fits --plot --plotdir './' ")
    #parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Run sextractor to create object catalogs')
    parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of default config files')
    parser.add_argument('--image1',dest = 'image1', default = None,  help = 'image used to define apertures (R-band)')
    parser.add_argument('--image2',dest = 'image2', default = None,  help = 'image used to for measuring phot based on image1 (typically this is the Halpha image)')
    parser.add_argument('--plot',dest = 'plot', default = False, action = 'store_true', help = 'make diagnostic plots')
    parser.add_argument('--plotdir',dest = 'plotdir', default = '.', help = 'directory for saving plots')

    args = parser.parse_args()

    # get input files
    #print 'cp ' +args.d + '/default.* .'
    os.system('cp ' +args.d + '/default.* .')
    #files = sorted(glob.glob(args.filestring))

    #nfiles = len(files)
    i = 1
    ZP1,zp1flag,ZP2,zp2flag = run_sextractor(args.image1, args.image2, default_se_dir=args.d)
    print()
    print("in filterratio main, zpflag = ",zp1flag,zp2flag,(zp1flag and zp2flag))
    if zp1flag and zp2flag:
        #print("got ZP ratio")
        zpargs = (ZP1,ZP2)
    if args.plot:
        make_plot(args.image1, args.image2, plotdir = args.plotdir, zps = zpargs)

    
