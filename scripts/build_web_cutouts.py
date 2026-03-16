#!/usr/bin/env python

'''
GOAL:
* create web page to inspect the cutouts

USAGE:
* run from cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)
* https://docs.astropy.org/en/stable/visualization/normalization.html#:~:text=The%20astropy.visualization%20module%20provides%20a%20framework%20for%20transforming,x%20represents%20the%20values%20in%20the%20original%20image%3A

TO DO: (I think these are all fixed!)
* DONE need to figure out how to handle repeated observations
  - don't overwrite directory

* fix how I combine unwise images when multiple images are returned
* same for galex

2024-03-10 : 
* add panel to show color-based continuum subtract (added this to Halpha panel, instead of showing CS with two stretches)
* stellar mass, sfr and ssfr images

'''

import os
import sys
import numpy as np
import glob

from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

from scipy.stats import scoreatpercentile
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch, PercentileInterval
from astropy.visualization import ImageNormalize
from astropy.visualization import LinearStretch,SinhStretch
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip

import warnings

import multiprocessing as mp

from hapy.imagetools.plotting import display_image
from hapy.geometry.adapters import pa_ccw_north_to_photutils_theta
from PIL import Image

homedir = os.getenv("HOME")

os.sys.path.append(homedir+'/github/virgowise/')
import rungalfit as rg #This code has galfit defined functions 

#from build_web_coadds import get_galaxies_fov, plot_vf_gals
from build_web_common import *
###########################################################
####  GLOBAL VARIABLES
###########################################################


VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'
VFEPHOT_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits'
haimaging_path = os.path.join(homedir,'github/HalphaImaging/python3/')
#sys.path.append(haimaging_path)

#vfmain = fits.getdata(VFMAIN_PATH)
residual_stretch = LinearStretch(slope=0.5, intercept=0.5) + SinhStretch() + \
    LinearStretch(slope=2, intercept=-1)
###########################################################
####  FUNCTIONS
###########################################################
def get_params_from_name(image_name):
    #print(t)
    tels = ['BOK','HDI','INT','MOS']
    for t in tels:
        if t in image_name:
            telescope = t
            break
    t = os.path.basename(image_name).split('-')
    for item in t:
        if item.startswith('20'):
            dateobs = item
            break
    pointing = t[-1]

    return telescope,dateobs,pointing

def buildone(subdir,outdir,flist):
    print(subdir)

    telescope,dateobs,pointing = get_params_from_name(subdir)
    run = dateobs+'-'+pointing
    #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
    if os.path.isdir(subdir):
        print('##########################################')
        print('##########################################')        
        print('WORKING ON DIRECTORY: ',subdir)
        print('##########################################')
        print('##########################################')
        
        #try:
        # move to subdirectory
        # adding the telescope and run so that we don't write over
        # images if galaxy was observed more than once
        gal_outdir = os.path.join(outdir,subdir+"")
        #print('out directory for this galaxy = ',gal_outdir)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(gal_outdir):
            os.mkdir(gal_outdir)

        p = cutout_dir(cutoutdir=subdir,outdir=gal_outdir)
        p.runall()

        i = flist.index(args.oneimage)            
        # define previous gal for html links
        if i > 0:
            previous = (flist[i-1])
            #print('previous = ',previous)
        else:
            previous = None
        if i < len(flist)-1:
            next = flist[i+1]
            #print('next = ',next)
        else:
            next = None
        h = build_html_cutout(p,gal_outdir,previous=previous,next=next,tel=telescope,run=run)
        h.build_html()
        #except:
        #    print('WARNING: problem building webpage for ',subdir)
    

# def display_image(image,percentile1=.5,percentile2=99.5,stretch='asinh',mask=None,sigclip=True,zoom=None):
#     if zoom is not None:
#         print("who's zoomin' who?")
#         # display central region of image
        
#         # get image dimensions and center
#         xmax,ymax = image.shape
#         xcenter = int(xmax/2)
#         ycenter = int(ymax/2)
        
#         # calculate new size to display based on zoom factor
#         new_xradius = int(xmax/2/(float(zoom)))
#         new_yradius = int(ymax/2/(float(zoom)))
        
#         # calculate pixels to keep based on zoom factor
#         x1 = xcenter - new_xradius
#         x2 = xcenter + new_xradius
#         y1 = ycenter - new_yradius
#         y2 = ycenter + new_yradius
        
#         # check to make sure limits are not outsize image dimensions
#         if (x1 < 1):
#             x1 = 1
#         if (y1 < 1):
#             y1 = 1
#         if (x2 > xmax):
#             x2 = xmax
#         if (y2 > ymax):
#             y2 = ymax

#         # cut images to new size
#         image = image[x1:x2,y1:y2]
#         if mask is not None:
#             mask = mask[x1:x2,y1:y2]
#     # use inner 80% of image
#     xdim,ydim = image.shape
#     xmin = int(.1*xdim)
#     xmax = int(.9*xdim)    
#     ymin = int(.1*ydim)
#     ymax = int(.9*ydim)
#     if mask is not None:
#         imdata = np.ma.array(image,mask=mask)
        
#     else:
#         imdata = image
#     v1 = scoreatpercentile(imdata,percentile1)    
#     v2 = scoreatpercentile(imdata,percentile2)
    
#     if mask is not None:
#         statim = image[~mask]
#     else:
#         statim = image

#     if sigclip:
#         if mask is not None:
#             clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax][~mask[xmin:xmax,ymin:ymax]],sigma_lower=1.5,sigma_upper=1.5,grow=10,stdfunc='mad_std')
#         else:
#             clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=1.5,sigma_upper=1.5,grow=10,stdfunc='mad_std')            
#     else:
#         clipped_data = image[xmin:xmax,ymin:ymax]

#     norm = simple_norm(clipped_data, stretch=stretch,max_percent=percentile2,min_percent=percentile1)

#     plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')#,vmin=v1,vmax=v2)
#     #plt.imshow(imdata, norm=norm,origin='lower')#,vmin=v1,vmax=v2)
    

def make_png(fitsimage,outname,mask=None,ellipseparams=None,zoom=None):
    imdata,imheader = fits.getdata(fitsimage,header=True)
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(projection=wcs.WCS(imheader))
    plt.subplots_adjust(top=.95,right=.95,left=.2,bottom=.15)
    display_image(imdata,sigclip=True,mask=mask,zoom=zoom)
    plt.xlabel('RA (deg)',fontsize=16)
    plt.ylabel('DEC (deg)',fontsize=16)
    if ellipseparams is not None:
        ax = plt.gca()
        ellipseparams[-1] += 90
        plot_ellipse(ax,ellipseparams)
    plt.savefig(outname)        
    plt.close(fig)

def make_mask_png(fitsimage,outname,ellipseparams=None,zoom=None):
    imdata,imheader = fits.getdata(fitsimage,header=True)
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(projection=wcs.WCS(imheader))
    plt.subplots_adjust(top=.95,right=.95,left=.2,bottom=.15)
    plt.imshow(imdata, cmap='viridis')
    plt.xlabel('RA (deg)',fontsize=16)
    plt.ylabel('DEC (deg)',fontsize=16)
    if ellipseparams is not None:
        ax = plt.gca()
        ellipseparams[-1] += 90
        plot_ellipse(ax,ellipseparams)
    plt.savefig(outname)        
    plt.close(fig)
    
def plot_ellipse(ax,ellipseparams):

    xc,yc,r,BA,PA = ellipseparams

    #print("just checking - adding ellipse drawing ",self.ellipseparams)
    
    # need to reset b to be consistent with galfit ellipticity
    BA = float(BA)
    PA = float(PA)
    #print('THETA inside phot wrapper',THETA, BA)
    b = BA*r
    eps = 1 - BA
    #print(self.b,self.eps,self.sma,BA)
    t = PA
    if t < 0:
        theta = (180. + t)
    else:
        theta = (t) # orientation in radians
    # EllipticalAperture gives rotation angle in radians from +x axis, CCW
    # matplotlib uses total width and height, not semi-major /minor axes
    ellipse = Ellipse((xc,yc), 2*r, 2*b, angle=theta,facecolor='None',edgecolor='r',lw=2)
    ax.add_patch(ellipse)

    
def display_galfit_model(galfile,prefix="",percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,cmap='viridis',zoom=None,outdir=None,mask=None,ellipseparams=None):
      '''
      ARGS:
      galfile = galfit output image (with image, model, residual)
      percentile1 = min percentile for stretch of image and model
      percentile2 = max percentile for stretch of image and model
      p1residual = min percentile for stretch of residual
      p2residual = max percentile for stretch of residual
      cmap = colormap, default is viridis
      mask = bad pixel mask, with bad values set to True
      '''
      # model name

      #print("inside display_galfit_model, mask = ",mask)
      image,h = fits.getdata(galfile,1,header=True)
      model = fits.getdata(galfile,2)
      residual = fits.getdata(galfile,3)

      if zoom is not None:
         print("who's zoomin' who?")
         # display central region of image

         # get image dimensions and center
         ymax,xmax = image.shape
         xcenter = int(xmax/2)
         ycenter = int(ymax/2)

         # calculate new size to display based on zoom factor
         new_xradius = int(xmax/2/(float(zoom)))
         new_yradius = int(ymax/2/(float(zoom)))

         # calculate pixels to keep based on zoom factor
         x1 = xcenter - new_xradius
         x2 = xcenter + new_xradius
         y1 = ycenter - new_yradius
         y2 = ycenter + new_yradius
         
         # check to make sure limits are not outsize image dimensions
         if (x1 < 1):
            x1 = 1
         if (y1 < 1):
            y1 = 1
         if (x2 > xmax):
            x2 = xmax
         if (y2 > ymax):
            y2 = ymax

         # cut images to new size
         # python is data[row,col]
         image = image[y1:y2,x1:x2]
         model = model[y1:y2,x1:x2]
         residual = residual[y1:y2,x1:x2]         
         pass
      imwcs = wcs.WCS(h)
      images = [image,model,residual]
      titles = ['image','model','residual']
      if mask is not None:
          #print("\nshape of image = ",image.shape)
          #print("shape of mask = ",mask.shape)
          #print()
          im = image[~mask]
          res = residual[~mask]
          norms = [simple_norm(im,'asinh',max_percent=percentile2),
                   simple_norm(im,'asinh',max_percent=percentile2),
                   simple_norm(res,'asinh',max_percent=percentile2,min_percent=20)]

      else:
          norms = [simple_norm(image,'asinh',max_percent=percentile2),
                   simple_norm(image,'asinh',max_percent=percentile2),
                   simple_norm(residual,'asinh',max_percent=percentile2)]

      outim = [prefix+imname for imname in ['galfit_image.png','galfit_model.png','galfit_residual.png']]
      if outdir is not None:
          outim = [os.path.join(outdir,f) for f in outim]
      for i,im in enumerate(images):
          fig = plt.figure(figsize=(6,6))          
          plt.subplot(1,1,1,projection=imwcs)
          plt.subplots_adjust(top=.95,right=.95,left=.2,bottom=.15)
          plt.imshow(im,origin='lower',cmap=cmap,norm=norms[i])
          #plt.colorbar(fraction=.08)
          plt.xlabel('RA (deg)',fontsize=16)
          plt.ylabel('DEC (deg)',fontsize=16)
          #plt.title(titles[i],fontsize=16)
          
          # TODO add ellipse to the residual image
          if (i == 2) and (ellipseparams is not None):

              # plot the ellipse
              plot_ellipse(plt.gca(),ellipseparams)
          plt.savefig(outim[i])
          #plt.close(fig)

###########################################################
####  CLASSES
###########################################################

    
    
class cutout_dir():

    def __init__(self,cutoutdir=None,outdir=None):
        '''
        INPUT:
        * directory containing cutouts
        * output directory for png images

        This creates the png images for different cutouts
        '''
        self.gname = os.path.basename(os.path.dirname(cutoutdir))
        if len(self.gname) < 1:
            self.gname = cutoutdir
        self.vfid = self.gname.split('-')[0]
        #print('inside cutoutdir, gname = ',self.gname)
        #print('cutoutdir = ',cutoutdir)
        #print('outdir = ',outdir)        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.outdir = outdir
        self.cutoutdir = cutoutdir

        # optional products that may or may not be created
        self.results = None
        self.results_file = None

        self.legacy_jpg = None

        self.galimage = None
        self.galmodel = None
        self.galresidual = None

        self.cgalimage = None
        self.cgalmodel = None
        self.cgalresidual = None

        self.r_phot = None
        self.cs_phot = None
        self.phot_tables_ok = False

        self.efluxsma_png = None
        self.emagsma_png = None
        self.sbfluxsma_png = None
        self.sbmagsma_png = None

        self.ellipseparams = None
        self.galfit_ellipseparams = None

        self.legacy_flag = False
        self.wise_flag = False
        self.galex_flag = False
        self.nuv_flag = False
        
    def get_results_table(self):
        """Read the per-galaxy HAPY results table."""
        from astropy.table import Table
        from pathlib import Path

        results_files = list(Path(self.cutoutdir).glob("*-results.ecsv"))
        if len(results_files) == 0:
            self.results = None
            self.results_file = None
            print(f"WARNING: no results.ecsv found in {self.cutoutdir}")
            return

        self.results_file = results_files[0]
        tab = Table.read(self.results_file, format="ascii.ecsv")
        self.results = tab[0] if len(tab) > 0 else None
        #print("DEBUG: found results.ecsv file", self.results_file)
        #print("DEBUD: results colnames:")
        #print(self.results.colnames)
        
    def runall(self):
        self.get_results_table()                
        self.get_halpha_names()

        self.get_ellipse_params()
        try:
            self.get_legacy_names()
            self.legacy_flag = True
            self.get_legacy_jpg()                    
        except IndexError:
            print('WARNING: problem with legacy images')
            self.legacy_flag = False
        try:
            self.get_wise_names()
            self.wise_flag = True            
        except:
            print('WARNING: problem with wise images')
            self.wise_flag = False
        try:
            self.get_galex_names()
            self.galex_flag = True            
        except:
            print('WARNING: problem with wise images')
            self.galex_flag = False
        self.make_png_plots()
        self.make_cs_png()
        #self.make_cs_png(gr=True)
        #self.make_cs_png(grauto=True)                
        #self.get_galfit_model()
        self.get_galfit_images()

        self.get_cgalfit_images()
        self.get_galfit_results_nc()
        self.get_galfit_results_cv()        
        

        self.read_phot_tables()
        if getattr(self, "phot_tables_ok", False):
            self.plot_phot_tables()
        else:
            print("WARNING: skipping phot profile plots because phot tables are missing")
        
        #try:
        #    self.get_phot_tables()
        #except FileNotFoundError:
        #    print('WARNING: no phot files found - check this out')
    def get_halpha_names(self):
        search_string = os.path.join(self.cutoutdir,self.gname+'*-R.fits')
        #print(search_string)
        t = glob.glob(search_string)
        #print(t)
        
        self.rimage = t[0]
        self.haimage = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-Ha.fits'))[0]
        self.csimage = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-CS-ZP.fits'))[0]
        #self.csgrimage = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-CS-gr.fits'))[0]
        #self.csgrimageauto = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-CS-gr-auto.fits'))[0]
        self.csgrimage = None
        self.csgrimageauto = None
        self.maskimage = self.rimage.replace('-R.fits','-mask.fits').replace('-r.fits','-mask.fits')

        try:
            self.conscale_auto = fits.getheader(self.csgrimageauto)['CONSCALE']
        except:
            self.conscale_auto = -99


    def get_ellipse_params(self):
        """
        Get ellipse parameters for plotting.

        Priority:
            1) values stored in results.ecsv
            2) values stored in mask FITS header
        """

        # --- first try results.ecsv ---
        if getattr(self, "results", None) is not None:

            try:
                xc = self.results["ELLIP_XCENTROID"]
                yc = self.results["ELLIP_YCENTROID"]
                r  = self.results["ELLIP_SMA_PIX"]
                ba = self.results["ELLIP_BA"]
                pa = np.degrees(self.results["ELLIP_THETA_RAD"])+90

                if np.isfinite(r) and r > 0:
                    self.ellipseparams = [xc, yc, r, ba, pa]
                    print("DEBUG: ellipseparams from results table = ",self.ellipseparams)
                    return

            except Exception:
                print("WARNING: could not get ellipse params!")
                pass

        # --- fallback to mask header ---
        try:

            header = fits.getheader(self.maskimage)

            xc = header["ELLIP_XC"]
            yc = header["ELLIP_YC"]
            r  = header["ELLIP_A"]
            ba = header["ELLIP_BA"]
            pa = header["ELLIP_PA"]

            self.ellipseparams = [xc, yc, r, ba, pa]

        except Exception:

            print("\nWARNING: could not determine ellipse parameters")
            print(self.maskimage)

            self.ellipseparams = None
        

        
    def get_legacy_names(self):
        ''' get names of legacy images  '''
        legdir = os.path.join(self.cutoutdir,'legacy')
        self.legacy_g = glob.glob(legdir+'/*-g.fits')[0]
        self.legacy_r = glob.glob(legdir+'/*-r.fits')[0]
        self.legacy_z = glob.glob(legdir+'/*-z.fits')[0]
    def get_legacy_jpg(self):
        ''' copy jpg to local directory '''
        legdir = os.path.join(self.cutoutdir,'legacy')        
        legacy_jpg = glob.glob(legdir+'/*.jpg')[0]
        jpeg_data = Image.open(legacy_jpg)
        fig = plt.figure(figsize=(6,6))


        #hdu = fits.open(filename)[0]
        header = fits.getheader(self.legacy_g)
        imwcs = wcs.WCS(header)
        plt.subplot(projection=imwcs)
        plt.subplots_adjust(left=.2,bottom=.15,top=.95,right=.95)        
        plt.imshow(jpeg_data, origin='lower')
        plt.xlabel('RA (deg)',fontsize=16)
        plt.ylabel('Dec (deg)',fontsize=16)
        self.legacy_jpg = os.path.join(self.outdir,os.path.basename(legacy_jpg))        
        plt.savefig(self.legacy_jpg)
        plt.close(fig)

    def get_wise_names(self):
        ''' get names of wise images  '''
        # need to fix this to check for combined unwise images
        wisedir = os.path.join(self.cutoutdir,'unwise')
        self.w1 = glob.glob(wisedir+'/*-w1-img-m.fits')[0]
        self.w2 = glob.glob(wisedir+'/*-w2-img-m.fits')[0]
        self.w3 = glob.glob(wisedir+'/*-w3-img-m.fits')[0]
        self.w4 = glob.glob(wisedir+'/*-w4-img-m.fits')[0]                        
    def get_galex_names(self):
        ''' get names of galex images  '''        
        galdir = os.path.join(self.cutoutdir,'galex')
        galexfiles = glob.glob(galdir+'/*nuv*.fits')
        if len(galexfiles) > 0:
            for f in galexfiles:
                if f.find('nuv') > -1:
                    self.nuv = f
                    self.nuv_flag = True
        else:
            self.nuv_flag = False

    def define_png_names(self):
        pass
    def make_png_plots(self, zoom=None):
        # fitsimages and pngimages should be dictionaries
        # so I am not relying on where they are in the list
        imnames = ['r','ha','cs','csgr','csgr_auto','legacy_g','legacy_r','legacy_z',\
                   'w1','w2','w3','w4',\
                   'mask','nuv']
        self.image_keys = imnames
        # build dictionaries to store fits and png images,
        # setting to None if image is not available
        self.fitsimages = {i:None for i in imnames}
        self.pngimages = {i:None for i in imnames}
        keys = imnames[0:4]
        imlist = [self.rimage,self.haimage,self.csimage,self.csgrimage,self.csgrimageauto]
        for i,k in enumerate(keys):
            self.fitsimages[k] = imlist[i]
        
        if self.legacy_flag:
            keys = imnames[5:8]
            imlist = [self.legacy_g,self.legacy_r,self.legacy_z]
            for i,k in enumerate(keys):
                self.fitsimages[k] = imlist[i]
            
        if self.wise_flag:
            keys = imnames[8:12]            
            imlist = [self.w1,self.w2,self.w3,self.w4]
            for i,k in enumerate(keys):
                self.fitsimages[k] = imlist[i]
            
        self.fitsimages['mask'] = self.maskimage
        if self.nuv_flag:
            self.fitsimages['nuv'] = self.nuv

        mask = fits.getdata(self.maskimage)
        mask = mask > 0
        # plt.figure()
        # plt.imshow(mask)
        # plt.savefig("mask.png")
        # plt.close("all")
        key_list = list(self.fitsimages.keys())
        for i,f in enumerate(self.fitsimages): # loop over keys

            try:
                pngfile = os.path.join(self.outdir,os.path.basename(self.fitsimages[f]).replace('.fits','.png'))
            except TypeError:
                print(f"WARNING: problem with {f}")
                continue

            try:
                if i < 4:
                    make_png(self.fitsimages[f],pngfile,mask=mask)

                elif key_list[i] == 'mask':
                    if self.ellipseparams is not None:
                        make_mask_png(self.fitsimages[f],pngfile,ellipseparams=self.ellipseparams,zoom=zoom)
                    else:
                        make_mask_png(self.fitsimages[f],pngfile,zoom=zoom)
                else:
                    make_png(self.fitsimages[f],pngfile,zoom=zoom)                    
                self.pngimages[f] = pngfile
            except FileNotFoundError:
                print('WARNING: can not find ',self.fitsimages[f])

            except TypeError:
                print('WARNING: problem making png for ',self.fitsimages[f])


    def make_cs_png(self,gr=False,grauto=False,zoom=None):
        if gr:
            csdata,csheader = fits.getdata(self.csgrimage,header=True)
            #imx,imy,keepflag = get_galaxies_fov(self.csimage,vfmain['RA'],vfmain['DEC'])
        elif grauto:
            csdata,csheader = fits.getdata(self.csgrimageauto,header=True)
            #imx,imy,keepflag = get_galaxies_fov(self.csimage,vfmain['RA'],vfmain['DEC'])
        else:
            csdata,csheader = fits.getdata(self.csimage,header=True)
            #imx,imy,keepflag = get_galaxies_fov(self.csimage,vfmain['RA'],vfmain['DEC'])
            
        mask = fits.getdata(self.maskimage)
        mask = mask > 0
        #galsize=60/(abs(csheader['CD1_1'])*3600)        
        p2 = [99.5,99.9]
        percentile1 = .5
        percentile2 = 99.5
        
        stretchs = ['asinh','linear']
        for i,s in enumerate(stretchs):
            fig = plt.figure(figsize=(6,6))
            plt.subplot(projection = wcs.WCS(csheader))
            plt.subplots_adjust(bottom=.15,left=.2,right=.95,top=.95)
            ax = plt.gca()
            
            #clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=1.5,sigma_upper=1.5,grow=10,stdfunc='mad_std')
            if s == "linear":
                display_image(csdata,percent=p2[i],mask=mask,zoom=zoom,lowrange=True)
            else:
                display_image(csdata,percent=p2[i],mask=mask,zoom=zoom)
            # mark VF galaxie
            #plot_vf_gals(imx,imy,keepflag,vfmain,ax,galsize=galsize)
            suffix = "-{}.png".format(p2[i])
            if gr:
                pngfile = os.path.join(self.outdir,os.path.basename(self.csgrimage).replace('.fits',suffix))
            if grauto:
                pngfile = os.path.join(self.outdir,os.path.basename(self.csgrimageauto).replace('.fits',suffix))
            else:
                pngfile = os.path.join(self.outdir,os.path.basename(self.csimage).replace('.fits',suffix))
            plt.xlabel('RA (deg)',fontsize=16)
            plt.ylabel('DEC (deg)',fontsize=16)        
            plt.savefig(pngfile)
            plt.close(fig)
            if i == 0:
                if gr:
                    self.csgr_png1 = pngfile
                elif grauto:
                    self.csgrauto_png1 = pngfile
                else:
                    self.cs_png1 = pngfile
            elif i == 1:
                if gr:
                    self.csgr_png2 = pngfile
                elif grauto:
                    self.csgrauto_png2 = pngfile
                else:
                    self.cs_png2 = pngfile 

 

    def get_galfit_results_nc(self):
        """
        Read GALFIT results from results.ecsv and store on the cutout object.
        """

        if getattr(self, "results", None) is None:
            print("WARNING: results.ecsv not loaded — cannot read GALFIT results")
            return

        r = self.results

        def getval(key):
            try:
                val = r[key]
                return val if np.isfinite(val) else np.nan
            except Exception:
                return np.nan

        self.xc = getval("GAL_XC")
        self.xc_err = getval("GAL_XC_ERR")

        self.yc = getval("GAL_YC")
        self.yc_err = getval("GAL_YC_ERR")

        self.mag = getval("GAL_MAG")
        self.mag_err = getval("GAL_MAG_ERR")

        self.re = getval("GAL_RE")
        self.re_err = getval("GAL_RE_ERR")

        self.nsersic = getval("GAL_N")
        self.nsersic_err = getval("GAL_N_ERR")

        self.BA = getval("GAL_BA")
        self.BA_err = getval("GAL_BA_ERR")

        self.PA = getval("GAL_PA")
        self.PA_err = getval("GAL_PA_ERR")

        self.sky = getval("GAL_SKY")
        self.sky_err = getval("GAL_SKY_ERR")

        self.error = int(getval("GAL_ERROR"))
        self.chi2nu = getval("GAL_CHISQ")

    def get_galfit_results_cv(self):
        """
        Read convolved GALFIT results from results.ecsv.
        """

        if getattr(self, "results", None) is None:
            print("WARNING: results.ecsv not loaded — cannot read GALFIT CV results")
            return

        r = self.results

        def getval(key):
            try:
                val = r[key]
                return val if np.isfinite(val) else np.nan
            except Exception:
                return np.nan

        self.cxc = getval("GAL_CXC")
        self.cxc_err = getval("GAL_CXC_ERR")
        self.cyc = getval("GAL_CYC")
        self.cyc_err = getval("GAL_CYC_ERR")
        self.cmag = getval("GAL_CMAG")
        self.cmag_err = getval("GAL_CMAG_ERR")
        self.cre = getval("GAL_CRE")
        self.cre_err = getval("GAL_CRE_ERR")
        self.cnsersic = getval("GAL_CN")
        self.cnsersic_err = getval("GAL_CN_ERR")
        self.cBA = getval("GAL_CBA")
        self.cBA_err = getval("GAL_CBA_ERR")
        self.cPA = getval("GAL_CPA")
        self.cPA_err = getval("GAL_CPA_ERR")
        self.csky = getval("GAL_CSKY")
        self.csky_err = getval("GAL_CSKY_ERR")
        self.cerror = getval("GAL_CERROR")
        self.cchi2nu = getval("GAL_CCHISQ")


    def get_galfit_ellipse_params(self):
        """Construct ellipse parameters from GALFIT results."""

        if getattr(self, "results", None) is None:
            self.galfit_ellipseparams = None
            return

        r = self.results

        try:
            if r["GAL_CV_OK"]:
                xc = r["GAL_CXC"]
                yc = r["GAL_CYC"]
                sma = r["GAL_CRE"]
                ba = r["GAL_CBA"]
                pa = pa_ccw_north_to_photutils_theta(r["GAL_CPA"])
            else:
                xc = r["GAL_XC"]
                yc = r["GAL_YC"]
                sma = r["GAL_RE"]
                ba = r["GAL_BA"]
                pa = pa_ccw_north_to_photutils_theta(r["GAL_PA"])

            self.galfit_ellipseparams = [xc, yc, sma, ba, pa]

        except Exception:
            self.galfit_ellipseparams = None
        

    def get_galfit_images(self):
        """Read in GALFIT model and make png."""
        self.galfit = self.rimage.replace('.fits', '-1Comp-galfit-out.fits')
        if not os.path.exists(self.galfit):
            self.galfit = self.rimage.replace('-R.fits', '-1Comp-galfit-out.fits')

        self.get_galfit_ellipse_params()

        if os.path.exists(self.galfit):
            mask = fits.getdata(self.maskimage)
            mask = mask > 0

            display_galfit_model(
                self.galfit,
                outdir=self.outdir,
                mask=mask,
                ellipseparams=self.galfit_ellipseparams
            )

            outim = ['galfit_image.png', 'galfit_model.png', 'galfit_residual.png']

            self.galimage = os.path.join(self.outdir, outim[0])
            self.galmodel = os.path.join(self.outdir, outim[1])
            self.galresidual = os.path.join(self.outdir, outim[2])

        else:
            warnings.warn(f"No image {self.galfit}")
            self.galimage = None
            self.galmodel = None
            self.galresidual = None
    



    def get_cgalfit_images(self):
        ''' read in galfit model and make png '''
        self.cgalfit = self.rimage.replace('.fits','-1Comp-galfit-out-conv.fits')
        if not os.path.exists(self.cgalfit):
            self.cgalfit = self.rimage.replace('-R.fits','-1Comp-galfit-out-conv.fits')
        self.get_galfit_ellipse_params()
        if os.path.exists(self.cgalfit):
            # store fit results

            mask = fits.getdata(self.maskimage)
            mask = mask > 0


            display_galfit_model(self.cgalfit,prefix='c',outdir=self.outdir,mask=mask,ellipseparams=self.galfit_ellipseparams)

            outim = ['cgalfit_image.png','cgalfit_model.png','cgalfit_residual.png']
        
            self.cgalimage = os.path.join(self.outdir,outim[0])
            self.cgalmodel = os.path.join(self.outdir,outim[1])
            self.cgalresidual = os.path.join(self.outdir,outim[2])
        else:
            warnings.warn(f"No image {self.cgalfit}")
            self.cgalimage = None
            self.cgalmodel = None
            self.cgalresidual = None

        print("self.cgalimage = ",self.cgalresidual)

    def read_phot_tables(self):
        """Read photutils photometry tables for R and Halpha images."""
        from astropy.table import Table

        self.r_phot_file = self.rimage.replace('.fits', '_phot.fits')
        self.cs_phot_file = self.csimage.replace('.fits', '_phot.fits')

        print("phot files = ", self.cs_phot_file)

        self.r_phot = None
        self.cs_phot = None
        self.phot_tables_ok = False

        missing = []

        if os.path.exists(self.r_phot_file):
            try:
                self.r_phot = Table.read(self.r_phot_file)
            except Exception as e:
                print(f"WARNING: could not read {self.r_phot_file}: {e}")
        else:
            missing.append(self.r_phot_file)

        if os.path.exists(self.cs_phot_file):
            try:
                self.cs_phot = Table.read(self.cs_phot_file)
            except Exception as e:
                print(f"WARNING: could not read {self.cs_phot_file}: {e}")
        else:
            missing.append(self.cs_phot_file)

        if missing:
            print("WARNING: missing phot tables:")
            for m in missing:
                print("   ", m)

        if (self.r_phot is not None) and (self.cs_phot is not None):
            self.phot_tables_ok = True



    def plot_phot_tables(self):
        """Plot flux, magnitude, and surface-brightness profiles from photutils tables."""

        if (self.r_phot is None) or (self.cs_phot is None):
            print("WARNING: phot tables not available; skipping profile plots")
            return
        
        tabs = [self.r_phot, self.cs_phot]
        labels = ['photutils r', 'photutils Halpha x100']
        alphas = [1.0, 0.4]
        mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        def get_plotflag(tab):
            if 'snr_per_pixel' in tab.colnames:
                return np.isfinite(tab['snr_per_pixel']) & (tab['snr_per_pixel'] > 2)
            elif 'sb_avg_snr' in tab.colnames:
                return np.isfinite(tab['sb_avg_snr']) & (tab['sb_avg_snr'] > 2)
            else:
                return np.ones(len(tab), dtype=bool)

        # --------------------------------------------------
        # enclosed flux
        # --------------------------------------------------
        fig = plt.figure(figsize=(6, 6))
        plt.subplots_adjust(left=.15, bottom=.1, right=.95, top=.95)

        plotflag = tabs[0]['snr_per_pixel'] > 2
        for i, t in enumerate(tabs):
            

            x = np.asarray(t['sma_arcsec'])[plotflag]
            y0 = np.asarray(t['flux_cgs'])[plotflag]
            yerr = np.asarray(t['flux_cgs_err'])[plotflag]

            y1 = y0 + yerr
            y2 = y0 - yerr

            if i == 1:
                y0 = y0 * 100
                y1 = y1 * 100
                y2 = y2 * 100

            good = np.isfinite(x) & np.isfinite(y0) & np.isfinite(y1) & np.isfinite(y2) & (y0 > 0)
            if np.any(good):
                plt.fill_between(x[good], y1[good], y2[good],
                                 label=labels[i], alpha=alphas[i], color=mycolors[i])
                plt.plot(x[good], y0[good], '-', lw=2, color=mycolors[i])

        plt.xlabel('SMA (arcsec)', fontsize=16)
        plt.ylabel('Flux (erg/s/cm$^2$)', fontsize=16)
        plt.gca().set_yscale('log')
        plt.legend(loc='lower right')
        self.efluxsma_png = os.path.join(self.outdir, self.gname + '-enclosed-flux.png')
        plt.savefig(self.efluxsma_png)
        plt.close(fig)

        # --------------------------------------------------
        # enclosed magnitude
        # --------------------------------------------------
        fig = plt.figure(figsize=(6, 6))
        plt.subplots_adjust(left=.15, bottom=.1, right=.95, top=.95)

        labels_mag = ['photutils r', 'photutils Halpha']

        for i, t in enumerate(tabs):
            #plotflag = get_plotflag(t)

            x = np.asarray(t['sma_arcsec'])[plotflag]
            y0 = np.asarray(t['mag_cum'])[plotflag]
            yerr = np.asarray(t['mag_cum_err'])[plotflag]

            y1 = y0 + yerr
            y2 = y0 - yerr

            good = np.isfinite(x) & np.isfinite(y0) & np.isfinite(y1) & np.isfinite(y2)
            if np.any(good):
                plt.fill_between(x[good], y1[good], y2[good],
                                 label=labels_mag[i], alpha=alphas[i], color=mycolors[i])
                plt.plot(x[good], y0[good], '-', lw=2, color=mycolors[i])

        plt.xlabel('SMA (arcsec)', fontsize=16)
        plt.ylabel('Magnitude (AB)', fontsize=16)
        plt.gca().invert_yaxis()
        plt.legend(loc='lower right')
        self.emagsma_png = os.path.join(self.outdir, self.gname + '-mag-sma.png')
        plt.savefig(self.emagsma_png)
        plt.close(fig)

        # --------------------------------------------------
        # surface brightness in cgs
        # --------------------------------------------------
        fig = plt.figure(figsize=(6, 6))
        plt.subplots_adjust(left=.15, bottom=.1, right=.95, top=.95)

        for i, t in enumerate(tabs):
            #plotflag = get_plotflag(t)

            x = np.asarray(t['sma_arcsec'])[plotflag]
            y0 = np.asarray(t['sb_cgs_arcsec2'])[plotflag]
            yerr = np.asarray(t['sb_cgs_arcsec2_err'])[plotflag]

            y1 = y0 + yerr
            y2 = y0 - yerr

            if i == 1:
                y0 = y0 * 100
                y1 = y1 * 100
                y2 = y2 * 100

            good = np.isfinite(x) & np.isfinite(y0) & np.isfinite(y1) & np.isfinite(y2) & (y0 > 0)
            if np.any(good):
                plt.fill_between(x[good], y1[good], y2[good],
                                 label=labels[i], alpha=alphas[i], color=mycolors[i])
                plt.plot(x[good], y0[good], '-', lw=2, color=mycolors[i])

        plt.xlabel('SMA (arcsec)', fontsize=16)
        plt.ylabel('SB (erg/s/cm$^2$/arcsec$^2$)', fontsize=16)
        plt.gca().set_yscale('log')
        plt.legend()
        self.sbfluxsma_png = os.path.join(self.outdir, self.gname + '-sb-sma.png')
        plt.savefig(self.sbfluxsma_png)
        plt.close(fig)

        # --------------------------------------------------
        # surface brightness in mag / arcsec^2
        # --------------------------------------------------
        fig = plt.figure(figsize=(6, 6))
        plt.subplots_adjust(left=.15, bottom=.1, right=.95, top=.95)

        labels_sbmag = ['photutils r', 'photutils Halpha']

        for i, t in enumerate(tabs):
            #plotflag = get_plotflag(t)

            x = np.asarray(t['sma_arcsec'])[plotflag]
            y0 = np.asarray(t['sb_mag_arcsec2'])[plotflag]
            yerr = np.asarray(t['sb_mag_arcsec2_err'])[plotflag]

            y1 = y0 + yerr
            y2 = y0 - yerr

            good = np.isfinite(x) & np.isfinite(y0) & np.isfinite(y1) & np.isfinite(y2)
            if np.any(good):
                plt.fill_between(x[good], y1[good], y2[good],
                                 label=labels_sbmag[i], alpha=alphas[i], color=mycolors[i])
                plt.plot(x[good], y0[good], '-', lw=2, color=mycolors[i])

        plt.xlabel('SMA (arcsec)', fontsize=16)
        plt.ylabel('Surface Brightness (mag/arcsec$^2$)', fontsize=16)
        plt.gca().invert_yaxis()
        plt.legend()
        self.sbmagsma_png = os.path.join(self.outdir, self.gname + '-sbmag-sma.png')
        plt.savefig(self.sbmagsma_png)
        plt.close(fig)
    
    def get_phot_tables(self):
        ''' read in phot tables and make plot of flux and sb vs sma '''


        # open data files
        cs_galfit_phot = self.csimage.replace('.fits','-GAL_phot.fits')
        cs_gphot = fits.getdata(cs_galfit_phot)
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_galfit_phot = self.rimage.replace('.fits','-GAL_phot.fits')
        r_gphot = fits.getdata(r_galfit_phot)

        # photutils flux
        cs_galfit_phot = self.csimage.replace('.fits','_phot.fits')
        cs_phot = fits.getdata(cs_galfit_phot)
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_galfit_phot = self.rimage.replace('.fits','_phot.fits')
        r_phot = fits.getdata(r_galfit_phot)
        

        # define colors - need this for plotting line and fill_between in the same color
        mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        # plot enclosed flux        
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
        alphas = [1,.4,.6,.4]
        for i,t in enumerate(tabs):
            y0 = t['flux_erg']            
            y1 = t['flux_erg']+t['flux_erg_err']
            y2 = t['flux_erg']-t['flux_erg_err']

            if (i == 1) + (i == 3):
                y0=y0*100
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])

        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend(loc='lower right')        
        self.efluxsma_png = os.path.join(self.outdir,self.gname+'-enclosed-flux.png')
        plt.savefig(self.efluxsma_png)
        plt.close(fig)
        
        # plot mag sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        #tabs = [r_gphot,cs_gphot]
        #labels = ['r','Halpha']
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halpha','photutil r','photutil Halpha']
        ncolor = 0
        for i,t in enumerate(tabs):
            y0 = t['mag']
            y1 = t['mag']+t['mag_err']
            y2 = t['mag']-t['mag_err']
            if (i == 1) + (i == 3):
                alpha=.4
            else:
                alpha=1

            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])

            
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('magnitude (AB)',fontsize=16)
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()        
        self.emagsma_png = os.path.join(self.outdir,self.gname+'-mag-sma.png')
        plt.legend(loc='lower right')        
        plt.savefig(self.emagsma_png)
        plt.close(fig)
        
        # plot sb erg vs sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']

        for i,t in enumerate(tabs):
            y0 = t['sb_erg_sqarcsec']
            y1 = t['sb_erg_sqarcsec']+t['sb_erg_sqarcsec_err']
            y2 = t['sb_erg_sqarcsec']-t['sb_erg_sqarcsec_err']
            if (i == 1) + (i == 3):
                y0 = y0*100
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])

                
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('SB (erg/s/cm^2/Hz/arcsec^2)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend()
        self.sbfluxsma_png = os.path.join(self.outdir,self.gname+'-sb-sma.png')
        plt.savefig(self.sbfluxsma_png)
        plt.close(fig)
                             
        # plot mag sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halpha','photutil r','photutil Halpha']
        
        for i,t in enumerate(tabs):
            y0 = t['sb_mag_sqarcsec']
            y1 = t['sb_mag_sqarcsec']+t['sb_mag_sqarcsec_err']
            y2 = t['sb_mag_sqarcsec']-t['sb_mag_sqarcsec_err']

            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])
            
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Surface Brightness (mag/arcsec^2)',fontsize=16)
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()
        plt.legend()        
        self.sbmagsma_png = os.path.join(self.outdir,self.gname+'-sbmag-sma.png')
        plt.savefig(self.sbmagsma_png)
        plt.close(fig)
    
class build_html_cutout():

    def __init__(self,cutoutdir,outdir,previous=None,next=None,tel=None,run=None):
        ''' pass in instance of cutout_dir class and output directory '''
        #print("in build_html_cutout!")
        self.cutout = cutoutdir


        outfile = os.path.join(outdir,self.cutout.gname+'.html')
        print("outfile = ",outfile)
        #vfindices = np.arange(len(vfmain))
        #self.vfindex = vfindices[vfmain['VFID'] == self.cutout.vfid]
        #print('inside build html')
        #print('coutdir = ',coutdir)
        #print('outfile = ',outfile)        
        self.html = open(outfile,'w')
        self.htmlhome = 'index.html'
        self.next = next
        self.previous = previous

        self.telescope = tel
        self.run=run
        
        # for reference, this is the order of the png images
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\
        #              self.legacy_g,self.legacy_r,self.legacy_z,\
        #              self.w1,self.w2,self.w3,self.w4]

        
        #self.build_html()

    def _get_result(self, key, default=np.nan):
        #t = self.cutout.results[key]
        if getattr(self.cutout, "results", None) is None:
            print(f"WARNING: no {key} in results table")
            return default
        try:
            return self.cutout.results[key]
        except Exception:
            print(f"WARNING: could not retrieve {key}")
            return default

    def _fmt_result(self, key, fmt="{:.2f}", default="--"):
        val = self._get_result(key, np.nan)
        try:
            if np.isfinite(val):
                return fmt.format(val)
        except Exception:
            if isinstance(val, str) and len(val) > 0:
                return val
        return default

    # def _status_cell(self, flag):
    #     print("in status cell, flag= ",flag," flag==True = ",flag==True)
    #     if flag == True:
    #         return '<span class="ok">OK</span>'
    #     elif flag == False:
    #         return '<span class="fail">FAIL</span>'
    #     else:
    #         return '<span class="warn">--</span>'
    

    def _status_cell(self, flag):
        if flag == True:
            return '<span style="color:green; font-size:30px;">&#9679;</span> OK'
        elif flag == False:
            return '<span style="color:red; font-size:30px;">&#9679;</span> FAIL'
        else:
            return '<span style="color:goldenrod; font-size:30px;">&#9679;</span> --'
    def build_html(self):
        print("building html ",self.html)
        self.write_header()
        self.write_navigation_links()
        # adding this here so we can inspect the masks quickly
        # can remove once we are done with masks

        self.write_image_stats()
        #print("Pipeline flags:")
        #for k in ["MASK_OK","PHOT_OK","PSF_OK","GAL_NC_OK","GAL_CV_OK"]:
        #    print(k, self._get_result(k))
        self.write_pipeline_status()
        #self.write_galfit_images()        
            

        #self.write_sfr_images()
        #if self.cutout.wise_flag:
        #    self.write_wise_images()
        self.write_halpha_images()
        
        #if self.cutout.legacy_flag:
        #    self.write_legacy_images()
        
        #if self.cutout.galimage is not None:
        #    self.write_galfit_images()
        #    self.write_galfit_table()
        if getattr(self.cutout, "phot_tables_ok", False):
            self.write_phot_profiles()
        else:
            self.html.write('<h2>Elliptical Photometry</h2>\n')
            self.html.write('<p>Photometry profile files not available.</p>\n')
    
        self.write_mag_table()
        self.write_morph_table()
        self.write_statmorph_table()
        self.write_galfit_images()
        self.write_galfit_table()
        self.write_cgalfit_images()        
        self.write_navigation_links()
        self.close_html()
    def write_header(self):
        # title
        # home link
        # previous
        # next
        self.html.write('<html><body>\n')
        self.html.write('<style type="text/css">\n')
        self.html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')

        self.html.write('.ok {background-color:#c7f5c7}\n')     # green
        self.html.write('.warn {background-color:#fff3b0}\n')   # yellow
        self.html.write('.fail {background-color:#f8c6c6}\n')   # red

        self.html.write('</style>\n')

    def write_navigation_links(self):
        # Top navigation menu--
        self.html.write('<h1>{}</h1>\n'.format(self.cutout.gname))

        self.html.write('<a href="../{}">Home</a>\n'.format(self.htmlhome))
        self.html.write('<br />\n')
        if self.previous is not None:
            previoushtml = "{}.html".format(self.previous)
            self.html.write('<a href="../{}/{}">Previous ({})</a>\n'.format(self.previous,previoushtml, self.previous))
            self.html.write('<br />\n')        
        if self.next is not None:
            nexthtml = "{}.html".format(self.next)
            self.html.write('<a href="../{}/{}">Next ({})</a>\n'.format(self.next,nexthtml,self.next))
            self.html.write('<br />\n')

    def write_image_stats(self):
        self.html.write('<h2>Image Statistics</h2>\n')

        labels = [
            'VFID',
            'Galaxy',
            'HAPY Version',
            'Run Date',
            'Telescope',
            'Run',
            'Pointing',
            'R FWHM <br> (arcsec)',
            'H&alpha; FWHM <br> (arcsec)',
            'Filter Ratio',
            'Filter Correction',
            'Status',
            'Stage',
        ]

        pointing = get_result("POINTING", "")
        if isinstance(pointing, str) and len(pointing) > 0:
            pointing_str = f'<a href="../../coadds/{pointing}/{pointing}.html">{pointing}</a>'
        else:
            pointing_str = "--"

        data = [
            get_result(self.cutout.results,"VFID", self.cutout.vfid),
            get_result(self.cutout.results,"GALNAME", self.cutout.gname),
            get_result(self.cutout.results,"HAPY_VERSION", ""),
            get_result(self.cutout.results,"RUN_DATE", ""),
            get_result(self.cutout.results,"TELESCOPE", self.telescope),
            self.run,
            pointing_str,
            fmt_result(self.cutout.results,"R_FWHM", "{:.2f}"),
            fmt_result(self.cutout.results,"H_FWHM", "{:.2f}"),
            fmt_result(self.cutout.results,"FILTER_RATIO", "{:.4f}"),
            fmt_result(self.cutout.results,"FILTER_CORRECTION", "{:.2f}"),
            get_result(self.cutout.results,"STATUS", ""),
            get_result(self.cutout.results,"STAGE", ""),
        ]

        write_text_table(self.html, labels, data)

    def write_pipeline_status(self):

        self.html.write('<h2>Pipeline Status</h2>\n')

        labels = [
            'Mask',
            'Photometry',
            'PSF',
            'R Profile',
            'Hα Profile',
            'R Statmorph',
            'H&alpha; Statmorph',             
            'GALFIT NC',
            'GALFIT CV',
            ]

        data = [
            status_cell(get_result(self.cutout.results,'MASK_OK')),
            status_cell(get_result(self.cutout.results,'PHOT_OK')),
            status_cell(get_result(self.cutout.results,'PSF_OK')),
            status_cell(get_result(self.cutout.results,'R_PROFILE_OK')),
            status_cell(get_result(self.cutout.results,'HA_PROFILE_OK')),
            status_cell(get_result(self.cutout.results,'R_SM_FLAG')),
            status_cell(get_result(self.cutout.results,'H_SM_FLAG')),            
            status_cell(get_result(self.cutout.results,'GAL_NC_OK')),
            status_cell(get_result(self.cutout.results,'GAL_CV_OK')),
            ]

        write_text_table(self.html, labels, data)
    

        
    def write_sfr_images(self):
        ''' legacy jpeg, galex nuv, halpha, w4 '''
        self.html.write('<h2>Star-Formation Indicators</h2>\n')
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\ # 0,1,2
        #                   self.legacy_g,self.legacy_r,self.legacy_z,\ # 3,4,5
        #                   self.w1,self.w2,self.w3,self.w4] # 6,7,8,9

        if self.cutout.nuv_flag:
            images = [self.cutout.pngimages['nuv'],\
                      self.cutout.cs_png1,\
                      self.cutout.pngimages['w3'],self.cutout.pngimages['w4']]
            labels = ['NUV','H&alpha;','W3','W4']
        else:
            try:
                images = [self.cutout.legacy_jpg,\
                          self.cutout.pngimages['cs'],\
                          self.cutout.pngimages['w3'],self.cutout.pngimages['w4']]                      

                labels = ['Legacy','H&alpha;','W3','W4']
            except IndexError:
                print("WARNING: problem plotting sfr images")
                return
            except AttributeError:
                print("WARNING: problem plotting sfr images, probably with legacy jpg")
                return
        newim = []
        for i in images:
            if i is not None:
                newim.append(os.path.basename(i))
            else:
                newim.append(i)
        images = newim
        #images = [os.path.basename(i) for i in images]            
        write_table(self.html,images=images,labels=labels)
        pass

    
    def write_halpha_images(self):
        '''  r, halpha, cs, and mask images '''
        self.html.write('<h2>Halpha Images</h2>\n')        
        images = [self.cutout.legacy_jpg,self.cutout.pngimages['r'],self.cutout.pngimages['ha'],self.cutout.cs_png1]#,self.cutout.cs_png2]
        #images = [self.cutout.pngimages['r'],self.cutout.pngimages['ha'],self.cutout.cs_png1,self.cutout.csgr_png1,self.cutout.csgrauto_png1]

        # removing r-band
        #images = [self.cutout.pngimages['ha'],self.cutout.cs_png1]#,self.cutout.csgr_png1,self.cutout.csgrauto_png1]
        # just changing order to see if halpha image is still the biggest in the table, re issue #15
        # the second was still the biggest
        # so what if we also change the label
        # seems to scale with label
        #images = [self.cutout.pngimages['ha'],self.cutout.pngimages['r'],self.cutout.cs_png1,self.cutout.cs_png2]        
        images = [os.path.basename(i) for i in images]

        labels = ['Legacy grz','R-band Image','H&alpha;+Cont','CS from ZP']#,'CS, stretch 2']

        #labels = ['R-band Image','H&alpha;+Cont','CS from ZP ratio','CS from ZP and g-r cor',f'CS g-r auto scale={self.cutout.conscale_auto:.2f}']
        #labels = ['H&alpha;+Cont','CS from ZP ratio','CS from ZP and g-r cor',f'CS g-r auto scale={self.cutout.conscale_auto:.2f}']        
        #labels = ['Halpha+Cont','R','CS, stretch 1','CS, stretch 2']        
        write_table(self.html,images=images,labels=labels)

    def write_mstar_sfr_images(self):
        '''  TODO : add panel for stellar mass, sfr and ssfr '''
        self.html.write('<h2>Halpha Images</h2>\n')        
        images = [self.cutout.pngimages['r'],self.cutout.pngimages['ha'],self.cutout.cs_png1,self.cutout.cs_png2]
        # just changing order to see if halpha image is still the biggest in the table, re issue #15
        # the second was still the biggest
        # so what if we also change the label
        # seems to scale with label
        #images = [self.cutout.pngimages['ha'],self.cutout.pngimages['r'],self.cutout.cs_png1,self.cutout.cs_png2]        
        images = [os.path.basename(i) for i in images]

        labels = ['R-band Image','H&alpha;+Cont','CS, stretch 1','CS, stretch 2']
        #labels = ['Halpha+Cont','R','CS, stretch 1','CS, stretch 2']        
        write_table(self.html,images=images,labels=labels)
        
    def write_legacy_images(self):
        ''' jpg, g,r,z legacy images '''
        self.html.write('<h2>Legacy Images</h2>\n')

        images = [self.cutout.pngimages['legacy_g'],self.cutout.pngimages['legacy_r'],self.cutout.pngimages['legacy_z']]
        images = [os.path.basename(i) for i in images]        
        images.insert(0,os.path.basename(self.cutout.legacy_jpg))        
        labels = ['grz','g','r','z']
        write_table(self.html,images=images,labels=labels)

    def write_wise_images(self):
        ''' w1 - w4 images '''
        self.html.write('<h2>WISE Images</h2>\n')
        pngimages = [self.cutout.pngimages['w1'],self.cutout.pngimages['w2'],\
                     self.cutout.pngimages['w3'],self.cutout.pngimages['w4']]
        wlabels = ['W1','W2','W3','W4']
        images=[]
        labels=[]
        for i,im in enumerate(pngimages):
            if im is not None:
                images.append(os.path.basename(im))
                labels.append(wlabels[i])

        write_table(self.html,images=images,labels=labels)
    
    def write_galex_images(self):
        ''' right now just nuv '''
        self.html.write('<h2>GALEX Images</h2>\n')                
        pass

    def write_galfit_images(self):
        ''' display galfit model and fit parameters for r-band image '''
        self.html.write('<h2>GALFIT r-band Modeling </h2>\n')                
        if self.cutout.galimage is not None:
            
            images = [self.cutout.galimage,self.cutout.galmodel,self.cutout.galresidual,\
                      self.cutout.pngimages['mask']]

            cimages = [self.cutout.galimage,self.cutout.galmodel,self.cutout.galresidual,\
                      self.cutout.pngimages['mask']]
            images = [os.path.basename(i) for i in images]        
            labels = ['Image', 'Model', 'Residual','Mask']
            write_table(self.html,images=images,labels=labels)
        else:
            print("no self.cutout.galimage")

    def write_cgalfit_images(self):
        ''' display galfit model and fit parameters for r-band image '''
        self.html.write('<h2>GALFIT r-band Modeling + Convolution </h2>\n')      
        if self.cutout.cgalimage is not None:
            
            images = [self.cutout.cgalimage,self.cutout.cgalmodel,self.cutout.cgalresidual,\
                      self.cutout.pngimages['mask']]

            images = [os.path.basename(i) for i in images]        
            labels = ['Image', 'Model+Conv', 'Residual','Mask']
            write_table(self.html,images=images,labels=labels)
        else:
            print("no self.cutout.cgalimage")


    def write_galfit_table(self):
        """Display GALFIT Sersic parameters from results table."""

        self.html.write('<h3>GALFIT Sersic Parameters</h3>\n')

        labels = [
            'Fit',
            'Status',
            'XC',
            'YC',
            'MAG',
            'RE',
            'N',
            'BA',
            'PA',
            'SKY',
            'CHI2NU',
            ]

        data = [
            'NC',
            status_cell(get_result(self.cutout.results,'GAL_NC_OK')),
            fmt_result(self.cutout.results,'GAL_XC', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_YC', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_MAG', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_RE', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_N', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_BA', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_PA', '{:.1f}'),
            fmt_result(self.cutout.results,'GAL_SKY', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CHISQ', '{:.2f}'),
            ]



        data2 = [
            'CV',
            status_cell(get_result(self.cutout.results,'GAL_CV_OK')),
            fmt_result(self.cutout.results,'GAL_CXC', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CYC', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CMAG', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CRE', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CN', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CBA', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CPA', '{:.1f}'),
            fmt_result(self.cutout.results,'GAL_CSKY', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_CCHISQ', '{:.2f}'),
            ]

        write_text_table(self.html, labels, data, data2=data2)
        #write_text_table(self.html, labels, data, data2=data2)

            
 
    def write_phot_profiles(self):
        ''' photometry table with galfit and photutil results '''
        self.html.write('<h2>Elliptical Photometry</h2>\n')
        #self.html.write('<p>using galfit and photutil geometry</p>\n')                        
        images = [self.cutout.efluxsma_png,self.cutout.emagsma_png,self.cutout.sbfluxsma_png,self.cutout.sbmagsma_png]
        images = [os.path.basename(i) for i in images]
        labels = ['Enclosed Flux','Enclosed Magnitude','Surface Brightness','Surface Brightness']
        write_table(self.html,images=images,labels=labels)

    def write_mag_table(self):
        self.html.write('<h2>R-band Magnitudes / Flux Summaries</h2>\n')

        labels = ['R24', 'R25 iso', 'R25.5', 'Petro', 'GALFIT', 'Segment']
        data = [
            fmt_result(self.cutout.results,'R24_MAG', '{:.2f}'),
            fmt_result(self.cutout.results,'R25_ISO_MAG', '{:.2f}'),
            fmt_result(self.cutout.results,'R25P5_MAG', '{:.2f}'),
            fmt_result(self.cutout.results,'R_PETRO_MAG', '{:.2f}'),
            fmt_result(self.cutout.results,'GAL_MAG', '{:.2f}'),
            fmt_result(self.cutout.results,'ELLIP_SEGMENT_MAG', '{:.2f}'),
        ]
        write_text_table(self.html, labels, data)

        self.html.write('<h2>Halpha Flux / Size Summaries</h2>\n')

        labels = [
            'H&alpha; Tot Flux',
            'H&alpha; Iso 5e-17',
            'H&alpha; Iso 1.7e-17',
            'H&alpha; R24 Flux',
            'H&alpha; C30(R24)',
            'H&alpha; Max Det Radius',
        ]
        data = [
            fmt_result(self.cutout.results,'HA_TOT_FLUX_CGS', '{:.2e}'),
            fmt_result(self.cutout.results,'HA_ISO5E17_FLUX_CGS', '{:.2e}'),
            fmt_result(self.cutout.results,'HA_ISO17E18_FLUX_CGS', '{:.2e}'),
            fmt_result(self.cutout.results,'HA_R24_FLUX_CGS', '{:.2e}'),
            fmt_result(self.cutout.results,'HA_C30_R24', '{:.2f}'),
            fmt_result(self.cutout.results,'HA_MAXDET_ARCSEC', '{:.1f}'),
        ]
        write_text_table(self.html, labels, data)
        
    def write_morph_table(self):
        self.html.write('<h2>Morphology / Profile Summary</h2>\n')

        labels = ['Band', 'Gini', 'M20', 'Asym', 'C30', 'Petro Con', 'Profile OK']

        data = [
            'r',
            fmt_result(self.cutout.results,'ELLIP_GINI_DET', '{:.2f}'),
            fmt_result(self.cutout.results,'R_M20', '{:.2f}'),
            fmt_result(self.cutout.results,'R_ASYM', '{:.2f}'),
            fmt_result(self.cutout.results,'R_C30', '{:.2f}'),
            fmt_result(self.cutout.results,'R_PETRO_CON', '{:.2f}'),
            str(bool(get_result(self.cutout.results,'R_PROFILE_OK', False))),
        ]

        data2 = [
            'Halpha',
            '--',
            fmt_result(self.cutout.results,'H_M20', '{:.2f}'),
            fmt_result(self.cutout.results,'H_ASYM', '{:.2f}'),
            fmt_result(self.cutout.results,'HA_C30_R24', '{:.2f}'),
            fmt_result(self.cutout.results,'HA_PETRO_CON', '{:.2f}'),
            str(bool(get_result(self.cutout.results,'HA_PROFILE_OK', False))),
        ]

        write_text_table(self.html, labels, data, data2=data2)

    def write_statmorph_table(self):
        self.html.write('<h2>Statmorph Parameters</h2>\n')

        labels = ['Band', 'XC', 'YC', 'Gini', 'M20', 'C', 'A', 'S', 'Rhalf']

        data = [
            'r',
            fmt_result(self.cutout.results,'R_SM_XCENTROID', '{:.2f}'),
            fmt_result(self.cutout.results,'R_SM_YCENTROID', '{:.2f}'),
            fmt_result(self.cutout.results,'R_SM_GINI', '{:.2f}'),
            fmt_result(self.cutout.results,'R_SM_M20', '{:.2f}'),
            fmt_result(self.cutout.results,'R_SM_C', '{:.2f}'),
            fmt_result(self.cutout.results,'R_SM_A', '{:.2f}'),
            fmt_result(self.cutout.results,'R_SM_S', '{:.2f}'),
            fmt_result(self.cutout.results,'R_SM_RHALF_ELLIP', '{:.2f}'),
        ]

        data2 = [
            'Halpha',
            fmt_result(self.cutout.results,'H_SM_XCENTROID', '{:.2f}'),
            fmt_result(self.cutout.results,'H_SM_YCENTROID', '{:.2f}'),
            fmt_result(self.cutout.results,'H_SM_GINI', '{:.2f}'),
            fmt_result(self.cutout.results,'H_SM_M20', '{:.2f}'),
            fmt_result(self.cutout.results,'H_SM_C', '{:.2f}'),
            fmt_result(self.cutout.results,'H_SM_A', '{:.2f}'),
            fmt_result(self.cutout.results,'H_SM_S', '{:.2f}'),
            fmt_result(self.cutout.results,'H_SM_RHALF_ELLIP', '{:.2f}'),
        ]

        write_text_table(self.html, labels, data, data2=data2)

    def close_html(self):
        self.html.close()
# wrap

if __name__ == '__main__':

    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description='Build HTML inspection pages for HAPY cutouts')

    parser.add_argument('--cutoutdir', dest='cutoutdir', default=None,
                        help='Directory containing cutout folders')

    parser.add_argument('--oneimage', dest='oneimage', default=None,
                        help='Process only one cutout directory')

    parser.add_argument('--outdir', dest='outdir',
                        default='/data-pool/Halpha/html_dev/cutouts/',
                        help='Output directory for HTML pages')

    parser.add_argument('--hacat', dest='hacat',
                        default='../halphagui-output-combined-2024-Jul-08.fits')

    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    # optional external catalog (may be phased out later)
    #fullha = fits.getdata(args.hacat)

    if args.cutoutdir is not None:
        os.chdir(args.cutoutdir)

    # find cutout directories
    rfiles = sorted([p.name for p in Path().iterdir() if p.is_dir()])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.oneimage is not None:

        if args.oneimage not in rfiles:
            print(f"Could not find {args.oneimage}")
            sys.exit()

        buildone(args.oneimage, outdir, rfiles)

    else:

        import multiprocessing as mp

        pool = mp.Pool(mp.cpu_count())

        results = [
            pool.apply_async(buildone, args=(subdir, outdir, rfiles))
            for subdir in rfiles
        ]

        pool.close()
        pool.join()
