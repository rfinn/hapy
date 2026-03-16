import sys
import wget
import tarfile
import os
import numpy as np
import glob

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from PIL import Image
from scipy.stats import scoreatpercentile

from urllib.parse import urlencode
from urllib.request import urlretrieve

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip
from astroquery.mast import Observations

#from .plotting import display_image

vmin = .5
vmax = 50.

UNWISE_PIXSCALE = 2.75
LEGACY_PIXSCALE = 1

#def get_legacy_images(ra,dec,galid='VFID0',pixscale=1,imsize='60',band='g',makeplots=False,subfolder=None,verbose=False):
def get_legacy_images(
    ra, dec, galid='VFID0', pixscale=0.262, imsize='60', band='grz',
    makeplots=False, subfolder=None, verbose=False, layer='ls-dr9'):
    """
    Download legacy image for a particular ra, dec
    
    Inputs:
    * ra
    * dec
    * galid = galaxy id (e.g. VFID0001); used for naming the image files
    * imsize = size of cutout in pixels
    * band = filter(s) for the fits images that will be returned.
+             e.g. 'g' or 'r' or 'z' or 'grz'
    * pixscale = pixel scale of cutout in arcsec; native is 0.262 for legacy
    * makeplots = boolean, generate plot of image
    * subfolder = default is None; you can specify a name of a subfolder to 
                  save the data in, e.g., subfolder='legacy-images'
    Returns:
    * fits_name = fits image name (single band) or dict of band:file names (multi-band)
    * jpeg_name = jpeg image name
    """
    imsize = int(imsize)

    # make output image names
    if subfolder is not None:
        # check if subfolder exists. if not, make it.
        if not os.path.exists(subfolder):
            os.mkdir(subfolder)
        rootname = subfolder+'/'+str(galid)+'-legacy-'+str(imsize)
    else:
        rootname = str(galid)+'-legacy-'+str(imsize)        
    jpeg_name = rootname+'.jpg'
    fits_name = rootname+'-'+band+'.fits'
    band_fits_names = {b: rootname+'-'+b+'.fits' for b in band}


    #print('legacy imsize = ',imsize)
    
    # check if images already exist
    # if not download images
    if not(os.path.exists(jpeg_name)):
        if verbose:
            print('retrieving ',jpeg_name)
        url = (
            'https://www.legacysurvey.org/viewer/jpeg-cutout?'
            + urlencode({
                'ra': ra,
                'dec': dec,
                'layer': layer,
                'size': imsize,
                'pixscale': pixscale,
                })
            )
        print("legacy download: \n\t",url)
        #url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr9&size='+str(imsize)+'&pixscale='+str(pixscale)
        urlretrieve(url, jpeg_name)
    else:
        if verbose:
            print('previously downloaded ',jpeg_name)

    need_fits_download = False
    if len(band) == 1:
        need_fits_download = not os.path.exists(fits_name)
    else:
        # download the combined MEF if it does not exist
        need_fits_download = not os.path.exists(fits_name)

    if need_fits_download:
    #if not(os.path.exists(fits_name)):
        if verbose:
            print('retrieving ',fits_name)
        url = (
            'https://www.legacysurvey.org/viewer/cutout.fits?'
            + urlencode({
                'ra': ra,
                'dec': dec,
                'layer': layer,
                'size': imsize,
                'pixscale': pixscale,
                'bands': band,
                })
            )
        #url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(pixscale)+'&bands='+band
        if verbose:
            print(url)
        urlretrieve(url, fits_name)
    else:
        if verbose:
            print('previously downloaded ',fits_name)

    # try to read the data in
    try:
        hdul = fits.open(fits_name)

    except Exception:
        print('problem accessing image')
        print(fits_name)
        print(url)
        return None

    # single-band case
    if len(band) == 1:
        try:
            if len(hdul) > 1 and hdul[1].data is not None:
                t = hdul[1].data
                h = hdul[1].header
            else:
                t = hdul[0].data
                h = hdul[0].header
        except Exception:
            print('problem accessing image data')
            print(fits_name)
            hdul.close()
            return None

        # trigger if image is outside footprint
        if t is None or np.all(np.asarray(t) == 0):
            hdul.close()
            return None

    # multi-band case: split extensions 1..N into separate files
    else:
        for i, b in enumerate(band, start=1):
            outname = band_fits_names[b]
            if verbose:
                print(f"working on extracting {outname} from MEF {fits_name}")

            if os.path.exists(outname):
                if verbose:
                    print('previously downloaded ', outname)
                continue

            if i >= len(hdul) or hdul[i].data is None:
                print(f'problem accessing extension {i} for band {b}')
                print(fits_name)
                hdul.close()
                return None

            data = hdul[i].data
            hdr = hdul[i].header

            # outside footprint
            if data is None or np.all(np.asarray(data) == 0):
                hdul.close()
                return None

            fits.writeto(outname, data, header=hdr, overwrite=True)
            if verbose:
                print('wrote ', outname)

        # use g-band as default image for plotting
        t = hdul[1].data
        h = hdul[1].header

    hdul.close()

    # plot the images
    if makeplots:
        if os.path.exists(jpeg_name):
            t = Image.open(jpeg_name)
            plt.imshow(t,origin='upper')
        else:
            norm = simple_norm(t,stretch='asinh',percent=99.5)            
            plt.imshow(t,origin='upper',cmap='gray_r', norm=norm)

    # return the name of the fits images and jpeg image
    if len(band) == 1:
        return fits_name, jpeg_name
    else:
        return band_fits_names, jpeg_name


def get_unwise_image(ra,dec,galid='VFID0',pixscale=2.75,imsize='60',bands='1234',makeplots=False,subfolder=None,verbose=False):
    """
    Download unwise image for a particular ra, dec
    
    Inputs:
    * ra
    * dec
    * galid = galaxy id (e.g. VFID0001)
    * imsize = size of cutout in pixels
    * pixscale = pixel scale of cutout in arcsec
      - native is 0.262 for legacy; 
      - 2.75 for wise
    """
    downloadwise = True
    # check if images already exist
    if subfolder is not None:
        image_names = glob.glob(subfolder+'/'+galid+'-unwise*img-m.fits')
    else:
        image_names = glob.glob(galid+'-unwise*img-m.fits')
    if len(image_names) > 3:
        if verbose:
            print('unwise images already downloaded')
            print(image_names)
        # should be only one *-img-m.fits image per band
        if len(image_names) > len(bands):
            multiframe=True
        else:
            multiframe = False
        weight_names = glob.glob(galid+'-unwise*std*.fits')
        if not multiframe:
            return image_names,weight_names,multiframe
        else:
            if verbose:
                print('going to try new stacking for wise multiframe')
            downloadwise = False
    if downloadwise:
        imsize = int(imsize)
        if verbose:
            print('wise image size = ',imsize)
        baseurl = 'http://unwise.me/cutout_fits?version=allwise'
        imurl = baseurl +'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(ra,dec,imsize,bands)
        if verbose:
            print('downloading unwise images')
            print(imurl)
        wisetar = wget.download(imurl)
        tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
        wnames = tartemp.getnames()

        #print(wnames)
        # check for multiple pointings - means galaxy is split between images
        multiframe = False
        if len(wnames) > 4*len(bands):
            multiframe = True
    
        wmembers = tartemp.getmembers()
        image_names = []
        weight_names = []
        tartemp.extractall()
        for fname in wnames:
            #print(fname)
            t = fname.split('-')
            if subfolder is not None:
                rename = subfolder+'/'+str(galid)+'-'+t[0]+'-'+t[1]+'-'+t[2]+'-'+t[3]+'-'+t[4]
            else:
                rename = str(galid)+'-'+t[0]+'-'+t[1]+'-'+t[2]+'-'+t[3]+'-'+t[4]
            #print('rename = ',rename)
            #print(rename.find('gz'))
            #if os.path.exists(rename): # this should only occur if multiple images are returned from wise
            #    os.remove(rename)
            os.rename(fname, rename)
            if rename.find('.gz') > -1:
                #print('hello????')
                os.system('gunzip '+rename)
                rename = rename.split('.gz')[0]
                #print('after gunzip, rename = ',rename)
            if rename.find('img') > -1:
                image_names.append(rename)
            if rename.find('std') > -1:
                # move ivar images to imagename.weight.fits
                outim = rename.replace('std-m','img-m.std')
                os.rename(rename,outim)
                weight_names.append(outim)
        os.remove(wisetar)
    
    # if multiframe
    # run swarp to create coadded image
    if multiframe:
        image_names=[]
        weight_names=[]
        for b in bands:
            if verbose:
                print('running swarp to combine multiple unwise images in band ',b)
            #########################################
            ## COMBINE THE IMAGE FRAMES USING AVERAGE
            #########################################        
            # create default.swarp
            os.system('swarp -d > default.swarp')
            # run swarp
            matchstring = "*w{}-img-m.fits".format(b)            
            if subfolder is not None:
                allfiles = glob.glob(subfolder+'/'+galid+matchstring)


                
            else:
                allfiles = glob.glob(galid+matchstring)

            all_images = " ".join(allfiles)
            output_image = str(galid)+'-'
            s = 'swarp '+all_images+' -COMBINE_TYPE AVERAGE -WEIGHT_SUFFIX .std.fits -SUBTRACT_BACK N'
            if verbose:
                print(s)
            os.system(s)

            # rename coadd.fits to the output image name
            outimage = str(galid)+'-unwise-w'+str(b)+'-coadd.fits'
            if subfolder is not None:
                os.rename('coadd.fits',os.path.join(subfolder,outimage))
                image_names.append(os.path.join(subfolder,outimage))                
            else:
                os.rename('coadd.fits',outimage)
                image_names.append(outimage)

            
            #os.rename('coadd.fits','unwise/'+outimage)

            #########################################
            ## COMBINE THE STD FRAMES USING SUM
            #########################################
            # just combine using average and then multiply by sqrt 2
            matchstring = "*w{}-img-m.std.fits".format(b)
            
            if subfolder is not None:
                allfiles = glob.glob(subfolder+'/'+galid+matchstring)
            else:
                allfiles = glob.glob(galid+matchstring)
            #print('allfiles with std images = ',allfiles)
            all_images = " ".join(allfiles)
            
            output_image = str(galid)+'-'
            s = 'swarp '+all_images+' -COMBINE_TYPE MEDIAN -WEIGHT_TYPE NONE -SUBTRACT_BACK N'
            #print(s)
            os.system(s)

            # now take sqrt of image values
            #im,h = fits.getdata('coadd.fits',header=True)
            #im = np.sqrt(2)*im
            # divide by number of images because we took average of data
            #im = im/len(all_images)
            
            weightname = str(galid)+'-unwise-w'+str(b)+'-coadd.std.fits'  
            #fits.writeto(weightname,im,header=h,overwrite=True)


            # just trying average!
            if subfolder is not None:
                os.rename('coadd.fits',os.path.join(subfolder,weightname))
                weight_names.append(os.path.join(subfolder,weightname))
            else:
                os.rename('coadd.fits',weightname)
                weight_names.append(weightname)
            
        
    if makeplots:
        ##### DISPLAY IMAGE
        im = fits.getdata(rename)
        norm = simple_norm(im, stretch='asinh',percent=99)
        plt.imshow(im, norm=norm,origin='upper')
        #plt.show()
    #print(image_names)
    #print(multiframe)


    return image_names,weight_names,multiframe

def get_galex_image(ra,dec,imsize):
    """

    get galex image of a galaxy
    
    Input:
    * ra in deg
    * dec in deg
    * imsize in arcsec
    
    Returns:
    * image
    """
    

    # following procedure outlined here:
    # https://astroquery.readthedocs.io/en/latest/mast/mast.html

    # get data products in region near ra,dec
    obs_table = Observations.query_region("%12.8f %12.8f"%(ra,dec),radius=.1*u.arcmin)
    # create a flag to select galex data
    galexFlag = obs_table['obs_collection'] == 'GALEX'

    # separate out galex data
    data_products = Observations.get_product_list(obs_table[galexFlag])

    # download the observations
    manifest = Observations.download_products(data_products,productType="SCIENCE")

    for m in manifest:
        # choose the first NUV image
        if m['Local Path'].find('nd-int') > -1:
            nuv_path = m['Local Path']
            break
        
    # should be able to construct path from the obs_id in data_products
    # this will let us check if the image is already downloaded
    #
    # DONE: I can also save the cutout in a GALEX folder, and
    # look for the image before calling this function
    #
    # but I should also look for the image, because if I change the image
    # size of the cutout, I don't need to download the big FOV
    
    nuv,nuv_header = fits.getdata(nuv_path,header=True)

    # this is a big image, so we need to get a cutout

    nuv_wcs = WCS(nuv_header)
    position = SkyCoord(ra,dec,unit="deg",frame='icrs')
    try:
        cutout = Cutout2D(nuv,position,(imsize*u.arcsec,imsize*u.arcsec),wcs=nuv_wcs)
    except:
        print('WARNING: problem getting galex cutout')
        cutout = None

    return cutout
    

