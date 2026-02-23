#!/usr/bin/env python


'''
GOAL:
- take output from Halpha gui
- for each galaxy, make a 3panel plot with
  - R cutout
  - Halpha cutout
  - superimposed profiles

PROCEDURE:
- search directory for cutouts
- for each galaxy, plot cutouts
- if photometry exists, plot profiles
- mark r24 and r17

OUTPUT:
- pdf files of results


'''

import glob
from astropy.visualization import simple_norm
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import os
from photutils import EllipticalAperture

prefix = 'NRGs27'
photfile = 'halpha-data-rfinn-2019-Sep-17.fits'
photfile = 'v17p3-data-rfinn-2020-Jun-26.fits'

rcutouts = glob.glob(prefix+'*-R.fits')

cat = fits.getdata(photfile)
pixel_scale = 0.43

def mark_ellipse(cat, nsaid):
    pass
    # ellipse showing R24
    index = np.where(cat.NEDname == nsaid)
    #print('start of mark_ellipse',nsaid, index)
    radii = ['GAL_R24','GAL_HR17']
    colors = ['cyan','magenta']
    apertures = []
    for i,sma in enumerate(radii):
        #print('in loop',sma,index)
        rad = cat[sma][index][0]
        #print(sma, rad)
        if rad > 0.1:
            t = cat['GAL_PA'][index][0]
            if t < 0:
                theta = np.radians(180. + t+90)
            else:
                theta = np.radians(t+90) # orientation in radians

            apertures.append(EllipticalAperture((cat['GAL_XC'][index][0],cat['GAL_YC'][index][0]),\
                                            cat[sma][index][0]/pixel_scale, \
                                            cat[sma][index][0]*cat['GAL_BA'][index][0]/pixel_scale,\
                                            theta = theta))
            
    for i,aperture in enumerate(apertures):
        aperture.plot(color=colors[i],lw=1.5)



def mark_radii(cat,nsaid):
    pass
    # vertical line  showing R24
    index = np.where(cat.NSAID == nsaid)
    apertures = ['GAL_R24','GAL_HR17']
    colors = ['cyan','magenta']
    for i,sma in enumerate(apertures):
        plt.axvline(x=cat[sma][index][0], color=colors[i], ls='--', label=sma)
        
    # vertical line showing R17


nperpage = 3
npages = len(rcutouts)/nperpage
ngal = 0
nplot = 0
#plt.figure(figsize=(12,8))
#plt.clf()
#plt.subplots_adjust(hspace=0,wspace=0)
row = 0
haoffset = -4.

rcutouts = glob.glob('*-R.fits')
for rc in rcutouts:
    ngal += 1
    row += 1
    #print(ngal,nperpage)
    if np.mod(ngal,nperpage) == 1:
        fig=plt.figure(figsize=(10,8))
        plt.clf()
        plt.subplots_adjust(hspace=.25,wspace=0.4)
        row = 1
    #print('row = ',row, ' plot # = ',nperpage*(row-1)+1)
    if rc.find('-R.fits') > -1:
        base_filename = rc.split('-R.fits')[0]
    elif rc.find('-r.fits') > -1:
        base_filename = rc.split('-r.fits')[0]
    mask = base_filename+'-R-mask.fits'
    rphot = 'GAL_'+base_filename+'-R_phot.fits'
    himage = base_filename+'-CS.fits'
    hphot = 'GAL_'+base_filename+'-CS_phot.fits'
    nsaid = (base_filename.split('-')[1])
    rdat, rheader = fits.getdata(rc, header=True)
    if os.path.exists(mask):
        mdat = fits.getdata(mask)
        masked_image = np.ma.array(rdat, mask=mdat)
        maskflag=True
    else:
        masked_image = rdat
        maskflag=False
    if os.path.exists(himage):
        haflag = True
        hdat, hheader = fits.getdata(himage, header=True)
        if maskflag:
            masked_himage = np.ma.array(hdat, mask=mdat)
        else:
            masked_himage = hdat
    galindex = np.where(cat.NEDname == nsaid)

    #r24 = cat['GAL_R24'][galindex][0]
    #r24pix = r24/pixel_scale
    #Halpha plus continuum
    ax1=plt.subplot2grid((nperpage,4),(row-1,0))
    norm = simple_norm(masked_image, stretch='asinh')#, percent=99.5)
    ax1.imshow(rdat, norm=norm, origin='lower', cmap='gray_r')
    ax1.set_title(str(nsaid)+': R')
    #mark_ellipse(cat, nsaid)
    #plt.gca().set_yticks(())
    if haflag:
        ax2=plt.subplot2grid((nperpage,4),(row-1,1))
        norm = simple_norm(masked_himage, 'power', percent=99.5)
        ax2.imshow(hdat,norm=norm, origin='lower', cmap='gray_r')
        ax2.set_title('Halpha')
        #plt.gca().set_yticks(())
        #mark_ellipse(cat, nsaid)
    #plt.gca().set_yticks(())
    #plt.xlabel('NSA ID '+str(nsaid),fontsize=14)
    #R
    ax3=plt.subplot2grid((nperpage,4),(row-1,2),colspan=2)
    # plot profiles
    if os.path.exists(rphot):
        rtab = fits.getdata(rphot)
        ax3.errorbar(rtab.sma_arcsec, rtab.sb_mag_sqarcsec, yerr=rtab.sb_mag_sqarcsec_err,fmt='b.', label='R',markersize=6)
        #mark_radii(cat,nsaid)

        if haflag:
            htab = fits.getdata(hphot)
            flag = htab.sb_erg_sqarcsec > 0.
            ax3.errorbar(htab.sma_arcsec[flag], htab.sb_mag_sqarcsec[flag]+haoffset, yerr=htab.sb_mag_sqarcsec_err[flag],fmt='c.', label=r'$H\alpha +$'+str(haoffset),markersize=6)

        #ax3.set_xlim(0,3*r24)
        ax3.set_ylim(16,28)        
        ax3.invert_yaxis()
        ax3.legend()
        ax3.set_ylabel(r'$\rm \mu_r, \ \mu_{H\alpha}+$'+str(haoffset),fontsize=12)

    else:
        print('could not find ',rphot)
    if np.mod(ngal, nperpage) == 0:
        # increment plot number
        nplot += 1
        # save fig
        plt.savefig('profiles-'+str(nplot)+'.png')
