from astropy.io import fits
import os
import numpy as np
from matplotlib import pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u

class image_pair():
    def __init__(self,image1,image2,prefix=None,default_se_dir = '/Users/rfinn/github/halphagui/astromatic'):
        self.image1 = image1
        self.image2 = image2
        self.sepath = default_se_dir
        base = os.path.basename(self.image1)
        self.froot1 = os.path.splitext(base)[0]
        base = os.path.basename(self.image2)
        self.froot2 = os.path.splitext(base)[0]
        if prefix is not None:
            self.prefix = prefix
        else:
            self.prefix = 'coadd'
    def link_files(self):
        # these are the sextractor files that we need
        # set up symbolic links from sextractor directory to the current working directory
        sextractor_files=['default.sex.HDI','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('ln -s '+self.sepath+'/'+file+' .')

    def run_se(self):
        # get magnitude zeropoint for image 1 and image 2
        header1 = fits.getheader(self.image1)
        header2 = fits.getheader(self.image2)
        try:
            ZP1 = header1['PHOTZP']
            print('got ZP = ',ZP1)
            zp1flag = True
        except KeyError:
            print('no PHOTZP found in image 1 header.  Too bad :(')
            print('did you run getzp.py?')
            zp1flag = False
        try:
            ZP2 = header2['PHOTZP']
            print('got ZP = ',ZP2)
            zp2flag = True
        except KeyError:
            print('no PHOTZP found in image 2 header.  Too bad :(')
            print('did you run getzp.py?')
            zp2flag = False

        print('RUNNING SEXTRACTOR')

        if zp1flag:
            os.system('sex ' + self.image1 + ' -c default.sex.hdi -CATALOG_NAME ' + self.froot1 + '.cat -MAG_ZEROPOINT '+str(ZP1))
        else:
            os.system('sex ' + self.image1 + ' -c default.sex.hdi -CATALOG_NAME ' + self.froot1 + '.cat')
        if zp2flag:
            os.system('sex ' + self.image2 + ' -c default.sex.hdi -CATALOG_NAME ' + self.froot2 + '.cat -MAG_ZEROPOINT '+str(ZP2))
        else:
            os.system('sex ' + self.image2 + ' -c default.sex.hdi -CATALOG_NAME ' + self.froot2 + '.cat')
    def read_catalogs(self):
        self.cat1 = fits.getdata(self.froot1 + '.cat',2)
        self.cat2 = fits.getdata(self.froot2 + '.cat',2)
        
    def match_positions(self):
        '''
        match positions in two catalogs
        '''
        c1 = SkyCoord(ra = self.cat1.ALPHA_J2000*u.degree, dec = self.cat1.DELTA_J2000*u.degree)
        c2 = SkyCoord(ra = self.cat2.ALPHA_J2000*u.degree, dec = self.cat2.DELTA_J2000*u.degree)
        self.idx, self.d2d, self.d3d = c2.match_to_catalog_sky(c1)
        pass
    
    def compare_positions(self):
        # plot positions in two filters to check for alignment errors
        plt.figure()
        #sep = np.sqrt((self.cat1.X_IMAGE[flag]-self.cat2.X_IMAGE[flag])**2+(self.cat1.Y_IMAGE[flag]-self.cat2.Y_IMAGE[flag])**2)
        plt.plot(self.cat1.X_IMAGE,self.cat1.Y_IMAGE,'r.',self.cat2.X_IMAGE,self.cat2.Y_IMAGE,'b.',markersize=2)

        #plt.scatter(self.cat1.X_IMAGE[flag],self.cat1.Y_IMAGE[flag],c=sep,s=10)
    
    def compare_match(self):
        plt.figure(figsize=(6,4))
        self.flag = (self.cat2.FLAGS == 0) & (self.cat2.CLASS_STAR > .98)
        flag = self.flag
        #plt.subplot(1,2,1)
        #plt.plot(self.cat1.X_IMAGE[self.idx][flag],self.cat1.Y_IMAGE[self.idx][flag],'r.',self.cat2.X_IMAGE[flag],self.cat2.Y_IMAGE[flag],'b.',markersize=2)
        #plt.subplot(1,2,2)
        plt.scatter(self.cat1.X_IMAGE[self.idx][flag],self.cat1.Y_IMAGE[self.idx][flag],c=self.d2d.to('arcsec')[flag],s=5,vmin=0.,vmax=.1)
        plt.colorbar(label=r'$\Delta \theta \ b/w \ R \ and \ Halpha \ (arcsec)$')
        plt.xlabel('X_IMAGE (pixels)')
        plt.ylabel('Y_IMAGE (pixels)')
        plt.title(self.prefix+': median offset = %.3f arcsec'%(np.median(self.d2d.to('arcsec')[self.flag]).value))
        plt.savefig(self.prefix+'-offsets.png')

if __name__ == '__main__':

    
    prefix = 'NRGs27'
    im1 = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_R.coadd.fits'
    im2 = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_ha16.coadd.fits'
    #prefix = 'MKW8'
    #im1 = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/MKW8_R.coadd.fits'
    #im2 = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/MKW8_ha16.coadd.fits'
    ip = image_pair(im1,im2,prefix=prefix)
    ip.link_files()
    #ip.run_se()
    ip.read_catalogs()
    ip.match_positions()
    ip.compare_match()
