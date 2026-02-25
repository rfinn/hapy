
from astropy.table import Table
from astropy.io import fits

from astropy.table import Table, Column


class OutputTable(output_table_view):
    """
    output table that stores all of the measured values for each galaxy in FOV
    """
    def initialize_results_table(self, prefix=None,virgo=False,uat=False,nogui=False):
        #print("in initialize_results_table, value of uat = ", uat)
        #print("in initialize_results_table, value of virgo = ", virgo)        
        self.nogui = nogui
        
        '''
        Data to store:
        - NSAID
        * AGC number
        - RA
        - DEC
        - filter_ratio
        - cutout_size
        - xmin:xmax,ymin:ymax from parent image
        - ha_flag -- boolean
        - ha_class -- category of Halpha emission
        - psf_fwhm
        - galfit re
        - galfit n
        - galfit BA
        - galfit PA
        - galfit (xc,yc)
        - galfit (RA,DEC) - translate the pixel coords to RA and Dec of galaxy center
        - galfit mag
        - galfit sky
        - ellipse PA
        - ellipse BA
        - ellipse Gini
        - ellipse skynoise
        - ellipse mag R
        - ellipse mag Ha
        - ellipse SFR Ha
        - profiles Re r
        - profiles Re Ha
        - becky inner ssfr
        - becky outer ssfr
        - becky C30
        - becky C70
        '''
        ## define fits table output name
        # get directory after Users - this should be username for 
        user = os.getenv('USER')
        today = date.today()
        str_date_today = today.strftime('%Y-%b-%d')
        if prefix is None:
            self.output_table = 'halpha-data-'+user+'-'+str_date_today+'.fits'
        else:
            self.output_table = prefix+'-data-'+user+'-'+str_date_today+'.fits'
        ## check for existing table
        ##
        ## load if it exists

        # why is this block at the beginning???
        #if virgo:
        #    self.create_table_virgo()
        #elif uat:
        #    print("just to check, I am running create_table_uat")
        #    self.create_table_uat()
        #else:
        #    self.create_table()


        
        if os.path.exists(self.output_table):
            if virgo:
                self.read_table_virgo()
            elif uat:
                self.read_table_uat()
            else:
                self.read_table()
                self.agc2 = fits.getdata(self.prefix+'-agc-matched.fits')
        ## if not, create table                
        else:
            if virgo:
                self.create_table_virgo()
            elif uat:
                self.create_table_uat()
            else:
                self.create_table()
            # call other methods to add columns to the table
            self.add_part1()
            # skipping for now b/c this will have to be different for virgo
            #self.add_nsa()
            self.add_flags()            
            self.add_cutout_info()
            self.add_galfit_r()
            #self.add_galfit_ha()            
            self.add_ellipse()
            self.add_profile_fit()
            self.add_photutils()
            self.add_statmorph()            
        if not self.nogui:
            self.update_gui_table()
                
    def read_table(self):
        ''' read in output from previous run, if it exists'''
        self.table = Table(fits.getdata(self.output_table))
        self.gredshift = self.table['REDSHIFT']
        self.ngalaxies = len(self.gredshift)
        self.ra = self.table['NSA_RA']*self.table['NSA_FLAG']+ self.table['AGC_RA']*(~self.table['NSA_FLAG'])
        self.dec = self.table['NSA_DEC']*self.table['NSA_FLAG']+ self.table['AGC_DEC']*(~self.table['NSA_FLAG'])
        self.gradius = self.table['SERSIC_TH50']*self.table['NSA_FLAG']/self.pixelscale + 100.*np.ones(self.ngalaxies)*(~self.table['NSA_FLAG'])
        self.gzdist = self.table['ZDIST']
        charar1 = np.chararray(self.ngalaxies)
        charar1[:] = 'N'
        charar2 = np.chararray(self.ngalaxies)
        charar2[:] = '-A'
        self.galid=np.zeros(self.ngalaxies, dtype='U15')
        for i in np.arange(self.ngalaxies):
            self.galid[i] = 'N'+str(self.table['NSAID'][i])+'-A'+str(self.table['AGCNUMBER'][i])
        # read in nsa2
        self.nsa2 = fits.getdata(self.prefix+'-nsa-matched.fits')
        # read in agc2
        self.agc2 = fits.getdata(self.prefix+'-agc-matched.fits')                                                                              
        ## if not, create table
    def read_table_virgo(self):
        self.table = Table(fits.getdata(self.output_table))
        self.gredshift = self.table['REDSHIFT']
        self.ngalaxies = len(self.gredshift)
        self.ra = self.table['RA']
        self.dec = self.table['DEC']
        self.gradius = self.table['radius']
        self.gzdist = self.table['ZDIST']
        self.galid=self.table['VFID']
        self.NEDname=self.table['NEDname']
        self.gprefix=self.table['prefix']                

    def read_table_uat(self):
        self.table = Table(fits.getdata(self.output_table))
        self.gredshift = self.table['vopt']/3.e5
        self.ngalaxies = len(self.gredshift)
        self.ra = self.table['RA']
        self.dec = self.table['DEC']
        self.gradius = self.table['a']
        self.gzdist = self.table['vopt']/3.e5
        self.galid= self.table['AGCnr']
        self.NEDname= self.table['AGCnr']
        self.haflag= self.table['AGCnr']        
        try:
            self.gprefix= self.table['prefix']
        except:
            print("WARNING: no prefix in gal table")
            self.gprefix = None
        
    def create_table_virgo(self):
        # updating this part for virgo filament survey 

        self.table = self.defcat.cat['VFID','RA','DEC','vr','radius','NEDname','prefix']
        self.table['VFID'].description = 'ID from Virgo Filament catalog'                
        self.table['RA'].unit = u.deg
        self.table['RA'].description = 'RA from VF catalog'        
        self.table['DEC'].unit = u.deg
        self.table['DEC'].description = 'DEC from VF catalog'                
        self.table['vr'].unit = u.km/u.s
        self.table['vr'].description = 'recession velocity from VF catalog'
        self.table['radius'].unit = u.arcsec
        self.table['radius'].description = 'radius from VF catalog'        
        self.ngalaxies = len(self.table)
        #print('number of galaxies = ',self.ngalaxies)
        self.haflag = np.zeros(self.ngalaxies,'bool')
        self.galid = self.table['VFID']
        self.NEDname = self.table['NEDname']                
        self.gredshift = self.defcat.cat['vr']/3.e5
        self.gzdist = self.defcat.cat['vr']/3.e5

        ##
        # update this to use the SMA_SB24
        ##
        #self.gradius = self.defcat.cat['radius']/self.pixelscale
        self.gradius = self.radius_arcsec/self.pixelscale

        
        self.ra = self.defcat.cat['RA']
        self.dec = self.defcat.cat['DEC']        
        c1 = Column(self.haflag, name='HAflag', description='Halpha flag')
        c2 = Column(self.gredshift, name='REDSHIFT', description='redshift')
        c3 = Column(self.gredshift, name='ZDIST', description='redshift')        
        self.table.add_columns([c1,c2,c3])

    def create_table_uat(self):
        # updating this part for uat Halpha Groups

        self.table = self.defcat.cat['AGCnr','RA','DEC','vopt','v21','a','b','hiflux',]
        self.table['AGCnr'].description = 'Number in AGC catalog'
        self.table['RA'].unit = u.deg
        self.table['RA'].description = 'RA from AGC catalog'        
        self.table['DEC'].unit = u.deg
        self.table['DEC'].description = 'DEC from AGC catalog'                
        self.table['vopt'].unit = u.km/u.s
        self.table['vopt'].description = 'optical recession velocity from AGC catalog'
        self.table['v21'].unit = u.km/u.s
        self.table['v21'].description = 'HI recession velocity from AGC catalog'
        self.table['a'].unit = u.arcmin
        self.table['a'].description = 'AGC semi-major axis'
        self.table['b'].unit = u.arcmin
        self.table['b'].description = 'AGC semi-minor axis'
        self.ngalaxies = len(self.table)
        #print('number of galaxies = ',self.ngalaxies)
        self.haflag = np.zeros(self.ngalaxies,'bool')
        self.galid = self.table['AGCnr']
        #self.galid = ['AGC'+str{i} for i in self.galid] # prepend 'AGC' to the AGCnr
        self.NEDname = None #self.table['NEDname']                
        self.gredshift = self.defcat.cat['vopt']/3.e5
        self.gzdist = self.defcat.cat['vopt']/3.e5

        ##
        # update this to use the SMA_SB24
        ##
        #self.gradius = self.defcat.cat['radius']/self.pixelscale
        self.gradius = self.radius_arcsec/self.pixelscale

        
        self.ra = self.defcat.cat['RA']
        self.dec = self.defcat.cat['DEC']        
        c1 = Column(self.haflag, name='HAflag', description='Halpha flag')
        c2 = Column(self.gredshift, name='REDSHIFT', description='redshift')
        c3 = Column(self.gredshift, name='ZDIST', description='redshift')        
        self.table.add_columns([c1,c2,c3])
        
    def create_table(self):
        # updating this part to make use of NSA and AGC catalogs
        # not going to make this backward compatible, meaning you need to enter both catalogs
        # probably just being lazy, but it's giving me a headache...

        # much better approach would be to match entire NSA and AGC catalogs
        # and then just use that,
        # but forging ahead for now.



        # this returns new nsa and agc catalogs, that are row matched to the joined table
        if self.agcflag:
            self.nsa2, self.nsa_matchflag, self.agc2, self.agc_matchflag = make_new_cats(self.nsa.cat, self.agc.cat)
            # write matched tables
            self.nsa2.write(self.prefix+'-nsa-matched.fits', format='fits', overwrite=True)
            self.agc2.write(self.prefix+'-agc-matched.fits', format='fits', overwrite=True)

            self.ngalaxies = len(self.nsa_matchflag)        
            # create arrays that we need for other parts of the programs, like ra, dec, size
            self.ra = self.nsa2['RA']*self.nsa_matchflag + self.agc2['RA']*(~self.nsa_matchflag)

            self.dec = self.nsa2['DEC']*self.nsa_matchflag + self.agc2['DEC']*(~self.nsa_matchflag)
 
            self.gradius = self.nsa2['SERSIC_TH50']*self.nsa_matchflag/self.pixelscale + 100.*np.ones(self.ngalaxies)*(~self.nsa_matchflag)
            charar1 = np.chararray(self.ngalaxies)
            charar1[:] = 'N'
            charar2 = np.chararray(self.ngalaxies)
            charar2[:] = '-A'
            self.galid=np.zeros(self.ngalaxies, dtype='U15')
            for i in np.arange(self.ngalaxies):
                self.galid[i] = 'N'+str(self.nsa2['NSAID'][i])+'-A'+str(self.agc2['AGCnr'][i])
            voptflag = self.agc2['vopt'] > 0.
            agcredshift = self.agc2['vopt']/3.e5*voptflag + self.agc2['v21']/3.e5*(~voptflag)
            self.gredshift = self.nsa2['Z']*self.nsa_matchflag + agcredshift*(~self.nsa_matchflag)
        else:
            self.nsa2 = self.nsa.cat
            self.nsa_matchflag = np.ones(len(self.nsa2))
            self.ra = self.nsa2['RA']

            self.dec = self.nsa2['DEC']

            self.gradius = self.nsa2['SERSIC_TH50']*self.nsa_matchflag/self.pixelscale 
            self.ngalaxies = len(self.nsa_matchflag)
            self.galid = self.nsa2['NSAID']*self.nsa_matchflag
            self.gredshift = self.nsa2['Z']*self.nsa_matchflag
            
        # number of galaxies in the joined table
        self.haflag = np.zeros(self.ngalaxies,'bool')
        c0 = Column(self.galid,name='ID')
        c1 = Column(self.gredshift,name='REDSHIFT')
        c2 = Column(self.nsa2['NSAID'], name='NSAID',dtype=np.int32, description='NSAID')
        c3 = Column(self.nsa_matchflag, name='NSA_FLAG',dtype='bool', description='NSA_FLAG')
        c4 = Column(self.nsa2['RA'], name='NSA_RA',dtype='f', unit=u.deg)
        c5 = Column(self.nsa2['DEC'], name='NSA_DEC',dtype='f', unit=u.deg)
        self.table = Table([c0,c1,c2,c3,c4,c5])
        if self.agcflag:
            c5 = Column(self.agc2['AGCnr'], name='AGCNUMBER',dtype=np.int32, description='AGC ID NUMBER')
            c6 = Column(self.agc_matchflag, name='AGC_FLAG',dtype='bool', description='AGC_FLAG')
            c7 = Column(self.agc2['RA'], name='AGC_RA',dtype='f', unit=u.deg)
            c8 = Column(self.agc2['DEC'], name='AGC_DEC',dtype='f', unit=u.deg)
            self.table.add_columns([c5,c6,c7,c8])
        
    def add_part1(self):
        g1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_RA', unit=u.deg,description='R-band center RA from galfit')
        g2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_DEC', unit=u.deg,description='R-band center DEC from galfit')
        g3 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HRA', unit=u.deg,description='HA center RA from galfit')
        g4 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HDEC', unit=u.deg,description='HA center DEC from galfit')
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_RA', unit=u.deg,description='R-band center RA from photutil centroid')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_DEC', unit=u.deg,description='R-band center DEC from photutil centroid')
        e1a = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HRA', unit=u.deg,description='Halpha center RA from photutil centroid')
        e2a = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HDEC', unit=u.deg,description='Halpha center DEC from photutil centroid')
        c9 = Column(self.haflag, name='HA_FLAG',description='shows HA emission')
        c10 = Column(np.ones(self.ngalaxies,'f'),name='FILT_COR',unit='', description='max filt trans/trans at gal z')
        c11 = Column(np.zeros(self.ngalaxies,'f'),name='R_FWHM',unit=u.arcsec, description='R FWHM in arcsec')
        c12 = Column(np.zeros(self.ngalaxies,'f'),name='H_FWHM',unit=u.arcsec, description='HA FWHM in arcsec')
        c13 = Column(np.zeros(self.ngalaxies,dtype='|S40'),name='POINTING', description='string specifying year and pointing')
        c14 = Column(np.zeros(self.ngalaxies,dtype='|S3'),name='TEL', description='telescope/instrument')
        c15 = Column(np.zeros(self.ngalaxies,dtype='i'),name='DATE-OBS', description='string specifying date of observation')                        

        self.table.add_columns([g1,g2,g3,g4,e1,e2,e1a,e2a,c9,c10,c11,c12,c13,c14,c15])
    def add_nsa(self):
        # add some useful info from NSA catalog (although matching to NSA could be done down the line)
        r = 22.5 - 2.5*np.log10(self.nsa2['NMGY'][:,4])
        c11 = Column(r,name='NSA_RMAG',unit=u.mag,description='NSA r mag')
        c12 = Column(self.nsa2['SERSIC_TH50'],name='SERSIC_TH50', unit=u.arcsec,description='NSA SERSIC_TH50')
        c13 = Column(self.nsa2['SERSIC_N'],name='SERSIC_N',description='NSA SERSIC index')
        c14 = Column(self.nsa2['SERSIC_BA'],name='SERSIC_BA',description='NSA SERSIC B/A')
        c15 = Column(self.nsa2['SERSIC_PHI'],name='SERSIC_PHI', unit=u.deg,description='NSA SERSIC PHI')
        self.gzdist = self.nsa2['ZDIST']*self.nsa_matchflag + self.gredshift*~self.nsa_matchflag
        c16 = Column(self.gzdist,name='ZDIST',description='NSA ZDIST')
        self.table.add_columns([c11,c12,c13,c14,c15,c16,])
    def add_cutout_info(self):
        # cutout region in coadded images
        c1 = Column(np.zeros(len(self.table),dtype='U22'), name='BBOX',description='location of galaxy cutout in mosaic')
        # R-band scale factor for making continuum-subtracted image
        c2 = Column(np.zeros(len(self.table),'f'), name='FILTER_RATIO',description='R/Ha ratio used in cont subtraction')
        # r-band ZP
        c3 = Column(np.zeros(len(self.table),'f'), name='RZP',description='R-band ZP')
        # Halpha ZP        
        c4 = Column(np.zeros(len(self.table),'f'), name='HZP',description='Halpha ZP')
        c5 = Column(np.zeros(len(self.table),'f'), name='PIXSCALE',description='Pixel scale')                

        self.table.add_columns([c1,c2,c3,c4,c5])
    def add_galfit_r(self):
        ##############################################3
        ### GALFIT R-BAND FITS
        ##############################################3

        fields = ['XC','YC','MAG','RE','N','BA','PA']
        units = ['pixel','pixel','mag','arcsec',None,'deg',None]
        descriptions = ['R-band center from galfit (pix)',\
                        'R-band center from galfit (pix)',\
                        'R-band mag from galfit',\
                        # this is currently written in pixels - need to write out in arcsec
                        #'R-band effective radius from galfit (pix)',\
                        'R-band effective radius from galfit (arcsec)',\
                        'R-band sersic index from galfit',\
                        'R-band axis ratio from galfit',\
                        'R-band position angle from galfit']
        i=0
        for f,unit in zip(fields,units):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR',description='err in '+descriptions[i])
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f, unit=unit,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR', unit=unit,description='err in '+descriptions[i])
            #print(c1)
            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_SKY',unit=u.adu,description='sky from galfit')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_CHISQ',description='chisq of galfit sersic model')
        #c3 = Column(np.zeros(self.ngalaxies,'f'), name='GAL_GINI')
        #c4 = Column(np.zeros(self.ngalaxies), name='GAL_GINI2')
        #c5 = Column(np.zeros(self.ngalaxies,'f'), name='GAL_ASYM')
        #c6 = Column(np.zeros(self.ngalaxies,'f'), name='GAL_ASYM2')
        self.table.add_columns([c1,c2])#,c3,c4,c5,c6])

    def add_galfit_2comp_r(self):
        # galfit sersic parameters from 2 comp fit
        c16 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_2SERSIC',description='galfit R-band 2comp fit')
        c17 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_2SERSIC_ERR',description='galfit R-band 2comp fit errors')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_2SERSIC_ERROR',description='galfit R-band 2comp fit num err flag')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_2SERSIC_CHISQ',description='galfit R-band 2comp chi sq')
        self.table.add_columns([c16,c17,c18,c19])

    def add_galfit_1comp_with_asymmetry_r(self):

        # galfit 1 comp with asymmetry
        c16 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_SERSASYM',description='galfit R-band 1comp sersic w/asymmetry')
        c17 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_SERSASYM_ERR')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_ERROR',description='galfit R-band 1comp sersic w/asymmetry num err flag')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_CHISQ',description='galfit R-band 1comp sersic w/asymmetry chi sq')
        c20 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_RA',unit='deg',description='RA from galfit R-band 1comp sersic w/asymmetry')
        c21 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_DEC',unit='deg',description='DEC from galfit R-band 1comp sersic w/asymmetry')
        self.table.add_columns([c16,c17,c18,c19,c20,c21])
    def add_galfit_ha(self):
        ##############################################
        ### GALFIT Halpha FITS
        ##############################################

        fields = ['XC','YC','MAG','RE','N','BA','PA']
        units = ['pixel','pixel','mag','arcsec',None,'deg',None]
        descriptions = ['HA center from galfit (pix)',\
                        'HA center from galfit (pix)',\
                        'HA mag from galfit',\
                        # currently written in pixels
                        # need to convert to arcsec
                        'HA effective radius from galfit (arcsec)',\
                        'HA sersic index from galfit',\
                        'HA axis ratio from galfit',\
                        'HA position angle from galfit']
        i=0
        for f,unit in zip(fields,units):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f+'_ERR',description='err in '+descriptions[i])
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f, unit=unit,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f+'_ERR', unit=unit,description='err in '+descriptions[i])

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HSKY',description='galfit HA sky')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HCHISQ',description='galfit chisq of HA model')
        self.table.add_columns([c1,c2])#,c3,c4,c5,c6])

    def add_galfit_2comp_ha(self):
        # galfit sersic parameters from 2 comp fit
        c16 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_H2SERSIC',description='galfit HA 2-comp fit')
        c17 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_H2SERSIC_ERR')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_H2SERSIC_ERROR',description='galfit HA 2-comp num error code')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_H2SERSIC_CHISQ',description='galfit HA 2-comp chisq')
        self.table.add_columns([c16,c17,c18,c19])

    def add_galfit_1comp_with_asymmetry_ha(self):
        # galfit 1 comp with asymmetry
        c16 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_HSERSASYM',description='galfit HA model w/asym')
        c17 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_HSERSASYM_ERR')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_ERROR',description='galfit HA asym num error code')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_CHISQ',description='galfit HA asym chisq')
        c20 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_RA',unit='deg')
        c21 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_DEC',unit='deg')
        self.table.add_columns([c16,c17,c18,c19,c20,c21])

    def add_ellipse(self):
        #####################################################################
        # ellipse output
        # xcentroid, ycentroid, eps, theta, gini, sky_centroid, area, background_mean, source_sum, source_sum_err
        #####################################################################
        e0 = Column(np.zeros(self.ngalaxies,'bool'), name='BADGAL',description='bad galaxy flag - maybe partial coverage')
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_XCENTROID', unit='pixel',description='xcentroid from ellipse')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_YCENTROID', unit='pixel',description='ycentroid from ellipse')
        e3 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_EPS',description='axis ratio from ellipse')
        e4 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_THETA', unit=u.degree,description='position angle from ellipse')
        e5 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_GINI',description='gini coeff from ellipse')
        e6 = Column(np.zeros(self.ngalaxies), name='ELLIP_HGINI',description='gini coeff method 2')
        e7 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_M20',description='M20 for r image')
        e8 = Column(np.zeros(self.ngalaxies), name='ELLIP_HM20',description='M20 for Halpha image ')
        e9 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_UNMASKED_AREA',description='unmasked source area from photutils')
        e9b = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_TOTAL_AREA',description='total source area from photutils')
        e10 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM', unit = u.erg/u.s/u.cm**2,description='total flux from ellipse')
        e11 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM_MAG', unit = u.mag,description='mag from ellipse')
        e12 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM',description='asym from ellipse')
        e13 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM_ERR')
        e14 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM', unit=u.erg/u.s/u.cm**2,description='HA flux from ellipse')
        e15 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM_MAG', unit=u.mag,description='HA mag from ellipse')
        e16 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM',description='HA asymmetry from ellipse')
        e17 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM_ERR')
        e18 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKYNOISE',description='R skynoise in 1E-17 erg/s/cm^2/arcsec^2')
        e19 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKYNOISE',description='HA skynoise in 1E-17  erg/s/cm^2/arcsec^2')
        e20 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKY',description='R sky level in ADU')
        e21 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKY',description='HA sky level in ADU')

        # photutils radii
        e22 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R30',description='photutils R flux frac 30')
        e23 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R50',description='photutils R flux frac 50')
        e24 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R90',description='photutils R flux frac 90')
        
        e25 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR30',description='photutils Halpha flux frac 30')
        e26 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR50',description='photutils Halpha flux frac 50')
        e27 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR90',description='photutils Halpha flux frac 90')        

        
        self.table.add_columns([e0,e1,e2,e3,e4,e5,e6,e7,e8, e9, e9b,e10, e11, e12, e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24,e25,e26,e27])
    def add_profile_fit(self):
        #####################################################################
        # profile fitting using galfit geometry
        #####################################################################
        #
        # r-band parameters
        # 
        self.fields_r = ['R24','R25','R26','R_F25','R24V','R25V','R_F50','R_F75','M24','M25','M26', 'F_30R24','F_R24','C30',\
                    'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        self.units_r = [u.arcsec,u.arcsec,u.arcsec,u.arcsec,u.arcsec,\
                   u.arcsec,u.arcsec,u.arcsec,\
                   u.mag, u.mag, u.mag, \
                   u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,'',\
                   u.arcsec,u.erg/u.s/u.cm**2,u.arcsec, u.arcsec,'',u.mag
                   ]
        self.descriptions= ['isophotal radius at 24mag/sqarc AB',\
                            'isophotal radius at 25mag/sqarc AB',\
                            'isophotal radius at 26mag/sqarc AB',\
                            'radius that encloses 25% of total flux',\
                            'isophotal radius at 24mag/sqarc Vega',\
                            'isophotal radius at 24mag/sqarc Vega',\
                            'radius that encloses 50% of total flux',\
                            'radius that encloses 75% of total flux',\
                            'isophotal mag within R24',\
                            'isophotal mag within R25',\
                            'isophotal mag within R26',\
                            'flux within 30% of R24',\
                            'flux within R24',\
                            'C30 = flux w/in 0.3 r24 / flux w/in r24',\
                            'petrosian radius: where sb is 0.2 times mean sb',\
                            'flux enclosed within 2xpetro radius',\
                            'radius enclosing 50% of petrosian flux',\
                            'radius enclosing 90% of petrosian flux',\
                            '90% petro radius / 50% petro radius',\
                            'magnitude of petrosian flux']
        i=0
        for f,unit in zip(self.fields_r,self.units_r):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f,description='galfit '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f, unit=unit,description='galfit '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        #
        # Halpha parameters
        # 
        self.fields_ha = ['R16','R17',\
                  'R_F25','R_F50','R_F75',\
                  'M16','M17', \
                  'F_30R24','F_R24','C30',\
                  'R_F95R24','F_TOT',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG'
                  ]
        self.units_ha = [u.arcsec,u.arcsec,\
                 u.arcsec,u.arcsec, u.arcsec, \
                 u.mag, u.mag, \
                 u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2, '',\
                 u.arcsec,u.erg/u.s/u.cm**2,\
                 u.arcsec,u.erg/u.s/u.cm**2,u.arcsec, u.arcsec,'',u.mag]
        self.descriptions_ha= ['HA isophotal radius at 16erg/s/cm^2',\
                            'HA isophotal radius at 17erg/s/cm^2',\
                            'HA radius that encloses 25% of total flux',\
                            'HA radius that encloses 50% of total flux',\
                            'HA radius that encloses 75% of total flux',\
                            'HA isophotal radius at 16erg/s/cm^s',\
                            'HA isophotal radius at 17erg/s/cm^2',\
                            'HA flux within 30% of R-band R24',\
                            'HA flux within R-band R24',\
                            'HA C30 = flux w/in 0.3 R-band r24 / flux w/in R-band r24',\
                            'HA flux within 30% of R-band R24',\
                            'HA total flux',\
                            'petrosian radius: where sb is 0.2 times mean sb',\
                            'flux enclosed within 2xpetro radius',\
                            'radius enclosing 50% of petrosian flux',\
                            'radius enclosing 90% of petrosian flux',\
                            '90% petro radius / 50% petro radius',\
                            'magnitude of petrosian flux']
        i=0
        for f,unit in zip(self.fields_ha,self.units_ha):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f,description='galfit '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f, unit=unit,description='galfit '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1           
        f='GAL_'+'LOG_SFR_HA'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=u.M_sun/u.yr,description='log10 of HA SFR in Msun/yr')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR',unit=u.M_sun/u.yr,description='error in log10 of HA SFR in Msun/yr')
        c3 = Column(np.zeros(self.ngalaxies,'bool'),name=f+'_FLAG')        
        self.table.add_column(c1)
        self.table.add_column(c2)
        self.table.add_column(c3)        
        
        f='GAL_'+'SSFR_IN'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(r) within 0.3 R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_column(c1)
        self.table.add_column(c2)
        f='GAL_'+'SSFR_OUT'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(r) within 0.3 R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_column(c1)
        self.table.add_column(c2)
    def add_photutils(self):
        #####################################################################
        # profile fitting using photutils geometry
        #####################################################################
        #
        # r-band parameters
        #
        i=0
        for f,unit in zip(self.fields_r,self.units_r):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='ellipse '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=unit,description='ellipse '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        #
        # Halpha parameters
        #
        i=0
        for f,unit in zip(self.fields_ha,self.units_ha):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f,description='ellipse '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f, unit=unit,description='ellipse '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        f='LOG_SFR_HA'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=u.M_sun/u.yr,description='log10 of HA SFR in Msun/yr')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR',unit=u.M_sun/u.yr)
        c3 = Column(np.zeros(self.ngalaxies,'bool'),name=f+'_FLAG')
        print('testing: colname = ',f+'_FLAG')
        self.table.add_columns([c1,c2,c3])

        ######################################################################
        ### LAST TWO QUANTITIES, I SWEAR!
        ######################################################################        
        
        f='SSFR_IN'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(r) within 0.3 R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_columns([c1,c2])

        f='SSFR_OUT'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(R) outside 0.3R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_columns([c1,c2])



        
        self.table.add_column(Column(np.zeros(self.ngalaxies,dtype='U50'), name='COMMENT'))
        #print(self.table)
        
    def add_statmorph(self):
        #####################################################################
        # statmorph output
        #####################################################################

        # rband area
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_XCENTROID', unit='pixel',description='xcentroid from ellipse')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_YCENTROID', unit='pixel',description='ycentroid from ellipse')
        e3 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RPETRO_CIRC', unit='arcsec',description='rpetro circ')
        e4 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RPETRO_ELLIP', unit='arcsec',description='rpetro ellip')
        e5 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RHALF_ELLIP', unit='arcsec',description='rhalf ellip')
        e6 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_R20', unit='arcsec',description='R20')
        e7 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_R80', unit='arcsec',description='R80')
        e8 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_GINI',description='statmorph gini')
        e8b = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_M20',description='statmorph M20')        
        e9 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_F_GM20',description='statmorph F(G,M20)')
        e10 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_S_GM20',description='statmorph S(G,M20)')
        e11 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_C',description='statmorph concentration')
        e12 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_A',description='statmorph asymmetry')
        e13 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_S',description='statmorph smoothness')
        e14 = Column(np.zeros(self.ngalaxies,'bool'), name='SMORPH_FLAG',description='statmorph flag')                 

        ## Halpha parameters
        h1 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HXCENTROID', unit='pixel',description='xcentroid from ellipse')
        h2 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HYCENTROID', unit='pixel',description='ycentroid from ellipse')
        h3 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRPETRO_CIRC', unit='arcsec',description='rpetro circ')
        h4 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRPETRO_ELLIP', unit='arcsec',description='rpetro ellip')
        h5 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRHALF_ELLIP', unit='arcsec',description='rhalf ellip')
        h6 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HR20', unit='arcsec',description='R20')
        h7 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HR80', unit='arcsec',description='R80')
        h8 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HGINI',description='statmorph gini')
        h8b = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HM20',description='statmorph M20')        
        h9 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HF_GM20',description='statmorph F(G,M20)')
        h10 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HS_GM20',description='statmorph S(G,M20)')
        h11 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HC',description='statmorph concentration')
        h12 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HA',description='statmorph asymmetry')
        h13 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HS',description='statmorph smoothness')
        h14 = Column(np.zeros(self.ngalaxies,'bool'), name='SMORPH_HFLAG',description='statmorph flag')                 
        
        
        self.table.add_columns([e1,e2,e3,e4,e5,e6,e7,e8,e8b, e9, e10, e11, e12, e13,e14,\
                                h1,h2,h3,h4,h5,h6,h7,h8,h8b, h9, h10, h11, h12, h13,h14])

    def add_flags(self):
        '''
        these are common comments that the user will be able to select
        '''
        names = ['CONTSUB_FLAG','MERGER_FLAG','SCATLIGHT_FLAG','ASYMR_FLAG','ASYMHA_FLAG','OVERSTAR_FLAG','OVERGAL_FLAG','PARTIAL_FLAG','EDGEON_FLAG','NUC_HA']
        descriptions =  ['Halpha Emission','Cont Sub Prob','merger/tidal','scattered light','asymmetric R-band', 'asymmetric Ha','foreground star', 'foreground gal','galaxy is edge-on','galaxy is only partially covered by mosaic','nuclear ha emission'] 
        for i,n in enumerate(names):
            #print(n)
            c = Column(np.zeros(self.ngalaxies,'bool'),name=n,description=descriptions[i])
            self.table.add_column(c)
            
    def append_column(self, variable, var_name, var_dtype=None, var_unit=None):
        if (var_dtype != None) & (var_unit != None):
            z = Column(variable, name=var_name, dtype = var_dype, unit  = var_unit)
        elif (var_dtype != None) & (var_unit == None):
            z = Column(variable, name=var_name, dtype = var_dype)
        elif (var_dtype == None) & (var_unit != None):
            z = Column(variable, name=var_name, unit  = var_unit)
        else:
            z = Column(variable, name=var_name)
        self.table.add_column(z)
    def write_fits_table(self):
        if (self.igal is not None) & (not self.auto):
            #print(self.ui.commentLineEdit.text())
            t = str(self.ui.commentLineEdit.text())
            if len(t) > 1:
                self.table['COMMENT'][self.igal] = t
                # don't call update_gui_table_cell here - keep the fits and gui table calls separate
                #self.update_gui_table_cell(self.igal, 'COMMENT',t)
        #fits.writeto('halpha-data-'+user+'-'+str_date_today+'.fits',self.table, overwrite=True)
        if self.prefix is not None:
            # this is not working when running gui - need to feed in the r-band image name
            try:
                telescope,dateobs,p = get_params_from_name(self.prefix)
            except UnboundLocalError:
                if self.uat:
                    telescope,dateobs,p = get_params_from_name_uat(self.rcoadd_fname)
                else:
                    telescope,dateobs,p = get_params_from_name(self.rcoadd_fname)
                #print(f"telescope={telescope}, dateobs={dateobs}, p={p}")
            for i in range(len(self.table)):
                self.table['POINTING'][i] = self.prefix
                self.table['TEL'][i] = telescope
                self.table['DATE-OBS'] = dateobs
        self.table.write(self.output_table, format='fits', overwrite=True)

    def set_galfit_r_row(self, i: int, res):
        """Write a GalfitResult into the per-galaxy row i (R-band)."""
        self.table["GAL_XC"][i] = res.comp1.xc
        self.table["GAL_XC_ERR"][i] = res.comp1.xc_err
        self.table["GAL_YC"][i] = res.comp1.yc
        self.table["GAL_YC_ERR"][i] = res.comp1.yc_err
        self.table["GAL_MAG"][i] = res.comp1.mag
        self.table["GAL_MAG_ERR"][i] = res.comp1.mag_err
        self.table["GAL_RE"][i] = res.comp1.re
        self.table["GAL_RE_ERR"][i] = res.comp1.re_err
        self.table["GAL_N"][i] = res.comp1.n
        self.table["GAL_N_ERR"][i] = res.comp1.n_err
        self.table["GAL_BA"][i] = res.comp1.ba
        self.table["GAL_BA_ERR"][i] = res.comp1.ba_err
        self.table["GAL_PA"][i] = res.comp1.pa
        self.table["GAL_PA_ERR"][i] = res.comp1.pa_err

        self.table["GAL_SKY"][i] = res.sky
        # If you want SKY_ERR, add a column and store res.sky_err
        self.table["GAL_CHISQ"][i] = res.chi2nu
