# WARNING - UNDER CONSTRUCTION

def run_photutil(self, snrcut=1.5,npixels=10):
    ''' 
    run photutils detect_sources to find objects in fov.  
    you can specify the snrcut, and only pixels above this value will be counted.

    this also measures the sky noise as the mean of the threshold image
    '''
    self.threshold = detect_threshold(self.image, nsigma=snrcut)
    segment_map = detect_sources(self.image, self.threshold, npixels=npixels)
    # deblind sources a la source extractor
    # tried this, and the deblending is REALLY slow
    # going back to source extractor
    self.segmentation = deblend_sources(self.image, segment_map,
                           npixels=10, nlevels=32, contrast=0.001)        
    self.maskdat = self.segmentation.data
    #self.cat = source_properties(self.image, self.segmentation)
    self.cat = SourceCatalog(self.image, self.segmentation)        
    # get average sky noise per pixel
    # threshold is the sky noise at the snrcut level, so need to divide by this
    self.sky_noise = np.mean(self.threshold)/snrcut
    #self.tbl = self.cat.to_table()

    if self.off_center_flag:
        print('setting center object to objid ',self.galaxy_id)
        self.center_object = self.galaxy_id
    else:
        distance = np.sqrt((self.cat.xcentroid - self.xc)**2 + (self.cat.ycentroid - self.yc)**2)
        # save object ID as the row in table with source that is closest to center
        objIndex = np.arange(len(distance))[(distance == min(distance))][0]
        # the value in shown in the segmentation image is called 'label'
        self.center_object = self.cat.label[objIndex]

    self.maskdat[self.maskdat == self.center_object] = 0
    self.update_mask()


def get_photutils_mask(self,galaxy_id = None):
    # TODO make an alternate function that creates segmentation image from photutils
    from astropy.stats import sigma_clipped_stats
    from photutils import make_source_mask

    # create mask to cut low SNR pixels based on SNR in SFR image
    mask = make_source_mask(imdat,nsigma=self.snr,npixels=self.minarea,dilate_size=5)
    masked_data = np.ma.array(imdat,mask=mask)


    self.catname = self.image_name.replace('.fits','.cat')
    self.segmentation = self.image_name.replace('.fits','-segmentation.fits')


    self.maskdat = fits.getdata(self.segmentation)
    # grow masked areas
    bool_array = np.array(self.maskdat.shape,'bool')
    #for i in range(len(self.xsex)):
    # check to see if the object is not centered in the cutout

