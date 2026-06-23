# HAPY Run History

Started this a bit late. Sorry about that!

## hapy-output-20260517 
Input coadds:
- coadds-v20260330

Major changes:
- 

Known issues:
- 

Reason for rerun:
- 


## hapy-output-20260517-pre2025coadds
Input coadds:
- coadds-pre2025-hapy

Major changes:
- 

Known issues:
- swarp alignment of coadds was not done correctly

Reason for rerun:
- to compare hapy results from new coadds with hapy results from prior reduction of the coadds (before swarp fixes)


## hapy-output-20260519

Date:
2026-05-19

Input coadds:
- coadds-v20260518

Major changes:
- First hybrid INT dataset
- Old INT Halpha + new INT r-band
- Manual masks transferred
- CS-gr images generated
- Simple INT weight maps

Known issues:
- INT r-band reprojection not aligned with Halpha
- ~24 Legacy cutouts missing
- Image2 aperture noise bug still present

Status:
SUPERSEDED

## hapy-output-20260518-hybrid

Date:
2026-05-18

Input coadds:
- coadds-v20260518

Major changes:
- The new reduction seems to introduce larger photometric errors in the INT halpha images.  
- This run uses the new reduction for BOK, HDI, MOS, and the INT r-band images, but the pre-2025 reduction for the INT halpha images.

Pipeline changes:
- Fixed image2 aperture-noise calculation
- FILTER_RATIO derived from PHOTZP
- CSGR schema initialization fixed

Reason for rerun:
The new reduction seems to introduce larger photometric errors in the INT halpha images.  

## hapy-output-20260519-hybrid

Date:
2026-05-19

Input coadds:
- coadds-v20260518

Major changes:
- not sure

Pipeline changes:
- 

Reason for rerun:


## hapy-output-20260609-hybrid

Date:
2026-06-09

Input coadds:
- coadds-v20260609

Major changes:
- Rebuilt INT hybrid coadds
- SCAMP+SWarp astrometric reprojection
- Halpha and r-band resampled to common grid
- Filter ratios recomputed

Pipeline changes:
- Fixed image2 aperture-noise calculation
- FILTER_RATIO derived from PHOTZP
- CSGR schema initialization fixed

Reason for rerun:
The first hybrid INT reduction (coadds-v20260518)  produced misaligned Hα and r-band images because the SWarp-only reprojection did not correctly handle the INT astrometric solution. The coadds were rebuilt using the SCAMP+SWarp workflow.


Status:
more recent hybrid


## hapy-output-20260612

Date:
2026-06-12

Input coadds:
- coadds-v20260330

Major changes:
- going back to all new coadds

Pipeline changes:
- Fixed image2 aperture-noise calculation
- FILTER_RATIO derived from PHOTZP
- CSGR schema initialization fixed
- additional parameters are propagated into results.ecsv tables

Reason for rerun:
The hybrid INT images are still not aligned well, and I am going to table the INT reduction issues until after Ascona conference.


Status:
- supercedes hapy-output-20260517
- CURRENT


## hapy-output-20260620

Date:
2026-06-20

Input coadds:
- coadds-v20260330

Major changes:
- implemented min cutout size
- downloaded legacy images at their native pixel scale of 0.262"/pixel
- updated a bunch of masks as well
- changed the background estimation threshold to 2 sigma, 10 pixels.
in imutils.estimate_and_subtract_sky and
calculate_background_photutils.  was getting sky estimates that were
too big for halpha, so systematically negative values after sky
subtraction.  lowering the threshold seemed to help.

Pipeline changes:
- implemented min cutout size in get_cutouts of 90 arcsec

Reason for rerun:
- implementing a min cutout size should fix problems with a number.
  Problems are due to cutout images that are too small, e.g. can't get
  a good sky determination so mask is not good, phot is not good,
  etc. 
