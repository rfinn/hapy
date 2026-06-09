# Overview

The duplicate validation and old/new reduction comparisons indicate that the new INT r-band reductions are robust and significantly improved relative to the older reductions, with very small spatially uniform residuals ((< 0.02) dex) and improved duplicate reproducibility. In contrast, the new INT H-alpha reductions exhibit substantially larger duplicate scatter and clear spatially coherent residual structure across the detector. Rebuilding the CS-ZP images locally from the cutouts reduced some of the spatial systematics, but did not fully recover the duplicate consistency seen in the older H-alpha reductions. Additional tests using polynomial flattening with two rounds of illumination correction also failed to reproduce the quality of the older reductions, suggesting that the issue is more fundamental to the newer INT H-alpha processing. Based on these results, the current plan is to adopt a hybrid INT reduction for the catalog paper: use the older INT H-alpha coadds together with the newer INT r-band coadds, reprojecting the r-band images onto the native H-alpha grid to preserve the H-alpha image quality and PSF.

To improve the reliability of the INT continuum-subtracted H-alpha measurements, we constructed a hybrid INT dataset that combines the original pre-2025 H-alpha reductions with the newer 2026 r-band reductions. The newer INT H-alpha reductions showed large spatially dependent residuals and increased scatter in duplicate measurements, likely associated with the illumination-correction and continuum-subtraction processing. We therefore retained the older INT H-alpha coadds, while reprojecting the newer INT r-band coadds and corresponding weight images onto the astrometric grid of the older H-alpha images using SWarp. New PSF models were generated for the reprojected r-band images, while the original H-alpha PSFs were preserved. Continuum-subtracted images (CS-ZP) were then rebuilt directly from the matched r and H-alpha cutouts after local sky subtraction. This hybrid approach substantially improved the repeatability of duplicate measurements while preserving the superior quality of the newer INT r-band reductions.


# Alignment Issues
 
 

# Coadd with no cutouts b/c of bad weight image
 
 
 
 ```
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-162.760+32.934-INT-20190205-p065-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5
 ```
 
``` 
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-162.760+32.934-INT-20190205-p065-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5
 ```
 
 ```
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-162.760+32.934-INT-20190205-p065-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite
 ```
 
 ```
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-208.804+05.187-INT-20190206-p120-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
 ```
 
 
 ```
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-162.760+32.934-INT-20190205-p065-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
 ```
 
 ```
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-221.102+01.782-INT-20190209-p149-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
 ```
 
 ```
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
 ```
 
 ```
 get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
 ```



```
 python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts
```

CUTOUT SUMMARY
--------------
Input coadds:              222
Cutout directories:        824
Coadds with no cutouts:    2
Cutout dirs missing R:     3
Cutout dirs missing CS:    3
Bad coadd names:           0
Bad cutout dir names:      0

First few coadds with no cutouts:
  /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits
  /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits
  /data-pool/Halpha/coadds-v20260518/VF-226.932+01.661-INT-20190602-p033-r.fits

First few cutout dirs missing R:
  cutouts/VFID2927-SDSSJ121943.54+293931.7-INT-20220505-VFID2935
  cutouts/VFID2933-NGC4274-INT-20220505-VFID2935
  cutouts/VFID2959-NGC4286-INT-20220505-VFID2935

First few cutout dirs missing CS:
  cutouts/VFID2927-SDSSJ121943.54+293931.7-INT-20220505-VFID2935
  cutouts/VFID2933-NGC4274-INT-20220505-VFID2935
  cutouts/VFID2959-NGC4286-INT-20220505-VFID2935


# Working on /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits


get this output from get_cutouts

```
get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
Writing cutouts to: /data-pool/Halpha/hapy-output-20260519/cutouts
2026-05-19 22:32:54,380 INFO [pid=3040091] cutouts.VF-184.950+29.460-INT-20220505-VFID2935-r: Starting get_cutouts with rimage=/data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits scheme=virgo
WARNING: did not find CS image!
Loaded image /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits with shape (6800, 6727) and pixel scale 0.331 arcsec/pix
WARNING: FITSFixedWarning: 'datfix' made the change 'Set DATE-OBS to '2022-05-05T21:17:57.520' from MJD-OBS'. [astropy.wcs.wcs]
2026-05-19 22:32:54,397 WARNING [pid=3040091] astropy: FITSFixedWarning: 'datfix' made the change 'Set DATE-OBS to '2022-05-05T21:17:57.520' from MJD-OBS'.
Loaded image /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-Halpha.fits with shape (6800, 6727) and pixel scale 0.331 arcsec/pix
testing, self.hafilter =  Halpha INT
filter file =  WFC-Ha-197.fits
number of galaxies based on keepflag  = 5
found 5 after RA/DEC cuts
testing, self.hafilter =  Halpha INT
number of galaxies in FOV = 5
Skipping VFID2927-SDSSJ121943.54+293931.7: invalid cutout region (r_invalid); ra=184.931555,dec=29.658845
Skipping VFID2933-NGC4274: invalid cutout region (r_invalid); ra=184.960745,dec=29.614339
Skipping VFID2959-NGC4286: invalid cutout region (r_invalid); ra=185.175327,dec=29.345862
Skipping VFID2963-NGC4283: invalid cutout region (r_invalid); ra=185.086565,dec=29.310803
Skipping VFID2966-NGC4278: invalid cutout region (r_invalid); ra=185.028264,dec=29.280640
Wrote summary table: cutouts_summary/cutouts_summary-VF-184.950+29.460-INT-20220505-VFID2935-rfinn-20260519.ecsv

```

### Try moving r-band weight to BAD_WEIGHTS

in coadd directory:
```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ mv VF-184.950+29.460-INT-20220505-VFID2935-r.weight.fits BAD_WEIGHTS/.
```

then rerun `get_cutouts`, and it works fine!


### Try making a simple weight:
```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ python ~/github/hapy/scripts/make_simple_weight_from_coadd.py VF-184.950+29.460-INT-20220505-VFID2935-r.fits
Wrote VF-184.950+29.460-INT-20220505-VFID2935-r.weight.fits
Good fraction: 0.7876
```

Then rerun `get_cutouts` in output directory:
```
get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-184.950+29.460-INT-20220505-VFID2935-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
```
and that worked fine too.  Much easier then trying to retrofit the way get_cutouts was identifying bad regions in the coadd images.



### One image skipped

```
Skipping VFID2933-NGC4274: invalid cutout region (ha_invalid); ra=184.960745,dec=29.614339
```

- will check the location using a region in ds9
- this is one of the main target galaxies, and it's on a good region on the halpha image.  
- going to make a simple weight image and rerun

```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ mv VF-184.950+29.460-INT-20220505-VFID2935-Halpha.weight.fits BAD_WEIGHTS/.
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ python ~/github/hapy/scripts/make_simple_weight_from_coadd.py VF-184.950+29.460-INT-20220505-VFID2935-Halpha.fits 
WARNING: VerifyWarning: Card is too long, comment will be truncated. [astropy.io.fits.card]
Wrote VF-184.950+29.460-INT-20220505-VFID2935-Halpha.weight.fits
Good fraction: 0.7882

```
# Working on /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits

In coadd directory:
```
get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
```

I'm getting a similar output:

```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
Writing cutouts to: /data-pool/Halpha/hapy-output-20260519/cutouts
2026-05-19 22:41:50,504 INFO [pid=3049055] cutouts.VF-211.560+06.016-INT-20220505-VFID5726-r: Starting get_cutouts with rimage=/data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits scheme=virgo
WARNING: did not find CS image!
Loaded image /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits with shape (6808, 6727) and pixel scale 0.331 arcsec/pix
WARNING: FITSFixedWarning: 'datfix' made the change 'Set DATE-OBS to '2022-05-05T01:26:33.521' from MJD-OBS'. [astropy.wcs.wcs]
2026-05-19 22:41:50,521 WARNING [pid=3049055] astropy: FITSFixedWarning: 'datfix' made the change 'Set DATE-OBS to '2022-05-05T01:26:33.521' from MJD-OBS'.
Loaded image /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-Halpha.fits with shape (6808, 6727) and pixel scale 0.331 arcsec/pix
testing, self.hafilter =  Halpha INT
filter file =  WFC-Ha-197.fits
number of galaxies based on keepflag  = 1
found 1 after RA/DEC cuts
testing, self.hafilter =  Halpha INT
number of galaxies in FOV = 1
Skipping VFID5709-NGC5470: invalid cutout region (ha_invalid); ra=211.632965,dec=6.029321
Wrote summary table: cutouts_summary/cutouts_summary-VF-211.560+06.016-INT-20220505-VFID5726-rfinn-20260519.ecsv

```

### Try moving coadd weight to BAD_WEIGHTS and make a simple weight

```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ mv VF-211.560+06.016-INT-20220505-VFID5726-r.weight.fits BAD_WEIGHTS/
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ python ~/github/hapy/scripts/make_simple_weight_from_coadd.py VF-211.560+06.016-INT-20220505-VFID5726-r.fits
Wrote VF-211.560+06.016-INT-20220505-VFID5726-r.weight.fits
Good fraction: 0.7891
```

The rerun `get_cutouts`. 
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
Writing cutouts to: /data-pool/Halpha/hapy-output-20260519/cutouts
2026-05-19 22:43:41,655 INFO [pid=3050572] cutouts.VF-211.560+06.016-INT-20220505-VFID5726-r: Starting get_cutouts with rimage=/data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits scheme=virgo
WARNING: did not find CS image!
Loaded image /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-r.fits with shape (6808, 6727) and pixel scale 0.331 arcsec/pix
WARNING: FITSFixedWarning: 'datfix' made the change 'Set DATE-OBS to '2022-05-05T01:26:33.521' from MJD-OBS'. [astropy.wcs.wcs]
2026-05-19 22:43:41,671 WARNING [pid=3050572] astropy: FITSFixedWarning: 'datfix' made the change 'Set DATE-OBS to '2022-05-05T01:26:33.521' from MJD-OBS'.
Loaded image /data-pool/Halpha/coadds-v20260518/VF-211.560+06.016-INT-20220505-VFID5726-Halpha.fits with shape (6808, 6727) and pixel scale 0.331 arcsec/pix
testing, self.hafilter =  Halpha INT
filter file =  WFC-Ha-197.fits
number of galaxies based on keepflag  = 1
found 1 after RA/DEC cuts
testing, self.hafilter =  Halpha INT
number of galaxies in FOV = 1
Skipping VFID5709-NGC5470: invalid cutout region (ha_invalid); ra=211.632965,dec=6.029321
Wrote summary table: cutouts_summary/cutouts_summary-VF-211.560+06.016-INT-20220505-VFID5726-rfinn-20260519.ecsv

```

Oops - didn't realize that the error was with the halpha image the first time.


```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ mv BAD_WEIGHTS/VF-211.560+06.016-INT-20220505-VFID5726-r.weight.fits .
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ mv VF-211.560+06.016-INT-20220505-VFID5726-Halpha.weight.fits BAD_WEIGHTS/
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ python ~/github/hapy/scripts/make_simple_weight_from_coadd.py VF-211.560+06.016-INT-20220505-VFID5726-Halpha.fits
WARNING: VerifyWarning: Card is too long, comment will be truncated. [astropy.io.fits.card]
Wrote VF-211.560+06.016-INT-20220505-VFID5726-Halpha.weight.fits
Good fraction: 0.7871
```

Then get_cutouts ran fine.


**Need to track these messages in the log file so we can look for problems.**

# Working on /data-pool/Halpha/coadds-v20260518/VF-226.932+01.661-INT-20190602-p033-r.fits

```
hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-226.932+01.661-INT-20190602-p033-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
Writing cutouts to: /data-pool/Halpha/hapy-output-20260519/cutouts
2026-05-19 22:51:19,829 INFO [pid=3058416] cutouts.VF-226.932+01.661-INT-20190602-p033-r: Starting get_cutouts with rimage=/data-pool/Halpha/coadds-v20260518/VF-226.932+01.661-INT-20190602-p033-r.fits scheme=virgo
WARNING: did not find CS image!
Loaded image /data-pool/Halpha/coadds-v20260518/VF-226.932+01.661-INT-20190602-p033-r.fits with shape (6798, 6617) and pixel scale 0.332 arcsec/pix
WARNING: FITSFixedWarning: 'datfix' made the change 'Invalid parameter values: MJD-OBS and DATE-OBS are inconsistent'. [astropy.wcs.wcs]
2026-05-19 22:51:19,845 WARNING [pid=3058416] astropy: FITSFixedWarning: 'datfix' made the change 'Invalid parameter values: MJD-OBS and DATE-OBS are inconsistent'.
Loaded image /data-pool/Halpha/coadds-v20260518/VF-226.932+01.661-INT-20190602-p033-Halpha.fits with shape (6798, 6617) and pixel scale 0.332 arcsec/pix
testing, self.hafilter =  Halpha INT
filter file =  WFC-Ha-197.fits
number of galaxies based on keepflag  = 10
found 10 after RA/DEC cuts
testing, self.hafilter =  Halpha INT
number of galaxies in FOV = 10
Skipping VFID6350-2MASXJ15084716+0154003: invalid cutout region (r_invalid); ra=227.196612,dec=1.900052
Skipping VFID6372-WISEAJ150822.72+014754.7: invalid cutout region (r_invalid); ra=227.094510,dec=1.798534
Skipping VFID6386-WISEAJ150825.51+014224.3: invalid cutout region (r_invalid); ra=227.106873,dec=1.706933
Skipping VFID6392-WISEAJ150819.90+014123.6: invalid cutout region (r_invalid); ra=227.082771,dec=1.689625
Skipping VFID6403-CGCG021-013: invalid cutout region (r_invalid); ra=227.023395,dec=1.651582
Skipping VFID6420-CGCG021-015: invalid cutout region (r_invalid); ra=227.038538,dec=1.608510
Skipping VFID6434-WISEAJ150634.27+013331.7: invalid cutout region (r_invalid); ra=226.642704,dec=1.558832
Skipping VFID6437-NGC5850: invalid cutout region (r_invalid); ra=226.781914,dec=1.544552
Skipping VFID6447-SDSSJ150812.35+012959.7: invalid cutout region (r_invalid); ra=227.051591,dec=1.499728
Skipping VFID6463-WISEAJ150809.13+012516.6: invalid cutout region (r_invalid); ra=227.038311,dec=1.420577
Wrote summary table: cutouts_summary/cutouts_summary-VF-226.932+01.661-INT-20190602-p033-rfinn-20260519.ecsv
```

another one where most of the weight image is zero.  move it to bad and make a simple weight.

```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ mv VF-226.932+01.661-INT-20190602-p033-r.weight.fits BAD_WEIGHTS/.
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ python ~/github/hapy/scripts/make_simple_weight_from_coadd.py VF-226.932+01.661-INT-20190602-p033-r.fits 
Wrote VF-226.932+01.661-INT-20190602-p033-r.weight.fits
Good fraction: 0.7989

```

the rerun `get_cutouts`

all cutouts saved except:

```
Skipping VFID6434-WISEAJ150634.27+013331.7: invalid cutout region (ha_invalid); ra=226.642704,dec=1.558832
```

and this is really off the edge of the image - excellent!


# Rechecking cutouts

```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              222
Cutout directories:        834
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0

```

# Summary

I was able to extra cutouts from 3 problematic coadd pairs by removing the original weight images and making a simple weight image.  the simple weight image has a value of 1 inside the ccds, and 0 outside.   It identifies bad regions as where image values are zero.


Now this doesn't mean that all of the cutouts were created ok.  there could others that were skipped and we just don't know it. 

We could:
- add logging info to get_cutouts so we can check logs for skipped galaxies
- convert all of the INT masks to simple masks - this ensures that only galaxies on the edges will be cut.
- or do both.


INT hybrid coadds: use simple science-footprint weights
Other coadds: keep existing weights


# Creating Simple Weights for all INT coadds
```
mv BAD_WEIGHTS  BAD_WEIGHTS_ORIG
```

```
mv *INT*weight.fits ORIGINAL_INT_WEIGHTS/.
```
with the exception of the 3 that already went into `BAD_WEIGHTS_ORIG`


create a list of all INT coadds:
```
ls VF*INT*r.fits VF*INT*Halpha.fits VF*INT*Ha6657.fits > INT_all_coadds.txt
```
or

```
find . -maxdepth 1 -type f \( -name "VF*INT*r.fits" -o -name "VF*INT*Halpha.fits" -o -name "VF*INT*Ha6657.fits" \) | sort > INT_all_coadds.txt
```

```
parallel --bar -j 16 python ~/github/hapy/scripts/make_simple_weight_from_coadd.py {} :::: INT_all_coadds.txt
```

```
parallel --bar -j 16 --joblog make_simple_weight_INT.joblog --results make_simple_weight_logs python ~/github/hapy/scripts/make_simple_weight_from_coadd.py {} :::: INT_all_coadds.txt
```

### Checking output

```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ awk 'NR==1 || $7 != 0 {print}' cutouts_parallel.log
Seq	Host	Starttime	JobRuntime	Send	Receive	Exitval	Signal	Command
69	:	1779233598.248	    12.708	0	1193	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-162.760+32.934-INT-20190205-p065-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
100	:	1779233610.030	    16.910	0	2902	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-179.740+50.875-INT-20220503-VFID1277-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
97	:	1779233608.816	    70.387	0	1233	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-178.160+52.290-INT-20220503-VFID1213-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
156	:	1779233664.870	    41.046	0	2833	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-208.804+05.187-INT-20190206-p120-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
179	:	1779233690.251	    84.259	0	1193	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-221.102+01.782-INT-20190209-p149-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata

```

five different coadds - woo hoo!


and two different cutouts:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              222
Cutout directories:        856
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    2
Bad coadd names:           0
Bad cutout dir names:      0

First few cutout dirs missing CS:
  cutouts/VFID1212-NGC3953-INT-20220503-VFID1213
  cutouts/VFID1297-UGC06992-INT-20220503-VFID1277
```


Let's go through coadds one at a time.


### undid slicing of halpha and rband cutouts

- so had to recreate the wheel on that one...
- the other 4 coadds then ran fine!


### Updated summary

```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              222
Cutout directories:        856
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0

```

### going to run `get_cutouts` one more time...

just so that all of the logs are updated.


```
cat fullpath_rcoadds_all.txt | parallel -j 16 --bar --joblog cutouts_parallel.log get_cutouts --rimage {} --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata
```

- no jobs failed
- got the same results - 856 cutouts.


```
grep -R "SKIP_INVALID_REGION" logs/*.cutouts.log
```

- found 33 skipped galaxies with a mix of INT (16), BOK (13), HDI (1), MOS (3)


#### spot checking skipped galaxies

```
python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts
```
  



```
grep -R "SKIP_INVALID_REGION" logs/*.cutouts.log | awk -F"status=" '{print $2}' | sort | uniq -c
      2 ha_invalid
     31 r_invalid
```

logs/VF-124.537+24.716-INT-20190207-lmp002-r.cutouts.log:2026-05-20 00:18:57,220 WARNING [pid=3807413] cutouts.VF-124.537+24.716-INT-20190207-lmp002-r: SKIP_INVALID_REGION galid=VFID3338-KUG0814+251 ra=124.33766550 dec=24.96266860 size_arcsec=84.664 status=r_invalid

- galaxy is in uncovered region in NW of coadd

logs/VF-162.760+32.934-INT-20190205-p065-r.cutouts.log:2026-05-20 00:19:30,991 WARNING [pid=3823704] cutouts.VF-162.760+32.934-INT-20190205-p065-r: SKIP_INVALID_REGION galid=VFID2701-2MASXiJ1052228+323813 ra=163.09537950 dec=32.63733260 size_arcsec=107.138 status=r_invalid

- off the bottom edge


# Summary
So the final state is now much cleaner:

get_cutouts logic remains simple and fast
precheck behavior preserved
problematic INT weight maps replaced upstream
true off-FOV systems still rejected
slice-sharing between Hα and r ensures identical cutout dimensions
skip logging lets you audit failures later


# Next Steps to get the hybrid sample through run_analysis

### Already done
- gaia catalogs are downloaded
- psf were rebuilt for INT
- cutouts are mad

##### rsync existing legacy images
- rsync legacy images for BOK, HDI, MOS from /data-pool/Halpha/hapy-output-v20260330/cutouts/ to /data-pool/Halpha/hapy-output-v20260519/cutouts/

Dry run:
```
rsync -avn --ignore-existing  --include="*BOK*" --include="*HDI*" --include="*MOS*" --exclude="*INT*"  --relative /data-pool/Halpha/hapy-output-20260417/cutouts/./*/legacy/* /data-pool/Halpha/hapy-output-20260519/cutouts/
```

```
rsync -av --ignore-existing  --include="*BOK*" --include="*HDI*" --include="*MOS*" --exclude="*INT*"  --relative /data-pool/Halpha/hapy-output-20260417/cutouts/./*/legacy/* /data-pool/Halpha/hapy-output-20260519/cutouts/
```

- DONE rsync manual masks for BOK, HDI, MOS from /data-pool/Halpha/hapy-output-20260417/cutouts


##### rsync or  REPROJECT existing manual masks


  - look in source directories: /data-pool/Halpha/hapy-output-20260417/cutouts/<tag>/*-mask-manual.fits
  - if INT in <tag>
	- reproject mask-manual.fits onto /data-pool/Halpha/hapy-output-20260519/cutouts/<tag>/*Ha.fits and save it in /data-pool/Halpha/hapy-output-20260519/cutouts/<tag>/
  - if BOK, MOS, HDI in <tag>:
    - rsync the *mask-manual.fits to /data-pool/Halpha/hapy-output-20260519/cutouts/<tag>/.

```
python ~/github/hapy/scripts/transfer_manual_masks.py --dry-run > transfer_manual_mask_output.txt
```

Check output for any unmatched coadds (needed to allow for different dateobs).  Then:
```
python ~/github/hapy/scripts/transfer_manual_masks.py
```

The end of the output is:
```
Summary
  copied: 60
  no_manual_mask: 723
  reprojected: 41
```

##### make cs-gr images
  - reproject legacy g and r to halpha
  - construct the CS-gr image
  

```
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID1934-NGC2799-INT-20190205-p026 --auto-contscale --auto-contscale-method ratio --auto-contscale-percentile 30 --overwrite
```
```
parallel --bar -j 16 --joblog cs_gr_auto.joblog --results \
cs_gr_auto_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {} \
--auto-contscale --auto-contscale-percentile 30 --overwrite :::: \
reproject_cutout_list.txt
```
### In Progress 
#####  download legacy images for INT coadds 
```
find /data-pool/Halpha/hapy-output-20260519/cutouts/ -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep INT |sort > INT_cutout_list.txt

ROOTDIR=/data-pool/Halpha/hapy-output-20260519/

parallel --bar -j 2 --joblog fetch_legacy.joblog --results fetch_legacy_logs python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir "$ROOTDIR/cutouts/{}"  :::: INT_cutout_list.txt
```


Checking results, 
```
awk 'NR==1 || $7 != 0 {print}' fetch_legacy.joblog  |wc -l
```
and 35 failed.


- I retested one, and it looks like it was a problem downloading the grz image all at once.
- I removed the *grz.fits, and then it worked ok.
- Unfortunately, the script looks for the grz.fits image and won't redownload.
- should look instead for the g.fits, r.fits and z.fits
- fixed 

This is taking forever.  I've tried rerunning the parallel command, and I still  have 24 without a download. I am going to remove asking legacy to do the reprojection and just download the native pixel scale.  I'

**Legacy website is down :(**


##### then run `run_analysis`

```
run_analysis --make-mask  --psf-dir /data-pool/Halpha/psf-images-v20260518/ --statmorph \
--galfit --convflag --log-to-console --gaia-dir \
/data-pool/Halpha/coadds-v20260518/gaia_catalogs/ --cutout-dir \
cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```


```
parallel --bar  -j 16  --memfree 30G --joblog run_analysis.joblog \
--results parallel-logs run_analysis --cutout-dir "{}" --make-mask \
--psf-dir /data-pool/Halpha/psf-images-v20260518/ --statmorph --galfit \
--convflag --gaia-dir /data-pool/Halpha/coadds-v20260518/gaia_catalogs/ :::: cutout_with_dir.txt
```

Started at 2:54 PM, 2026-05-20

5:02 PM, @99% 853:3

5:16 PM @99% 854:2
5:42 PM @99% 855:1=11s
##### Merge Results
ran this while `run_analysis` was still at 99%
```
merge_results --indir cutouts --mode run_analysis --review-csv review_sample_20260514.csv
```

```
	validated 603/818 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260519/merged_results_virgo_20260520.fits
Done.
Final table rows: 603
Final table columns: 494
```


##### Select Best Duplicate

```
python ~/github/hapy/scripts/inspect_cs_images.py make-table merged_results_virgo_20260520.fits --outdir cs_image_inspection --min-dups 1
```
### Still to do 






- check duplicates 
  - create visualizations
  - run validation notebook

other issues:
- the old INT halpha images used a different color transformation when solving for the ZP with `getzp`, and the same solution is used for Halpha and Ha6657 filters.  Do we want to rerun getzp for these?  Let's look at results first...

- create a new set of coadd webpages

what else am I missing? oh yeah, write the paper!


# Working on INT 2019 Halpha issue

Made test directory:
```
/data-pool/Halpha/hapy-test-INT-2019-halpha-gain-adjustment
```

set the gain to 1800:
```
sethead GAIN=1800
cutouts/VFID2313-NGC3294-INT-20190211-p059/VFID2313-NGC3294-INT-20190211-p059-Ha.fits
```


```
sethead GAIN=1800
cutouts/VFID2313-NGC3294-INT-20190211-p059/VFID2313-NGC3294-INT-20190211-p059-CS-ZP.fits
```

```
sethead GAIN=1800
cutouts/VFID2313-NGC3294-INT-20190211-p059/VFID2313-NGC3294-INT-20190211-p059-CS-gr.fits
```

Then run `run_analysis`
```
run_analysis --make-mask --psf-dir /data-pool/Halpha/psf-images-v20260518/ --log-to-console --gaia-dir \
/data-pool/Halpha/coadds-v20260518/gaia_catalogs/ --cutout-dir \
cutouts/VFID2313-NGC3294-INT-20190211-p059
```

### Found problem
- photometry.py was using one function to get the noise in the
  aperture for both image1 and image2.  so it was using the wrong gain
  AND the wrong sky noise. 
-  I updated the class function, but also made a regular function not attached to the class that explicitly takes in the sky noise and gain.  I am using this for calculating the noise in image2 apertures.
- I tested on VFID2313, and we now get values for the Halpha profiles.
- also, to speed up testing, I added a --csgr flag to run_analysis so we can skip running on the CS-gr images if so desired.
- I also made run_analysis get filter_ratio from the PHOTZP keywords in the r_fits and cs_fits images


This probably warrants a HAPY version update


### Fixing INT Data
- [x] in 2019 INT Halpha coadds, set GAIN = 1800
```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ sethead GAIN=1800 VF*INT*2019*Halpha.fits
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260518$ sethead GAIN=1800 VF*INT*2019*Ha6657.fits
```
- [x] in hapy-output-20260519-hybrid, for 2019 Halpha cutouts, set GAIN = 1800
  - [x] cutouts/*INT*2019*/*Ha.fits
  - [x] cutouts/*INT*2019*/*CS-gr.fits
  - [x] cutouts/*INT*2019*/*CS-ZP.fits
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519-hybrid/cutouts$ sethead GAIN=1800 *INT*2019*/*Ha.fits
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519-hybrid/cutouts$ sethead GAIN=1800 *INT*2019*/*CS-ZP.fits
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519-hybrid/cutouts$ sethead GAIN=1800 *INT*2019*/*CS-gr.fits
```

- [x] test `run_analysis` on VFID2313 and make sure `H_PROFILE` values
  are populated in the results.ecsv table.
```
run_analysis --make-mask --psf-dir /data-pool/Halpha/psf-images-v20260518/ --log-to-console --gaia-dir \
/data-pool/Halpha/coadds-v20260518/gaia_catalogs/ --cutout-dir \
cutouts/VFID2313-NGC3294-INT-20190211-p059
```
- [x] rerun run_analysis on 2019 INT cutouts
  - [x] create input list
  ```
  cat cutout_with_dir.txt | grep INT |grep 2019 >
  INT_2019_cutouts_with_dir.txt
  ```
  - [x] run parallel just on this list, and with `--csgr` flag
```
parallel --bar  -j 16  --memfree 30G --joblog run_analysis_INT_2019.joblog \
--results parallel-logs run_analysis --cutout-dir "{}" --make-mask \
--psf-dir /data-pool/Halpha/psf-images-v20260518/ --statmorph --galfit \
--convflag --gaia-dir \
/data-pool/Halpha/coadds-v20260518/gaia_catalogs/ --csgr :::: INT_2019_cutouts_with_dir.txt
```
- [ ] then rerun on full sample b/c error in photometry.py will affect
  the noise estimates in all of the halpha phot profiles
- [ ] Should consider rerunning on the pre2025 data as well
  - [x] in hapy-output-20260517-pre2025coadds, for 2019 Halpha cutouts, set GAIN = 1800
	  - [x] cutouts/*INT*2019*/*Ha.fits
	  - [ ] cutouts/*INT*2019*/*CS-gr.fits (N/A)
	  - [x] cutouts/*INT*2019*/*CS-ZP.fits

```
parallel --bar  -j 16  --memfree 30G --joblog run_analysis.joblog \
--results parallel-logs run_analysis --cutout-dir "{}" --make-mask \
--psf-dir /data-pool/Halpha/psf-images-v20260518/ --statmorph --galfit \
--convflag --gaia-dir \
/data-pool/Halpha/coadds-pre2025-hapy/gaia_catalogs/  :::: cutout_run_analysis_list.txt
```
