# Rename the files to remove shifted and create a dedicated directory
```
python ~/github/hapy/scripts/copy_clean_coadds.py /data-pool/Halpha/coadds-pre2025/all-virgo-coadds /data-pool/Halpha/coadds-pre2025-hapy
```

# Download the gaia files
Move to directory containing coadded images

```bash
cd /data-pool/Halpha/coadds-pre2025-hapy
```

Then start the download.

```bash
python ~/github/hapy/scripts/download_gaia_coadd_catalogs.py
```

# Add H-alpha image to headers


# Prepare the output directory
### make a directory for analysis run

```bash
cd /data-pool/Halpha/
```

```bash
mkdir hapy-output-pre2025coadds-20260516
```

```bash
cd hapy-output-pre2025coadds-20260516
```


### make a list of coadds to analyze

```bash
find /data-pool/Halpha/coadds-pre2025-hapy/ -maxdepth 1 -type f \( -name "VF*r.fits" -o -name "VF*R.fits" \) | sort > fullpath_rcoadds_all.txt
```


# Make Cutouts

### Test on One Image

### 

```bash
get_cutouts --rimage
/data-pool/Halpha/coadds-pre2025-hapy/VF-126.297+27.994-HDI-20180313-p004-R.fits
--catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits
--scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-pre2025/
```
this should create two folders in the cutouts/ directory.

get_cutouts --rimage /data-pool/Halpha/coadds-pre2025-hapy/VF-266.477+58.350-BOK-20220423-VFID0783-r.fits --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-pre2025/



get_cutouts --rimage /data-pool/Halpha/coadds-pre2025-hapy/VF-147.027+33.416-INT-20190207-p039-r.fits --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-pre2025/



### Run on full Virgo sample


Or run on the full list of coadds.
```bash
cat fullpath_rcoadds_all.txt | parallel -j 16 --bar --joblog
cutouts_parallel.log get_cutouts --rimage {} --catalog
~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo
--maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-pre2025/ --overwrite
```

### Check output
```bash
python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts
```
Output:
```
CUTOUT SUMMARY
--------------
Input coadds:              225
Cutout directories:        747
Coadds with no cutouts:    8
Cutout dirs missing R:     0
Cutout dirs missing CS:    1
Bad coadd names:           0
Bad cutout dir names:      0

First few coadds with no cutouts:
  /data-pool/Halpha/coadds-pre2025-hapy/VF-139.836+25.869-BOK-20210316-VFID3299-r.fits
  /data-pool/Halpha/coadds-pre2025-hapy/VF-211.560+06.016-INT-20220505-VFID5726-r.fits
  /data-pool/Halpha/coadds-pre2025-hapy/VF-217.470+3.421-MOS-20120422-NGC5846_07-R.fits
  /data-pool/Halpha/coadds-pre2025-hapy/VF-222.729+2.735-MOS-20120423-NGC5846_06-R.fits
  /data-pool/Halpha/coadds-pre2025-hapy/VF-223.132+3.623-MOS-20130415-NGC5846_05-R.fits
  /data-pool/Halpha/coadds-pre2025-hapy/VF-225.207+1.904-MOS-20110404-NGC5846_02-R.fits
  /data-pool/Halpha/coadds-pre2025-hapy/VF-225.738+0.831-MOS-20120418-NGC5846_04-R.fits
  /data-pool/Halpha/coadds-pre2025-hapy/VF-254.754+23.221-BOK-20210418-VFID3459-r.fits

First few cutout dirs missing CS:
  cutouts/VFID2057-NGC5371-HDI-20170522-p016
```

Moving forward anyway - who knows why these are missing...

## Merge get_cutouts tables
```bash
merge_results --mode get_cutouts --indir cutouts_summary --out merged_cutouts_results.fits
```

Example output:
```
Searching for files  cutouts_summary*.ecsv
Found 215 result files.
Merging 215 result files.
Reading tables...
Validating schema...
WAIT!!! Problem with table cutouts_summary/cutouts_summary-VF-211.560+06.016-INT-20220505-VFID5726-rfinn-20260516.ecsv!!!
Schema mismatch detected.

Missing columns:
['cutout_ell0_missing_frac_h', 'cutout_ell0_missing_frac_max', 'cutout_ell0_missing_frac_r', 'cutout_ell0_npix_good_h', 'cutout_ell0_npix_good_r', 'cutout_ell0_npix_onimage_h', 'cutout_ell0_npix_onimage_r', 'cutout_ell0_npix_total_h', 'cutout_ell0_npix_total_r', 'cutout_root', 'dateobs', 'dec', 'filter_correction', 'filter_warning', 'hafilter', 'objid', 'parent_haimage', 'parent_rimage', 'pointing', 'ra', 'scheme', 'size_arcsec', 'tag', 'telescope', 'valid_region', 'valid_status', 'x_parent', 'y_parent']

WAIT!!! Problem with table cutouts_summary/cutouts_summary-VF-217.470+3.421-MOS-20120422-NGC5846_07-rfinn-20260516.ecsv!!!
Schema mismatch detected.

Missing columns:
['cutout_ell0_missing_frac_h', 'cutout_ell0_missing_frac_max', 'cutout_ell0_missing_frac_r', 'cutout_ell0_npix_good_h', 'cutout_ell0_npix_good_r', 'cutout_ell0_npix_onimage_h', 'cutout_ell0_npix_onimage_r', 'cutout_ell0_npix_total_h', 'cutout_ell0_npix_total_r', 'cutout_root', 'dateobs', 'dec', 'filter_correction', 'filter_warning', 'hafilter', 'objid', 'parent_haimage', 'parent_rimage', 'pointing', 'ra', 'scheme', 'size_arcsec', 'tag', 'telescope', 'valid_region', 'valid_status', 'x_parent', 'y_parent']

WAIT!!! Problem with table cutouts_summary/cutouts_summary-VF-222.729+2.735-MOS-20120423-NGC5846_06-rfinn-20260516.ecsv!!!
Schema mismatch detected.

Missing columns:
['cutout_ell0_missing_frac_h', 'cutout_ell0_missing_frac_max', 'cutout_ell0_missing_frac_r', 'cutout_ell0_npix_good_h', 'cutout_ell0_npix_good_r', 'cutout_ell0_npix_onimage_h', 'cutout_ell0_npix_onimage_r', 'cutout_ell0_npix_total_h', 'cutout_ell0_npix_total_r', 'cutout_root', 'dateobs', 'dec', 'filter_correction', 'filter_warning', 'hafilter', 'objid', 'parent_haimage', 'parent_rimage', 'pointing', 'ra', 'scheme', 'size_arcsec', 'tag', 'telescope', 'valid_region', 'valid_status', 'x_parent', 'y_parent']

WAIT!!! Problem with table cutouts_summary/cutouts_summary-VF-223.132+3.623-MOS-20130415-NGC5846_05-rfinn-20260516.ecsv!!!
Schema mismatch detected.

Missing columns:
['cutout_ell0_missing_frac_h', 'cutout_ell0_missing_frac_max', 'cutout_ell0_missing_frac_r', 'cutout_ell0_npix_good_h', 'cutout_ell0_npix_good_r', 'cutout_ell0_npix_onimage_h', 'cutout_ell0_npix_onimage_r', 'cutout_ell0_npix_total_h', 'cutout_ell0_npix_total_r', 'cutout_root', 'dateobs', 'dec', 'filter_correction', 'filter_warning', 'hafilter', 'objid', 'parent_haimage', 'parent_rimage', 'pointing', 'ra', 'scheme', 'size_arcsec', 'tag', 'telescope', 'valid_region', 'valid_status', 'x_parent', 'y_parent']

WAIT!!! Problem with table cutouts_summary/cutouts_summary-VF-225.207+1.904-MOS-20110404-NGC5846_02-rfinn-20260516.ecsv!!!
Schema mismatch detected.

Missing columns:
['cutout_ell0_missing_frac_h', 'cutout_ell0_missing_frac_max', 'cutout_ell0_missing_frac_r', 'cutout_ell0_npix_good_h', 'cutout_ell0_npix_good_r', 'cutout_ell0_npix_onimage_h', 'cutout_ell0_npix_onimage_r', 'cutout_ell0_npix_total_h', 'cutout_ell0_npix_total_r', 'cutout_root', 'dateobs', 'dec', 'filter_correction', 'filter_warning', 'hafilter', 'objid', 'parent_haimage', 'parent_rimage', 'pointing', 'ra', 'scheme', 'size_arcsec', 'tag', 'telescope', 'valid_region', 'valid_status', 'x_parent', 'y_parent']

WAIT!!! Problem with table cutouts_summary/cutouts_summary-VF-225.738+0.831-MOS-20120418-NGC5846_04-rfinn-20260516.ecsv!!!
Schema mismatch detected.

Missing columns:
['cutout_ell0_missing_frac_h', 'cutout_ell0_missing_frac_max', 'cutout_ell0_missing_frac_r', 'cutout_ell0_npix_good_h', 'cutout_ell0_npix_good_r', 'cutout_ell0_npix_onimage_h', 'cutout_ell0_npix_onimage_r', 'cutout_ell0_npix_total_h', 'cutout_ell0_npix_total_r', 'cutout_root', 'dateobs', 'dec', 'filter_correction', 'filter_warning', 'hafilter', 'objid', 'parent_haimage', 'parent_rimage', 'pointing', 'ra', 'scheme', 'size_arcsec', 'tag', 'telescope', 'valid_region', 'valid_status', 'x_parent', 'y_parent']

	validated 209/215 tables
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-pre2025coadds-20260516/merged_cutouts_results.fits
Done.
Final table rows: 706
Final table columns: 28

```


# Run Analysis 

## Run on One Cutout


```bash
run_analysis --make-mask  --psf-dir /data-pool/Halpha/psf-images-pre2025/ --statmorph --galfit --convflag --log-to-console --gaia-dir /data-pool/Halpha/coadds-pre2025-hapy/gaia_catalogs/ --cutout-dir cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```

woo hoo!

## Running on a larger sample

### Create a list of cutouts:

```bash
find cutouts/ -mindepth 1 -maxdepth 1 -type d ! -name "cutouts_summary" | sort > cutout_with_dir.txt
```

check that the file contains, e.g., `cutouts/VFID2943-NGC2604-INT-20190204-p010`
```bash
head cutout_list.txt
wc -l cutout_list.txt
```

### Test smaller samples
Test on 5 galaxies:
```bash
head -5 cutout_list.txt | parallel --bar -j 2 --joblog run_analysis.joblog \
  --results parallel-logs \
  run_analysis --cutout-dir "{}" --make-mask \
  --psf-dir /data-pool/Halpha/psf-images/ \
  --statmorph --galfit --convflag \
  --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/
```


### Run in Parallel


Copy results from manual inspection that will exclude bad targets:
```bash
cp ../hapy-output-20260417/review_sample_20260514.csv .
```



Create an input list that removes objects that have `CATALOG_USE==EXCLUDE`:
```
python ~/github/hapy/scripts/make_run_analysis_list.py --cutout-dir
cutouts --review review_sample_20260514.csv --outfile
cutout_run_analysis_list.txt
```



If manual masking has been done, do not force rebuilding of masks:
```bash
parallel --bar  -j 16  --memfree 30G --joblog run_analysis.joblog
--results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir /data-pool/Halpha/psf-images-pre2025/ --statmorph --galfit
--convflag --gaia-dir /data-pool/Halpha/coadds-pre2025-hapy/gaia_catalogs/ :::: cutout_run_analysis_list.txt
```




## Merge Results From `run_analysis`

```
cp ../hapy-output-20260417/review_sample_20260514.csv .
```

```bash
merge_results --indir cutouts --mode run_analysis --review-csv review_sample_20260514.csv
```

Example output on 2026-05-16:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-pre2025coadds-20260516$ merge_results --indir cutouts --mode run_analysis --review-csv review_sample_20260514.csv
Searching for files  *results.ecsv
Found 716 result files.
Skipping 0 files with CATALOG_USE == EXCLUDE
Merging 716 result files.
Reading tables...
Validating schema...
	validated 716/716 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-pre2025coadds-20260516/merged_results_virgo_20260517.fits
Done.
Final table rows: 716
Final table columns: 386

```



Example output from new sample as of 2026-05-14:
```
Searching for files  *results.ecsv
Found 821 result files.
Skipping 39 files with CATALOG_USE == EXCLUDE
Merging 782 result files.
Reading tables...
Validating schema...
	validated 782/782 tables
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260417/merged_results_virgo_20260514.fits
Done.
Final table rows: 782
Final table columns: 358
```

### Summarize statistics in  `merged_results.py`
```
python ~/github/hapy/scripts/summarize_run.py --infile merged_results.fits --scheme virgo
```

Results from pre2025 coadds on 2026-04-17:
```
HAPY RUN SUMMARY
----------------
Total galaxies: 682
Unique galaxies: 571

Pipeline completion
-------------------
PSF_OK            :  682 OK  |    0 FAIL  (100.0%)
MASK_OK           :  682 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  430 OK  |  252 FAIL  ( 63.0%)
HAPY_MORPH_OK     :  374 OK  |  308 FAIL  ( 54.8%)
R_PROFILE_OK      :  428 OK  |  254 FAIL  ( 62.8%)
H_PROFILE_OK      :  343 OK  |  339 FAIL  ( 50.3%)
R_SM_OK           :  357 OK  |  325 FAIL  ( 52.3%)
H_SM_OK           :  357 OK  |  325 FAIL  ( 52.3%)
GAL_NC_OK         :  259 OK  |  423 FAIL  ( 38.0%)
GAL_CV_OK         :  258 OK  |  424 FAIL  ( 37.8%)
R_PETRO_OK        :  424 OK  |  258 FAIL  ( 62.2%)
R_EXPFIT_OK       :  428 OK  |  254 FAIL  ( 62.8%)
R_LOGFIT_OK       :  428 OK  |  254 FAIL  ( 62.8%)
H_PETRO_OK        :  309 OK  |  373 FAIL  ( 45.3%)
H_EXPFIT_OK       :  299 OK  |  383 FAIL  ( 43.8%)
H_LOGFIT_OK       :  299 OK  |  383 FAIL  ( 43.8%)
BRIGHT_STAR_FLAG  :    0 OK  |  682 FAIL  (  0.0%)
HAPY_MORPH_FLAG   :   64 OK  |  618 FAIL  (  9.4%)
R_SM_FLAG         :  383 OK  |  299 FAIL  ( 43.8%)
R_SM_SERSIC_FLAG  :  329 OK  |  353 FAIL  ( 51.8%)
H_SM_FLAG         :  467 OK  |  215 FAIL  ( 31.5%)
H_SM_SERSIC_FLAG  :  363 OK  |  319 FAIL  ( 46.8%)
CSGR_PHOT_OK      :    0 OK  |  682 FAIL  (  0.0%)
CSGR_HAPY_MORPH_OK:    0 OK  |  682 FAIL  (  0.0%)
CSGR_HAPY_MORPH_FLAG:    0 OK  |  682 FAIL  (  0.0%)


PROFILES_BOTH     :  343 ( 50.3%)

STATMORPH_BOTH    :   38

STATUS counts
-------------
ok          : 424
running     : 258

STAGE counts
------------
done        : 424
init        : 31
mask        : 221
phot        : 6

Runtime medians (sec)
---------------------
MASK_SEC    : 0.08
PHOT_SEC    : 3.67
SM_SEC      : 3.25
TOTAL_SEC   : 12.88

Number with bad phot = 252

```


Results from new coadds on 2026-04-17:

```
HAPY RUN SUMMARY
----------------
Total galaxies: 823
Unique galaxies: 659

Pipeline completion
-------------------
PSF_OK            :  823 OK  |    0 FAIL  (100.0%)
MASK_OK           :  823 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  809 OK  |   14 FAIL  ( 98.3%)
HAPY_MORPH_OK     :  771 OK  |   52 FAIL  ( 93.7%)
R_PROFILE_OK      :  808 OK  |   15 FAIL  ( 98.2%)
H_PROFILE_OK      :  652 OK  |  171 FAIL  ( 79.2%)
R_SM_OK           :  724 OK  |   99 FAIL  ( 88.0%)
H_SM_OK           :  724 OK  |   99 FAIL  ( 88.0%)
GAL_NC_OK         :  657 OK  |  166 FAIL  ( 79.8%)
GAL_CV_OK         :  653 OK  |  170 FAIL  ( 79.3%)
R_PETRO_OK        :  775 OK  |   48 FAIL  ( 94.2%)
R_EXPFIT_OK       :  808 OK  |   15 FAIL  ( 98.2%)
R_LOGFIT_OK       :  808 OK  |   15 FAIL  ( 98.2%)
H_PETRO_OK        :  492 OK  |  331 FAIL  ( 59.8%)
H_EXPFIT_OK       :  495 OK  |  328 FAIL  ( 60.1%)
H_LOGFIT_OK       :  495 OK  |  328 FAIL  ( 60.1%)
BRIGHT_STAR_FLAG  :    9 OK  |  814 FAIL  (  1.1%)
HAPY_MORPH_FLAG   :   38 OK  |  785 FAIL  (  4.6%)
R_SM_FLAG         :  215 OK  |  608 FAIL  ( 73.9%)
R_SM_SERSIC_FLAG  :  103 OK  |  720 FAIL  ( 87.5%)
H_SM_FLAG         :  380 OK  |  443 FAIL  ( 53.8%)
H_SM_SERSIC_FLAG  :  175 OK  |  648 FAIL  ( 78.7%)


PROFILES_BOTH     :  652 ( 79.2%)

STATMORPH_BOTH    :   87

STATUS counts
-------------
ok          : 809
running     : 14

STAGE counts
------------
done        : 809
mask        : 14

Runtime medians (sec)
---------------------
MASK_SEC    : 0.17
PHOT_SEC    : 8.27
SM_SEC      : 4.38
TOTAL_SEC   : 21.30

Number with bad phot = 14

```

After fixing p012-p013 craziness, one less gal with bad phot!:
```
HAPY RUN SUMMARY
----------------
Total galaxies: 821
Unique galaxies: 659

Pipeline completion
-------------------
PSF_OK            :  821 OK  |    0 FAIL  (100.0%)
MASK_OK           :  821 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  809 OK  |   12 FAIL  ( 98.5%)
HAPY_MORPH_OK     :  773 OK  |   48 FAIL  ( 94.2%)
R_PROFILE_OK      :  809 OK  |   12 FAIL  ( 98.5%)
H_PROFILE_OK      :  654 OK  |  167 FAIL  ( 79.7%)
R_SM_OK           :  729 OK  |   92 FAIL  ( 88.8%)
H_SM_OK           :  729 OK  |   92 FAIL  ( 88.8%)
GAL_NC_OK         :  653 OK  |  168 FAIL  ( 79.5%)
GAL_CV_OK         :  650 OK  |  171 FAIL  ( 79.2%)
R_PETRO_OK        :  773 OK  |   48 FAIL  ( 94.2%)
R_EXPFIT_OK       :  809 OK  |   12 FAIL  ( 98.5%)
R_LOGFIT_OK       :  809 OK  |   12 FAIL  ( 98.5%)
H_PETRO_OK        :  493 OK  |  328 FAIL  ( 60.0%)
H_EXPFIT_OK       :  492 OK  |  329 FAIL  ( 59.9%)
H_LOGFIT_OK       :  492 OK  |  329 FAIL  ( 59.9%)
BRIGHT_STAR_FLAG  :    9 OK  |  812 FAIL  (  1.1%)
HAPY_MORPH_FLAG   :   36 OK  |  785 FAIL  (  4.4%)
R_SM_FLAG         :  210 OK  |  611 FAIL  ( 74.4%)
R_SM_SERSIC_FLAG  :   99 OK  |  722 FAIL  ( 87.9%)
H_SM_FLAG         :  381 OK  |  440 FAIL  ( 53.6%)
H_SM_SERSIC_FLAG  :  179 OK  |  642 FAIL  ( 78.2%)


PROFILES_BOTH     :  654 ( 79.7%)

STATMORPH_BOTH    :   92

STATUS counts
-------------
ok          : 809
running     : 12

STAGE counts
------------
done        : 809
mask        : 12

Runtime medians (sec)
---------------------
MASK_SEC    : 0.16
PHOT_SEC    : 8.73
SM_SEC      : 4.39
TOTAL_SEC   : 22.71

Number with bad phot = 12

```
## Make qc plots

Create some basic qc plots:
```
python ~/github/hapy/scripts/qc_results.py merged_results.fits --scheme virgo
```

Output from 2026-04-17:
```
Read 823 rows from merged_results.fits
REVIEW_PRIORITY SUMMARY
{np.str_('high'): np.int64(141), np.str_('low'): np.int64(256), np.str_('medium'): np.int64(426)}
ELL_MISMATCH 269
FILTER_WARNING 79
WARN_MASK 29
BRIGHT_STAR_FLAG 9
WARN_WEAK_HA 297
Unable to revert mtime: /usr/local/share/fonts
/home/siena.edu/rfinn/github/hapy/scripts/qc_results.py:564: UserWarning: Warning: converting a masked element to nan.
  out[i] = float(v)
Wrote QC products to qc
Number of high priority in qc/tables/review = 141
writing  qc/tables/review/review_sample.csv
```

```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260417$ python ~/github/hapy/scripts/qc_results.py merged_results_virgo_20260421.fits --scheme virgo
Read 821 rows from merged_results_virgo_20260421.fits
REVIEW_PRIORITY SUMMARY
{np.str_('high'): np.int64(141), np.str_('low'): np.int64(278), np.str_('medium'): np.int64(402)}
ELL_MISMATCH 265
FILTER_WARNING 79
WARN_MASK 17
BRIGHT_STAR_FLAG 9
WARN_WEAK_HA 296
Unable to revert mtime: /usr/local/share/fonts
/home/siena.edu/rfinn/github/hapy/scripts/qc_results.py:564: UserWarning: Warning: converting a masked element to nan.
  out[i] = float(v)
Wrote QC products to qc
Number of high priority in qc/tables/review = 141
writing  qc/tables/review/review_sample.csv


REVIEW PRIORITY SUMMARY
{np.str_('high'): np.int64(141), np.str_('low'): np.int64(278), np.str_('medium'): np.int64(402)}

HIGH PRIORITY DRIVERS
NOT_PHOT_OK                 : total=  12  in_high=  12
NOT_HAPY_MORPH_OK           : total=  48  in_high=  48
BRIGHT_STAR_FLAG            : total=   9  in_high=   9
WARN_MASK                   : total=  17  in_high=  17
SEVERE_CEN_ANY              : total=  89  in_high=  89
WARN_CUTOUT_MISSING_SHAPE   : total=   0  in_high=   0

MEDIUM PRIORITY DRIVERS
ELL_MISMATCH                : total= 265  in_medium= 153
WARN_WEAK_HA                : total= 296  in_medium= 230
WARN_CEN_ANY                : total= 217  in_medium= 118
WARN_R_PROFILE_PEAK         : total= 142  in_medium=  81
WARN_CUTOUT_MISSING         : total=   0  in_medium=   0

OVERLAP AMONG HIGH DRIVERS
NOT_PHOT_OK                  & NOT_HAPY_MORPH_OK           :   12
NOT_PHOT_OK                  & WARN_MASK                   :    1
NOT_HAPY_MORPH_OK            & BRIGHT_STAR_FLAG            :    1
NOT_HAPY_MORPH_OK            & WARN_MASK                   :    4
NOT_HAPY_MORPH_OK            & SEVERE_CEN_ANY              :    6
BRIGHT_STAR_FLAG             & WARN_MASK                   :    4
BRIGHT_STAR_FLAG             & SEVERE_CEN_ANY              :    5
WARN_MASK                    & SEVERE_CEN_ANY              :    9

```



Inspecting duplicate observations:
```
usage: qc_duplicates.py [-h] [--outdir OUTDIR] [--max-ha-filter-correction MAX_HA_FILTER_CORRECTION] table

QC analysis for duplicate HAPY observations.

positional arguments:
  table                 Merged HAPY results table (e.g. merged_results.fits)

options:
  -h, --help            show this help message and exit
  --outdir OUTDIR       Output directory
  --max-ha-filter-correction MAX_HA_FILTER_CORRECTION
                        Maximum Halpha filter correction for 'good' Halpha duplicate comparisons

```

To run:
```
python ~/github/hapy/scripts/qc_duplicates.py merged_results.fits --scheme virgo
```

```
python ~/github/hapy/scripts/validate_dashboards.py merged_results.fits --sample ALL
```


# Build Webpages to Review Cutouts

## Download Legacy Images
### To copy legacy images from a prior run

Run this command from the directory that contains
e.g. `hapy-output-20260313` and `hapy-output-20260319`.
```
rsync -av hapy-output-20260330/cutouts/ hapy-output-20260417/cutouts/
--include '*/' --include 'legacy/***' --exclude '*' --exclude '*logs*'
--ignore-existing --prune-empty-dirs
```

```
rsync -av hapy-output-20260330/cutouts/ hapy-output-20260417/cutouts/
--include '*/' --include 'legacy/***' --exclude '*' --exclude '*logs*'
--ignore-existing --prune-empty-dirs
```

### To download images...
To test on one cutout:
```bash
python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir cutouts/VFID2891-UGC04559-HDI-20200225-p004/
```

### To run on full sample

```bash
find /data-pool/Halpha/hapy-output-pre2025coadds-20260516/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list.txt
```

```bash
ROOTDIR=/data-pool/Halpha/hapy-output-pre2025coadds-20260516
```
then
```bash
parallel --bar -j 2 --joblog fetch_legacy.joblog --results fetch_legacy_logs python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir "$ROOTDIR/cutouts/{}"  :::: cutout_list.txt
```
  
### to resume any failed jobs
```
parallel --resume-failed --joblog fetch_legacy.joblog \
  python ~/github/hapy/hapy/imagetools/fetch_legacy_cutouts.py \
    --cutout-dir "$ROOTDIR/cutouts/{}" \
    --layer ls-dr9 \
  :::: cutout_list.txt
```



## Build Cutout Webpages
Create a list of the cutout images:
```bash
find /data-pool/Halpha/hapy-output-20260417/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
```

Test on one directory:
```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260417
```
```bash
python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir
$ROOTDIR/cutouts --outdir $ROOTDIR/html/cutouts --oneimage VFID2891-UGC04559-HDI-20200225-p004
```


```bash 
parallel --bar -j 20 --memfree 60G --joblog build_web_cutouts.joblog --results build_web_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir $ROOTDIR/cutouts --oneimage "{}" --outdir $ROOTDIR/html/cutouts :::: cutout_list_buildwebpages.txt
```

## Build cutout index

```
python ~/github/hapy/scripts/build_cutout_index.py --help
```

```bash
python ~/github/hapy/scripts/build_cutout_index.py --runroot /data-pool/Halpha/hapy-output-20260417/ --results-table /data-pool/Halpha/hapy-output-20260417/merged_results_virgo_20260421.fits
```
```
python ~/github/hapy/scripts/build_cutout_index.py --runroot
/data-pool/Halpha/hapy-output-20260417/ --results-table
/data-pool/Halpha/hapy-output-20260417/merged_results.fits


```

## Rsync files

from html directory:
```
rsync -avz cutouts fitsxfr.siena.edu:/var/www/html/fits/virgo/.
```

# Commands for Testing a Virgo Subsample
Testing directories
```
find /data-pool/Halpha/hapy-test1/ -maxdepth 1 -type f \( -name "VF*r.fits" -o -name "VF*R.fits" \) | sort > fullpath_rcoadds_all.txt
```
```bash
cat fullpath_rcoadds_all.txt | parallel -j 32 --bar --joblog
cutouts_parallel.log get_cutouts --rimage {} --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits
--scheme virgo --psfdir /data-pool/Halpha/psf-images/
```

```
python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts
```

```
find cutouts/ -mindepth 1 -maxdepth 1 -type d ! -name "cutouts_summary" | sort > cutout_list.txt
```

```
parallel --bar  -j 16  --memfree 60G --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/ :::: cutout_list.txt
```
```
merge_results --indir cutouts --mode run_analysis
```

```
python ~/github/hapy/scripts/summarize_run.py merged_results.fits --scheme virgo
```


# Testing on laptop

```
run_analysis --make-mask  --psf-dir /Users/rfinn/research/Virgo/hatools_test/coadds-virgo-test3. --statmorph
--galfit --convflag --log-to-console --gaia-dir
/Users/rfinn/research/Virgo/hatools_test/coadds-virgo-test3/ --cutout-dir cutouts/VFID1588-NGC5169-HDI-20170523-p023
```

```
run_analysis --make-mask  --psf-dir /Users/rfinn/research/Virgo/hatools_test/coadds-virgo-test3 --gaia-dir
/Users/rfinn/research/Virgo/hatools_test/coadds-virgo-test3/ --cutout-dir cutouts/VFID1588-NGC5169-HDI-20170523-p023
```


# After completing manual masking

- I identified galaxies that needed manual masking and then created
  the masks.
- I have to rerun `run_analysis` on these galaxies

## Identify the galaxies that have a manual mask:
```bash
find cutouts -mindepth 2 -maxdepth 2 -name "*-mask-manual.fits" \
  | xargs -n1 dirname \
  | sort -u > cutouts_with_manual_masks.txt
```

Check the number:
```
wc -l cutouts_with_manual_masks.txt
head cutouts_with_manual_masks.txt
```


> [!NOTE] 
> In `virgo_qc_results`, 87 have `MASK_FIXED == YES`.  This is the same number in the file.  Check!

## 2. Rerun `run_analysis` only on these

```
parallel --bar -j 16 --memfree 60G --joblog
  run_analysis_manual_masks.joblog --results
  parallel-logs-manual-masks run_analysis --cutout-dir {} --make-mask
  --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit
  --convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/ 
  :::: cutouts_with_manual_masks.txt
```

## 3. Check Completion

```bash
awk 'NR==1 || $7 != 0' run_analysis_manual_masks.joblog
```

Find any missing:
```bash
while read d; do
  tag=$(basename "$d")
  test -f "$d/${tag}-results.ecsv" || echo "$d"
done < cutouts_with_manual_masks.txt > missing_manual_reruns.txt

wc -l missing_manual_reruns.txt
cat missing_manual_reruns.txt
```

## 4. Regenerate merged table

```bash
merge_results --indir cutouts --mode run_analysis
```

Update files used to build cutout index (I think):

```bash
python ~/github/hapy/scripts/qc_results.py merged_results.fits --scheme virgo
```

After mask updates:
```
Read 823 rows from merged_results.fits
UPDATE: adding VFINDEX
REVIEW_PRIORITY SUMMARY
{np.str_('high'): np.int64(143), np.str_('low'): np.int64(278), np.str_('medium'): np.int64(402)}
ELL_MISMATCH 267
FILTER_WARNING 79
WARN_MASK 18
BRIGHT_STAR_FLAG 9
WARN_WEAK_HA 297
Unable to revert mtime: /usr/local/share/fonts
/home/siena.edu/rfinn/github/hapy/scripts/qc_results.py:564: UserWarning: Warning: converting a masked element to nan.
  out[i] = float(v)
Wrote QC products to qc
Number of high priority in qc/tables/review = 143
writing  qc/tables/review/review_sample.csv


REVIEW PRIORITY SUMMARY
{np.str_('high'): np.int64(143), np.str_('low'): np.int64(278), np.str_('medium'): np.int64(402)}

HIGH PRIORITY DRIVERS
NOT_PHOT_OK                 : total=  13  in_high=  13
NOT_HAPY_MORPH_OK           : total=  49  in_high=  49
BRIGHT_STAR_FLAG            : total=   9  in_high=   9
WARN_MASK                   : total=  18  in_high=  18
SEVERE_CEN_ANY              : total=  89  in_high=  89
WARN_CUTOUT_MISSING_SHAPE   : total=   0  in_high=   0

MEDIUM PRIORITY DRIVERS
ELL_MISMATCH                : total= 267  in_medium= 154
WARN_WEAK_HA                : total= 297  in_medium= 230
WARN_CEN_ANY                : total= 217  in_medium= 117
WARN_R_PROFILE_PEAK         : total= 142  in_medium=  81
WARN_CUTOUT_MISSING         : total=   0  in_medium=   0

OVERLAP AMONG HIGH DRIVERS
NOT_PHOT_OK                  & NOT_HAPY_MORPH_OK           :   13
NOT_PHOT_OK                  & WARN_MASK                   :    1
NOT_HAPY_MORPH_OK            & BRIGHT_STAR_FLAG            :    1
NOT_HAPY_MORPH_OK            & WARN_MASK                   :    4
NOT_HAPY_MORPH_OK            & SEVERE_CEN_ANY              :    6
BRIGHT_STAR_FLAG             & WARN_MASK                   :    4
BRIGHT_STAR_FLAG             & SEVERE_CEN_ANY              :    5
WARN_MASK                    & SEVERE_CEN_ANY              :    9

```
## 5. Rebuild webpages

Build the file list:

```bash
find cutouts -mindepth 2 -maxdepth 2 -name "*-mask-manual.fits" | xargs -n1 dirname | xargs -n1 basename | sort -u > tags_with_manual_masks.txt
```

Rebulid webpages:
```bash 
parallel --bar -j 20 --memfree 60G --joblog build_web_cutouts_manual_masks.joblog --results build_web_manual_masks_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir $ROOTDIR/cutouts --oneimage "{}" --outdir $ROOTDIR/html/cutouts :::: tags_with_manual_masks.txt
```


Rebuild cutout index:
```bash
python ~/github/hapy/scripts/build_cutout_index.py --runroot /data-pool/Halpha/hapy-output-20260417/ --results-table /data-pool/Halpha/hapy-output-20260417/merged_results_virgo_20260501.fits
```


# Create cs-gr images


## Reproject Legacy images

```bash
find cutouts -mindepth 1 -maxdepth 1 -type d -name 'VFID*' | sort > reproject_cutout_list.txt
```


```bash
parallel --bar -j 20 --results legacy_reproject_logs python ~/github/hapy/scripts/make_legacy_reprojections.py "{}" :::: reproject_cutout_list.txt
```

## Then make CS-gr images

Test one image:
```bash
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```


```bash
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID0481-NGC6307-INT-20190602-p010
```

```
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID0481-NGC6307-INT-20190602-p010 --auto-contscale
--auto-contscale-percentile 30 --overwrite
```

To solve by fitting the low end of the tail
```
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID1934-NGC2799-INT-20190205-p026 --auto-contscale --auto-contscale-method negtail --overwrite
```

To solve using percentile of ratio of R/Halpha flux:
```
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID1934-NGC2799-INT-20190205-p026 --auto-contscale --auto-contscale-method ratio --auto-contscale-percentile 30 --overwrite
```

```bash
parallel --bar -j 16 --joblog csgr.joblog --results csgr_logs python
~/github/hapy/hapy/scripts/make_cs_gr.py "{}"  :::: reproject_cutout_list.txt
```

```
parallel --bar -j 16 --joblog cs_gr_auto.joblog --results
cs_gr_auto_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {}
--auto-contscale --auto-contscale-percentile 30 --overwrite ::::
reproject_cutout_list.txt
```


Needed to rerun on the INT images (filter lookup issue)

```
find cutouts -mindepth 1 -maxdepth 1 -type d -name 'VFID*INT*' | sort > int_cutout_list.txt
```

```
parallel --bar -j 16 --joblog cs_gr_int.joblog --results cs_gr_int_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {} :::: int_cutout_list.txt
```

## Check failures

```bash
grep -R "FAILED\|Traceback\|ERROR\|can't find\|problem getting" logs_legacy_reproject logs_cs_gr
```


```bash
parallel --bar  -j 16  --memfree 60G --joblog csgr.joblog --results csgr-logs ~/github/hapy/hapy/scripts/make_cs_gr.fits "{}" :::: cutout_with_dir.txt
```

## Visualize CS Images and Select Best Duplicate

To just write the duplicates table:
```
python ~/github/hapy/scripts/inspect_cs_images.py make-table merged_results.fits --outdir cs_image_inspection --min-dups 1
```

To create an input list for running in parallel:
```
python ~/github/hapy/scripts/inspect_cs_images.py list-groups cs_image_inspection/cs_image_inspection_groups.ecsv > cs_group_list.txt
```
To test one:
```
python ~/github/hapy/scripts/inspect_cs_images.py plot-one cs_image_inspection/cs_image_inspection_groups.ecsv {} --cutout-dir
cutouts --outdir cs_image_inspection :::: cs_group_list.txt
```

To build the plots in parallel:
```
parallel --bar -j 16 --joblog cs_image_plot.joblog --results
cs_image_plot_logs python ~/github/hapy/scripts/inspect_cs_images.py plot-one
cs_image_inspection/cs_image_inspection_groups.ecsv {} --cutout-dir
cutouts --outdir cs_image_inspection :::: cs_group_list.txt
```


