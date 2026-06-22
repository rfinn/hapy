
# Creating Simple Weights for all INT coadds
you only need to do this once for the coadd directory
```
mkdir ORIGINAL_INT_WEIGHTS
```

```
mv *INT*weight.fits ORIGINAL_INT_WEIGHTS/.
```

create a list of all INT coadds:
```
ls VF*INT*r.fits VF*INT*Halpha.fits VF*INT*Ha6657.fits > INT_all_coadds.txt
```


```
parallel --bar -j 16 --joblog make_simple_weight_INT.joblog --results make_simple_weight_logs python ~/github/hapy/scripts/make_simple_weight_from_coadd.py {} :::: INT_all_coadds.txt
```


# Download the gaia files
Move to directory containing coadded images

```bash
cd /data-pool/Halpha/coadds-v20260330
```

Then start the download, if catalogs are needed.

```bash
python ~/github/hapy/scripts/download_gaia_coadd_catalogs.py 
```

# Prepare the output directory
### make a directory for analysis run

```bash
cd /data-pool/Halpha/
```

```bash
mkdir hapy-output-20260612
```

```bash
cd hapy-output-20260612
```


### make a list of coadds to analyze

```bash
find /data-pool/Halpha/coadds-v20260330/ -maxdepth 1 -type f \( -name "VF*r.fits" -o -name "VF*R.fits" \) | sort > fullpath_rcoadds_all.txt
```
This contains 226 coadds as of 2026-06-12.


For hybrid coadds:
```bash
find /data-pool/Halpha/coadds-v20260609/ -maxdepth 1 -type f \( -name "VF*r.fits" -o -name "VF*R.fits" \) | sort > fullpath_rcoadds_all.txt
```
This contains 223 coadds as of 2026-06-09.



If you have coadds that are still under review, make a copy of the
full list and remove any coadds that are not ready.
```bash
cp fullpath_rcoadds_all.txt fullpath_rcoadds_hapy_ready.txt
```
then remove lines

As of 2026-Mar-09:
```bash
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260309$ wc -l *.txt
  227 fullpath_rcoadds_all.txt
  211 fullpath_rcoadds_hapy_ready.txt
  438 total
```

As of 2026-Mar-30:
- all coadds are ready!
# Make Cutouts

### Test on One Image

#### New Coadds (post Dec 2025)
```bash
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits \
--scheme virgo --maxcorrection 5 --psfdir \
/data-pool/Halpha/psf-images-v20260330/ --rimage \
/data-pool/Halpha/coadds-v20260330/VF-126.291+27.988-HDI-20180313-p004-R.fits
```
This should create two folders in the cutouts/ directory.

```
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits \
--scheme virgo --maxcorrection 5 --psfdir \
/data-pool/Halpha/psf-images-v20260518/ --rimage \
/data-pool/Halpha/coadds-v20260330/VF-177.143+56.083-INT-20220502-VFID0957-r.fits
```
This says there are 6 galaxies in FOV, but one is skipped b/c it is
outside usable area:
```
Skipping VFID0974-PGC12807144NED001: invalid cutout region (r_invalid); ra=177.526207,dec=55.782624
```

#### Hybrid coadds
For hybrid coadds (I think):
```bash
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits \
--scheme virgo --maxcorrection 5 --psfdir \
/data-pool/Halpha/psf-images-v20260518/ --rimage \
/data-pool/Halpha/coadds-v20260609/VF-126.291+27.988-HDI-20180313-p004-R.fits
```

```
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits \
--scheme virgo --maxcorrection 5 --psfdir \
/data-pool/Halpha/psf-images-v20260518/ --rimage \
/data-pool/Halpha/coadds-v20260609/VF-177.200+56.055-INT-20220502-VFID0957-r.fits
```

### Run on full Virgo sample

If not all coadds are ready for `hapy`,
```bash
cat fullpath_rcoadds_hapy_ready.txt | parallel -j 16 --bar --joblog cutouts_parallel.log get_cutouts --rimage {} --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5
```

Or run on the full list of coadds.
```bash
parallel -j 16 --bar --joblog \
cutouts_parallel.log get_cutouts --rimage {} --catalog \
~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo \
--maxcorrection 5 --overwrite --overwrite_metadata :::: fullpath_rcoadds_all.txt
```

### First check for any failures

```
awk 'NR==1 || $7 != 0 {print}' cutouts_parallel.log
```

Example output for hybrid coadd set:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260519$ awk 'NR==1 || $7 != 0 {print}' cutouts_parallel.log
Seq	Host	Starttime	JobRuntime	Send	Receive	Exitval	Signal	Command
69	:	1779219600.267	    12.449	0	1193	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-162.760+32.934-INT-20190205-p065-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5
156	:	1779219651.982	    31.741	0	2833	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-208.804+05.187-INT-20190206-p120-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5
179	:	1779219671.326	    45.526	0	1193	1	0	get_cutouts --rimage /data-pool/Halpha/coadds-v20260518/VF-221.102+01.782-INT-20190209-p149-r.fits --catalog /home/siena.edu/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5

```

### Check output

```bash
python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts
```

Example output:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260310$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_hapy_ready.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              211
Cutout directories:        768
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0

```

From 2026-03-30:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260330$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              227
Cutout directories:        823
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0

```

From 2025-05-17:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260517$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              226
Cutout directories:        853
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

From 2026-06-09, first pass:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260609-hybrid$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              223
Cutout directories:        856
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```
I love how the numbers don't match the previous run...


Lastest run on new coadds 2026-06-12:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260612$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              226
Cutout directories:        853
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

2026-06-20 after implementing min size in `get_cutouts`:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260620$ python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts

CUTOUT SUMMARY
--------------
Input coadds:              226
Cutout directories:        849
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0

```
To test one:

```
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images/ --rimage /data-pool/Halpha/coadds-v20260609/VF-178.160+52.290-INT-20220503-VFID1213-r.fits
```

## Merge get_cutouts tables
```bash
merge_results --mode get_cutouts --indir cutouts_summary --out merged_cutouts_results.fits
```

Example output:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260310$ merge_results --mode get_cutouts --indir cutouts_summary --out merged_cutouts_results.fits
Searching for files  cutouts_summary*.ecsv
Found 211 result files.
Reading tables...
Validating schema...
	validated 211/211 tables
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260310/merged_cutouts_results.fits
Done.
Final table rows: 768
Final table columns: 19

```

And after 20260330:
```
Searching for files  cutouts_summary*.ecsv
Found 227 result files.
Reading tables...
Validating schema...
	validated 227/227 tables
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260417/merged_cutouts_results.fits
Done.
Final table rows: 831
Final table columns: 28
```

And 20260517 (same for 2026-06-12):
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260517$ merge_results --mode get_cutouts --indir cutouts_summary --out merged_cutouts_results.fits
Searching for files  cutouts_summary*.ecsv
Found 226 result files.
Merging 226 result files.
Reading tables...
Validating schema...
	validated 226/226 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260517/merged_cutouts_results.fits
Done.
Final table rows: 853
Final table columns: 28
```

And after 20260609:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260609-hybrid$ merge_results --mode get_cutouts --indir cutouts_summary --out merged_cutouts_results.fits
Searching for files  cutouts_summary*.ecsv
Found 223 result files.
Merging 223 result files.
Reading tables...
Validating schema...
	validated 223/223 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260609-hybrid/merged_cutouts_results.fits
Done.
Final table rows: 860
Final table columns: 28
```

# Create cs-gr images

## Download Legacy Images
### To copy legacy images from a prior run

Run this command from the directory that contains the hapy output
directories (this is most likely `/data-pool/Halpha`).
```
cd /data-pool/Halpha/
```

Then do a test run with rsync:
```
rsync -av hapy-output-20260417/cutouts/ \
hapy-output-20260612/cutouts/ --include '*/' --include 'legacy/***' \
--exclude '*' --exclude '*logs*' --ignore-existing --prune-empty-dirs --dry-run
```

```
rsync -av hapy-output-20260612/cutouts/ \
hapy-output-20260620/cutouts/ --include '*/' --include 'legacy/***' \
--exclude '*' --exclude '*logs*' --ignore-existing --prune-empty-dirs --dry-run
```

#### Hybrid coadds
```
rsync -av hapy-output-20260519-hybrid/cutouts/ \
hapy-output-20260609-hybrid/cutouts/ --include '*/' --include 'legacy/***' \
--exclude '*' --exclude '*logs*' --ignore-existing --prune-empty-dirs --dry-run
```

If all looks good, remove the `--dry-run` flag and run again.

#### Older versions

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
find /data-pool/Halpha/hapy-output-20260310/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list.txt
```

```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260310
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


### Find cutouts with missing

```
python ~/github/hapy/scripts/find_missing_legacy_cutouts.py
```


For 2026-06-20, got this (same cutouts that are missing - this is an
issue with the parent coadd)

```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260620$ python ~/github/hapy/scripts/find_missing_legacy_cutouts.py
MISSING legacy files: VFID3297-SDSSJ091525.83+252511.1-INT-20190208-p024 (g=False, r=False, z=False)
MISSING legacy files: VFID6386-WISEAJ150825.51+014224.3-INT-20190601-p033 (g=False, r=False, z=False)
MISSING legacy files: VFID6392-WISEAJ150819.90+014123.6-INT-20190601-p033 (g=False, r=False, z=False)
MISSING legacy files: VFID6447-SDSSJ150812.35+012959.7-INT-20190601-p033 (g=False, r=False, z=False)
MISSING legacy files: VFID6463-WISEAJ150809.13+012516.6-INT-20190601-p033 (g=False, r=False, z=False)
```

```
parallel --bar -j 2 \
  --joblog fetch_legacy_retry.joblog \
  python ~/github/hapy/scripts/fetch_legacy_cutouts.py \
  --cutout-dir {} \
  :::: missing_legacy_cutouts.txt
```

## Reproject Legacy images

```bash
find cutouts -mindepth 1 -maxdepth 1 -type d -name 'VFID*' | sort > reproject_cutout_list.txt
```

test one:
```
python ~/github/hapy/scripts/make_legacy_reprojections.py cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```


To run all:
```bash
parallel --bar -j 20 --results legacy_reproject_logs python ~/github/hapy/scripts/make_legacy_reprojections.py "{}" :::: reproject_cutout_list.txt
```

Or just to run on INT images
```bash
parallel --bar -j 20 --results legacy_reproject_logs python \
~/github/hapy/scripts/make_legacy_reprojections.py "{}" --overwrite :::: reproject_cutout_list_INT.txt
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

```bash
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID0569-NGC5989-BOK-20220424-VFID0607/ --auto-contscale --auto-contscale-method ratio --auto-contscale-percentile 30 --overwrite
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
parallel --bar -j 16 --joblog cs_gr_auto.joblog --results \
cs_gr_auto_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {} \
--auto-contscale --auto-contscale-percentile 30 --overwrite :::: \
reproject_cutout_list.txt
```

To run without overwriting an existing cs-gr images:
```
parallel --bar -j 16 --joblog cs_gr_auto.joblog --results \
cs_gr_auto_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {} \
--auto-contscale --auto-contscale-percentile 30 :::: \
reproject_cutout_list.txt
```

#### Hybrid coadds
Needed to rerun on the INT images (filter lookup issue)

```
find cutouts -mindepth 1 -maxdepth 1 -type d -name 'VFID*INT*' | sort > int_cutout_list.txt
```

```
parallel --bar -j 16 --joblog cs_gr_int.joblog --results cs_gr_int_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {} :::: int_cutout_list.txt
```

## Check failures

```
python ~/github/hapy/scripts/check_for_missing_csgr.py
```

then rerun `run_analysis`
```
parallel --bar -j 16 --joblog rerun_csgr.joblog --results rerun_csgr_logs run_analysis --cutout-dir {} --csgr :::: missing_csgr_or_phot.txt
```

Or 
```bash
grep -R "FAILED\|Traceback\|ERROR\|can't find\|problem getting" logs_legacy_reproject logs_cs_gr
```


```bash
parallel --bar  -j 16  --memfree 60G --joblog csgr.joblog --results csgr-logs ~/github/hapy/hapy/scripts/make_cs_gr.fits "{}" :::: cutout_with_dir.txt
```

# Make Mstar Images

```
python ~/github/hapy/hapy/scripts/make_mstar_map.py cutouts/VFIDxxxx-... --scheme virgo --overwrite
```

To run all:

```
find cutouts -mindepth 1 -maxdepth 1 -type d | sort > cutout_list.txt
```


```
parallel --bar -j 16 python ~/github/hapy/hapy/scripts/make_mstar_map.py {} --overwrite :::: cutout_list.txt
```

# Make SFR Images
```
python ~/github/hapy/hapy/scripts/make_sfr_map.py cutouts/VFIDxxxx-... --scheme virgo --overwrite
```

```
find cutouts -mindepth 1 -maxdepth 1 -type d | sort > cutout_list.txt
```

```
parallel --bar -j 16 --joblog make_sfr_map.joblog python ~/github/hapy/hapy/scripts/make_sfr_map.py {} --scheme virgo --overwrite :::: cutout_list.txt
```

# Copy manual masks

```
rsync -avzn \
  --include='*/' \
  --include='*-mask-manual.fits' \
  --exclude='*' \
  /data-pool/Halpha/hapy-output-20260417/cutouts/ \
  /data-pool/Halpha/hapy-output-20260612/cutouts/
```

If that looks good, then remove `n` to exit `DRY_RUN` mode:

```
rsync -avz \
  --include='*/' \
  --include='*-mask-manual.fits' \
  --exclude='*' \
  /data-pool/Halpha/hapy-output-20260417/cutouts/ \
  /data-pool/Halpha/hapy-output-20260612/cutouts/
```

For 2026-06-20, I also finished making adjustments to any galaxies I
flagged as needing more editing on mask.

# Run Analysis 

### Run on One Cutout


```bash
run_analysis --make-mask  --psf-dir \
/data-pool/Halpha/psf-images-v20260518/ --statmorph \
--galfit --convflag --log-to-console --gaia-dir \
/data-pool/Halpha/coadds-v20260330/gaia_catalogs/ --cutout-dir \
cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```


```bash
run_analysis --make-mask  --psf-dir /data-pool/Halpha/psf-images-v20260330/ --statmorph \
--galfit --convflag --csgr --log-to-console --gaia-dir \
/data-pool/Halpha/coadds-v20260330/gaia_catalogs/ --cutout-dir \
cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```

### Hybrid Coadd sample
```
run_analysis --make-mask  --psf-dir /data-pool/Halpha/psf-images-v20260518/ --statmorph \
--galfit --convflag --log-to-console --gaia-dir \
/data-pool/Halpha/coadds-v20260518/gaia_catalogs/ --cutout-dir \
cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```



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

Test on 20 galaxies:
```bash
head -20 cutout_list.txt | parallel --bar -j 4 --joblog run_analysis.joblog \
  --results parallel-logs \
  run_analysis --cutout-dir "{}" --make-mask \
  --psf-dir /data-pool/Halpha/psf-images/ \
  --statmorph --galfit --convflag \
  --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/
```

Run the next 20 galaxies:
```bash
sed -n '21,40p' cutout_list.txt | parallel --bar -j 8 --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/
```

### Run in Parallel

Copy results from manual inspection that will exclude bad targets:
```bash
cp ../hapy-output-20260417/review_sample_20260514.csv .
```

Create an input list that removes objects that have `CATALOG_USE==EXCLUDE`:
```
python ~/github/hapy/scripts/make_run_analysis_list.py --cutout-dir \
cutouts --review review_sample_20260514.csv --outfile \
cutout_run_analysis_list.txt
```


If manual masking has been done, don't force rebuilding of masks:
```bash
parallel --bar  -j 16  --memfree 30G --joblog run_analysis.joblog \
--results parallel-logs run_analysis --cutout-dir "{}" --make-mask \
--psf-dir /data-pool/Halpha/psf-images-v20260330/ --statmorph --galfit \
--convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/ :::: cutout_with_dir.txt
```


#### Run with cs-gr enabled and manual masking

```
run_analysis --make-mask  --psf-dir /data-pool/Halpha/psf-images-v20260330/ --statmorph \
--galfit --convflag --csgr --log-to-console --gaia-dir \
/data-pool/Halpha/coadds-v20260330/gaia_catalogs/ --cutout-dir \
cutouts/VFID0377-IC1210-BOK-20210414-VFID0422
```

```bash
parallel --bar  -j 16  --memfree 30G --joblog run_analysis.joblog \
--results parallel-logs run_analysis --cutout-dir "{}" --make-mask \
--psf-dir /data-pool/Halpha/psf-images-v20260330/ --statmorph --galfit \
--csgr \
--convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/ :::: cutout_with_dir.txt
```


#### Rerunning INT cutouts

To rerun INT cutouts only (because I gave the wrong psf dir :( )

```
find cutouts -mindepth 1 -maxdepth 1 -type d -name "*INT*" | sort > INT_cutouts.txt
```

```
parallel --bar -j 16 --joblog run_analysis_INT_rerun.joblog \
--memfree 30G --results run_analysis_INT_logs \
run_analysis --cutout-dir "{}" --make-mask \
--psf-dir /data-pool/Halpha/psf-images-v20260330/ --statmorph --galfit \
--convflag --gaia-dir \
/data-pool/Halpha/coadds-v20260330/gaia_catalogs/ :::: INT_cutouts.txt
```


#### Other special cases

```bash
parallel --bar  -j 16  --memfree 30G --joblog \
run_analysis.missing_csgr.joblog \
--results parallel-logs-csgr-missing run_analysis --cutout-dir "{}" \
--make-mask --csgr --psf-dir /data-pool/Halpha/psf-images/ --statmorph \
--galfit --convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/ :::: cutout_run_analysis_list.txt
```

For rerunning individual  galaxies:

```
run_analysis --make-mask --csgr --psf-dir \
/data-pool/Halpha/psf-images/ --statmorph --galfit \
--convflag --gaia-dir \
/data-pool/Halpha/coadds-v20260330/gaia_catalogs/ --cutout-dir cutouts/VFID6463-WISEAJ150809.13+012516.6-INT-20190601-p033
```

and for :
```
cutouts/VFID6463-WISEAJ150809.13+012516.6-INT-20190601-p033
cutouts/VFID2766-WISEAJ133704.60+315337.9-HDI-20170522-p006
cutouts/VFID2745-UGC08602-HDI-20170522-p006
```

#### Tracking progress
You can monitor memory usage with 
```
htop
```


## Merge Results From `run_analysis`

```bash
merge_results --indir cutouts --mode run_analysis --review-csv review_sample_20260612.csv
```

Example output as of 2026-04-17:
```
Searching for files  *results.ecsv
Found 823 result files.
Reading tables...
Validating schema...
	validated 823/823 tables
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260417/merged_results.fits
Done.
Final table rows: 823
Final table columns: 344
```

Example output as of 2026-05-14:
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

Example output as of 2026-05-16
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260517$ merge_results --indir cutouts --mode run_analysis --review-csv review_sample_20260514.csv
Searching for files  *results.ecsv
Found 814 result files.
Skipping 0 files with CATALOG_USE == EXCLUDE
Merging 814 result files.
Reading tables...
Validating schema...
	validated 814/814 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260517/merged_results_virgo_20260518.fits
Done.
Final table rows: 814
Final table columns: 390
```

Hybrid sample as of 2026-05-16
```
100% 767:0=0s cutouts/VFID6634-WISEAJ143246.75+000131.3-INT-20220503-VFID6620                                                                                                                          
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260518$ merge_results --indir cutouts --mode run_analysis --review-csv review_sample_20260514.csv
Searching for files  *results.ecsv
Found 764 result files.
Skipping 0 files with CATALOG_USE == EXCLUDE
Merging 764 result files.
Reading tables...
Validating schema...
	validated 764/764 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260518/merged_results_virgo_20260519.fits
Done.
Final table rows: 764
Final table columns: 390
```

New sample as of 2026-06-12:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260612$ merge_results --indir cutouts --mode run_analysis
Searching for files  *results.ecsv
Found 853 result files.
Merging 853 result files.
Reading tables...
Validating schema...
/data-pool/github/hapy/hapy/scripts/merge_results.py:201: RuntimeWarning: invalid value encountered in cast
  tab[col] = np.array(tab[col], dtype=ref_dtype)
WARNING: adding 1 missing columns to cutouts/VFID3714-GALEXASCJ163239.29+194733.3-BOK-20220425-VFID3714/VFID3714-GALEXASCJ163239.29+194733.3-BOK-20220425-VFID3714-results.ecsv
WARNING: adding 1 missing columns to cutouts/VFID3719-GALEXASCJ163419.44+194311.5-BOK-20220425-VFID3714/VFID3719-GALEXASCJ163419.44+194311.5-BOK-20220425-VFID3714-results.ecsv
WARNING: adding 1 missing columns to cutouts/VFID4063-SDSSJ153524.36+161944.7-INT-20220502-VFID4037/VFID4063-SDSSJ153524.36+161944.7-INT-20220502-VFID4037-results.ecsv
WARNING: adding 1 missing columns to cutouts/VFID5868-SDSSJ135615.94+050954.0-HDI-20180313-p056/VFID5868-SDSSJ135615.94+050954.0-HDI-20180313-p056-results.ecsv
WARNING: adding 1 missing columns to cutouts/VFID6454-LEDA4609875-BOK-20210418-VFID6406/VFID6454-LEDA4609875-BOK-20210418-VFID6406-results.ecsv
WARNING: adding 1 missing columns to cutouts/VFID6454-LEDA4609875-MOS-20110404-NGC5846_02/VFID6454-LEDA4609875-MOS-20110404-NGC5846_02-results.ecsv
	validated 853/853 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260612/merged_results_virgo_20260613.fits
Done.
Final table rows: 853
Final table columns: 491

```

Sample as of 2026-06-20:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260620$ merge_results --indir cutouts --mode run_analysis --review-csv review_sample_20260612.csv 
Searching for files  *results.ecsv
Found 849 result files.
Skipping 52 files with CATALOG_USE == EXCLUDE
Merging 797 result files.
Reading tables...
Validating schema...
/data-pool/github/hapy/hapy/scripts/merge_results.py:201: RuntimeWarning: invalid value encountered in cast
  tab[col] = np.array(tab[col], dtype=ref_dtype)
	validated 797/797 tables
Normalizing string columns...
Stacking tables...
Writing merged table → /data-pool/Halpha/hapy-output-20260620/merged_results_virgo_20260622.fits
Done.
Final table rows: 797
Final table columns: 499
```
854 cutout directories, so 4 didn't run.  probably the same INT
cutouts that didn't actually have cutouts b/c of weird zeros and
clipping in the coadded images.



### Summarize statistics in  `merged_results.py`
```
python ~/github/hapy/scripts/summarize_run.py --infile merged_results_virgo_20260519.fits --scheme virgo
```


Hybrid sample:
```
HAPY RUN SUMMARY
----------------
Total galaxies: 764
Unique galaxies: 616

Pipeline completion
-------------------
PSF_OK            :  481 OK  |  283 FAIL  ( 63.0%)
MASK_OK           :  764 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  760 OK  |    4 FAIL  ( 99.5%)
HAPY_MORPH_OK     :  760 OK  |    4 FAIL  ( 99.5%)
R_PROFILE_OK      :  760 OK  |    4 FAIL  ( 99.5%)
H_PROFILE_OK      :  418 OK  |  346 FAIL  ( 54.7%)
R_SM_OK           :  671 OK  |   93 FAIL  ( 87.8%)
H_SM_OK           :  671 OK  |   93 FAIL  ( 87.8%)
GAL_NC_OK         :  585 OK  |  179 FAIL  ( 76.6%)
GAL_CV_OK         :  384 OK  |  380 FAIL  ( 50.3%)
R_PETRO_OK        :  733 OK  |   31 FAIL  ( 95.9%)
R_EXPFIT_OK       :  760 OK  |    4 FAIL  ( 99.5%)
R_LOGFIT_OK       :  760 OK  |    4 FAIL  ( 99.5%)
H_PETRO_OK        :  380 OK  |  384 FAIL  ( 49.7%)
H_EXPFIT_OK       :  362 OK  |  402 FAIL  ( 47.4%)
H_LOGFIT_OK       :  362 OK  |  402 FAIL  ( 47.4%)
BRIGHT_STAR_FLAG  :    6 OK  |  758 FAIL  (  0.8%)
HAPY_MORPH_FLAG   :   35 OK  |  729 FAIL  (  4.6%)
R_SM_FLAG         :  222 OK  |  542 FAIL  ( 70.9%)
R_SM_SERSIC_FLAG  :   95 OK  |  669 FAIL  ( 87.6%)
H_SM_FLAG         :  388 OK  |  376 FAIL  ( 49.2%)
H_SM_SERSIC_FLAG  :  162 OK  |  602 FAIL  ( 78.8%)
CSGR_PHOT_OK      :    0 OK  |  764 FAIL  (  0.0%)
CSGR_HAPY_MORPH_OK:    0 OK  |  764 FAIL  (  0.0%)
CSGR_HAPY_MORPH_FLAG:    0 OK  |  764 FAIL  (  0.0%)


PROFILES_BOTH     :  418 ( 54.7%)

STATMORPH_BOTH    :  100

STATUS counts
-------------
ok          : 760
running     : 4

STAGE counts
------------
done        : 760
mask        : 4

Runtime medians (sec)
---------------------
MASK_SEC    : 0.50
PHOT_SEC    : 9.53
SM_SEC      : 4.83
TOTAL_SEC   : 23.17

```

Results from 2026-04-17:

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

#### 2026-06-12 Run:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260612$ python ~/github/hapy/scripts/summarize_run.py --infile merged_results_virgo_20260613.fits --scheme virgo

HAPY RUN SUMMARY
----------------
Total galaxies: 853
Unique galaxies: 675

Pipeline completion
-------------------
PSF_OK            :  853 OK  |    0 FAIL  (100.0%)
MASK_OK           :  853 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  847 OK  |    6 FAIL  ( 99.3%)
HAPY_MORPH_OK     :  844 OK  |    9 FAIL  ( 98.9%)
R_PROFILE_OK      :  847 OK  |    6 FAIL  ( 99.3%)
H_PROFILE_OK      :  780 OK  |   73 FAIL  ( 91.4%)
R_SM_OK           :  757 OK  |   96 FAIL  ( 88.7%)
H_SM_OK           :  757 OK  |   96 FAIL  ( 88.7%)
GAL_NC_OK         :  264 OK  |  589 FAIL  ( 30.9%)
GAL_CV_OK         :  261 OK  |  592 FAIL  ( 30.6%)
R_PETRO_OK        :  799 OK  |   54 FAIL  ( 93.7%)
R_EXPFIT_OK       :  847 OK  |    6 FAIL  ( 99.3%)
R_LOGFIT_OK       :  847 OK  |    6 FAIL  ( 99.3%)
H_PETRO_OK        :  577 OK  |  276 FAIL  ( 67.6%)
H_EXPFIT_OK       :  572 OK  |  281 FAIL  ( 67.1%)
H_LOGFIT_OK       :  572 OK  |  281 FAIL  ( 67.1%)
BRIGHT_STAR_FLAG  :    0 OK  |  853 FAIL  (  0.0%)
HAPY_MORPH_FLAG   :   89 OK  |  764 FAIL  ( 10.4%)
R_SM_FLAG         :  299 OK  |  554 FAIL  ( 64.9%)
R_SM_SERSIC_FLAG  :  107 OK  |  746 FAIL  ( 87.5%)
H_SM_FLAG         :  371 OK  |  482 FAIL  ( 56.5%)
H_SM_SERSIC_FLAG  :  188 OK  |  665 FAIL  ( 78.0%)
CSGR_PHOT_OK      :  847 OK  |    6 FAIL  ( 99.3%)
CSGR_HAPY_MORPH_OK:  846 OK  |    7 FAIL  ( 99.2%)
CSGR_HAPY_MORPH_FLAG:   27 OK  |  826 FAIL  (  3.2%)
CSGR_H_EXPFIT_OK  :  525 OK  |  328 FAIL  ( 61.5%)
CSGR_H_LOGFIT_OK  :  525 OK  |  328 FAIL  ( 61.5%)
CSGR_H_PETRO_OK   :  535 OK  |  318 FAIL  ( 62.7%)
CSGR_H_PROFILE_OK :  810 OK  |   43 FAIL  ( 95.0%)
CSGR_R_EXPFIT_OK  :  847 OK  |    6 FAIL  ( 99.3%)
CSGR_R_LOGFIT_OK  :  847 OK  |    6 FAIL  ( 99.3%)
CSGR_R_PETRO_OK   :  799 OK  |   54 FAIL  ( 93.7%)
CSGR_R_PROFILE_OK :  847 OK  |    6 FAIL  ( 99.3%)


PROFILES_BOTH     :  780 ( 91.4%)

STATMORPH_BOTH    :  224

STATUS counts
-------------
running     : 513
ok          : 340

STAGE counts
------------
phot        : 507
done        : 340
init        : 6

Runtime medians (sec)
---------------------
MASK_SEC    : 0.00
PHOT_SEC    : 9.04
SM_SEC      : 4.16
TOTAL_SEC   : 0.00

Number with bad phot = 6
```

A lot of galfit failures.  Need to revisit this.


#### 2026-06-20 Run:
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260620$ python ~/github/hapy/scripts/summarize_run.py --infile merged_results_virgo_20260622.fits --scheme virgo

HAPY RUN SUMMARY
----------------
Total galaxies: 797
Unique galaxies: 641

Pipeline completion
-------------------
PSF_OK            :  797 OK  |    0 FAIL  (100.0%)
MASK_OK           :  797 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  797 OK  |    0 FAIL  (100.0%)
HAPY_MORPH_OK     :  795 OK  |    2 FAIL  ( 99.7%)
R_PROFILE_OK      :  797 OK  |    0 FAIL  (100.0%)
H_PROFILE_OK      :  753 OK  |   44 FAIL  ( 94.5%)
R_SM_OK           :  727 OK  |   70 FAIL  ( 91.2%)
H_SM_OK           :  727 OK  |   70 FAIL  ( 91.2%)
GAL_NC_OK         :  615 OK  |  182 FAIL  ( 77.2%)
GAL_CV_OK         :  610 OK  |  187 FAIL  ( 76.5%)
R_PETRO_OK        :  781 OK  |   16 FAIL  ( 98.0%)
R_EXPFIT_OK       :  797 OK  |    0 FAIL  (100.0%)
R_LOGFIT_OK       :  797 OK  |    0 FAIL  (100.0%)
H_PETRO_OK        :  517 OK  |  280 FAIL  ( 64.9%)
H_EXPFIT_OK       :  514 OK  |  283 FAIL  ( 64.5%)
H_LOGFIT_OK       :  514 OK  |  283 FAIL  ( 64.5%)
BRIGHT_STAR_FLAG  :    0 OK  |  797 FAIL  (  0.0%)
HAPY_MORPH_FLAG   :   67 OK  |  730 FAIL  (  8.4%)
R_SM_FLAG         :  275 OK  |  522 FAIL  ( 65.5%)
R_SM_SERSIC_FLAG  :   81 OK  |  716 FAIL  ( 89.8%)
H_SM_FLAG         :  353 OK  |  444 FAIL  ( 55.7%)
H_SM_SERSIC_FLAG  :  145 OK  |  652 FAIL  ( 81.8%)
CSGR_PHOT_OK      :  797 OK  |    0 FAIL  (100.0%)
CSGR_HAPY_MORPH_OK:  796 OK  |    1 FAIL  ( 99.9%)
CSGR_HAPY_MORPH_FLAG:   12 OK  |  785 FAIL  (  1.5%)
CSGR_H_EXPFIT_OK  :  566 OK  |  231 FAIL  ( 71.0%)
CSGR_H_LOGFIT_OK  :  566 OK  |  231 FAIL  ( 71.0%)
CSGR_H_PETRO_OK   :  558 OK  |  239 FAIL  ( 70.0%)
CSGR_H_PROFILE_OK :  789 OK  |    8 FAIL  ( 99.0%)
CSGR_R_EXPFIT_OK  :  797 OK  |    0 FAIL  (100.0%)
CSGR_R_LOGFIT_OK  :  797 OK  |    0 FAIL  (100.0%)
CSGR_R_PETRO_OK   :  781 OK  |   16 FAIL  ( 98.0%)
CSGR_R_PROFILE_OK :  797 OK  |    0 FAIL  (100.0%)


PROFILES_BOTH     :  753 ( 94.5%)

STATMORPH_BOTH    :  196

STATUS counts
-------------
ok          : 797

STAGE counts
------------
done        : 797

Runtime medians (sec)
---------------------
MASK_SEC    : 0.11
PHOT_SEC    : 9.19
SM_SEC      : 4.18
TOTAL_SEC   : 33.43

Number with bad phot = 0

```


## Make qc plots

Create some basic qc plots:
```
python ~/github/hapy/scripts/qc_results.py --scheme virgo merged_results.fits 
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

#### 2026-06-12 Run
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260612$ python ~/github/hapy/scripts/qc_results.py --scheme virgo merged_results_virgo_20260613.fits 
Read 853 rows from merged_results_virgo_20260613.fits
UPDATE: adding VFINDEX
REVIEW_PRIORITY SUMMARY
{np.str_('high'): np.int64(101), np.str_('low'): np.int64(345), np.str_('medium'): np.int64(407)}
ELL_MISMATCH 280
FILTER_WARNING 81
WARN_MASK 31
BRIGHT_STAR_FLAG 0
WARN_WEAK_HA 236
/home/siena.edu/rfinn/github/hapy/scripts/qc_results.py:564: UserWarning: Warning: converting a masked element to nan.
  out[i] = float(v)
Wrote QC products to qc
Number of high priority in qc/tables/review = 101
writing  qc/tables/review/review_sample.csv


REVIEW PRIORITY SUMMARY
{np.str_('high'): np.int64(101), np.str_('low'): np.int64(345), np.str_('medium'): np.int64(407)}

HIGH PRIORITY DRIVERS
NOT_PHOT_OK                 : total=   6  in_high=   6
NOT_HAPY_MORPH_OK           : total=   9  in_high=   9
BRIGHT_STAR_FLAG            : total=   0  in_high=   0
WARN_MASK                   : total=  31  in_high=  31
SEVERE_CEN_ANY              : total=  78  in_high=  78
WARN_CUTOUT_MISSING_SHAPE   : total=   0  in_high=   0

MEDIUM PRIORITY DRIVERS
ELL_MISMATCH                : total= 280  in_medium= 190
WARN_WEAK_HA                : total= 236  in_medium= 208
WARN_CEN_ANY                : total= 206  in_medium= 119
WARN_R_PROFILE_PEAK         : total= 137  in_medium=  90
WARN_CUTOUT_MISSING         : total=   0  in_medium=   0

OVERLAP AMONG HIGH DRIVERS
NOT_PHOT_OK                  & NOT_HAPY_MORPH_OK           :    6
NOT_PHOT_OK                  & WARN_MASK                   :    1
NOT_HAPY_MORPH_OK            & WARN_MASK                   :    3
NOT_HAPY_MORPH_OK            & SEVERE_CEN_ANY              :    2
WARN_MASK                    & SEVERE_CEN_ANY              :   14

```

#### 2026-06-20 run
```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260620$ python ~/github/hapy/scripts/qc_results.py --scheme virgo merged_results_virgo_20260622.fits 
Read 797 rows from merged_results_virgo_20260622.fits
UPDATE: adding VFINDEX
REVIEW_PRIORITY SUMMARY
{np.str_('high'): np.int64(79), np.str_('low'): np.int64(322), np.str_('medium'): np.int64(396)}
ELL_MISMATCH 241
FILTER_WARNING 76
WARN_MASK 11
BRIGHT_STAR_FLAG 0
WARN_WEAK_HA 184
/home/siena.edu/rfinn/github/hapy/scripts/qc_results.py:564: UserWarning: Warning: converting a masked element to nan.
  out[i] = float(v)
Wrote QC products to qc
Number of high priority in qc/tables/review = 79
writing  qc/tables/review/review_sample.csv


REVIEW PRIORITY SUMMARY
{np.str_('high'): np.int64(79), np.str_('low'): np.int64(322), np.str_('medium'): np.int64(396)}

HIGH PRIORITY DRIVERS
NOT_PHOT_OK                 : total=   0  in_high=   0
NOT_HAPY_MORPH_OK           : total=   2  in_high=   2
BRIGHT_STAR_FLAG            : total=   0  in_high=   0
WARN_MASK                   : total=  11  in_high=  11
SEVERE_CEN_ANY              : total=  70  in_high=  70
WARN_CUTOUT_MISSING_SHAPE   : total=   0  in_high=   0

MEDIUM PRIORITY DRIVERS
ELL_MISMATCH                : total= 241  in_medium= 169
WARN_WEAK_HA                : total= 184  in_medium= 174
WARN_CEN_ANY                : total= 208  in_medium= 134
WARN_R_PROFILE_PEAK         : total= 147  in_medium= 105
WARN_CUTOUT_MISSING         : total=   0  in_medium=   0

OVERLAP AMONG HIGH DRIVERS
WARN_MASK                    & SEVERE_CEN_ANY              :    4

```

### Duplicates

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



## Build Cutout Webpages
Create a list of the cutout images:
```bash
find /data-pool/Halpha/hapy-output-20260417/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
```
```
find /data-pool/Halpha/hapy-output-20260612/cutouts -mindepth 1 \
-maxdepth 1 -type d -printf "%f\n" | sort > \
cutout_list_buildwebpages.txt
```
Test on one directory:
```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260417
```
```bash
python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir \
$ROOTDIR/cutouts --outdir $ROOTDIR/html/cutouts --oneimage VFID2891-UGC04559-HDI-20200225-p004
```


```bash 
parallel --bar -j 20 --joblog build_web_cutouts.joblog --results build_web_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir $ROOTDIR/cutouts --oneimage {} --outdir $ROOTDIR/html/cutouts :::: cutout_list_buildwebpages.txt
```

## Build cutout index

```
python ~/github/hapy/scripts/build_cutout_index.py --help
```

```bash
python ~/github/hapy/scripts/build_cutout_index.py --runroot /data-pool/Halpha/hapy-output-20260417/ --results-table /data-pool/Halpha/hapy-output-20260417/merged_results_virgo_20260421.fits
```


```
python ~/github/hapy/scripts/build_cutout_index.py --runroot \
/data-pool/Halpha/hapy-output-20260519-hybrid/ --results-table \
/data-pool/Halpha/hapy-output-20260519-hybrid/merged_results_virgo_20260521_with_best_duplicate.fits
```


```
python ~/github/hapy/scripts/build_cutout_index.py --runroot \
/data-pool/Halpha/hapy-output-20260620/ --results-table \
/data-pool/Halpha/hapy-output-20260620/merged_results_virgo_20260622.fits
```

## Rsync files

from ROOTDIR:
```
rsync -avz html/cutouts fitsxfr.siena.edu:/var/www/html/fits/virgo/.
```


# Visualize CS Images and Select Best Duplicate

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
python ~/github/hapy/scripts/inspect_cs_images.py plot-one \
cs_image_inspection/cs_image_inspection_groups.ecsv VFID0481 \
--cutout-dir cutouts --outdir cs_image_inspection 
```

To build the plots in parallel:
```
parallel --bar -j 16 --joblog cs_image_plot.joblog --results \
cs_image_plot_logs python ~/github/hapy/scripts/inspect_cs_images.py \
plot-one cs_image_inspection/cs_image_inspection_groups.ecsv {} \
--cutout-dir cutouts --outdir cs_image_inspection :::: cs_group_list.txt
```

# Run photometry on Halpha image with continuum for qc/validation

To run on one galaxy:
```
python ~/github/hapy/scripts/measure_ha_with_continuum_profiles.py --cutout-dir cutouts/VFID1934-NGC2799-INT-20190205-p026 --overwrite
```

```
find cutouts -mindepth 1 -maxdepth 1 -type d | sort > cutout_with_dir.txt
```


```
parallel --bar -j 16 --joblog measure_ha_continuum.joblog --results measure_ha_continuum_logs python ~/github/hapy/scripts/measure_ha_with_continuum_profiles.py --cutout-dir {} --overwrite :::: cutout_with_dir.txt
```


```bash
merge_results --mode ha_continuum --indir cutouts --review-csv review_sample_20260514.csv 
```

Running this on the following directories:
- /data-pool/Halpha/hapy-output-20260517-pre2025coadds
- /data-pool/Halpha/hapy-output-20260517
- /data-pool/Halpha/hapy-output-20260609-hybrid


### To transfer data tables to my laptop:
```
rsync -avz
draco:/data-pool/Halpha/hapy-output-20260517-pre2025coadds/merged"*".fits
hapy-output-20260517-pre2025coadds/.
```

```
rsync -avz
draco:/data-pool/Halpha/hapy-output-20260609-hybrid/merged"*".fits
hapy-output-20260609-hybrid/.
```


```
rsync -avz draco:/data-pool/Halpha/hapy-output-20260517/merged"*".fits hapy-output-20260517/.
```


# After completing manual masking

- I identified galaxies that needed manual masking and then created
  the masks using `run_maskgui`
- I have to rerun `run_analysis` on these galaxies

## 1. Identify the galaxies that have a manual mask:
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


# Create output table that is row-matched to Virgo Filament Tables

### 1. Merge all HAPY cutout results
```
merge_results --indir cutouts --mode run_analysis \
--review-csv review_sample_20260514.csv 
```
This will add the columns from the `review_sample.csv` file.  The
script currently eliminates cutouts with `VIS

### 2. Make duplicate-inspection table
```
python ~/github/hapy/scripts/inspect_cs_images.py make-table \
    merged_results.fits \
    --outdir cs_image_inspection \
    --min-dups 1
```

This add information about duplicates, including `BEST_DUPLICATE`.

### 3. Make final row-matched VFS HAPY table


To create an output table that is row-matched to Virgo Filament Tables run:
```
python ~/github/hapy/scripts/make_vfs_hapy_rowmatched.py merged_results_virgo_20260615_with_best_duplicate.fits
```

Then copy the resulting table (`vf_v2_halpha_YYYYMMDD.fits`) to the
vfs table directory, and update `hapypost/hapypost/io/vfs_tables.py`
to read the new table.
