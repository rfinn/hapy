

# Download the gaia files
Move to directory containing coadded images

```bash
cd /data-pool/Halpha/coadds-2025DEC
```

Then start the download.

```bash
python ~/github/hapy/scripts/download_gaia_coadd_catalogs.py 
```

# Prepare the output directory
### make a directory for analysis run

```bash
cd /data-pool/Halpha/
```

```bash
mkdir hapy-output-20260309
```

```bash
cd hapy-output-20260309
```


### make a list of coadds to analyze

```bash
find /data-pool/Halpha/coadds-2025DEC/ -maxdepth 1 -type f \( -name "VF*r.fits" -o -name "VF*R.fits" \) | sort > fullpath_rcoadds_all.txt
```

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
# Make Cutouts

### Test on One Image

### 

```bash
get_cutouts --rimage
/data-pool/Halpha/coadds-2025DEC/VF-126.291+27.988-HDI-20180313-p004-R.fits
--catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits
--scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images/
```


### Run on full Virgo sample

```bash
cat fullpath_rcoadds_hapy_ready.txt | parallel -j 16 --bar --joblog cutouts_parallel.log get_cutouts --rimage {} --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5
```

### Check output
```bash
python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_hapy_ready.txt cutouts
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


# Run Analysis 

## Run on One Cutout


```bash
run_analysis --make-mask  --psf-dir /data-pool/Halpha/psf-images/ --statmorph
--galfit --convflag --log-to-console --gaia-dir
/data-pool/Halpha/coadds-2025DEC/gaia_catalogs/ --cutout-dir cutouts/VFID3084-NGC3512-HDI-20200226-p012
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
  --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/
```

Test on 20 galaxies:
```bash
head -5 cutout_list.txt | parallel --bar -j 4 --joblog run_analysis.joblog \
  --results parallel-logs \
  run_analysis --cutout-dir "{}" --make-mask \
  --psf-dir /data-pool/Halpha/psf-images/ \
  --statmorph --galfit --convflag \
  --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/
```

Run the next 20 galaxies:
```bash
sed -n '21,40p' cutout_list.txt | parallel --bar -j 8 --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/
```

### Run in Parallel

```bash
parallel --bar  -j 16  --memfree 60G --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/ :::: cutouts_with_dir.txt
```

```bash
parallel --eta -j 4 \
  --joblog run_analysis.joblog \
  --results parallel-logs \
  run_analysis --cutout-dir "{}" --make-mask \
  --psf-dir /data-pool/Halpha/psf-images/ \
  --statmorph --galfit --convflag \
  --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/ \
  :::: cutout_list.txt
```

You can monitor memory usage with 
```
htop
```


## Merge Results From `run_analysis`

```bash
merge_results --indir cutouts --mode run_analysis
```

### Summarize statistics in  `merged_results.py`
```
python ~/github/hapy/scripts/summarize_run.py --infile merged_results.fits --scheme virgo
```

## Make qc plots

Create some basic qc plots:
```
python ~/github/hapy/scripts/qc_results.py merge_results.fits --scheme virgo
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
python ~/github/hapy/scripts/qc_duplicates.py merged_results.fits
```



# Build Webpages to Review Cutouts

## Download Legacy Images
### To copy legacy images from a prior run

Run this command from the directory that contains
e.g. `hapy-output-20260313` and `hapy-output-20260319`.
```
rsync -av hapy-output-20260313/cutouts/ hapy-output-20260319/cutouts/
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



## Build Cutout Webpages
Create a list of the cutout images:
```bash
find /data-pool/Halpha/hapy-output-20260310/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
```

Test on one directory:
```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260313
```
```bash
python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir
$ROOTDIR/cutouts --outdir $ROOTDIR/html/cutouts --oneimage VFID2891-UGC04559-HDI-20200225-p004
```

```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260313
```

```bash 
parallel --bar -j 16 --memfree 60G --joblog build_web_cutouts.joblog --results build_web_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir $ROOTDIR/cutouts --oneimage "{}" --outdir $ROOTDIR/html/cutouts :::: cutout_list_buildwebpages.txt
```

## Build cutout index

```
python ~/github/hapy/scripts/build_cutout_index.py --help
```

```
python ~/github/hapy/scripts/build_cutout_index.py --runroot /data-pool/Halpha/hapy-output-20260310/
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
