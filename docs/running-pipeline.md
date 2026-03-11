

# Download the gaia files
Move to directory containing coadded images

```
cd /data-pool/Halpha/coadds-2025DEC
```

Then start the download.

```
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

If you have coadds that are still under review, kake a copy of the
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

```bash
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260309$ get_cutouts --help
usage: get_cutouts [-h] [--rimage RIMAGE] [--outdir OUTDIR] [--psfdir PSFDIR] [--cutout_scale CUTOUT_SCALE] [--overwrite_metadata] [--no-skysub]
                   [--catalog CATALOG] [--scheme {generic,virgo,agc}] [--maxcorrection MAXCORRECTION]

create psf image from image that contains stars

options:
  -h, --help            show this help message and exit
  --rimage RIMAGE       r-band image name
  --outdir OUTDIR       Directory where cutouts/ will be created (default: current working directory).
  --psfdir PSFDIR       set to coadd directory
  --cutout_scale CUTOUT_SCALE
                        multiplicative scale factor for increasing the size of cutout images
  --overwrite_metadata  Set this to overwrite metadata.json. Will store a *.bak file.
  --no-skysub           Disable local sky subtraction in cutouts (default: sky is subtracted).
  --catalog CATALOG     full path to galaxy catalog to use for cutouts.
  --scheme {generic,virgo,agc}
                        Filename parsing scheme for coadd images.
  --maxcorrection MAXCORRECTION
                        maximum filter correction for galaxies in FOV. default is 3, so galaxies whose redshift falls where filter transmission < 33 percent will
                        be skipped.

```

### Test on One Image

### 

```bash
get_cutouts --rimage
/data-pool/Halpha/coadds-2025DEC/VF-126.291+27.988-HDI-20180313-p004-R.fits
--catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits
--scheme virgo --maxcorrection 5
```

### Run on full Virgo sample

```bash
cat fullpath_rcoadds_hapy_ready.txt | parallel -j 32 --bar --joblog
cutouts_parallel.log get_cutouts --rimage {} --catalog
~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5
```

### Check output
```
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
```
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


```
run_analysis --cutout-dir cutouts/VFID3084-NGC3512-HDI-20200226-p012
--make-mask  --psf-dir /data-pool/Halpha/psf-images/ --statmorph
--galfit --convflag --log-to-console --gaia-dir
/data-pool/Halpha/coadds-2025DEC/gaia_catalogs/
```



## Running on a larger sample

### Create a list of cutouts:

```bash
find cutouts/ -mindepth 1 -maxdepth 1 -type d ! -name "cutouts_summary" | sort > cutout_list.txt
```

check that the file contains, e.g., `cutouts/VFID2943-NGC2604-INT-20190204-p010`
```
head cutout_list.txt
wc -l cutout_list.txt
```

### Test smaller samples
Test on 5 galaxies:
```
head -5 cutout_list.txt | parallel --bar -j 2 --joblog run_analysis.joblog \
  --results parallel-logs \
  run_analysis --cutout-dir "{}" --make-mask \
  --psf-dir /data-pool/Halpha/psf-images/ \
  --statmorph --galfit --convflag \
  --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/
```

Test on 20 galaxies:
```
head -5 cutout_list.txt | parallel --bar -j 4 --joblog run_analysis.joblog \
  --results parallel-logs \
  run_analysis --cutout-dir "{}" --make-mask \
  --psf-dir /data-pool/Halpha/psf-images/ \
  --statmorph --galfit --convflag \
  --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/
```

Run the next 20 galaxies:
```
sed -n '21,40p' cutout_list.txt | parallel --bar -j 8 --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/
```

### Run in Parallel

```bash
parallel --bar  -j 16  --memfree 60G --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-2025DEC/gaia_catalogs/ :::: cutout_list.txt
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


# Merge Results

```bash
merge_results --indir cutouts --mode run_analysis
```

## Summarize statistics in merged_results.py

```
python ~/github/hapy/scripts/summarize_run.py  merged_results.fits
```

# Build Webpages to review Cutouts

## Download Legacy Images
```
python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir cutouts/VFID2891-UGC04559-HDI-20200225-p004/
```

```
parallel --bar -j 8 --joblog fetch_legacy.joblog --results fetch_legacy_logs python ~/github/hapy/hapy/imagetools/fetch_legacy_cutouts.py --cutout-dir "$RUNROOT/cutouts/{}" --layer ls-dr10 : cutout_list.txt
```
  
## Build webpage

```
find /data-pool/Halpha/hapy-output-20260310/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list.txt
```

```
```
```
python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir
/data-pool/Halpha/hapy-output-20260310/cutouts --oneimage
VFID2891-UGC04559-HDI-20200225-p004 --outdir
/data-pool/Halpha/hapy-output-20260310/html/cutouts
```

```
parallel --bar -j 8 --joblog build_web_cutouts.joblog --results build_web_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir /data-pool/Halpha/hapy-output-20260310/cutouts --oneimage "{}" --outdir /data-pool/Halpha/hapy-output-20260310/html/cutouts :::: cutout_list.txt
```

## Build cutout index

```
python ~/github/hapy/scripts/build_cutout_index.py --help
```

```
python ~/github/hapy/scripts/build_cutout_index.py --runroot /data-pool/Halpha/hapy-output-20260310/
```
