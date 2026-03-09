

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
/data-pool/Halpha/coadds-2025DEC/VF-126.291+27.988-HDI-20180313-p004-R.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo
```


### Running on VFS-Halpha coadds
```
get_cutouts --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo
```


The run parallel.  Testing on macbook:

```
parallel --eta -j 0 get_cutouts --catalog
/Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme
virgo --overwrite_metadata --rimage :::: coadd_list
```


# Run Analysis 

### Run on One Cutout


```
run_analysis --cutout-dir cutouts/VFID3084-NGC3512-HDI-20200226-p012 --make-mask  --psf-image VF-165.869+28.044-HDI-20200226-p012-r-psf.fits --statmorph --galfit --convflag --diagnostic-plots
```



### Running on a larger sample


#### Run in Parallel

```bash
parallel --eta -j 0 --joblog run_analysis_joblog.tsv  run_analysis --cutout-dir {} --make-mask --psf-dir
.  --log-to-console :::: cutout_list.txt 2>&1 | tee screen_dump.out
```



# Merge Results

```python
python ~/github/hapy/scripts/merge_results.py --indir cutouts/
```


