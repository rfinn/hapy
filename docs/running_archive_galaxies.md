# Overview
 
We are incorporating archival data for Virgo Cluster (KKY2001) and Isolated galaxies (KK2006) into our analysis of halpha/SFR properties as a function of environment.

# Sort Programs


# Update Headers
- Move directory that contains subdirectories for each galaxy.
- Clear any prior versions of fits files with updated headers:
```
rm */h*.fits
```
- Then update headers:
```
python ~/github/hapy/scripts/update_archive_headers.py
```
# Rename Directories and Images
```
python ~/github/hapy/scripts/build_metadata_archive.py --archive-root ~/research/Virgo/koopmann-images/Virgo/raw --output-root ~/research/Virgo/koopmann-images/Virgo/cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --image-list ~/research/Virgo/koopmann-images/virgo.list --shape-catalog ~/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits --dry-run
```


```
python ~/github/hapy/scripts/build_metadata_archive.py --archive-root ~/research/Virgo/koopmann-images/Virgo/raw --output-root ~/research/Virgo/koopmann-images/Virgo/cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --image-list ~/research/Virgo/koopmann-images/virgo.list --shape-catalog ~/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits --dry-run
```

# Run `run_analysis`
I then `rsync` the files to draco:

```
cd /Users/rfinn/research/Virgo/koopmann-images/Virgo/cutouts
```

```
rsync -avz * draco:/data-pool/HalphaArchive/virgo_cluster/cutouts/.
```

The remaining commands are executed on draco.

```
cd /data-pool/HalphaArchive/virgo_cluster
```

Build a directory list:
```
find cutouts/ -mindepth 1 -maxdepth 1 -type d ! -name "cutouts_summary" | sort > cutout_with_dir.txt
```

Test one galaxy:
```bash
run_analysis --cutout-dir cutouts/VFID3757-NGC4561-TI2-19990101-n4561 --make-mask --statmorph --galfit --convflag --no-gaia
```

Run sample
```bash
head -5 cutout_with_dir.txt | parallel --bar -j 2 --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir cutouts/VFID3757-NGC4561-TI2-19990101-n4561 --make-mask --statmorph --galfit --convflag --no-gaia
```

```bash
parallel --bar  -j 16  --memfree 60G --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --statmorph --galfit --convflag --no-gaia :::: cutout_with_dir.txt
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
find /data-pool/HalphaArchive/virgo_cluster/hapy-output-20260326/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_legacy.txt
```

```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260310
```
then
```bash
parallel --bar -j 2 --joblog fetch_legacy.joblog --results fetch_legacy_logs python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir "$ROOTDIR/cutouts/{}"  :::: cutout_list_legacy.txt
```
  
## Build Cutout Webpages
Create a list of the cutout images:
```bash
find /data-pool/HalphaArchive/virgo_cluster/hapy-output-20260326/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
```

Test on one directory:
```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260313
```
```bash
python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir $ROOTDIR/cutouts --outdir $ROOTDIR/html/cutouts --oneimage VFID3757-NGC4561-TI2-19990101-n4561 
```


```bash 
parallel --bar -j 16 --memfree 60G --joblog build_web_cutouts.joblog --results build_web_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir $ROOTDIR/cutouts --oneimage "{}" --outdir $ROOTDIR/html/cutouts :::: cutout_list_buildwebpages.txt
```

## Build cutout index

```
python ~/github/hapy/scripts/build_cutout_index.py --help
```

```
python ~/github/hapy/scripts/build_cutout_index.py --runroot /data-pool/HalphaArchive/virgo_cluster/hapy-output-20260326/
```

## Transfer to fitsxfr

```
rsync -avz html/cutouts fitsxfr.siena.edu:/var/www/html/fits/archive/.
```
