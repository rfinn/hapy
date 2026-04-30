# Overview
 
We are incorporating archival data for Virgo Cluster (KKY2001) and Isolated galaxies (KK2006) into our analysis of halpha/SFR properties as a function of environment.

I am doing a lot of the pre-pipeline setup on my macbook and then
moving the pipeline-ready files to draco to run `run_analysis`.
# Sort Images into Directories

TBA


# Update Headers
- Move to the directory that contains subdirectories for each galaxy
  (e.g. on macbook, this is
 `/Users/rfinn/research/Virgo/koopmann-images/Virgo/raw/` or 
  `/Users/rfinn/research/Virgo/koopmann-images/KK06-isolated/raw/`
- Clear any prior versions of fits files with updated headers:
```
rm */h*.fits
```
- Then update headers:
```
python ~/github/hapy/scripts/update_archive_headers.py --center-file ~/research/Virgo/koopmann-images/offcenter_virgo_centers.csv
```
# Rename Directories and Images

This will create a parallel set of directories that are named
according to the `Virgo Filament Survey` naming convention.  This
allows for easier incorporation into the analysis downstream.
```
python ~/github/hapy/scripts/build_metadata_archive.py --archive-root ~/research/Virgo/koopmann-images/Virgo/raw --output-root ~/research/Virgo/koopmann-images/Virgo/cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --image-list ~/research/Virgo/koopmann-images/virgo.list --shape-catalog ~/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits --dry-run
```

When the output looks correct, remove the `--dry-run` flag and rerun.

```
python ~/github/hapy/scripts/build_metadata_archive.py --archive-root ~/research/Virgo/koopmann-images/Virgo/raw --output-root ~/research/Virgo/koopmann-images/Virgo/cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --image-list ~/research/Virgo/koopmann-images/virgo.list --shape-catalog ~/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits 
```

## For the isolated galaxies
```
python ~/github/hapy/scripts/build_metadata_archive.py --archive-root ~/research/Virgo/koopmann-images/KK06-isolated/raw --output-root ~/research/Virgo/koopmann-images/KK06-isolated/cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --image-list ~/research/Virgo/koopmann-images/isolated3.list --shape-catalog ~/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits --dry-run
```

Only 11 of the isolated galaxies match the VFS catalog - the others
are outside the VFS footprint.  

When the output looks correct, remove the `--dry-run` flag and rerun
for the isolated sample.

```
python ~/github/hapy/scripts/build_metadata_archive.py --archive-root ~/research/Virgo/koopmann-images/KK06-isolated/raw --output-root ~/research/Virgo/koopmann-images/KK06-isolated/cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --image-list ~/research/Virgo/koopmann-images/isolated3.list --shape-catalog ~/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits 
```

# Run `run_analysis`
Make a directory on draco,
e.g. `/data-pool/HalphaArchive/virgo_cluster/hapy-output-20260401`

I then `rsync` the files to draco:

```
cd /Users/rfinn/research/Virgo/koopmann-images/Virgo/
```


```
rsync -avz cutouts draco:/data-pool/HalphaArchive/virgo_cluster/hapy-output-20260401/.
```

The remaining commands are executed on draco and mirror the steps
used for the VFS and UAT samples.

```
cd /data-pool/HalphaArchive/virgo_cluster/hapy-output-20260401
```

Build a directory list:
```
find cutouts/ -mindepth 1 -maxdepth 1 -type d ! -name "cutouts_summary" | sort > cutouts_with_dir.txt
```

Test one galaxy:
```bash
run_analysis --cutout-dir cutouts/VFID3757-NGC4561-TI2-19990101-n4561 --make-mask --statmorph --galfit --convflag --no-gaia
```

```bash
run_analysis --cutout-dir cutouts/VFID0155-NGC2787-T2KA-19950203-n2787 --make-mask --statmorph --galfit --convflag --no-gaia
```


Run sample
```bash
head -5 cutout_with_dir.txt | parallel --bar -j 2 --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir cutouts/VFID3757-NGC4561-TI2-19990101-n4561 --make-mask --statmorph --galfit --convflag --no-gaia
```

```bash
parallel --bar  -j 16  --memfree 60G --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --statmorph --galfit --no-gaia :::: cutouts_with_dir.txt
```

```
parallel --bar  -j 4 --memfree 60G --joblog run_analysis.joblog
--results parallel-logs run_analysis --cutout-dir "{}" --make-mask
--statmorph --galfit --no-gaia :::: cutouts_with_dir.txt
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




# Build Webpages to Review Cutouts

## Download Legacy Images
### To copy legacy images from a prior run

Run this command from the directory that contains
e.g. `hapy-output-20260401` and `hapy-output-20260429`.
```
rsync -av hapy-output-20260401/cutouts/ hapy-output-20260429/cutouts/
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

Isolated sample:
```
find /data-pool/HalphaArchive/isolated/hapy-output-20260419/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_legacy.txt
```


```bash
ROOTDIR=/data-pool/HalphaArchive/virgo_cluster/hapy-output-20260401/
```
then
```bash
parallel --bar -j 2 --joblog fetch_legacy.joblog --results fetch_legacy_logs python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir "$ROOTDIR/cutouts/{}"  :::: cutout_list_legacy.txt
```
  
## Build Cutout Webpages
Create a list of the cutout images:
```bash
find /data-pool/HalphaArchive/virgo_cluster/hapy-output-20260429/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
```


or 

```bash
find /data-pool/HalphaArchive/isolated/hapy-output-20260419/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
```


Set the ROOTDIR

```
ROOTDIR=/data-pool/HalphaArchive/virgo_cluster/hapy-output-20260429/
```

Test on one directory:

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

For cluster sample:
```
python ~/github/hapy/scripts/build_cutout_index.py --runroot
/data-pool/HalphaArchive/virgo_cluster/hapy-output-20260401/
--results-table merged_results.fits 
```

For isolated sample:
```
python ~/github/hapy/scripts/build_cutout_index.py --runroot
/data-pool/HalphaArchive/isolated/hapy-output-20260419/
--results-table merged_results.fits 
```

## Transfer to fitsxfr

```
rsync -avz html/cutouts fitsxfr.siena.edu:/var/www/html/fits/archive/cluster/.
```

From `/data-pool/HalphaArchive/isolated/hapy-output-20260419`:
```
rsync -avz html/cutouts
fitsxfr.siena.edu:/var/www/html/fits/archive/isolated/.
```
