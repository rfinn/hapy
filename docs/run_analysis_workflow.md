> This supercedes running-pipeline.md

# HAPY running-analysis rerun workflow

This file is the streamlined recipe for rerunning the Virgo HAPY analysis. Keep this file mostly command-focused. Put dated outputs, oddities, and interpretation in `running-analysis-results-log.md`.

## 0. Set run variables

Edit these at the top of each new run.

```bash
RUNDATE=20260626
ROOTDIR=/data-pool/Halpha/hapy-output-${RUNDATE}
COADDDIR=/data-pool/Halpha/coadds-v20260330
PSFDIR=/data-pool/Halpha/psf-images-v20260330
GAIADIR=${COADDDIR}/gaia_catalogs
CATALOG=~/research/Virgo/tables-north/v2/vf_v2_main.fits
REVIEW_CSV=review_sample_20260612.csv
NCPU=16
```

For a hybrid-coadd run, set `COADDDIR` and `PSFDIR` to the hybrid versions before creating lists or cutouts.

## 1. One-time preparation for INT simple weights

Only do this once for a given coadd directory, and only when replacing the original INT weights with simple good/bad masks.

```bash
cd ${COADDDIR}
mkdir -p ORIGINAL_INT_WEIGHTS
mv *INT*weight.fits ORIGINAL_INT_WEIGHTS/.
ls VF*INT*r.fits VF*INT*Halpha.fits VF*INT*Ha6657.fits > INT_all_coadds.txt
parallel --bar -j ${NCPU} --joblog make_simple_weight_INT.joblog --results make_simple_weight_logs python ~/github/hapy/scripts/make_simple_weight_from_coadd.py {} :::: INT_all_coadds.txt
```

## 2. Download Gaia catalogs for coadds

```bash
cd ${COADDDIR}
python ~/github/hapy/scripts/download_gaia_coadd_catalogs.py
```

## 3. Prepare the output directory

```bash
cd /data-pool/Halpha
mkdir -p ${ROOTDIR}
cd ${ROOTDIR}
```

Create the coadd list:

```bash
find ${COADDDIR}/ -maxdepth 1 -type f \( -name "VF*r.fits" -o -name "VF*R.fits" \) | sort > fullpath_rcoadds_all.txt
wc -l fullpath_rcoadds_all.txt
```

If some coadds are not ready, make a curated list:

```bash
cp fullpath_rcoadds_all.txt fullpath_rcoadds_hapy_ready.txt
```

Then manually remove any coadds that are still under review.

## 4. Make cutouts

Test on one coadd first:

```bash
get_cutouts --catalog ${CATALOG} --scheme virgo --maxcorrection 5 --psfdir ${PSFDIR} --rimage $(head -1 fullpath_rcoadds_all.txt)
```

Run the full list:

```bash
parallel -j ${NCPU} --bar --joblog cutouts_parallel.log get_cutouts --rimage {} --catalog ${CATALOG} --scheme virgo --maxcorrection 5 --overwrite --overwrite_metadata :::: fullpath_rcoadds_all.txt
```

Check for failed jobs:

```bash
awk 'NR==1 || $7 != 0 {print}' cutouts_parallel.log
```

Check cutout products:

```bash
python ~/github/hapy/scripts/check_cutouts.py fullpath_rcoadds_all.txt cutouts
```

Merge the get_cutouts summary tables:

```bash
merge_results --mode get_cutouts --indir cutouts_summary --out merged_cutouts_results.fits
```

## 5. Get Legacy images to use when making CS-gr images

### 5.1 Copy Legacy cutouts from a previous run, when possible

Run this from `/data-pool/Halpha`, not from inside the run directory.

```bash
cd /data-pool/Halpha
rsync -av hapy-output-20260620/cutouts/ hapy-output-${RUNDATE}/cutouts/ --include '*/' --include 'legacy/***' --exclude '*' --exclude '*logs*' --ignore-existing --prune-empty-dirs --dry-run
```

If the dry run looks good, repeat without `--dry-run`.

### 5.2 Download missing Legacy cutouts

From `${ROOTDIR}`:

```bash
find cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list.txt
parallel --bar -j 2 --joblog fetch_legacy.joblog --results fetch_legacy_logs python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir "${ROOTDIR}/cutouts/{}" :::: cutout_list.txt
python ~/github/hapy/scripts/find_missing_legacy_cutouts.py
```

Retry missing ones, if needed:

```bash
parallel --bar -j 2 --joblog fetch_legacy_retry.joblog python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir {} :::: missing_legacy_cutouts.txt
```

### 5.3 Reproject Legacy images

```bash
find cutouts -mindepth 1 -maxdepth 1 -type d -name 'VFID*' | sort > reproject_cutout_list.txt
python ~/github/hapy/scripts/make_legacy_reprojections.py $(head -1 reproject_cutout_list.txt)
parallel --bar -j 20 --results legacy_reproject_logs python ~/github/hapy/scripts/make_legacy_reprojections.py "{}" --overwrite :::: reproject_cutout_list.txt
```

### 5.4 OUTDATED - skip to step 6
> We now build the cs-gr images within `run_analysis`, so continue to Step 6.
#### Build CS-gr images

Test one:

```bash
python ~/github/hapy/hapy/scripts/make_cs_gr.py $(head -1 reproject_cutout_list.txt) --auto-contscale --auto-contscale-percentile 30 --overwrite
```

Run all:

```bash
parallel --bar -j ${NCPU} --joblog cs_gr_auto.joblog --results cs_gr_auto_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {} --auto-contscale --auto-contscale-percentile 30 --overwrite :::: reproject_cutout_list.txt
```

Check for missing CS-gr products:

```bash
python ~/github/hapy/scripts/check_for_missing_csgr.py
```

Search logs for failures:

```bash
grep -R "FAILED\|Traceback\|ERROR\|can't find\|problem getting" legacy_reproject_logs cs_gr_auto_logs
```

## 6. Make stellar-mass and SFR maps

```bash
find cutouts -mindepth 1 -maxdepth 1 -type d | sort > cutout_list.txt
parallel --bar -j ${NCPU} python ~/github/hapy/hapy/scripts/make_mstar_map.py {} --scheme virgo --overwrite :::: cutout_list.txt
parallel --bar -j ${NCPU} --joblog make_sfr_map.joblog python ~/github/hapy/hapy/scripts/make_sfr_map.py {} --scheme virgo --overwrite :::: cutout_list.txt
```

## 7. Copy manual masks from the previous reviewed run

Dry run first:

```bash
rsync -avzn --include='*/' --include='*-mask-manual.fits' --exclude='*' /data-pool/Halpha/hapy-output-20260/cutouts/ ${ROOTDIR}/cutouts/
```

Then copy:

```bash
rsync -avz --include='*/' --include='*-mask-manual.fits' --exclude='*' /data-pool/Halpha/hapy-output-20260417/cutouts/ ${ROOTDIR}/cutouts/
```

## 8. Run analysis

Test one cutout:

```bash
run_analysis --make-mask --psf-dir ${PSFDIR}/ --csgr --clumps --statmorph --galfit --convflag --log-to-console --gaia-dir ${GAIADIR}/ --cutout-dir $(find cutouts -mindepth 1 -maxdepth 1 -type d | sort | head -1)
```

Create a run list, using the review file to skip excluded targets:

```bash
python ~/github/hapy/scripts/make_run_analysis_list.py --cutout-dir cutouts --review ${REVIEW_CSV} --outfile cutout_run_analysis_list.txt
wc -l cutout_run_analysis_list.txt
```

Create list with all cutouts:
```bash
find cutouts/ -mindepth 1 -maxdepth 1 -type d ! -name "cutouts_summary" | sort > cutout_with_dir.txt
```
Run the sample:

```bash
parallel --bar -j ${NCPU} --memfree 30G --joblog run_analysis.joblog --results parallel-logs run_analysis --cutout-dir "{}" --make-mask --psf-dir ${PSFDIR}/ --csgr --clumps --statmorph --galfit --convflag --gaia-dir ${GAIADIR}/ :::: cutout_with_dir.txt
```

Monitor progress and failures:

```bash
htop
awk 'NR==1 || $7 != 0 {print}' run_analysis.joblog
grep -R "FAILED\|Traceback\|ERROR" parallel-logs
```

## 9. Merge run_analysis results

```bash
merge_results --indir cutouts --mode run_analysis --review-csv ${REVIEW_CSV} --out merged_results_virgo_${RUNDATE}.fits
```

## 10. Summarize and QC the run

```bash
python ~/github/hapy/scripts/summarize_run.py --infile merged_results_virgo_${RUNDATE}.fits --scheme virgo
python ~/github/hapy/scripts/qc_results.py merged_results_virgo_${RUNDATE}.fits --scheme virgo
python ~/github/hapy/scripts/validate_dashboards.py merged_results_virgo_${RUNDATE}.fits --sample ALL
```

## 11. Build webpages for review

```bash
find ${ROOTDIR}/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir ${ROOTDIR}/cutouts --outdir ${ROOTDIR}/html/cutouts --oneimage $(head -1 cutout_list_buildwebpages.txt)
parallel --bar -j 20 --joblog build_web_cutouts.joblog --results build_web_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir ${ROOTDIR}/cutouts --oneimage {} --outdir ${ROOTDIR}/html/cutouts :::: cutout_list_buildwebpages.txt
python ~/github/hapy/scripts/build_cutout_index.py --runroot /data-pool/Halpha/hapy-output-20260626/ --results-table merged_results_virgo_20260628.fits
```

Sync webpages to the server:

```bash
cd ${ROOTDIR}
rsync -avz html/cutouts fitsxfr.siena.edu:/var/www/html/fits/virgo/.
```

## 12. CS image inspection and best duplicate selection

```bash
python ~/github/hapy/scripts/inspect_cs_images.py make-table merged_results_virgo_20260628.fits --outdir cs_image_inspection --min-dups 1
python ~/github/hapy/scripts/inspect_cs_images.py list-groups
cs_image_inspection/cs_image_inspection_groups.ecsv >
cs_group_list.txt
python ~/github/hapy/scripts/inspect_cs_images.py plot-one cs_image_inspection/cs_image_inspection_groups.ecsv VFID0377 --cutout-dir cutouts --outdir cs_image_inspection
```

```
python ~/github/hapy/scripts/inspect_cs_images.py plot-one cs_image_inspection/cs_image_inspection_groups.ecsv $(head -1 cs_group_list.txt) --cutout-dir cutouts --outdir cs_image_inspection 
parallel --bar -j ${NCPU} --joblog cs_image_plot.joblog --results cs_image_plot_logs python ~/github/hapy/scripts/inspect_cs_images.py plot-one cs_image_inspection/cs_image_inspection_groups.ecsv {} --cutout-dir cutouts --outdir cs_image_inspection :::: cs_group_list.txt
```
