# HAPY running-analysis special cases and troubleshooting

Use this file for reruns, partial fixes, and command fragments that are useful but should not clutter the main rerun recipe.

## Rerun INT cutouts only

Useful after changing the INT PSF directory or INT-specific calibration behavior.

```bash
find cutouts -mindepth 1 -maxdepth 1 -type d -name "*INT*" | sort > INT_cutouts.txt
parallel --bar -j 16 --joblog run_analysis_INT_rerun.joblog --memfree 30G --results run_analysis_INT_logs run_analysis --cutout-dir "{}" --make-mask --csgr --psf-dir /data-pool/Halpha/psf-images-v20260330/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/ :::: INT_cutouts.txt
```

## Rerun cutouts with missing CS-gr products

```bash
python ~/github/hapy/scripts/check_for_missing_csgr.py
parallel --bar -j 16 --joblog rerun_csgr.joblog --results rerun_csgr_logs run_analysis --cutout-dir {} --make-mask --csgr --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/ :::: missing_csgr_or_phot.txt
```

## Rerun individual galaxies

```bash
run_analysis --make-mask --csgr --psf-dir /data-pool/Halpha/psf-images/ --statmorph --galfit --convflag --gaia-dir /data-pool/Halpha/coadds-v20260330/gaia_catalogs/ --cutout-dir cutouts/VFID6463-WISEAJ150809.13+012516.6-INT-20190601-p033
```

Other individual rerun candidates from the old notes:

```text
cutouts/VFID6463-WISEAJ150809.13+012516.6-INT-20190601-p033
cutouts/VFID2766-WISEAJ133704.60+315337.9-HDI-20170522-p006
cutouts/VFID2745-UGC08602-HDI-20170522-p006
```

## Retry missing Legacy cutouts

```bash
python ~/github/hapy/scripts/find_missing_legacy_cutouts.py
parallel --bar -j 2 --joblog fetch_legacy_retry.joblog python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir {} :::: missing_legacy_cutouts.txt
```

## Rebuild CS-gr images for INT cutouts only

Useful after a filter-lookup fix.

```bash
find cutouts -mindepth 1 -maxdepth 1 -type d -name 'VFID*INT*' | sort > int_cutout_list.txt
parallel --bar -j 16 --joblog cs_gr_int.joblog --results cs_gr_int_logs python ~/github/hapy/hapy/scripts/make_cs_gr.py {} :::: int_cutout_list.txt
```

## Try different CS-gr continuum-scale methods on one cutout

Default auto-contscale percentile method:

```bash
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID0481-NGC6307-INT-20190602-p010 --auto-contscale --auto-contscale-percentile 30 --overwrite
```

Ratio method:

```bash
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID0569-NGC5989-BOK-20220424-VFID0607/ --auto-contscale --auto-contscale-method ratio --auto-contscale-percentile 30 --overwrite
```

Negative-tail method:

```bash
python ~/github/hapy/hapy/scripts/make_cs_gr.py cutouts/VFID1934-NGC2799-INT-20190205-p026 --auto-contscale --auto-contscale-method negtail --overwrite
```

## Quick log checks

Cutout failures:

```bash
awk 'NR==1 || $7 != 0 {print}' cutouts_parallel.log
```

Run-analysis failures:

```bash
awk 'NR==1 || $7 != 0 {print}' run_analysis.joblog
grep -R "FAILED\|Traceback\|ERROR" parallel-logs
```

Legacy / CS-gr failures:

```bash
grep -R "FAILED\|Traceback\|ERROR\|can't find\|problem getting" legacy_reproject_logs cs_gr_auto_logs fetch_legacy_logs
```

Monitor memory usage:

```bash
htop
```

## Webpage rebuild only

```bash
ROOTDIR=/data-pool/Halpha/hapy-output-20260612
find ${ROOTDIR}/cutouts -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort > cutout_list_buildwebpages.txt
parallel --bar -j 20 --joblog build_web_cutouts.joblog --results build_web_logs python ~/github/hapy/scripts/build_web_cutouts.py --cutoutdir ${ROOTDIR}/cutouts --oneimage {} --outdir ${ROOTDIR}/html/cutouts :::: cutout_list_buildwebpages.txt
python ~/github/hapy/scripts/build_cutout_index.py --runroot ${ROOTDIR}/ --results-table ${ROOTDIR}/merged_results_virgo_20260613.fits
```
