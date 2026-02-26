## Make Cutouts

```
python ~/github/hapy/scripts/get_cutouts.py --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo
```

## Run Analysis on Cutouts

### Run on One Cutout

```
python ~/github/hapy/scripts/run_analysis.py --root cutouts/VFID3084-NGC3512-HDI-20200226-p012/VFID3084-NGC3512-HDI-20200226-p012 --make-mask --psf-image VF-165.869+28.044-HDI-20200226-p012-r-psf.fits --statmorph --image2-filter 4 --galfit
```

To make use of `--cutout-dir` instead of `--root`:

```
python ~/github/hapy/scripts/run_analysis.py --cutout-dir cutouts/VFID3084-NGC3512-HDI-20200226-p012 --make-mask --convflag 0 --psf-image VF-165.869+28.044-HDI-20200226-p012-r-psf.fits --statmorph --image2-filter 4 --galfit
```

### Running on a larger sample

```
mkdir -p logs
```

```
find cutouts -mindepth 1 -maxdepth 1 -type d > cutout_list.txt
```

```
parallel -j 0 python ~/github/hapy/scripts/run_analysis.py --cutout-dir {} :::: cutout_list.txt
```

```
parallel -j 0 python ~/github/hapy/scripts/run_analysis.py --cutout-dir {}   --make-mask --psf-dir . --statmorph --image2-filter 4 --galfit :::: cutout_list.txt
```

### With Logging and Failure Continuation

```
parallel -0 -j 0 --joblog logs/joblog.tsv \
  python ~/github/hapy/scripts/run_analysis.py --cutout-dir {} \
  :::: cutout_list.txt
```

# Merge Results
