## Make Cutouts

### Running on VFS-Halpha coadds
```
python ~/github/hapy/scripts/get_cutouts.py --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo
```

### Running with the AGC for UAT Halpha Groups
If using the AGC scheme, the code will assume that coadded images are named like (someday I will make this more robust than just splitting on '`'):
``` 
VF-202.260+46.474-HDI-20170523-p023-h01-R.fits
```

I copied one virgo coadd to test with both the virgo and agc schemes.  I first needed to alter the name to conform with the UAT naming convention.

First, add the extra '-h01' to the pointing name, as we do for the UAT groups.
```bash
for f in *p023*; do\n    mv "$f" "${f/p023/p023-h01}"\ndone
```

Second, update the name of HAIMAGE in the r-band header accordingly:
```bash
sethead VF-202.260+46.474-HDI-20170523-p023-h01-R.fits HAIMAGE=VF-202.260+46.474-HDI-20170523-p023-h01-ha4.fits
```


Then create cutouts:
```bash
python ~/github/hapy/scripts/get_cutouts.py --rimage VF-202.260+46.474-HDI-20170523-p023-h01-R.fits --catalog /Users/rfinn/research/AGC/agcnorthminus1.full200617.fits --scheme agc
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

To create a plot with the r-band image, the mask, the input ellipse, and the photutils ellipse, add the `--diagnostic-plots` argument.

### Make diagnostic plot for a mask

Here is one example of how to create the mask diagnostic plot (we will hard code this into `run_analysis.py` soon...)
```python
python ~/github/hapy/scripts/plot_mask_diagnostic.py --root cutouts/VFID1588-NGC5169-HDI-20170523-p023/VFID1588-NGC5169-HDI-20170523-p023
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
parallel -eta -j 0 python ~/github/hapy/scripts/run_analysis.py --cutout-dir {}   --make-mask --psf-dir . --statmorph --image2-filter 4 --galfit :::: cutout_list.txt
```

### With Logging and Failure Continuation

```
parallel -0 -j 0 --joblog logs/joblog.tsv \
  python ~/github/hapy/scripts/run_analysis.py --cutout-dir {} \
  :::: cutout_list.txt
```

# Merge Results

```python
python ~/github/hapy/scripts/merge_results.py --indir cutouts/
```
