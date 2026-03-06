## Make Cutouts

### Running on VFS-Halpha coadds
```
python ~/github/hapy/scripts/get_cutouts.py --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo
```

To run in parallel, you need a list of the r-band coadds. For example:
```
ls VF*R.fits > coadd_list
ls VF*r.fits >> coadd_list
```

The run parallel.  Testing on macbook:

```
parallel --eta -j 0 get_cutouts --catalog
/Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme
virgo --overwrite_metadata --rimage :::: coadd_list
```


## Run Analysis on Cutouts

### Run on One Cutout


```
run_analysis --cutout-dir cutouts/VFID3084-NGC3512-HDI-20200226-p012 --make-mask  --psf-image VF-165.869+28.044-HDI-20200226-p012-r-psf.fits --statmorph --galfit --convflag --diagnostic-plots
```



### Running on a larger sample


#### Make List of Cutouts
```
find cutouts -mindepth 1 -maxdepth 1 -type d > cutout_list.txt
```

To save the full path:

```
find cutouts -mindepth 1 -maxdepth 1 -type d \
  -exec test -f "{}/$(basename {}).results.ecsv" \; -print > cutout_list.txt
```
#### Run in Parallel

```bash
parallel --eta -j 0 --joblog run_analysis_joblog.tsv  run_analysis --cutout-dir {} --make-mask --psf-dir
.  --log-to-console :::: cutout_list.txt 2>&1 | tee screen_dump.out
```



# Merge Results

```python
python ~/github/hapy/scripts/merge_results.py --indir cutouts/
```


