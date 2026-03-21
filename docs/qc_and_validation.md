## QC 
```
python ~/github/hapy/scripts/qc_results.py merged_results.fits
```

## Using duplicates to assess
```
python ~/github/hapy/scripts/qc_duplicates.py merged_results.fits
```

## Validation

```
python ~/github/hapy/scripts/validate_measurements.py merged_results.fits --sample ALL
```
then examine plots in `validation/` subdirectory.


## Quick-Look Science

```
python ~/github/hapy/scripts/science_firstlook.py merged_results.fits --sample ALL
```
