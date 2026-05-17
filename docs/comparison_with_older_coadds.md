
# The newer reduction is clearly much better overall

Especially for:

* continuum structural quantities
* r-band radii
* same-telescope consistency
* morphology stability


| Quantity | Units | Old NMAD | New NMAD | Improvement Factor |
|---|---|---:|---:|---:|
| `R25_ARCSEC` | dex | 0.0571 | 0.0165 | 3.5× |
| `R50_ARCSEC` | dex | 0.0667 | 0.0174 | 3.8× |
| `R75_ARCSEC` | dex | 0.0858 | 0.0181 | 4.7× |
| `R24_FLUX_CGS` | dex | 0.0680 | 0.0728 | 0.9× |
| `R24_MAG` | mag | 0.0702 | 0.0803 | 0.9× |
| `R_C30` | value | 0.0127 | 0.0119 | 1.1× |
| `R_SM_GINI` | value | 0.0182 | 0.0117 | 1.6× |
| `R_HAPY_GINI` | value | 0.0814 | 0.0849 | 1.0× |
| `H50_ARCSEC` | dex | 0.1525 | 0.1195 | 1.3× |
| `H75_ARCSEC` | dex | 0.1923 | 0.1535 | 1.3× |
| `H_R24_FLUX_CGS` | dex | 0.1429 | 0.1715 | 0.8× |
| `H_TOT_FLUX_CGS` | dex | 0.1971 | 0.2352 | 0.8× |
| `H_ISO17E18_ARCSEC` | dex | 0.0613 | 0.0751 | 0.8× |
| `H_ISO17E18_FLUX_CGS` | dex | 0.1073 | 0.1827 | 0.6× |
| `H_HAPY_GINI` | value | 0.0327 | 0.0342 | 1.0× |
| `DELTA_GINI` | value | 0.0816 | 0.0705 | 1.2× |

Key takeaways:
- The newer reduction dramatically improves the reproducibility of continuum structural measurements, especially the radial size measurements (`R25`, `R50`, `R75`), with factors of ~4 improvement in NMAD.
- The new reductions also eliminate several catastrophic failures seen in the older reductions, especially for HDI duplicate pairs.
- Continuum morphology measurements (`R_SM_GINI`) also improve significantly.
- The integrated continuum flux measurements (`R24_FLUX_CGS`, `R24_MAG`) are similar between reductions, likely because cross-instrument filter differences dominate the scatter.
- Hα size measurements improve moderately in the newer reduction.
- Hα flux measurements do not improve overall and in some cases exhibit larger scatter, likely reflecting sensitivity to continuum subtraction and filter-transmission corrections.
- INT now stands out more clearly as the remaining noisier dataset, suggesting that many global pipeline failures have been removed and that the remaining scatter is likely instrument-specific.



# Most important result

## R-band structural quantities improved dramatically

### Old reduction

* `R50_ARCSEC` NMAD = 0.0667 dex 
* `R75_ARCSEC` NMAD = 0.0858 dex 

### New reduction

* `R50_ARCSEC` NMAD = 0.0174 dex 
* `R75_ARCSEC` NMAD = 0.0181 dex 

That is a factor of:

* ~4 improvement for R50
* ~5 improvement for R75


# Hα results improved for BOK and HDI but worse for INT



### `H_R24_FLUX_CGS` Duplicate NMAD Comparison

| Pair type | Old reduction | New reduction |
|---|---:|---:|
| BOK-BOK | 0.0545 | 0.0231 |
| HDI-HDI | 0.3949 | 0.0453 |
| BOK-HDI | 0.1499 | 0.0805 |
| INT-INT | 0.0480 | 0.2180 |
| BOK-INT | 0.1578 | 0.2511 |
| INT-HDI | 0.1422 | 0.2915 |
| INT-MOS | — | 0.3261 |

### `H_TOT_FLUX_CGS` Duplicate NMAD Comparison

| Pair type | Old reduction | New reduction |
|---|---:|---:|
| BOK-BOK | 0.0933 | 0.0364 |
| HDI-HDI | 0.4671 | 0.0632 |
| BOK-HDI | 0.1919 | 0.1304 |
| INT-INT | 0.0307 | 0.2925 |
| BOK-INT | 0.1390 | 0.2780 |
| INT-HDI | 0.2152 | 0.3019 |
| INT-MOS | — | 0.4016 |


Key point:
- The newer reduction dramatically stabilizes the BOK and HDI Hα total flux measurements.
- The INT Hα flux scatter got dramatically worse
- The increase in overall Hα flux scatter is dominated by duplicate pairs involving INT observations.


Using matched pairs between old and new reductions, I get the following statistics.

### `H_R24_FLUX_CGS` Duplicate NMAD Comparison

| Pair type | Old reduction | New reduction |
|---|---:|---:|
| BOK-BOK | 0.0993 | 0.0231 |
| HDI-HDI | 0.3949 | 0.0453 |
| BOK-HDI | 0.1002 | 0.0677 |
| INT-INT | 0.0480 | 0.2607 |

Key point:
- Non-INT duplicate consistency improves dramatically in the new reduction.
- INT-INT Hα flux reproducibility becomes substantially worse.

### `H_TOT_FLUX_CGS` Duplicate NMAD Comparison

| Pair type | Old reduction | New reduction |
|---|---:|---:|
| BOK-BOK | 0.1264 | 0.0364 |
| HDI-HDI | 0.4671 | 0.0632 |
| BOK-HDI | 0.1093 | 0.1108 |
| INT-INT | 0.0307 | 0.1404 |

Key point:
- BOK and HDI total Hα flux measurements become much more reproducible in the new reduction.
- INT-INT total Hα flux reproducibility worsens substantially.

### `H50_ARCSEC` Duplicate NMAD Comparison

| Pair type | Old reduction | New reduction |
|---|---:|---:|
| BOK-BOK | 0.2250 | 0.0075 |
| HDI-HDI | 0.7774 | 0.1979 |
| BOK-HDI | 0.1383 | 0.0411 |
| INT-INT | 0.0152 | 0.1029 |

Key point:
- The new reduction dramatically improves the reproducibility of Hα half-light radii for BOK and HDI observations.
- INT-INT Hα radius measurements become noisier in the new reduction, though not as severely as the integrated Hα fluxes.

# Focus on INT Data


## INT-INT R-band measurements improved 

## Example: R50

### Old:

* INT-INT = 0.0518 dex 

### New:

* INT-INT = 0.0109 dex 

That is  nearly a factor of 5 improvement.

## Scatter in R-band magnitudes is similar

* scatter is similarly large in old and new reductin.


Magnitudes remain dominated by filter-system differences,
while sizes became dramatically more reproducible.

## INT Halpha measurements got worse!

like a lot worse!

For matched/filter-cut duplicates:

|Quantity |	Old CS    |	New CS-ZP	| New CS-gr|
|---|---|---|---|
|`H_TOT_FLUX_CGS` | overall |	0.1264 |	0.1175 |	0.1070|
|INT-INT	|0.0307	|0.1404	|0.0900|

