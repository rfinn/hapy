# HAPY running-analysis results log

This file records dated outputs, diagnostics, and interpretation from HAPY reruns. Keep the reusable commands in `running-analysis-rerun-workflow.md`.

## Summary table

| Date / run | Coadds | Cutouts / rows | Key output | Notes |
|---|---:|---:|---|---|
| 2026-03-09 | 227 total; 211 ready | 768 cutout dirs | `merged_cutouts_results.fits` had 768 rows, 19 columns | Not all coadds were ready. |
| 2026-03-30 | 227 | 823 cutout dirs | All coadds ready | Later merged cutouts had 831 rows, 28 columns. |
| 2026-04-17 | — | 823 run rows | `PHOT_OK=809`, `HAPY_MORPH_OK=771`, `PROFILES_BOTH=652` | 14 bad-phot objects. QC high-priority count was 141. |
| 2026-04-21 | — | 821 run rows | QC: high=141, medium=402, low=278 | High-priority drivers included bad photometry, failed morphology, bright-star flag, mask warnings, and severe centering warnings. |
| 2026-05-17 / 2026-05-18 | 226 | 853 cutout dirs; 814 merged analysis rows | `merged_results_virgo_20260518.fits`, 390 columns | Same cutout count as 2026-06-12 before the min-size change. |
| 2026-05-19 hybrid | 223 | 764 merged analysis rows | `merged_results_virgo_20260519.fits`, 390 columns | Hybrid sample. |
| 2026-06-09 hybrid | 223 | 856 cutout dirs | First pass hybrid run | Cutout count differed from the previous run. Worth tracking why. |
| 2026-06-12 | 226 | 853 cutout dirs / 853 analysis rows | `PHOT_OK=847`, `HAPY_MORPH_OK=844`, `H_PROFILE_OK=780` | Big improvement in H profiles. GALFIT success was low and needs follow-up. |
| 2026-06-20 | 226 | 849 cutout dirs | After implementing minimum cutout size in `get_cutouts` | Five Legacy cutouts still missing; same parent-coadd issue. |

## 2026-03-09

Coadd list counts:

```text
227 fullpath_rcoadds_all.txt
211 fullpath_rcoadds_hapy_ready.txt
```

Cutout check:

```text
Input coadds:              211
Cutout directories:        768
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

Merged cutout results:

```text
Found 211 result files.
Final table rows: 768
Final table columns: 19
```

## 2026-03-30

All coadds were ready.

Cutout check:

```text
Input coadds:              227
Cutout directories:        823
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

Merged cutout results after 2026-03-30:

```text
Found 227 result files.
Final table rows: 831
Final table columns: 28
```

## 2026-04-17

Run summary:

```text
Total galaxies: 823
Unique galaxies: 659
PSF_OK            :  823 OK  |    0 FAIL  (100.0%)
MASK_OK           :  823 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  809 OK  |   14 FAIL  ( 98.3%)
HAPY_MORPH_OK     :  771 OK  |   52 FAIL  ( 93.7%)
R_PROFILE_OK      :  808 OK  |   15 FAIL  ( 98.2%)
H_PROFILE_OK      :  652 OK  |  171 FAIL  ( 79.2%)
PROFILES_BOTH     :  652 ( 79.2%)
STATUS counts: ok=809, running=14
STAGE counts: done=809, mask=14
Number with bad phot = 14
```

QC output:

```text
Read 823 rows from merged_results.fits
REVIEW_PRIORITY SUMMARY
high=141, medium=426, low=256
ELL_MISMATCH 269
FILTER_WARNING 79
WARN_MASK 29
BRIGHT_STAR_FLAG 9
WARN_WEAK_HA 297
```

Note: after fixing the p012-p013 issue, the bad-phot count decreased by one in a later summary.

## 2026-04-21

QC output:

```text
Read 821 rows from merged_results_virgo_20260421.fits
REVIEW_PRIORITY SUMMARY
high=141, medium=402, low=278
ELL_MISMATCH 265
FILTER_WARNING 79
WARN_MASK 17
BRIGHT_STAR_FLAG 9
WARN_WEAK_HA 296
```

High-priority drivers:

```text
NOT_PHOT_OK                 : total=  12  in_high=  12
NOT_HAPY_MORPH_OK           : total=  48  in_high=  48
BRIGHT_STAR_FLAG            : total=   9  in_high=   9
WARN_MASK                   : total=  17  in_high=  17
SEVERE_CEN_ANY              : total=  89  in_high=  89
WARN_CUTOUT_MISSING_SHAPE   : total=   0  in_high=   0
```

Medium-priority drivers:

```text
ELL_MISMATCH                : total= 265  in_medium= 153
WARN_WEAK_HA                : total= 296  in_medium= 230
WARN_CEN_ANY                : total= 217  in_medium= 118
WARN_R_PROFILE_PEAK         : total= 142  in_medium=  81
WARN_CUTOUT_MISSING         : total=   0  in_medium=   0
```

## 2026-05-17 / 2026-05-18

Cutout check:

```text
Input coadds:              226
Cutout directories:        853
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

Merged cutout results:

```text
Found 226 result files.
Merging 226 result files.
Final table rows: 853
Final table columns: 28
```

Merged run-analysis results:

```text
Found 814 result files.
Skipping 0 files with CATALOG_USE == EXCLUDE
Merging 814 result files.
Final table rows: 814
Final table columns: 390
```

## 2026-05-19 hybrid sample

Merged run-analysis results:

```text
Found 764 result files.
Skipping 0 files with CATALOG_USE == EXCLUDE
Merging 764 result files.
Final table rows: 764
Final table columns: 390
```

## 2026-06-09 hybrid first pass

Cutout check:

```text
Input coadds:              223
Cutout directories:        856
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

Note: the number of cutouts did not match the previous run. Track whether this is caused by coadd coverage, cutout validity checks, catalog changes, or minimum-size logic.

## 2026-06-12

Cutout check:

```text
Input coadds:              226
Cutout directories:        853
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

Run summary:

```text
Total galaxies: 853
Unique galaxies: 675
PSF_OK            :  853 OK  |    0 FAIL  (100.0%)
MASK_OK           :  853 OK  |    0 FAIL  (100.0%)
PHOT_OK           :  847 OK  |    6 FAIL  ( 99.3%)
HAPY_MORPH_OK     :  844 OK  |    9 FAIL  ( 98.9%)
R_PROFILE_OK      :  847 OK  |    6 FAIL  ( 99.3%)
H_PROFILE_OK      :  780 OK  |   73 FAIL  ( 91.4%)
R_SM_OK           :  757 OK  |   96 FAIL  ( 88.7%)
H_SM_OK           :  757 OK  |   96 FAIL  ( 88.7%)
GAL_NC_OK         :  264 OK  |  589 FAIL  ( 30.9%)
GAL_CV_OK         :  261 OK  |  592 FAIL  ( 30.6%)
R_PETRO_OK        :  799 OK  |   54 FAIL  ( 93.7%)
R_EXPFIT_OK       :  847 OK  |    6 FAIL  ( 99.3%)
CSGR_H_PROFILE_OK :  810 OK  |   43 FAIL  ( 95.0%)
CSGR_R_PROFILE_OK :  847 OK  |    6 FAIL  ( 99.3%)
PROFILES_BOTH     :  780 ( 91.4%)
STATMORPH_BOTH    :  224
```

Status snapshot:

```text
STATUS counts: running=513, ok=340
STAGE counts: phot=507, done=340, init=6
Runtime medians: PHOT_SEC=9.04, SM_SEC=4.16
Number with bad phot = 6
```

Interpretation / TODO:

- H-profile success improved substantially relative to April.
- GALFIT failures are numerous and should be revisited.
- The status/stage counts suggest some jobs may not have completed cleanly or result states need interpretation.

## 2026-06-20

After implementing a minimum cutout size in `get_cutouts`:

```text
Input coadds:              226
Cutout directories:        849
Coadds with no cutouts:    0
Cutout dirs missing R:     0
Cutout dirs missing CS:    0
Bad coadd names:           0
Bad cutout dir names:      0
```

Missing Legacy files:

```text
VFID3297-SDSSJ091525.83+252511.1-INT-20190208-p024
VFID6386-WISEAJ150825.51+014224.3-INT-20190601-p033
VFID6392-WISEAJ150819.90+014123.6-INT-20190601-p033
VFID6447-SDSSJ150812.35+012959.7-INT-20190601-p033
VFID6463-WISEAJ150809.13+012516.6-INT-20190601-p033
```

Note: these are the same missing cutouts as before and appear to be related to the parent coadd.

Manual-mask note: mask adjustments were finished for galaxies that were flagged as needing additional editing.

## Open issues / follow-up

- Resolve why hybrid cutout counts differ between runs.
- Revisit GALFIT failures, especially in the 2026-06-12 run.
- Confirm whether `STATUS=running` with `STAGE=phot` reflects incomplete jobs, stale result fields, or expected intermediate output.
- Keep tracking missing Legacy files tied to problematic parent coadds.
