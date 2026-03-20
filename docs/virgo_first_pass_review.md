
# HAPY Virgo First-Pass Review

This document records issues identified after the **first full run of `run_analysis` on the Virgo sample**.

The goals of the fixes below are:

1. Ensure **scientific correctness of measurements**
2. Improve **pipeline diagnostics and robustness**
3. Improve **cutout webpage usability for inspection**

---

# Priority 1 — Science-Critical Issues

These issues may affect **measured quantities** and should be addressed **before the next full pipeline run**.

## `get_cutouts.py`

### Exposure Time Handling

GALFIT magnitudes appear offset due to non-unit exposure times in the coadd headers.

**Fix**

* Save original exposure time as:

```
EXPTCOAD
```

* Set:

```
EXPTIME = 1
```

in the cutout headers.

This ensures:

* GALFIT magnitudes are computed correctly
* photometry uses flux-calibrated units.

---

### Gain Handling

Verify whether gain propagation should be:

```
newgain = oldgain * exptime
```

Confirm correct gain definition for:

* photutils
* GALFIT
* statmorph

---

## `run_analysis.py`

### Statmorph Flag Behavior

Observed issue:

* Some Hα statmorph runs produce **valid output tables**
* but `H_SM_FLAG` is **False**

Investigate:

* whether the flag logic is too strict
* whether a later failure resets the flag
* whether some output values should be treated as unreliable.

---

### GALFIT Success Flag

Clarify how the GALFIT success flag is derived.

Determine whether it differs from:

```
GAL_NUMERR
```

If so:

* propagate the **numerical error flag** into the results table
* expose both flags in the pipeline outputs.

---

### Masked Fraction of Initial Ellipse

Current pipeline stores masked fraction of the **photutils-derived ellipse**, but this can fail when a galaxy is close to a bright star.

Add measurement of:

```
masked_fraction_initial_ellipse
```

using the **input ellipse from metadata.json**.

This should be stored in the results table and propagated to webpages.

---

### Write Results After Each Phase

Ensure that results are written after each major stage:

```
mask
photometry
profiles
statmorph
galfit NC
galfit CV
```

This prevents loss of diagnostic information when later stages fail.

---

# Priority 2 — Pipeline Diagnostics

These improvements will make debugging and quality control significantly easier.

## `run_analysis.py`

### Photutils Logging

Add detailed logging for:

* ellipse initialization
* failure conditions
* photometry termination conditions.

---

### Statmorph Logging

Add logs for:

* statmorph input parameters
* exceptions
* quality flags returned by statmorph.

Clarify when results exist but flags are false.

---

### Bright Star Warning

Add a diagnostic flag for galaxies affected by nearby bright stars.

Possible implementation:

* use Gaia stars
* apply a **magnitude–radius relation**
* flag stars within *N × radius* of the galaxy.

Minimum requirement:

* flag stars with magnitude ≤ 10 near the galaxy.

---

# Priority 3 — Webpage Issues

These issues affect **inspection and visualization** but do not impact measured values.

---

# Cutout Webpages (`build_web_cutouts.py`)

### Image Display

* Remove zoom on legacy and other images
* Fix display stretch (currently saturated)
* Ensure all images use the standard `display_image()` routine

---

### GALFIT Visualization

* Correct PA conversion between:

  * photutils
  * GALFIT

* Remove zoom in GALFIT mask window

* Verify ellipse position is correct.

---

### Image Statistics Panel

Replace:

```
run date
telescope
run
pointing
```

with:

```
parent coadd image
```

and link to the **coadd webpage**.

Add:

```
fraction of masked pixels
```

---

### Mask Fraction Flag

Add a warning flag when:

```
masked_fraction_initial_ellipse > 0.5
```

Use the **input ellipse**, not the photutils ellipse.

---

### Pipeline Status Table

Fix rendering issues for the **Hα symbol** in:

```
Halpha profile label
```

Note: the symbol renders correctly in the Image Statistics table.

---

# Cutout Index (`build_cutout_index.py`)

### FWHM Display Patch

Temporary fix until pipeline rerun:

* if `R_FWHM == 0`, use value from `R_FHWM`
* same for Hα.

---

### Column Changes

Replace:

```
TAG
```

with

```
VFID
```

---

### Column Order

Reorder pipeline columns to:

```
PSF
R FWHM
Hα FWHM
Mask
Phot
Profiles
Statmorph
GALFIT
```

---

### Split Band Flags

Display separate flags for:

```
PSF (R / Hα)
Phot (R / Hα)
Profiles (R / Hα)
Statmorph (R / Hα)
```

under shared headings.

---

### Split GALFIT Status

Replace single flag with:

```
GALFIT NC
GALFIT CV
```

---

### Hα Rendering

Fix Hα symbol rendering in the **cutout index table**.

---

# Legacy Image Fetching (`fetch_legacy_images`)

Some cutout pages are missing Legacy images.

Investigate:

* missing downloads
* incorrect paths
* incomplete batch runs.

---

# `summarize_results.py`

Reorder output columns to match pipeline execution order:

```
mask
phot
profiles
statmorph
galfit NC
galfit CV
```

This will make pipeline summaries easier to interpret.

---

# Recommended Fix Order

To minimize rework:

### Phase 1 — Pipeline correctness

1. `get_cutouts.py`
2. `run_analysis.py`
3. `summarize_results.py`

---

### Phase 2 — Diagnostics

4. logging improvements
5. bright star flag
6. masked ellipse fraction

---

### Phase 3 — Webpages

7. `build_web_cutouts.py`
8. `build_cutout_index.py`

---

### Phase 4 — Ancillary tools

9. `fetch_legacy_images`

---

# Goal Before Next Virgo Re-run

Before running the pipeline again, confirm:

* exposure time and gain handling are correct
* statmorph and GALFIT flags behave consistently
* masked fraction of the initial ellipse is tracked
* pipeline logs are sufficiently detailed for debugging

