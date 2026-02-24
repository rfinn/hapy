# Ellipse Photometry in HAPY

## Overview

The `hapy.ellipse` package implements elliptical/isophotal photometry for R-band and continuum-subtracted images, optionally using masks produced by `hapy.masktools`. It is intended to run either:

- headless via the CLI (`python -m hapy.scripts.run_photometry ...`)
- from higher-level workflows (e.g., a future `HaSession` pipeline)
- with optional structural/morphology measurements

The design goal is to keep photometry logic separate from:
- image download/plotting utilities (`hapy.imagetools`)
- mask generation (`hapy.masktools`)
- GALFIT modeling (`hapy.galfittools`)

---

## Main entry point

```python
from hapy.ellipse.photometry import run_ellipse_photometry

You’re right — I accidentally made that whole section look like something you should paste verbatim, and the “end of the md” didn’t clearly close out as a single copyable block.

Here’s a clean version: **one copy/paste-ready `docs/ellipse.md`**, and everything after that is just normal chat commentary (not part of the file).

````markdown
# Ellipse Photometry in HAPY

## Overview

The `hapy.ellipse` package implements elliptical/isophotal photometry for R-band and continuum-subtracted images, optionally using masks produced by `hapy.masktools`. It is intended to run either:

- headless via the CLI (`python -m hapy.scripts.run_photometry ...`)
- from higher-level workflows (e.g., a future `HaSession` pipeline)
- with optional structural/morphology measurements

The design goal is to keep photometry logic separate from:
- image download/plotting utilities (`hapy.imagetools`)
- mask generation (`hapy.masktools`)
- GALFIT modeling (`hapy.galfittools`)

---

## Main entry point

```python
from hapy.ellipse.photometry import run_ellipse_photometry
```

This is the function called by `hapy/scripts/run_photometry.py`.

---

## Inputs

Typical inputs:

* `r_fits`: R-band FITS image (required)
* `cs_fits`: continuum-subtracted FITS image (optional)
* `mask_fits`: FITS mask (optional)
* `objra`, `objdec`: sky position (optional; can be used to select the target in catalogs)
* `fixcenter`: if set, hold the center fixed during geometry/photometry
* `run_statmorph`: run morphology measurements (optional)
* `write_prefix`: prefix tag for output tables

---

## High-level workflow

1. Read image(s) and WCS/header metadata
2. Load or build a mask (if provided)
3. Detect sources / build segmentation (if needed)
4. Select target object (based on RA/Dec or “central object” heuristic)
5. Estimate initial geometry (center, PA, axis ratio, size)
6. Perform elliptical photometry / profile extraction
7. Compute derived quantities (enclosed flux, radii, etc.)
8. Write output tables + optional diagnostic plots

---

## Outputs

Common products (names may vary with `write_prefix`):

* `*_phot.fits` : FITS table of radial profile / aperture photometry
* `*_phot.dat`  : ASCII summary
* `*_phot_apertures.png` : diagnostic apertures plot
* `*_enclosed_flux_fancy.png` : curve-of-growth plot

---


