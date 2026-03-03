# hapy

**hapy** (pronounced *happy*) is a Python package for astronomical image processing and survey-scale galaxy analysis.

Originally developed for Hα imaging, hapy now provides a full cutout-based analysis pipeline including masking, photometry, structural measurements, and GALFIT modeling.

---

# Overview

The typical survey workflow is:

1. Generate galaxy cutouts from coadded images
2. Run automated analysis (masking, photometry, statmorph, GALFIT)
3. Merge per-galaxy results into a survey-level table

Each galaxy is processed independently, enabling parallel execution.

---

# Core Features

## Cutout Generation

* Create galaxy cutouts from coadded survey images
* Store per-object metadata (`metadata.json`)
* Support multiple catalog schemes (Virgo, AGC, generic)

```bash
python scripts/get_cutouts.py \
    --rimage <R_COADD.fits> \
    --catalog <catalog.fits> \
    --scheme virgo
```

Each galaxy gets its own directory:

```
cutouts/<galaxy_tag>/
    metadata.json
    <tag>-R.fits
    <tag>-mask.fits
```

---

## Automated Analysis

Run masking, ellipse photometry, statmorph, and GALFIT:

```bash
python scripts/run_analysis.py \
    --root cutouts/<galaxy_tag>/<galaxy_tag> \
    --make-mask \
    --statmorph \
    --galfit
```

Outputs include:

* Mask FITS image
* Per-galaxy results table (`*-results.ecsv`)
* Optional diagnostic plot showing:

  * Input mask ellipse
  * Photutils-derived ellipse

GALFIT runs in two stages:

* NC (no convolution)
* CV (PSF convolution, optional)

Failures in the CV stage do not stop processing.

---

## Merge Results

Combine all galaxy results into a single table:

```bash
python scripts/merge_results.py --indir cutouts/
```

Each row corresponds to one independent observation.

---

# Masking & GUI Tools

Interactive mask editing is available:

```bash
python scripts/run_maskgui.py
```

Features:

* Segmentation mask creation
* Gaia star masking
* Mask growth tools
* Interactive editing (Qt-based)

See `docs/masking.md` for details.

---

# Conventions

## Pixel Coordinate System

* +x axis = West
* +y axis = North

## Position Angles

Internal convention:

* `PA_DEG` = degrees CCW from North

Conversion to photutils theta:

```python
theta_deg = (90 + PA_DEG) % 180
```

All ellipse comparisons are 180° periodic.

---

# Project Structure

```
hapy/
│
├── scripts/
│   ├── get_cutouts.py
│   ├── run_analysis.py
│   ├── merge_results.py
│   └── run_maskgui.py
│
├── hapy/
│   ├── hatools/
│   ├── catools/
│   ├── imagetools/
│   ├── masktools/
│   ├── maskgui/
│   ├── galfittools/
│   └── astromatic/
│
├── docs/
└── tests/
```

---

# Installation

```bash
git clone https://github.com/rfinn/hapy.git
cd hapy
pip install -e .
```

External dependencies:

* Source Extractor
* GALFIT (for modeling)

---

# Development Notes

* Each cutout directory contains a `metadata.json` file describing the object and input parameters.
* Output tables use a stable schema to support survey-level merging.
* Diagnostic plots can be enabled from `run_analysis.py`.
* Logging support is being added for batch-scale production runs.

---

# Philosophy

hapy separates:

* Core engines (masking, photometry, modeling)
* GUI tools
* Survey pipeline orchestration

This enables:

* Fully automated batch processing
* Interactive workflows
* Reproducible survey analysis

---
