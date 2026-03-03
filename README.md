

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

# 🚀 Quickstart

```bash
# 1️⃣ Generate galaxy cutouts
python scripts/get_cutouts.py \
    --rimage /path/to/R_coadd.fits \
    --catalog /path/to/catalog.fits \
    --scheme virgo \
    --outdir survey_run

# 2️⃣ Run automated analysis (repeat in parallel for all cutouts)
python scripts/run_analysis.py \
    --root survey_run/cutouts/<galaxy_tag>/<galaxy_tag> \
    --make-mask --statmorph --galfit

# 3️⃣ Merge all galaxy results into one table
python scripts/merge_results.py \
    --indir survey_run/cutouts \
    --out merged_results.fits \
    --outdir survey_run
```

After this, your survey directory will contain:

```
survey_run/
    cutouts/
        <galaxy_tag>/
            metadata.json
            <tag>-results.ecsv
            <tag>-diagnostic.png
    merged_results.fits
```

---

# ⚡ Parallel Processing Example

After generating cutouts, you can run analysis on all galaxies in parallel:

```bash id="n1hsu3"
find survey_run/cutouts -mindepth 1 -maxdepth 1 -type d | \
parallel -j 8 \
'python scripts/run_analysis.py \
    --root {}/$(basename {}) \
    --make-mask --statmorph --galfit'
```

Explanation:

* Each galaxy has its own directory inside `cutouts/`
* The cutout root is `<dir>/<basename>`
* `-j 8` runs 8 galaxies simultaneously (adjust for your machine)

---

If you prefer creating a list first:

```bash id="4wzqet"
find survey_run/cutouts -mindepth 1 -maxdepth 1 -type d > cutout_list.txt

parallel -j 8 \
'python scripts/run_analysis.py \
    --root {}/$(basename {}) \
    --make-mask --statmorph --galfit' \
:::: cutout_list.txt
```

---

This keeps the workflow:

1. Cutouts
2. Parallel analysis
3. Merge

clean and reproducible.

---


# Installation

```bash
git clone https://github.com/rfinn/hapy.git
cd hapy
pip install -e .
```

External dependencies:

* Source Extractor
* GALFIT
* PyQt5 (for GUI tools)

---

# Survey Workflow

---

## 1️⃣ Generate Cutouts

Create per-galaxy cutouts from a coadded image and catalog:

```bash
python scripts/get_cutouts.py \
    --rimage /path/to/R_coadd.fits \
    --catalog /path/to/catalog.fits \
    --scheme virgo \
    --outdir /path/to/survey_run
```

If `--outdir` is not provided, cutouts are written to the current working directory.

Output structure:

```
survey_run/
    cutouts/
        <galaxy_tag>/
            metadata.json
            <tag>-R.fits
            <tag>-mask.fits
```

Each galaxy directory contains a `metadata.json` file describing:

* Object ID
* Sky coordinates
* Initial ellipse parameters
* Parent image information
* PSF and image metadata

---

## 2️⃣ Run Automated Analysis

Run masking, ellipse photometry, statmorph, and GALFIT:

```bash
python scripts/run_analysis.py \
    --root survey_run/cutouts/<galaxy_tag>/<galaxy_tag> \
    --make-mask \
    --statmorph \
    --galfit
```

Outputs include:

* Updated mask image
* Per-galaxy results table (`*-results.ecsv`)
* Optional diagnostic plot showing:

  * Input mask ellipse
  * Photutils-derived ellipse

GALFIT runs in two stages:

* NC (no convolution)
* CV (PSF convolution, optional)

If convolution fails, processing continues and the failure is recorded in the results table.

---

## 3️⃣ Merge Results

Merge all per-galaxy results into a single table:

```bash
python scripts/merge_results.py \
    --indir survey_run/cutouts \
    --out merged_results.fits \
    --outdir survey_run
```

* `--indir` specifies where to search for `*-results.ecsv` files.
* `--outdir` specifies where the merged FITS table will be written.
* If `--outdir` is not provided, the merged table is written to the current directory.

Each row corresponds to one independent observation.

---

# Masking & GUI Tools

Interactive mask editing:

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

Ellipse angles are 180° periodic.

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

# Philosophy

hapy separates:

* Core engines (masking, photometry, modeling)
* GUI tools
* Survey orchestration

This enables:

* Fully automated batch processing
* Interactive workflows
* Reproducible survey analysis

---
