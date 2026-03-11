# hapy

**hapy** (pronounced *happy*) is a Python package for astronomical image processing and survey-scale galaxy analysis.

Originally developed for **Hα imaging surveys**, hapy now provides a **cutout-based analysis pipeline** including:

* automated masking
* elliptical photometry
* structural measurements
* statmorph morphology metrics
* GALFIT modeling

Each galaxy is processed independently, enabling **robust parallel processing of large surveys**.

---

# Overview

The typical survey workflow is:

1. Generate galaxy cutouts from coadded images
2. Run automated analysis for each galaxy
3. Merge per-galaxy results into a survey-level table
4. Inspect results using automatically generated webpages

Because every galaxy is analyzed independently, the pipeline scales well using **GNU Parallel or cluster schedulers**.

---

# 🚀 Quickstart

```bash
# 1️⃣ Generate galaxy cutouts
python scripts/get_cutouts.py \
    --rimage /path/to/R_coadd.fits \
    --catalog /path/to/catalog.fits \
    --scheme virgo \
    --outdir survey_run

# 2️⃣ Run automated analysis (parallel over galaxies)
python scripts/run_analysis.py \
    --root survey_run/cutouts/<galaxy_tag>/<galaxy_tag> \
    --make-mask --statmorph --galfit --convflag

# 3️⃣ Merge per-galaxy results
python scripts/merge_results.py \
    --indir survey_run/cutouts \
    --out merged_results.fits \
    --outdir survey_run
```

After this step you will have a **survey-level results table** ready for analysis.

---

# ⚡ Parallel Processing Example

After generating cutouts, run analysis on all galaxies:

```bash
find survey_run/cutouts -mindepth 1 -maxdepth 1 -type d | \
parallel -j 8 \
'python scripts/run_analysis.py \
    --root {}/$(basename {}) \
    --make-mask --statmorph --galfit --convflag'
```

This runs **8 galaxies simultaneously** (adjust to match available CPU and RAM).

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

Optional:

* statmorph
* GNU Parallel (recommended for large surveys)

---

# Survey Workflow

---

## 1️⃣ Generate Cutouts

Create per-galaxy cutouts from a coadded image and catalog.

```bash
python scripts/get_cutouts.py \
    --rimage /path/to/R_coadd.fits \
    --catalog /path/to/catalog.fits \
    --scheme virgo \
    --outdir survey_run
```

Output structure:

```
survey_run/
    cutouts/
        <tag>/
            metadata.json
            <tag>-R.fits
            <tag>-Ha.fits
            <tag>-CS-ZP.fits
```

Each galaxy directory contains:

* image cutouts
* metadata describing the observation and initial ellipse parameters

---

## 2️⃣ Run Automated Analysis

Run masking, photometry, statmorph, and GALFIT.

```bash
python scripts/run_analysis.py \
    --root survey_run/cutouts/<tag>/<tag> \
    --make-mask \
    --statmorph \
    --galfit --convflag
```

Outputs include:

```
cutouts/<tag>/
    metadata.json
    <tag>-results.ecsv
    <tag>-mask.fits
    <tag>-diagnostic.png
```

The `results.ecsv` file contains all measurements and status flags for that galaxy.

---

## 3️⃣ Merge Results

Merge all per-galaxy results into a survey table.

```bash
python scripts/merge_results.py \
    --indir survey_run/cutouts \
    --mode run_analysis \
    --out merged_results.fits
```

Each row corresponds to **one observation of one galaxy**.

---

# Result Inspection Webpages  ⭐ 

hapy can generate inspection webpages for each cutout.

These pages show:

* Legacy Survey color images
* R and Hα cutouts
* photometry profiles
* statmorph results
* GALFIT model images
* pipeline status flags

Typical workflow:

### Download Legacy images

```bash
python scripts/fetch_legacy_cutouts.py \
    --cutout-root survey_run/cutouts
```

### Build webpages

```bash
parallel -j 8 \
python scripts/build_web_cutouts.py \
    --cutoutdir survey_run/cutouts \
    --oneimage "{}" \
    --outdir survey_run/html/cutouts \
:::: cutout_list.txt
```

### Build index page

```bash
python scripts/build_cutout_index.py \
    --runroot survey_run
```

Final structure:

```
survey_run/
    cutouts/
        <tag>/
    html/
        cutouts/
            <tag>/
                <tag>.html
            index.html
```

The `index.html` page provides a **dashboard view of the entire survey**.

---

# Pipeline Status Flags  

Each galaxy includes boolean flags describing pipeline success.

Important flags include:

| Flag            | Meaning                                |
| --------------- | -------------------------------------- |
| `MASK_OK`       | segmentation mask successfully created |
| `PHOT_OK`       | elliptical photometry succeeded        |
| `PSF_OK`        | PSF successfully built                 |
| `R_PROFILE_OK`  | R-band profile computed                |
| `HA_PROFILE_OK` | Hα profile computed                    |
| `R_SM_FLAG`     | statmorph measurement succeeded (R)    |
| `H_SM_FLAG`     | statmorph measurement succeeded (Hα)   |
| `GAL_NC_OK`     | GALFIT non-convolved fit succeeded     |
| `GAL_CV_OK`     | GALFIT PSF-convolved fit succeeded     |

These flags are used for **quality control and sample selection**.

---

# Masking & GUI Tools

Interactive mask editing:

```bash
python scripts/run_maskgui.py
```

Features:

* segmentation mask creation
* Gaia star masking
* mask growth tools
* interactive editing (Qt-based)

See `docs/masking.md` for details.

---

# Conventions

## Pixel Coordinate System

* +x axis = West
* +y axis = North

## Position Angles

Internal convention:

* `PA_DEG` = degrees CCW from North

Conversion to photutils angle:

```python
theta_deg = (90 + PA_DEG) % 180
```

Ellipse angles are **180° periodic**.

---

# Project Structure
```
hapy/
│
├── scripts/
│   ├── get_cutouts.py
│   ├── run_analysis.py
│   ├── merge_results.py
│   ├── build_web_cutouts.py
│   ├── build_cutout_index.py
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

* **core analysis engines**
* **GUI tools**
* **survey orchestration scripts**

This enables:

* automated survey processing
* interactive inspection
* reproducible science pipelines

---

# Roadmap   

Planned improvements include:

* improved documentation
* enhanced QC and analysis scripts
* automated survey summary plots
* improved pipeline architecture and module interfaces

---
