## Meaning of Version Increments

#### Patch-level examples (+0.01)
- fix a bad filename
- fix one flag not being written correctly
- improve logging
- rerun after a narrow implementation bug

#### Minor-version examples (+0.10)
- fix units on statmorph radii
- add derived columns
- change how profiles are measured
- add corrected radii
- change masking behavior that affects science outputs

### Major-version examples (+1.00)
- redesign results table
- replace major parts of morphology/profile pipeline
- break backward compatibility

## v0.2 (2026-03-23)
- fixed pixel vs arcsec inconsistency in radial measurements; column names for radial sizes now explicitly have _ARCSEC
- added vr and redshift to metadata in `get_cutouts` and propagated through to `run_analysis`
- validated telescope-dependent offsets as a unit-conversion issue
- updated QC + validation framework
- flagged need for inclination correction (not yet applied)


## hapy v0.2.1 — Archive Image Support Update (2026-03-24)

### New Features
- Added support for **archival Hα imaging datasets** with non-standard formats
- Implemented `build_metadata_archive.py` to:
  - Ingest messy archive directories
  - Select correct **R / CS / mask** images from an authoritative file list
  - Standardize directory and filename structure to `<tag>` convention
  - Generate minimal, pipeline-compatible `metadata.json`

### Metadata Improvements
- Added support for **explicit file pointers** in metadata:
  - `r_fits`, `cs_fits`, `mask_fits`
- Included **VFID** as a primary identifier (alongside `objid`)
- Added **GALNAME** for easier downstream analysis
- Added **FWHM and filter metadata** from FITS headers:
  - `rimage_fwhm_psf_arcsec`, `himage_fwhm_psf_arcsec`
  - `rimage_filter`, `himage_filter`
- Added **shape provenance tracking**:
  - `shape_flag` (1 = valid, 0 = default/fallback)
  - `shape_source` (e.g., `shape_catalog`, `default`)

### Robust Handling of Missing Data
- Standardized missing shape parameters to:
  - `sma_arcsec = 0`, `ba = 0`, `pa_deg = 0`
  - Allows `run_analysis` to trigger fallback ellipse estimation
- Updated `run_analysis.py` to:
  - Safely handle missing metadata fields without crashing
  - Avoid `float(None)` errors via conditional assignments

### Archive-Specific Behavior
- Introduced `scheme = "archive"` to control behavior
- Automatically disables **Gaia star masking** for archive data
- Skips archive-inapplicable metadata fields (e.g., filter ratio, cutout scale)
- Supports operation without reliable WCS

### File & Directory Handling
- Added support for **separate raw and processed directories**:
  - `--archive-root` (input)
  - `--output-root` (standardized output)
- Switched to **copy/symlink workflow** (no destructive moves)
  - Enables safe reprocessing and debugging
- Improved filename parsing:
  - Handles inconsistent naming (e.g., leading `h` prefixes, mixed case `R/r`)

### Workflow Integration
- Archive datasets now fully compatible with existing pipeline:
  - `run_analysis.py` unchanged in core logic
  - GNU parallel workflow preserved
- Archive and survey data now share a **common metadata + directory interface**

### Overall
- Extended hapy from survey-only pipeline to a **general-purpose Hα analysis framework**
- Maintained backward compatibility with Virgo/AGC workflows


## hapy v0.2.2 — run_maskgui and added functionality (2026-03-27)

### Masking & GUI Updates

- Added interactive mask overlay modes (off / contour / filled with true mask footprint)
- Implemented synchronized pan (auto) and zoom (manual via `z`) across panels
- Added pixel inspection tool (`v`) showing r-band, Hα, and mask values at cursor
- Introduced CLI tools: `build_mask` (headless) and `run_maskgui` (interactive)
- Refactored masking architecture (separation of `masktools` and `maskgui`)
- Standardized masking arguments across pipeline and CLI
- Fixed overlay rendering (removed large bounding boxes; now pixel-accurate)
- Fixed GUI quit behavior and keybinding conflicts


## hapy v0.2.3 — improved centering and center-offset QC diagnostics (2026-04-17)

### Centering & Photometry Updates

- Reworked central-object centering in ellipse photometry to improve robustness
- Added flux-based centering options within the selected photutils segmentation region
- Added logic to compare peak-anchored and guess-anchored flux centers and adopt a best center
- Improved behavior for centrally concentrated galaxies while reducing some large center offsets
- Updated downstream geometry to use the improved adopted center where appropriate

### QC / Analysis Table Updates

- Added merged-table columns tracking center offsets between:
  - input coordinates and photutils center
  - input coordinates and GALFIT centers
  - photutils center and GALFIT centers
- Added center-offset warning and severe warning QC flags
- Incorporated center-offset diagnostics into QC / review-priority logic
- Added segmentation diagnostic plots showing:
  - r-band image
  - masked r-band image
  - SExtractor segmentation
  - mask image
  - photutils segmentation
- Fixed bug with gaia coverage when downloading catalogs for coadd
  images.  I was not accounting for cos(dec) term, so RA width of
  returned catalog was too small.  Changed code to use astropy's
  calc_footprint, which gives coordinates of the image corners.
- adding tracking of fraction of initial ellipse that has valid image
  data in both the r and halpha images
