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
