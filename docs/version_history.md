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
  - Updated masktools.engine to grow mask BEFORE removing central
    object to avoid neighboring objects bleeding into galaxy
- updated the `HAPY_MORPH_OK` flag so that this can still be true if the halpha flux is zero.  No halpha is still a valid measurement and not a failure of the morphology code.



## HAPY v0.2.4 notes

### Major fixes

- Fixed aperture-noise calculation for image2 / Hα profiles.
  - `photometry.py` had been using image1 sky noise and gain for both images.
  - Added an explicit standalone aperture-noise function that takes `sky_noise` and `gain`.
  - Image2/Hα aperture errors now use the correct Hα sky noise and gain.
  - Verified on `VFID2313`; Hα profiles now produce valid values.
  - this bug was identified b/c most of 2019 INT Halpha phot profiles
    were not fit - calculated snr was much lower than correct value.

- Updated `run_analysis` to derive `FILTER_RATIO` from `PHOTZP` keywords.
  - Preferred source is now the `PHOTZP` ratio from `r_fits` and `cs_fits`.
  - Metadata `filter_ratio` is now fallback behavior.

### CS-gr support

- Added optional `CS-gr` processing in `run_analysis`.
  - Detects `*-CS-gr.fits` when present.
  - Runs ellipse photometry and HAPY morphology on CS-gr.
  - Writes `CSGR_*` output columns.
  - Stores `CSGR_CONTSCL` from the CS-gr FITS header.
  - Added `--csgr` flag to allow skipping CS-gr processing for faster testing/reruns.

- Added `fileid` support to avoid overwriting diagnostics/products from the primary CS-ZP run.
  - CS-gr HAPY morphology diagnostics use a separate namespace.

### Cutout / coadd handling

- Updated `get_cutouts` / `coadd_images.py` so CS-ZP cutouts are generated from matched R and Hα cutouts, not from parent CS-ZP coadds.
- Added Hα-defined slice sharing for hybrid coadds.
  - Hα cutout defines the pixel slice.
  - R cutout uses the same parent-image slice.
  - Prevents one-pixel shape mismatches when R and Hα WCS/pixel scales differ slightly.
- Added or improved skip logging for invalid cutout regions.
- Built simple footprint-style weight images for problematic INT hybrid coadds.
  - `weight=1` for finite nonzero science pixels.
  - `weight=0` for zero/invalid/off-CCD pixels.

### Manual masks

- Added/used `transfer_manual_masks.py`.
  - Copies BOK/HDI/MOS manual masks directly.
  - Reprojects INT manual masks onto the new Hα cutout grid.
  - Preserves mask labels/object IDs rather than collapsing all masked regions to 1.
  - Handles strict ±1 day tag/date shifts for INT cutouts.

### Segmentation / CS-gr prerequisites

- Added `hapy/hatools/segmentation.py`.
  - Includes `make_simple_photutils_segmentation`.
  - Used by `make_cs_gr.py` when `*-R-phot-segmentation.fits` is missing.

### Legacy image download robustness

- Updated `get_legacy_images` logic.
  - For combined `grz` requests, checks for final `g/r/z` products rather than trusting stale `grz.fits`.
  - Removes stale combined files before retrying.
- Added/considered retry/backoff handling for transient Legacy Survey HTTP errors.

### Merge/results schema

- Updated `run_analysis` initialization for optional CS-gr columns.
  - Ensures rows without CS-gr still include expected `CSGR_*` columns.
- Updated `merge_results` schema handling.
  - Can fill missing columns using safe/default values.
  - Better support for partial result tables from before CS-gr schema was finalized.







# HAPY v0.3.0 (2026-06-12)

### Major changes:
- Morphology segmentation now uses a physical SB noise floor.
- Fixed image2 aperture-noise calculations.
- FILTER_RATIO derived from PHOTZP.
- Improved CS-gr support and schema handling.
- Added robust fallback photutils segmentation generation.

### Impact:
- Morphology measurements (Gini, M20, asymmetry, fill fraction)
  may differ from earlier releases.
- Hα profile uncertainties are more accurate.

### HAPY Morphology Updates 

- Morphology segmentation thresholds are now based on a minimum physical
  surface-brightness noise floor rather than purely the measured image noise.

- For each image, the effective noise used to build the morphology
  segmentation is:

      `sigma_eff = max(measured_sigma_SB, sigma_floor_SB)`

  where sigma_SB is the sky noise in physical surface-brightness units.

- Initial floor values:
      `Hα: 5e-17 erg s^-1 cm^-2 arcsec^-2`
      `R : 4e-16 erg s^-1 cm^-2 arcsec^-2`

- The adopted physical threshold is converted back to ADU before creating
  the segmentation image, preserving compatibility with the existing
  HAPY morphology code.

- This prevents very deep images (especially BOK) from using excessively
  low segmentation thresholds while still allowing noisier images (e.g.
  some HDI observations) to use their measured sky noise.

- Goal: reduce instrument-dependent systematics in Gini, M20, asymmetry,
  fill fraction, and segmentation area by measuring morphology at a more
  consistent physical surface-brightness limit across the survey.



# HAPY v0.3.1 (2026-06-16) 
- Added a minimum cutout size of 90 arcsec. For galaxies whose scaled cutout
would be smaller than this, the cutout is enlarged to ensure enough sky area
for background estimation, segmentation, morphology, and visual review.
Larger galaxies retain the historical cutout_scale × diameter sizing.
- Defaults: min_cutout_size = 75 arcsec, cutout_buffer = 75 arcsec.
- Fixed R-band weight cutout creation when using Hα-defined slices.
- Added CS-ZP weight cutouts from min(Ha weight, R weight).


## HAPY v0.4.0 (2026-06-27)

### Added
- Added `hapy/ellipse/clumps.py` for H-alpha clump detection and measurement.
- Added segmentation-based H-alpha clump catalogs within the central R-band segmentation footprint.
- Added optional local peak finding for H-alpha clump structure diagnostics.
- Added H-alpha clump summary quantities for inclusion in HAPY result catalogs.
- Added saved clump products:
  - `*-halpha-clumps.ecsv`
  - `*-halpha-clumps-segm.fits`
  - `*-halpha-clump-peaks.ecsv`
  - `*-halpha-clumps-diagnostic.png`
- Added `EllipsePhotometry.measure_halpha_clumps()`.
- Added `run_analysis.py --clumps` and related clump-analysis configuration options.

### Changed
- HAPY pipeline can now optionally quantify internal H-alpha star-forming structure.

### Notes
- Clump analysis is optional and off by default.
- Default clump detection parameters should be treated as provisional pending validation on test galaxies.


## HAPY v0.4.1 (2026-06-28)

### Changed
- Integrate CS-gr image creation into run_analysis after mask generation.
- Refactor `make_cs_gr.py` into callable `make_cs_gr_image()`.
- Use HAPY masks during CS-gr continuum scaling while preserving continuous output images.
- Add graceful handling for missing Legacy g/r reprojections.
- Improve robustness of CS-gr continuum-scale estimation.
