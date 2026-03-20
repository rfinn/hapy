# HAPY Results Table

## Overview

The HAPY results table contains one row per processed galaxy cutout.

Each row records:
- object identity and provenance
- input file paths and filter metadata
- pipeline stage completion flags
- masking and contamination diagnostics
- geometry and ellipse quantities
- image-quality and sky measurements
- custom HAPY morphology measurements
- statmorph outputs
- radial-profile and photometric quantities
- Hα extent and flux measurements
- GALFIT structural fits
- derived QC and science columns

This table is the primary per-object summary product of the HAPY pipeline.

---

## Row definition

Each row corresponds to:

> one processed observation of one galaxy cutout

This is not necessarily one unique galaxy.  
A galaxy may appear in multiple rows if it has duplicate observations.

---

## Naming conventions

The column names follow a few common conventions.

### Band prefixes
- `R_` : r-band / stellar-continuum quantity
- `H_` : Hα / continuum-subtracted quantity

### Status and diagnostics
- `_OK` : boolean success flag for a pipeline stage or measurement family
- `_FLAG` : integer diagnostic or bitmask flag
- `_WARN` : boolean warning flag

### Units
- `_ARCSEC` : angular size in arcseconds
- `_PIX` : quantity in pixels
- `_MAG` : magnitude
- `_CGS` : flux in CGS units
- `_ERR` : uncertainty or formal fit error

### Measurement families
- `_SM_` : statmorph output
- `_PETRO_` : Petrosian quantity
- `_EXPFIT_` : exponential-fit quantity from radial profile
- `_LOGFIT_` : log-linear fit quantity from radial profile

---

## General notes

- Angular quantities are in **arcseconds** unless otherwise specified.
- Pixel-based quantities use the cutout image pixel scale.
- Missing or invalid values are typically stored as `NaN`.
- Boolean `_OK` columns indicate whether a measurement completed successfully, not necessarily whether it is scientifically trustworthy.
- QC and science scripts may define additional derived columns such as `QC_TIER`, `WARN_*`, and size ratios.

---

# 1. Identity and provenance

## Core identifiers

### `VFID`
- **Meaning:** Virgo Filament Survey identifier
- **Type:** string

### `GALNAME`
- **Meaning:** Common or catalog galaxy name
- **Type:** string

### `OBJID`
- **Meaning:** Survey/object identifier used in upstream catalogs
- **Type:** string

### `RA`, `DEC`
- **Meaning:** Adopted sky position of the galaxy center
- **Units:** degrees
- **Type:** float

## Run provenance

### `HAPY_VERSION`
- **Meaning:** HAPY pipeline version used to produce the row
- **Type:** string

### `RUN_DATE`
- **Meaning:** Date when the pipeline run was executed
- **Type:** string

### `TELESCOPE`
- **Meaning:** Telescope or instrument associated with the observation
- **Type:** string

### `DATEOBS`
- **Meaning:** Observation date
- **Type:** string

### `POINTING`
- **Meaning:** Pointing or field identifier
- **Type:** string

### `SCHEME`
- **Meaning:** Processing scheme or sample type (e.g., `virgo`, `agc`)
- **Type:** string

### `TAG`
- **Meaning:** Unique tag used to identify the cutout/run products
- **Type:** string

### `CUTDIR`
- **Meaning:** Directory containing per-object output products
- **Type:** string

---

# 2. Input files and filter metadata

## Parent images

### `PARENT_RIMAGE`
- **Meaning:** Original parent r-band image used to generate the cutout
- **Type:** string

### `PARENT_HIMAGE`
- **Meaning:** Original parent Hα image used to generate the cutout
- **Type:** string

## Processed image products

### `MASK_FITS`
- **Meaning:** Final mask FITS file used in analysis
- **Type:** string

### `PSF_FITS`
- **Meaning:** PSF image used for modeling and/or diagnostics
- **Type:** string

### `R_FITS`
- **Meaning:** r-band cutout FITS file used in the pipeline
- **Type:** string

### `CS_FITS`
- **Meaning:** continuum-subtracted Hα FITS file
- **Type:** string

### `SIGMA_FITS`
- **Meaning:** sigma/noise image or images used by later analysis steps
- **Type:** string / JSON-like serialized content

## Filter metadata

### `HFILTER`
- **Meaning:** Hα filter name or label
- **Type:** string

### `RFILTER_FILENAME`
- **Meaning:** filename of the r-band filter transmission curve
- **Type:** string

### `RFILTER_CENTER`
- **Meaning:** effective/central wavelength of the r-band filter
- **Units:** Angstrom
- **Type:** float

### `RFILTER_WIDTH`
- **Meaning:** effective width of the r-band filter
- **Units:** Angstrom
- **Type:** float

### `HFILTER_FILENAME`
- **Meaning:** filename of the Hα filter transmission curve
- **Type:** string

### `HFILTER_CENTER`
- **Meaning:** effective/central wavelength of the Hα filter
- **Units:** Angstrom
- **Type:** float

### `HFILTER_WIDTH`
- **Meaning:** effective width of the Hα filter
- **Units:** Angstrom
- **Type:** float

### `PSF_SOURCE`
- **Meaning:** source of the PSF model/image
- **Type:** string

---

# 3. Pipeline stage status and timing

## Global stage/status fields

### `STAGE`
- **Meaning:** last completed pipeline stage or run state
- **Type:** string

### `STATUS`
- **Meaning:** overall run status
- **Type:** string
- **Typical values:** `ok`, failure states, or partial-completion states

## Runtime columns

### `MASK_SEC`
- **Meaning:** elapsed runtime for masking stage
- **Units:** seconds

### `PHOT_SEC`
- **Meaning:** elapsed runtime for photometry/profile stage
- **Units:** seconds

### `SM_SEC`
- **Meaning:** elapsed runtime for statmorph stage
- **Units:** seconds

### `GAL_NC_SEC`
- **Meaning:** elapsed runtime for GALFIT non-convolved fit
- **Units:** seconds

### `GAL_CV_SEC`
- **Meaning:** elapsed runtime for GALFIT convolved fit
- **Units:** seconds

### `TOTAL_SEC`
- **Meaning:** total runtime for this object
- **Units:** seconds

## Success flags

### `PSF_OK`
- **Meaning:** PSF-related setup completed successfully

### `MASK_OK`
- **Meaning:** mask generation completed successfully

### `PHOT_OK`
- **Meaning:** photometry stage completed successfully

### `HAPY_MORPH_OK`
- **Meaning:** custom HAPY morphology (e.g., custom Gini/M20 workflow) completed successfully

### `R_PROFILE_OK`, `H_PROFILE_OK`
- **Meaning:** radial-profile analysis succeeded in r-band / Hα

### `R_SM_OK`, `H_SM_OK`
- **Meaning:** statmorph analysis succeeded in r-band / Hα

### `GAL_NC_OK`, `GAL_CV_OK`
- **Meaning:** GALFIT non-convolved / convolved fit succeeded

### `GAL_CV_INIT_FROM_NC`
- **Meaning:** GALFIT convolved fit was initialized from the non-convolved solution

### `GAL_NC_RERUN_FIXEDN`, `GAL_CV_RERUN_FIXEDN`
- **Meaning:** GALFIT was rerun with a fixed Sérsic index after an initial problem or instability

### `ELL_MISMATCH`
- **Meaning:** warning that adopted ellipse/geometry solutions are inconsistent
- **Notes:** often useful as a QC warning

### `R_PETRO_OK`, `R_EXPFIT_OK`, `R_LOGFIT_OK`
- **Meaning:** specific r-band profile-derived quantities succeeded

### `H_PETRO_OK`, `H_EXPFIT_OK`, `H_LOGFIT_OK`
- **Meaning:** specific Hα profile-derived quantities succeeded

---

# 4. Masking and contamination diagnostics

## Bright star diagnostics

### `BRIGHT_STAR_FLAG`
- **Meaning:** a bright star is close enough to potentially affect the cutout or measurements

### `BRIGHT_STAR_DIST`
- **Meaning:** distance to nearest bright star
- **Units:** likely arcseconds

### `BRIGHT_STAR_MASKRAD_ARCSEC`
- **Meaning:** adopted mask radius for the nearest bright star
- **Units:** arcseconds

### `BRIGHT_STAR_MAG`
- **Meaning:** brightness of the bright star used in contamination logic
- **Units:** magnitude

## Central ellipse masking

### `ELL0_MASKFRAC`
- **Meaning:** fraction of masked pixels within the adopted central ellipse

### `ELL0_MASK_WARN`
- **Meaning:** boolean warning that masking in the central ellipse is large enough to be concerning

### `ELL0_NMASKPIX`
- **Meaning:** number of masked pixels within the central ellipse

### `ELL0_NTOTPIX`
- **Meaning:** total number of pixels within the central ellipse

## Ellipse-guess masking

### `AREA_GUESS_ELLIPSE_PIX`
- **Meaning:** geometric area of the initial guessed ellipse
- **Units:** pixels

### `AREA_GUESS_ELLIPSE_UNMASKED_PIX`
- **Meaning:** unmasked area of the guessed ellipse
- **Units:** pixels

### `MASKFRAC_GUESS_ELLIPSE`
- **Meaning:** masked fraction within the guessed ellipse

---

# 5. Cutout geometry and adopted ellipse quantities

## Cutout size

### `CUTOUT_SCALE`
- **Meaning:** cutout pixel scale
- **Units:** arcsec/pixel

### `CUTOUT_XSIZE`, `CUTOUT_YSIZE`
- **Meaning:** cutout dimensions
- **Units:** pixels

## Filter scaling

### `FILTER_CORRECTION`
- **Meaning:** correction factor associated with Hα filter throughput / line placement
- **Notes:** large values may indicate less reliable Hα flux interpretation

### `FILTER_RATIO`
- **Meaning:** ratio used in filter correction / scaling
- **Notes:** exact interpretation depends on implementation details

## Adopted ellipse (`ELL0_*`)

### `ELL0_SMA_ARCSEC`
- **Meaning:** adopted semi-major axis of the reference ellipse
- **Units:** arcseconds

### `ELL0_BA`
- **Meaning:** adopted axis ratio \(b/a\)

### `ELL0_PA_DEG`
- **Meaning:** adopted position angle
- **Units:** degrees

### `ELL0_XC`, `ELL0_YC`
- **Meaning:** adopted ellipse center in image coordinates
- **Units:** pixels

### `ELL0_SOURCE`
- **Meaning:** source of the adopted ellipse parameters
- **Examples:** metadata, upstream catalog, pipeline-derived estimate

## Measured ellipse quantities

### `ELLIP_XCENTROID`, `ELLIP_YCENTROID`
- **Meaning:** centroid measured from ellipse/segmentation workflow
- **Units:** pixels

### `ELLIP_SMA_PIX`
- **Meaning:** measured semi-major axis
- **Units:** pixels

### `ELLIP_B_PIX`
- **Meaning:** measured semi-minor axis
- **Units:** pixels

### `ELLIP_EPS`
- **Meaning:** ellipticity, usually \(1 - b/a\)

### `ELLIP_THETA_RAD`
- **Meaning:** position angle / orientation angle
- **Units:** radians

### `ELLIP_BA`
- **Meaning:** measured axis ratio from ellipse solution

### `ELLIP_SOURCE_SUM`
- **Meaning:** summed source flux over the ellipse/segment used for the ellipse calculation

### `ELLIP_SEGMENT_FLUX`
- **Meaning:** integrated flux in the segmentation region used by the ellipse workflow

### `ELLIP_SEGMENT_MAG`
- **Meaning:** magnitude corresponding to the ellipse segmentation flux

## Ellipse-comparison diagnostics

### `ELL_DC_PX`
- **Meaning:** center offset between two ellipse solutions
- **Units:** pixels

### `ELL_DBA`
- **Meaning:** difference in axis ratio between two ellipse solutions

### `ELL_DPA_DEG`
- **Meaning:** difference in position angle between two ellipse solutions
- **Units:** degrees

### `ELL_SMA_RATIO`
- **Meaning:** ratio of two semi-major-axis estimates

---

# 6. Image quality, sky, and scaling quantities

## Seeing / PSF diagnostics

### `R_FWHM_PSF`, `H_FWHM_PSF`
- **Meaning:** PSF full width at half maximum in r-band / Hα
- **Units:** likely arcseconds

### `R_FWHM_SE`, `H_FWHM_SE`
- **Meaning:** Source Extractor or alternative seeing estimate
- **Units:** likely arcseconds

## Sky measurements

### `R_SKYSTD_ADU`, `H_SKYSTD_ADU`
- **Meaning:** sky standard deviation in image units

### `R_SKYMED_ADU`, `H_SKYMED_ADU`
- **Meaning:** sky median in image units

### `R_SKYSTD_PHYS`, `H_SKYSTD_PHYS`
- **Meaning:** sky noise expressed in physical flux units

## Flux scaling

### `R_SCALE_ADU_CGS`, `H_SCALE_ADU_CGS`
- **Meaning:** conversion factor from image units to physical CGS flux units

---

# 7. Custom HAPY morphology

These are the custom HAPY morphology measurements and associated diagnostics.

## Detection/segmentation support

### `R_HAPY_NPIX`, `H_HAPY_NPIX`
- **Meaning:** number of pixels used in the custom morphology measurement

### `H_HAPY_FILLFRAC`
- **Meaning:** filling fraction of the Hα detection mask/region

### `R_HAPY_SNP_ALL`
- **Meaning:** signal-to-noise-like quantity for all pixels used in r-band morphology

### `H_HAPY_SNP_ALL`
- **Meaning:** signal-to-noise-like quantity for all Hα morphology pixels

### `H_HAPY_SNP_DET`
- **Meaning:** signal-to-noise-like quantity for detected Hα pixels

### `H_GINI_THRESHOLD`
- **Meaning:** Hα threshold used in custom Gini workflow

## Morphology values

### `R_HAPY_GINI`, `H_HAPY_GINI`
- **Meaning:** custom HAPY Gini coefficient in r-band / Hα

### `R_HAPY_M20`, `H_HAPY_M20`
- **Meaning:** custom HAPY M20 in r-band / Hα

## Flags

### `HAPY_MORPH_FLAG`
- **Meaning:** integer diagnostic flag or bitmask for custom morphology workflow
- **Notes:** decode using helper logic if available

## Legacy/simple morphology quantities

### `R_M20`, `H_M20`
- **Meaning:** M20 values from another morphology path or legacy implementation

### `R_ASYM`, `H_ASYM`
- **Meaning:** asymmetry values

### `R_ASYM_ERR`, `H_ASYM_ERR`
- **Meaning:** uncertainties on asymmetry

### `R_ELLIP_GINI`, `H_ELLIP_GINI`
- **Meaning:** Gini coefficient measured within ellipse-defined regions

---

# 8. Statmorph outputs

These columns summarize measurements from the `statmorph` package.

## r-band statmorph

### `R_SM_FLAG`
- **Meaning:** integer statmorph diagnostic flag for r-band measurement

### `R_SM_XCENTROID`, `R_SM_YCENTROID`
- **Meaning:** centroid used/measured by statmorph
- **Units:** pixels

### `R_SM_GINI`, `R_SM_M20`
- **Meaning:** statmorph Gini and M20 in r-band

### `R_SM_C`, `R_SM_A`, `R_SM_S`
- **Meaning:** concentration, asymmetry, smoothness in r-band

### `R_SM_RPETRO_ELLIP`
- **Meaning:** elliptical Petrosian radius from statmorph

### `R_SM_RHALF_ELLIP`
- **Meaning:** elliptical half-light radius from statmorph

### `R_SM_R20`, `R_SM_R50`, `R_SM_R80`
- **Meaning:** radii enclosing 20%, 50%, and 80% of the light

### `R_SM_RMAX_CIRCLE`, `R_SM_RMAX_ELLIP`
- **Meaning:** maximum radius in circular / elliptical definitions

### `R_SM_SERSIC_*`
- **Meaning:** Sérsic-model parameters returned by statmorph
- **Key columns:**
  - `R_SM_SERSIC_AMP`
  - `R_SM_SERSIC_RHALF`
  - `R_SM_SERSIC_N`
  - `R_SM_SERSIC_XC`
  - `R_SM_SERSIC_YC`
  - `R_SM_SERSIC_ELLIP`
  - `R_SM_SERSIC_THETA`
  - `R_SM_SERSIC_CHISQ_DOF`
  - `R_SM_SERSIC_FLAG`

### `R_SM_SN_PER_PIXEL`
- **Meaning:** signal-to-noise per pixel used in r-band statmorph

### `R_SM_SKY_MEAN`, `R_SM_SKY_MEDIAN`, `R_SM_SKY_SIGMA`
- **Meaning:** local sky statistics used by statmorph

## Hα statmorph

The Hα statmorph columns mirror the r-band set.

### `H_SM_FLAG`
- **Meaning:** integer statmorph diagnostic flag for Hα

### `H_SM_XCENTROID`, `H_SM_YCENTROID`
- **Meaning:** Hα centroid from statmorph

### `H_SM_GINI`, `H_SM_M20`
- **Meaning:** Hα statmorph Gini and M20

### `H_SM_C`, `H_SM_A`, `H_SM_S`
- **Meaning:** concentration, asymmetry, smoothness in Hα

### `H_SM_RPETRO_ELLIP`, `H_SM_RHALF_ELLIP`
- **Meaning:** Petrosian and half-light radii from Hα statmorph

### `H_SM_R20`, `H_SM_R50`, `H_SM_R80`
- **Meaning:** enclosed-light radii in Hα

### `H_SM_RMAX_CIRCLE`, `H_SM_RMAX_ELLIP`
- **Meaning:** maximum circular / elliptical Hα radii in statmorph

### `H_SM_SERSIC_*`
- **Meaning:** Sérsic parameters returned by Hα statmorph fit

### `H_SM_SN_PER_PIXEL`
- **Meaning:** Hα signal-to-noise per pixel used in statmorph

### `H_SM_SKY_MEAN`, `H_SM_SKY_MEDIAN`, `H_SM_SKY_SIGMA`
- **Meaning:** Hα local sky statistics used in statmorph

### `STATMORPH_SEC`
- **Meaning:** elapsed runtime for statmorph processing
- **Units:** seconds

---

# 9. r-band radial-profile and photometric quantities

## Basic profile sampling

### `R_PROFILE_NGOOD`
- **Meaning:** number of good radial bins used in the r-band profile

### `R_PROFILE_MASKFRAC_MAX`
- **Meaning:** maximum masked fraction in any r-band radial bin

## Enclosed-light radii

### `R25_ARCSEC`, `R50_ARCSEC`, `R75_ARCSEC`
- **Meaning:** radii enclosing 25%, 50%, and 75% of the r-band light

### `R25_PIX`, `R50_PIX`, `R75_PIX`
- **Meaning:** same radii in pixels

## Isophotal quantities

### `R24_ARCSEC`
- **Meaning:** radius of the 24 mag arcsec\(^{-2}\) isophote in the adopted magnitude system

### `R24_ARCSEC_ERR`
- **Meaning:** uncertainty on `R24_ARCSEC`

### `R24_MAG`
- **Meaning:** integrated magnitude within `R24_ARCSEC`

### `R24_MAG_ERR`
- **Meaning:** uncertainty on `R24_MAG`

### `R25_ISO_ARCSEC`, `R25_ISO_MAG`
- **Meaning:** radius and enclosed magnitude at the 25 mag arcsec\(^{-2}\) isophote

### `R25P5_ARCSEC`, `R25P5_MAG`
- **Meaning:** radius and enclosed magnitude at the 25.5 mag arcsec\(^{-2}\) isophote

### `R24_VEGA_*`, `R25_VEGA_*`
- **Meaning:** Vega-system versions of the corresponding isophotal quantities

## Fluxes and concentration

### `R30R24_FLUX_CGS`
- **Meaning:** flux within 0.3 × `R24`
- **Notes:** naming reflects the traditional concentration-style aperture choice

### `R24_FLUX_CGS`
- **Meaning:** total r-band flux within `R24` aperture

### `R_C30`
- **Meaning:** concentration-like quantity formed from flux ratio within inner vs `R24` aperture
- **Notes:** useful as a simple central-concentration measure

## Petrosian quantities

### `R_PETRO_RAD_ARCSEC`
- **Meaning:** Petrosian radius in r-band

### `R_PETRO_FLUX`, `R_PETRO_FLUX_CGS`
- **Meaning:** Petrosian flux in instrumental / physical units

### `R_PETRO_MAG`
- **Meaning:** Petrosian magnitude

### `R_PETRO_R50_ARCSEC`, `R_PETRO_R90_ARCSEC`
- **Meaning:** radii containing 50% and 90% of Petrosian flux

### `R_PETRO_CON`
- **Meaning:** Petrosian concentration measure

## Profile-fit quantities

### `R_EXPFIT_I0`
- **Meaning:** exponential-fit central intensity or normalization

### `R_EXPFIT_K`
- **Meaning:** exponential profile slope/scale parameter

### `R_EXPFIT_RE_ARCSEC`
- **Meaning:** effective radius implied by the exponential fit

### `R_LOGFIT_SLOPE`, `R_LOGFIT_INTERCEPT`
- **Meaning:** coefficients of a log-linear fit to the profile

### `R_LOGFIT_RE_ARCSEC`
- **Meaning:** effective radius implied by the log-fit model

## Total-S/N style quantities

### `R_TOT_MAG_SNR`
- **Meaning:** signal-to-noise proxy for total magnitude measurement

### `R_TOT_FLUX_CGS_SNR`
- **Meaning:** signal-to-noise proxy for total flux in physical units

### `R_TOT_FLUX_CGS_SNR_ERR`
- **Meaning:** uncertainty on the total flux SNR-like estimate

### `R_SNR_TRUNC_ARCSEC`
- **Meaning:** radius where the r-band profile becomes S/N-limited or truncated

---

# 10. Hα radial-profile, extent, and flux quantities

## Profile sampling

### `H_PROFILE_NGOOD`
- **Meaning:** number of good radial bins used in the Hα profile

### `H_PROFILE_LONGRUN`
- **Meaning:** length/extent of contiguous Hα detections in profile logic
- **Notes:** exact implementation meaning should follow code definition

### `H_NDET_RUNS`
- **Meaning:** number of detected radial runs or segments in Hα profile logic

### `H_PROFILE_MASKFRAC_MAX`
- **Meaning:** maximum masked fraction in any Hα profile bin

## Basic extent quantities

### `H_MAXDET_ARCSEC`
- **Meaning:** maximum radius with detected Hα emission

### `H_MAXDET_PIX`
- **Meaning:** same quantity in pixels

### `H25_ARCSEC`, `H50_ARCSEC`, `H75_ARCSEC`
- **Meaning:** radii enclosing 25%, 50%, and 75% of total Hα flux

### `H25_PIX`, `H50_PIX`, `H75_PIX`
- **Meaning:** same quantities in pixels

### `H_SNR_TRUNC_ARCSEC`
- **Meaning:** radius where the Hα profile becomes S/N-limited or truncated

## Integrated Hα fluxes

### `H_TOT_FLUX_CGS`
- **Meaning:** total integrated Hα flux

### `H_TOT_FLUX_CGS_ERR`
- **Meaning:** uncertainty on total Hα flux

### `H_R24_FLUX_CGS`
- **Meaning:** Hα flux enclosed within the r-band `R24` aperture

### `H_R24_FLUX_CGS_ERR`
- **Meaning:** uncertainty on `H_R24_FLUX_CGS`

### `H30R24_FLUX_CGS`
- **Meaning:** Hα flux within inner 0.3 × `R24` aperture

### `H30R24_FLUX_CGS_ERR`
- **Meaning:** uncertainty on `H30R24_FLUX_CGS`

### `H_C30_R24`
- **Meaning:** Hα concentration measure formed from inner-to-`R24` flux ratio

### `H_C30_R24_ERR`
- **Meaning:** uncertainty on `H_C30_R24`

## Hα isophotal quantities

### `H_ISO5E17_ARCSEC`
- **Meaning:** Hα isophotal radius at threshold \(5 \times 10^{-17}\) in pipeline physical units

### `H_ISO5E17_FLUX_CGS`
- **Meaning:** Hα flux within that isophotal threshold

### `H_ISO17E18_ARCSEC`
- **Meaning:** Hα isophotal radius at threshold \(1.7 \times 10^{-18}\) in pipeline physical units

### `H_ISO17E18_FLUX_CGS`
- **Meaning:** Hα flux within that isophotal threshold

### Corresponding `_ERR` columns
- **Meaning:** uncertainties on the isophotal radius and flux quantities

## Additional extent diagnostics

### `H_R95_R24_ARCSEC`
- **Meaning:** radius containing 95% of the Hα flux measured within the `R24` framework
- **Notes:** useful as an alternative outer-extent measure

## Hα Petrosian quantities

### `H_PETRO_RAD_ARCSEC`
- **Meaning:** Petrosian radius in Hα

### `H_PETRO_FLUX`, `H_PETRO_FLUX_CGS`
- **Meaning:** Hα Petrosian flux in instrumental / physical units

### `H_PETRO_MAG`
- **Meaning:** Hα Petrosian magnitude-like quantity

### `H_PETRO_R50_ARCSEC`, `H_PETRO_R90_ARCSEC`
- **Meaning:** radii containing 50% and 90% of Hα Petrosian flux

### `H_PETRO_CON`
- **Meaning:** Hα Petrosian concentration

## Hα profile fits

### `H_EXPFIT_I0`, `H_EXPFIT_K`, `H_EXPFIT_RE_ARCSEC`
- **Meaning:** exponential-fit parameters for the Hα radial profile

### `H_LOGFIT_SLOPE`, `H_LOGFIT_INTERCEPT`, `H_LOGFIT_RE_ARCSEC`
- **Meaning:** log-fit parameters and implied effective radius for Hα profile

---

# 11. GALFIT outputs

The table includes two GALFIT solution families:
- an unconstrained/non-convolved solution (`GAL_*`)
- a constrained or convolved solution (`GAL_C*`)

## Unconstrained/non-convolved fit

### `GAL_XC`, `GAL_YC`
- **Meaning:** best-fit center from GALFIT
- **Units:** pixels

### `GAL_XC_ERR`, `GAL_YC_ERR`
- **Meaning:** formal uncertainties on GALFIT center

### `GAL_MAG`
- **Meaning:** best-fit total magnitude from GALFIT

### `GAL_MAG_ERR`
- **Meaning:** formal uncertainty on `GAL_MAG`

### `GAL_RE`
- **Meaning:** effective radius from GALFIT
- **Units:** likely pixels unless explicitly converted elsewhere
- **Notes:** check implementation before mixing directly with `_ARCSEC` quantities

### `GAL_RE_ERR`
- **Meaning:** formal uncertainty on `GAL_RE`

### `GAL_N`
- **Meaning:** Sérsic index

### `GAL_N_ERR`
- **Meaning:** formal uncertainty on `GAL_N`

### `GAL_BA`
- **Meaning:** axis ratio \(b/a\)

### `GAL_BA_ERR`
- **Meaning:** formal uncertainty on `GAL_BA`

### `GAL_PA`
- **Meaning:** position angle
- **Units:** degrees

### `GAL_PA_ERR`
- **Meaning:** formal uncertainty on `GAL_PA`

### `GAL_SKY`
- **Meaning:** fitted sky level in GALFIT

### `GAL_SKY_ERR`
- **Meaning:** formal uncertainty on fitted sky

### `GAL_CHISQ`
- **Meaning:** fit goodness statistic

### `GAL_NUMERR`
- **Meaning:** GALFIT numerical error code/count

### `GAL_ERROR`
- **Meaning:** GALFIT summary error quantity

## Constrained/convolved fit

The constrained/convolved fit columns mirror the unconstrained set.

### `GAL_CXC`, `GAL_CYC`
- constrained-fit center

### `GAL_CMAG`
- constrained-fit magnitude

### `GAL_CRE`
- constrained-fit effective radius

### `GAL_CN`
- constrained-fit Sérsic index

### `GAL_CBA`
- constrained-fit axis ratio

### `GAL_CPA`
- constrained-fit position angle

### `GAL_CSKY`
- constrained-fit sky level

### `GAL_CCHISQ`
- constrained-fit goodness statistic

### `GAL_CNUMERR`
- constrained-fit numerical error code/count

### `GAL_CERROR`
- constrained-fit summary error quantity

### Corresponding `_ERR` columns
- **Meaning:** formal fit uncertainties for the constrained solution

---

# 12. Derived QC and science columns

These columns may be added downstream by QC or science scripts rather than by `run_analysis.py` itself.

## Example QC columns

### `WARN_MASK`
- **Meaning:** masking is large enough to potentially bias measurements

### `WARN_BRIGHT_STAR`
- **Meaning:** bright-star contamination may affect results

### `WARN_FILTER`
- **Meaning:** Hα filter correction is large enough to be concerning

### `WARN_ELLIPSE`
- **Meaning:** ellipse mismatch or geometry inconsistency warning

### `WARN_WEAK_HA`
- **Meaning:** Hα morphology/extent may be unreliable because the detection is weak

## Example usability columns

### `USE_R_STRUCTURE`
- **Meaning:** usability flag for r-band structural quantities
- **Values:** `0 = bad`, `1 = caution`, `2 = good`

### `USE_HA_EXTENT`
- **Meaning:** usability flag for Hα extent quantities

### `USE_HA_MORPH`
- **Meaning:** usability flag for Hα morphology quantities

### `USE_GALFIT`
- **Meaning:** usability flag for GALFIT-based structure

## Example summary QC tier

### `QC_TIER`
- **Meaning:** overall per-row QC tier summarizing scientific reliability
- **Typical interpretation:**
  - `A` : high-quality, science-ready
  - `B` : usable with minor limitations
  - `C` : warning-limited
  - `D` : weak / incomplete
  - `F` : failed

## Example derived science quantities

### `H50_R50_RATIO`
- **Meaning:** ratio of Hα to stellar half-light radius

### `H75_R75_RATIO`
- **Meaning:** ratio of Hα to stellar 75%-light radius

### `H_MAXDET_R25_RATIO`
- **Meaning:** ratio of maximum detected Hα extent to stellar `R25`

### `H_PETRO_R50_RATIO`
- **Meaning:** ratio of Hα to stellar Petrosian half-light radius

### `DELTA_GINI`
- **Meaning:** `H_HAPY_GINI - R_HAPY_GINI`
- **Interpretation:** positive values indicate Hα is more unequally distributed / clumpier or more centrally concentrated than the stellar light

### `DELTA_M20`
- **Meaning:** `H_HAPY_M20 - R_HAPY_M20`
- **Interpretation:** useful for comparing the spatial concentration of bright Hα-emitting regions relative to the stellar distribution

---

# 13. Suggested default science quantities

For first-pass science analysis, the following columns are often the most useful:

## Stellar structure
- `R50_ARCSEC`
- `R25_ARCSEC`
- `R24_ARCSEC`
- `R_C30`
- `R_PETRO_R50_ARCSEC`

## Hα extent and concentration
- `H50_ARCSEC`
- `H_MAXDET_ARCSEC`
- `H_R24_FLUX_CGS`
- `H_TOT_FLUX_CGS`
- `H_C30_R24`
- `H_R95_R24_ARCSEC`

## Morphology
- `R_HAPY_GINI`, `H_HAPY_GINI`
- `R_HAPY_M20`, `H_HAPY_M20`
- `DELTA_GINI`, `DELTA_M20`

## QC support
- `QC_TIER`
- `USE_HA_EXTENT`
- `USE_HA_MORPH`
- `WARN_WEAK_HA`
- `ELL0_MASK_WARN`
- `BRIGHT_STAR_FLAG`

---

# 14. Future improvements to this document

This file is intended to evolve alongside the pipeline. Useful future additions include:

- explicit unit confirmation for every column
- bitmask decoding tables for `_FLAG` columns
- examples of recommended science cuts
- cross-links to the code that produces each measurement
- a smaller “minimal science table” description for downstream papers
