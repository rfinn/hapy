# HAPY Hα Clump Analysis

HAPY includes an optional Hα clump-analysis module designed to quantify the internal structure of star formation within the stellar extent of each galaxy.

The clump analysis is implemented in:

```text
hapy/ellipse/clumps.py
```

and is run from `run_analysis.py` with:

```bash
run_analysis --cutout-dir <cutout_dir> --clumps
```

## Overview

The goal of the clump analysis is to identify and measure discrete Hα star-forming regions inside the R-band footprint of the central galaxy.

The current workflow is:

1. Use the existing R-band footprint/segmentation from `EllipsePhotometry`.
2. Restrict the continuum-subtracted Hα image to that R-band footprint.
3. Estimate a local Hα background and detection threshold using `hapy.imagetools.imutils.calculate_background_photutils`.
4. Run `photutils.segmentation.detect_sources` on the restricted Hα image.
5. Optionally run `photutils.segmentation.deblend_sources` to separate connected emission into individual clumps.
6. Measure the clumps with `photutils.segmentation.SourceCatalog`.
7. Save the clump catalog, segmentation map, local peak catalog, and diagnostic plot.
8. Store summary quantities in the main HAPY results table.

The primary clump definition is segmentation-based. Local peaks are also measured as a secondary diagnostic, but they are not treated as clumps unless they are part of the segmented clump catalog.

## Default input image

The first implementation runs on the standard CS-ZP continuum-subtracted Hα image.

In the output table, this is recorded as:

```text
HCL_INPUT_IMAGE = "CS-ZP"
```

Future runs on the CS-gr Hα image should use a different prefix, for example:

```text
CSGR_HCL
```

so that CS-ZP and CS-gr clump measurements can coexist in the same results table.

## Default parameters

The current recommended default configuration is:

```python
ClumpDetectionConfig(
    nsigma=3.0,
    npixels=5,
    deblend=True,
    nlevels=64,
    contrast=0.001,
    mode="linear",
    find_peaks=True,
)
```

The detection threshold is:

```text
threshold = background_median + nsigma * background_rms
```

where the background quantities are estimated using `calculate_background_photutils`.

## Saved output products

For a cutout tag `<tag>`, clump analysis saves the following products in the main cutout directory:

```text
<tag>-halpha-clumps.ecsv
<tag>-halpha-clumps-segm.fits
<tag>-halpha-clump-peaks.ecsv
<tag>-halpha-clumps-diagnostic.png
```

Optional point-source-like detections may also be saved if that mode is enabled:

```text
<tag>-halpha-clump-pointsrc.ecsv
```

## Main results-table columns

The examples below use the default prefix `HCL_`. If the clump analysis is run on CS-gr images, the same columns can be written with a different prefix such as `CSGR_HCL_`.

All clump offset and area summary quantities in the main HAPY results table are stored in angular units. The pixel-to-arcsec conversion uses the existing image pixel-scale column already stored in the main results table.

### Status and provenance

| Column            | Meaning                                         |
| ----------------- | ----------------------------------------------- |
| `HCL_OK`          | `True` if clump analysis completed successfully |
| `HCL_STATUS`      | "not_run", "ok", or "failed"                 |
| `HCL_INPUT_IMAGE` | Input Hα image used, e.g. `CS-ZP` or `CS-gr`    |

### Detection counts

| Column                | Meaning                                                                                   |
| --------------------- | ----------------------------------------------------------------------------------------- |
| `HCL_NCLUMP`          | Number of Hα clumps from `detect_sources` plus `deblend_sources`                          |
| `HCL_NCLUMP_GOOD`     | Number of clumps passing optional `min_flux` and `min_area` cuts                          |
| `HCL_NPEAK`           | Number of local Hα peaks found by `find_peaks`                                            |
| `HCL_NPEAK_IN_CLUMPS` | Number of local peaks that fall inside segmented clump regions                            |
| `HCL_NPOINTSRC`       | Number of optional point-source-like detections; usually zero unless this mode is enabled |

### Threshold and background

| Column           | Meaning                                                  |
| ---------------- | -------------------------------------------------------- |
| `HCL_THRESHOLD`  | Final Hα clump detection threshold                       |
| `HCL_BKG_MEAN`   | Background mean from `calculate_background_photutils`    |
| `HCL_BKG_MEDIAN` | Background median from `calculate_background_photutils`  |
| `HCL_BKG_RMS`    | Background RMS/std from `calculate_background_photutils` |

### Flux and area summaries

| Column                   | Meaning                                                       |
| ------------------------ | ------------------------------------------------------------- |
| `HCL_FLUX_FOOTPRINT`     | Total Hα flux inside the central R-band footprint             |
| `HCL_FLUX_SUM`           | Sum of fluxes in detected good Hα clumps                      |
| `HCL_FLUX_FRAC`          | `HCL_FLUX_SUM / HCL_FLUX_FOOTPRINT`                           |
| `HCL_BRIGHT_FLUX`        | Flux of the brightest Hα clump                                |
| `HCL_BRIGHT_FRAC`        | Brightest clump flux divided by total clump flux              |
| `HCL_TOP2_FRAC`          | Sum of the two brightest clumps divided by total clump flux   |
| `HCL_TOP3_FRAC`          | Sum of the three brightest clumps divided by total clump flux |
| `HCL_AREA_ARCSEC2`           | Number of pixels  (arcsec2) assigned to detected good clumps             |
| `HCL_FOOTPRINT_AREA_ARECSEC2` | Number of valid pixels (arcsec2) in the central R-band footprint        |
| `HCL_AREA_FRAC`          | `HCL_AREA_PIX / HCL_FOOTPRINT_AREA_PIX`                       |

### Peak and optional point-source summaries

| Column                  | Meaning                                                                                                      |
| ----------------------- | ------------------------------------------------------------------------------------------------------------ |
| `HCL_PEAK_MAX`          | Maximum local peak value from `find_peaks`                                                                   |
| `HCL_PEAK_SUM`          | Sum of local peak values from `find_peaks`                                                                   |
| `HCL_POINTSRC_FLUX_SUM` | Sum of flux proxy from optional point-source finder; blank or NaN when point-source detection is not enabled |

### Clump centroid summaries

| Column            | Meaning                                      |
| ----------------- | -------------------------------------------- |
| `HCL_XCEN_FLUXWT` | Flux-weighted mean x centroid of good clumps |
| `HCL_YCEN_FLUXWT` | Flux-weighted mean y centroid of good clumps |
| `HCL_XCEN_BRIGHT` | x centroid of the brightest clump            |
| `HCL_YCEN_BRIGHT` | y centroid of the brightest clump            |

### Nuclear Clumps
A nuclear clump is defined as a good Hα clump whose centroid lies within 1.5 times the measured Hα FWHM of the adopted galaxy center. The factor is recorded in `HCL_PARAM_NUCLEAR_RADIUS_FWHM`.


| Column | Meaning |
|---|---|
| `HCL_HAS_NUCLEAR` | `True` if at least one good clump lies within one Hα seeing FWHM of the adopted galaxy center |
| `HCL_NNUCLEAR` | Number of good clumps whose centroids lie within one Hα seeing FWHM |
| `HCL_NUCLEAR_FLUX_FRAC` | Fraction of total good-clump flux contained in nuclear clumps |
| `HCL_NUCLEAR_RADIUS_ARCSEC` | Nuclear radius used for the classification, equal to the measured Hα seeing FWHM |

This definition is intentionally tied to the image resolution rather than a fixed fraction of galaxy size, so the nuclear flag identifies clumps that are unresolved or marginally resolved from the adopted center.
### Clump centroid offsets from the galaxy center

These quantities use the clump centroids relative to the adopted galaxy center.

| Column              | Meaning                                                |
| ------------------- | ------------------------------------------------------ |
| `HCL_RMIN_ARCSEC`      | Minimum clump-centroid offset from the galaxy center   |
| `HCL_RMAX_ARCSEC`      | Maximum clump-centroid offset from the galaxy center   |
| `HCL_RMED_ARCSEC`      | Median clump-centroid offset                           |
| `HCL_RMEAN_ARCSEC`     | Mean clump-centroid offset                             |
| `HCL_RFLUXWT_ARCSEC`   | Flux-weighted mean clump-centroid offset               |
| `HCL_RBRIGHT_ARCSEC`   | Offset of the brightest clump from the galaxy center   |
| `HCL_DX_BRIGHT_ARCSEC` | x offset of the brightest clump from the galaxy center |
| `HCL_DY_BRIGHT_ARCSEC` | y offset of the brightest clump from the galaxy center |

### Saved product paths

| Column         | Meaning                                        |
| -------------- | ---------------------------------------------- |
| `HCL_CATALOG`  | Path to saved clump catalog ECSV               |
| `HCL_SEGMAP`   | Path to saved Hα clump segmentation FITS image |
| `HCL_PEAKS`    | Path to saved local peak table ECSV            |
| `HCL_POINTSRC` | Path to optional point-source table ECSV       |
| `HCL_DIAG`     | Path to clump diagnostic plot                  |

### Important configuration parameters

| Column                    | Meaning                                                 |
| ------------------------- | ------------------------------------------------------- |
| `HCL_PARAM_NSIGMA`        | Detection threshold in units of local background RMS    |
| `HCL_PARAM_NARCSECELS`       | Minimum connected pixels required for a clump detection |
| `HCL_PARAM_DEBLEND`       | Whether deblending was enabled                          |
| `HCL_PARAM_NLEVELS`       | Number of deblending levels                             |
| `HCL_PARAM_CONTRAST`      | Deblend contrast                                        |
| `HCL_PARAM_MODE`          | Deblend mode, e.g. `linear`, `exponential`, or `sinh`   |
| `HCL_PARAM_FIND_PEAKS`    | Whether local peak finding was enabled                  |
| `HCL_PARAM_PEAK_BOX_SIZE` | Box size used for `find_peaks`                          |
| `HCL_PARAM_PEAK_MIN_SEP`  | Minimum peak separation if used; `-1` means not set     |

## Interpretation notes

`HCL_NCLUMP` is the primary clump count. It measures robust segmented Hα clumps inside the R-band footprint.

`HCL_NPEAK` is a secondary diagnostic. It measures local maxima and can be sensitive to noise, especially in low surface-brightness regions.

A useful distinction is:

```text
HCL_NCLUMP = number of segmented Hα clumps
HCL_NPEAK = number of local Hα peaks
HCL_NPEAK_IN_CLUMPS = number of peaks inside segmented clumps
```

For example, a galaxy with one large Hα complex but many local maxima may have:

```text
HCL_NCLUMP = 1
HCL_NPEAK_IN_CLUMPS > 1
```

This indicates a multi-peaked star-forming complex rather than several fully separated clumps.

`HCL_BRIGHT_FRAC` and `HCL_TOP2_FRAC` are useful for distinguishing galaxies dominated by one central Hα region from galaxies dominated by two or more major clumps.

`HCL_RBRIGHT_ARCSEC`, `HCL_RFLUXWT_ARCSEC`, and `HCL_RMAX_ARCSEC` summarize how centrally concentrated or off-center the detected clumps are relative to the adopted galaxy center.

## Known caveats

The clump catalog is sensitive to detection threshold, minimum connected-pixel area, and deblending mode.

Faint visual features may not be segmented as clumps if they do not have enough connected pixels above the detection threshold.

The peak catalog can identify local maxima that are not robust enough to satisfy the segmentation-based clump definition.

For the first science runs, the clump analysis should be treated as a quantitative morphology diagnostic rather than a complete HII-region catalog.
