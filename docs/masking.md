# Masking Workflow

The `hapy.masktools` package provides tools for constructing and editing segmentation masks for optical imaging data. Masks are typically used to:

* Remove foreground stars
* Remove background galaxies
* Exclude artifacts
* Isolate a target galaxy for photometric analysis

The masking system is structured to clearly separate:

* **MaskEngine** — core masking logic and state management
* **maskops** — pixel-level geometry and mask utilities
* **maskgui** — interactive Qt-based GUI for visualization and editing

---

# Overview of the Masking Pipeline

When building an initial mask, the following steps are performed:

1. **Run Source Extractor**
   Sources are detected and a segmentation map is generated.

2. **Remove central object (optional)**
   If a target galaxy is specified, its segmentation region can be removed or preserved.

3. **Grow masked regions**
   Masked areas can be expanded to ensure full coverage of detected objects.

4. **Add Gaia star masks (optional)**
   Bright stars may be masked using Gaia catalog information.

5. **Write mask to disk**
   The output mask is written as:

   ```
   <image>-mask.fits
   ```

The resulting mask is a 2D FITS image with integer values:

* `0` → unmasked pixel
* `> 0` → masked object ID

Each object in the segmentation map receives a unique integer identifier.

---

# Running the Mask GUI

To launch the interactive GUI:

```
run_maskgui --image VFID2507-UGC08693-HDI-20170523-p007-R.fits --haimage VFID2507-UGC08693-HDI-20170523-p007-CS.fits
```


When launched, the GUI:

* Builds the initial mask using `MaskEngine`
* Displays three synchronized panels:

  * r-band image
  * Hα image (optional)
  * Mask
* Allows interactive editing of the mask


Here’s a concise description you can drop into `docs/masking.md`:

---

## Running `maskgui` with Metadata

You can launch the interactive masking GUI using the ellipse geometry and target information stored in a cutout’s `metadata.json`.

Use the helper script:

```bash
run_maskgui_from_metadata --cutout-dir cutouts/<tag>
```

This will:

* read `metadata.json` in the cutout directory
* load the r-band image (and Hα/CS image if present)
* pass the target galaxy parameters to `run_maskgui`:

  * RA / Dec
  * semi-major axis (`sma_arcsec`)
  * axis ratio (`b/a`)
  * position angle (`pa_deg`)

This ensures the GUI is initialized with the same geometry used by `run_analysis`.

### Optional arguments

```bash
--gaia-dir <path>    # use precomputed Gaia catalogs
--no-gaia            # disable Gaia star masking
```

Any additional arguments are passed directly to `run_maskgui`.

### Output

Mask files are written automatically to the cutout directory as:

```
<tag>-mask-manual.fits
<tag>-inv-mask-manual.fits
```

No need to change directories before running — the full image path is used to determine the output location.

---

# Interactive Editing

Click on any image panel, then use the following keyboard shortcuts:

| Key | Action                      |
| --- | --------------------------- |
| `r` | Remove object under cursor  |
| `c` | Add circular mask at cursor |
| `b` | Add square mask at cursor   |
| `g` | Grow masked regions         |
| `v` | print pixel values at cursor position|
| `m` | toggle mask view on current image (outline, filled, none) |
| `w` | Write mask to FITS file     |
| `h` | Print help menu             |
| `q` | Quit the GUI                |

Edits are applied immediately in memory.
Press `w` to write the updated mask to disk.

---

# Programmatic Usage (No GUI)

The masking engine can be used independently of the GUI:

```python
from hapy.masktools.api import MaskEngine

engine = MaskEngine(
    image_fits="image.fits",
    sepath="sex",
    config="default.sex.HDI.mask",
    threshold=0.005,
    snr=10,
    minarea=5,
)

mask = engine.build_initial_mask()
engine.write_mask("image-mask.fits")
```

This allows batch processing or integration into larger analysis pipelines.

---

# Core Components

## MaskEngine

`MaskEngine` is responsible for:

* Running Source Extractor
* Managing mask state (`maskdat`)
* Adding and removing objects
* Growing masked regions
* Writing FITS output

The engine owns all mask data and should be treated as the authoritative source of mask state.

---

## maskops

The `maskops` module contains pure functions for pixel-level mask operations, such as:

* Generating circular masks
* Applying box masks
* Growing masks
* Geometry utilities

These functions do not depend on Qt or GUI components.

---

## maskgui

The GUI provides:

* Image display using Ginga
* Interactive editing
* Keyboard-driven mask operations
* Visualization of the central galaxy ellipse (if provided)

The GUI does not implement masking logic directly. All mask modifications are delegated to `MaskEngine`.

---

# Architecture Philosophy

The masking system follows a simple Model–View–Controller pattern:

```
MaskWindow (GUI)
    ↓
MaskEngine (mask logic + state)
    ↓
maskops (pure geometry functions)
```

This separation allows:

* Batch masking without the GUI
* Interactive editing
* Easier testing of core logic
* Cleaner long-term maintenance

---

# Dependencies

The masking tools require:

* PyQt5
* Ginga
* Astropy
* NumPy
* Source Extractor (installed separately)

---




# Mask selection logic (run_analysis)

The pipeline selects a mask using the following priority order:

## 1. Command-line override (highest priority)

```text
--mask_fits
```

* If provided, this mask is always used.
* `MASK_SOURCE = "cli"`

---

## 2. Metadata / parameter override

```text
params["mask_fits"]
```

* Used if defined in metadata (e.g., `metadata.json`).
* `MASK_SOURCE = "params"`

---

## 3. Manual mask (GUI-edited)

```text
<tag>-mask-manual.fits
```

* Used if present in the cutout directory.
* Represents hand-edited mask from mask GUI.
* Overrides any auto-generated mask.
* `MASK_SOURCE = "manual"`

---

## 4. Existing auto mask

```text
<tag>-mask.fits
```

* Used if present and no higher-priority mask is found.
* Generated by pipeline.
* `MASK_SOURCE = "auto"`

---

## 5. Build new mask (if requested)

Triggered when:

```text
--make-mask
```

* Builds a new mask using pipeline masking logic.
* Writes `<tag>-mask.fits`
* `MASK_SOURCE = "auto"`

---

## 6. Force rebuild

```text
--force-mask
```

* Always rebuilds the auto mask (`<tag>-mask.fits`)
* Ignores existing masks (but does NOT overwrite manual masks)
* Archives previous auto mask before rebuilding
* `MASK_SOURCE = "auto"`

---

## 7. No mask

If none of the above apply:

* No mask is used
* `MASK_OK = False`
* `MASK_SOURCE = "none"`

---

## Key design principles

* Manual masks **always override** auto masks
* Auto masks are never overwritten by manual edits
* Mask provenance is tracked via `MASK_SOURCE`
* Mask naming is standardized:

  * auto: `<tag>-mask.fits`
  * manual: `<tag>-mask-manual.fits`

---

## Notes

* Mask quality directly affects all downstream measurements (photometry, morphology, sizes)
* Manual mask editing is expected for a subset of galaxies (e.g., bright stars, neighbors)
* `--force-mask` rebuilds only the auto mask and does not affect manual masks

---

If you want, I can also give you a short **one-line version** for commit messages or docstrings.
