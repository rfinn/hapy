# hapy

**hapy** (pronounced *happy*) is a Python package for astronomical image processing and analysis.
It provides tools for working with Hα imaging, segmentation masks, PSF construction, catalog handling, and general image utilities.

While originally developed for Hα projects, many modules are general-purpose and can be used in other astronomical workflows.

---

# Features

### Hα Image Processing

* Load, calibrate, and manipulate Hα and associated continuum images
* Handle image sets and coadds
* Perform photometric preparation steps

### Mask Construction & Editing

* Build segmentation masks using Source Extractor
* Grow masked regions
* Add Gaia star masks
* Interactive Qt-based GUI for mask editing
* Programmatic mask generation via `MaskEngine`

### Catalog Utilities

* Match and join astronomical catalogs
* Select and filter objects
* Cross-match with survey data

### PSF Tools

* Build PSFs from images
* Prepare PSFs for fitting software such as GALFIT

### GALFIT Tools

* Wrapper functions for running GALFIT
* GUI utilities for fitting workflows

### Image Utilities

* Download images from survey databases
* Flexible image display and visualization helpers

---

# Installation

Clone and install in editable mode:

```bash
git clone https://github.com/rfinn/hapy.git
cd hapy
pip install -e .
```

This installs the package and makes command-line scripts available.

---

# Usage

## Python API

Example usage in Python:

```python
from hapy.hatools import HalphaImageSet

image_set = HalphaImageSet(rimage, himage, psfdir=psfdir)
image_set.load_coadds()
```

---

## Masking (Interactive GUI)

An example script is provided:

```bash
python scripts/run_maskgui.py
```

Edit the script to point to your FITS image.

The GUI will:

* Build an initial segmentation mask
* Display r-band, Hα (optional), and mask panels
* Allow interactive mask editing
* Save output as:

```
<image>-mask.fits
```

Keyboard shortcuts inside the GUI:

| Key | Action                     |
| --- | -------------------------- |
| `r` | Remove object under cursor |
| `c` | Add circular mask          |
| `b` | Add square mask            |
| `g` | Grow mask                  |
| `w` | Write mask to disk         |
| `h` | Print help                 |
| `q` | Quit                       |

See `docs/masking.md` for a full description of the masking workflow.

---

## Command-Line Scripts

Example: create cutouts from an image and a catalog:

```bash
get_cutouts --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits \
            --catalog /path/to/catalog.fits \
            --virgo
```

The script processes the image, matches objects from the catalog, and generates cutouts.

---

# Project Structure

```
hapy/
│
├── README.md
├── docs/                 # Project documentation
│   └── masking.md
│
├── scripts/              # Runnable example scripts
│   └── run_maskgui.py
│
├── hapy/                 # Python package
│   ├── hatools/          # Hα-specific tools
│   ├── catools/          # Catalog utilities
│   ├── imagetools/       # Image utilities
│   ├── masktools/        # MaskEngine + mask operations
│   ├── maskgui/          # Qt GUI for mask editing
│   ├── galfittools/      # GALFIT helpers
│   └── astromatic/       # Configuration files for Astromatic tools
│
├── tests/                # Unit and integration tests
└── pyproject.toml
```

---

# Architecture Philosophy

Major components follow a clean separation of concerns:

* **Engine layers** contain core logic and data state
* **GUI layers** handle visualization and interaction
* **Utility modules** provide pure functions without UI dependencies

This allows:

* Batch processing without GUI
* Interactive workflows
* Easier testing and maintenance

---

# Dependencies

Core dependencies include:

* NumPy
* Astropy
* PyQt5 (for GUI tools)
* Ginga (image display in GUI)
* Source Extractor (external executable)

---

