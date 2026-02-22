# hapy

**hapy** (pronounced *happy*) is a Python package for astronomical image processing and analysis. It provides tools for working with Hα images, constructing PSFs, handling catalogs, and general image utilities. While originally developed for Hα projects, many modules are general-purpose and can be used in other astronomical workflows.

## Features

- **Hα image processing:** Load, calibrate, and manipulate Hα and associated images.  
- **Catalog handling:** Match and join astronomical catalogs, select objects, and filter based on criteria.  
- **PSF tools:** Build point-spread functions from images, for use in fitting software like GALFIT.
- **GALFIT tools:** Wrapper functions and gui for running GALFIT.  
- **Image utilities:** Download images from databases, display images with flexible visualization tools.  
- **Modular scripts:** Command-line scripts for common tasks such as making cutouts or analyzing images.  

## Installation

You can install **hapy** in editable mode for development:

```bash
git clone [https://github.com/your-username/hapy.git](https://github.com/rfinn/hapy.git)
source venv/bin/activate
pip install -e .
```

This will install the package and make scripts like get_cutouts available as commands.

## Usage
### Python API

Import modules in Python:

```python
from hapy.hatools import GalaxyCatalog, CoaddImage, HalphaImageSet, FilterTrace

# Load an image set
image_set = HalphaImageSet(args.rimage, himage, psfdir=args.psfdir)
image_set.load_coadds()
```

### Command-Line Scripts

Example: create cutouts from an image and a catalog:

```python
get_cutouts --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits \
            --catalog /path/to/catalog.fits \
            --virgo
```

The script will process the image, match objects from the catalog, and generate cutouts for analysis.

## Project Structure

```
hapy/
├── hapy/
│   ├── hatools/        # Hα-specific tools
│   ├── catools/        # Catalog utilities (general-purpose)
│   ├── imagetools/     # General image utilities
│   ├── galfittools/    # GALFIT helper functions
│   ├── scripts/        # Command-line scripts
│   └── astromatic/     # Configuration files for Astromatic tools
├── docs/               # Documentation, notebooks, wiki assets
├── tests/              # Unit and integration tests
├── pyproject.toml
└── README.md
```

