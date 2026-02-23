#import importlib.resources as pkg_resources
import os
import numpy as np

from astropy.wcs import WCS

# utils.py
from pathlib import Path
#import pkg_resources  # legacy, still works for pip-installed packages
#import importlib.resources as resources  # modern alternative
from hapy import hatools


from importlib import resources

def get_filter_file(name: str) -> str:
    """
    Return the full path to a filter file in hatools/filter_traces.

    Works both for development installs (editable) and normal installs.
    """
    try:
        p = resources.files("hapy.hatools.filter_traces").joinpath(name)
        if p.is_file():
            return str(p)
    except Exception:
        pass

    # editable/dev fallback
    base_dir = Path(hatools.__file__).parent
    fpath = base_dir / "filter_traces" / name
    if fpath.exists():
        return str(fpath)

    raise FileNotFoundError(f"Filter file {name} not found in hatools/filter_traces/")


