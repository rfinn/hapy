from pathlib import Path
import hapy

def astromatic_dir() -> Path:
    """
    Return the directory containing astromatic config files inside the installed hapy package.
    """
    return Path(hapy.__file__).resolve().parent / "astromatic"
