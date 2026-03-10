from pathlib import Path
import shutil
import subprocess
import numpy as np
from astropy.io import fits
from importlib.resources import files


ASTROMATIC_DIR = files("hapy") / "astromatic"

DEFAULT_SUPPORT_FILES = [
    "default.param",
    "default.conv",
    "default.nnw",
]


def resolve_config_path(config):
    """Return a Path to the requested SExtractor config file."""
    config_path = Path(config)
    if config_path.is_file():
        return config_path

    packaged = Path(ASTROMATIC_DIR) / config
    if packaged.exists():
        return packaged

    raise FileNotFoundError(f"SExtractor config not found: {config}")


def copy_sextractor_files(config, workdir="."):
    """
    Copy the requested config file plus standard support files
    into the working directory.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    config_path = resolve_config_path(config)
    files_to_copy = [config_path] + [Path(ASTROMATIC_DIR) / f for f in DEFAULT_SUPPORT_FILES]

    for src in files_to_copy:
        dst = workdir / src.name
        if not dst.exists():
            if not src.exists():
                raise FileNotFoundError(f"Missing SExtractor file: {src}")
            shutil.copy2(src, dst)


def catalog_path(image, outdir="SEcats_psf", passnum=1):
    """Return output catalog path for an image."""
    image = Path(image)
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    return outdir / f"{image.stem}-se{passnum}.cat"


def weight_path(image):
    """Return weight filename corresponding to image."""
    return Path(image).with_suffix(".weight.fits")


def read_ldac_catalog(catname, ext=2):
    """Read a FITS_LDAC catalog."""
    return fits.getdata(catname, ext)


def estimate_fwhm_arcsec(catname, pixelscale, ext=2):
    """Estimate median FWHM in arcsec from a SExtractor catalog."""
    secat = read_ldac_catalog(catname, ext=ext)
    return np.median(secat["FWHM_IMAGE"]) * pixelscale


def build_sextractor_command(
    image,
    config,
    catalog_name,
    saturate=None,
    seeing_fwhm=None,
    weight_image=None,
    extra_args=None,
):
    """Build a SExtractor command list."""
    image = Path(image)
    config_path = resolve_config_path(config)

    cmd = [
        "sex",
        str(image),
        "-c", str(config_path),
        "-CATALOG_NAME", str(catalog_name),
    ]

    if saturate is not None:
        cmd += ["-SATUR_LEVEL", str(saturate)]
    else:
        cmd += ["-SATUR_LEVEL", "40000.0"]

    if seeing_fwhm is not None:
        cmd += ["-SEEING_FWHM", f"{seeing_fwhm:.1f}"]

    if weight_image is not None and Path(weight_image).exists():
        cmd += ["-WEIGHT_IMAGE", str(weight_image), "-WEIGHT_TYPE", "MAP_WEIGHT"]

    if extra_args:
        cmd += list(extra_args)

    return cmd


def run_sextractor_once(
    image,
    config,
    catalog_name,
    saturate=None,
    seeing_fwhm=None,
    weight_image=None,
    extra_args=None,
    #outdir="SEcats",
):
    """Run SExtractor once."""
    cmd = build_sextractor_command(
        image=image,
        config=config,
        catalog_name=catalog_name,
        saturate=saturate,
        seeing_fwhm=seeing_fwhm,
        weight_image=weight_image,
        extra_args=extra_args,
    )
    copy_sextractor_files(config, workdir=Path(catalog_name).parent)
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_sextractor_two_pass(
    image,
    config,
    pixelscale,
    saturate=None,
    overwrite=False,
    outdir="SEcats_psf",
    ext=2,
    extra_args=None,
):
    """
    Run SExtractor twice:
    pass 1 for initial FWHM estimate,
    pass 2 with SEEING_FWHM set from pass 1.
    """
    image = Path(image)
    se1 = catalog_path(image, outdir=outdir, passnum=1)
    se2 = catalog_path(image, outdir=outdir, passnum=2)
    wgt = weight_path(image)

    if overwrite or not se1.exists():
        run_sextractor_once(
            image=image,
            config=config,
            catalog_name=se1,
            saturate=saturate,
            weight_image=wgt,
            extra_args=extra_args,
        )
    else:
        print(f"using existing catalog: {se1}")

    fwhm = estimate_fwhm_arcsec(se1, pixelscale, ext=ext)

    if overwrite or not se2.exists():
        run_sextractor_once(
            image=image,
            config=config,
            catalog_name=se2,
            saturate=saturate,
            seeing_fwhm=fwhm,
            weight_image=wgt,
            extra_args=extra_args,
        )
    else:
        print(f"using existing catalog: {se2}")

    return {
        "se1": str(se1),
        "se2": str(se2),
        "fwhm_arcsec": fwhm,
    }
