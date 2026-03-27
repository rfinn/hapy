from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils as wcs_utils

from hapy.geometry.adapters import pa_ccw_north_to_photutils_theta
from hapy.masktools.api import EllipseParams, MaskEngine
from hapy.scripts.common_args import add_masking_args, add_object_args, add_weight_args


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="build_mask",
        description="Build a mask automatically for a FITS image.",
    )

    parser.add_argument(
        "--image",
        required=True,
        help="Primary FITS image, usually the r-band image.",
    )
    parser.add_argument(
        "--haimage",
        default=None,
        help="Optional second image, typically the continuum-subtracted CS image.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output mask FITS path. Default: <image stem>-mask.fits",
    )
    parser.add_argument(
        "--inv-output",
        default=None,
        help="Optional inverse-mask FITS output path",
    )
    parser.add_argument(
        "--no-remove-center",
        action="store_true",
        help="Do not remove the central galaxy object",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress messages",
    )

    add_masking_args(parser, include_force_mask=False, include_gaia_dir=True)
    add_weight_args(parser)
    add_object_args(parser)

    return parser


def default_mask_path(image_path: str | Path) -> Path:
    image_path = Path(image_path)
    stem = image_path.name[:-5] if image_path.name.lower().endswith(".fits") else image_path.stem
    return image_path.with_name(f"{stem}-mask.fits")


def parse_object_args(args: argparse.Namespace) -> Optional[tuple[float, float, float, float, float]]:
    vals = (args.objra, args.objdec, args.objsma, args.objba, args.objpa)

    if all(v is None for v in vals):
        return None

    if not all(v is not None for v in vals):
        raise ValueError(
            "If any object parameter is provided, all of "
            "--objra --objdec --objsma --objba --objpa must be provided."
        )

    return (
        float(args.objra),
        float(args.objdec),
        float(args.objsma),
        float(args.objba),
        float(args.objpa),
    )


def build_ellipse_from_sky(
    image_fits: str | Path,
    ra_deg: float,
    dec_deg: float,
    sma_arcsec: float,
    ba: float,
    pa_deg: float,
) -> EllipseParams:
    _, header = fits.getdata(image_fits, header=True)
    w = WCS(header)

    xc, yc = w.wcs_world2pix(ra_deg, dec_deg, 0)

    pixel_scales_deg = wcs_utils.proj_plane_pixel_scales(w)
    pixscale_deg = float(np.mean(pixel_scales_deg))
    sma_pix = float(sma_arcsec) / (pixscale_deg * 3600.0)

    theta_deg = pa_ccw_north_to_photutils_theta(pa_deg)

    return EllipseParams(
        xc=float(xc),
        yc=float(yc),
        sma_pix=float(sma_pix),
        ba=float(ba),
        theta_deg=float(theta_deg),
    )


def progress_cb(*, stage: str, fraction: float, message: str | None = None) -> None:
    msg = message or ""
    print(f"[{stage:10s}] {fraction:5.1%} {msg}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    image_path = Path(args.image)
    if not image_path.exists():
        parser.error(f"Image not found: {image_path}")

    if args.haimage is not None and not Path(args.haimage).exists():
        parser.error(f"Second image not found: {args.haimage}")

    try:
        objparams = parse_object_args(args)
    except ValueError as e:
        parser.error(str(e))

    galaxy_ellipse = None
    if objparams is not None:
        ra_deg, dec_deg, sma_arcsec, ba, pa_deg = objparams
        galaxy_ellipse = build_ellipse_from_sky(
            image_fits=image_path,
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            sma_arcsec=sma_arcsec,
            ba=ba,
            pa_deg=pa_deg,
        )

    engine = MaskEngine(
        image_fits=str(image_path),
        ha_image_fits=args.haimage,
        sepath=args.sepath,
        config=args.seconfig,
        threshold=args.sethreshold,
        snr=args.sesnr,
        minarea=args.seminarea,
        add_gaia_stars=not args.no_gaia,
        verbose=args.verbose,
        logger=None,
    )

    mask = engine.build_initial_mask(
        weightim=args.weightim,
        weight_threshold=float(args.weight_threshold),
        remove_center_object=not args.no_remove_center,
        galaxy_ellipse=galaxy_ellipse,
        grow_size=int(args.grow_size),
        grow_iterations=int(args.grow_iterations),
        progress_callback=progress_cb if args.verbose else None,
    )

    outpath = Path(args.output) if args.output else default_mask_path(image_path)
    engine.write_mask(
        str(outpath),
        inv_mask_fits=args.inv_output,
        overwrite=True,
    )

    if args.verbose:
        nmask = int(np.sum(mask > 0))
        print(f"Wrote mask to {outpath}")
        print(f"Masked pixels: {nmask}")
    else:
        print(outpath)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
