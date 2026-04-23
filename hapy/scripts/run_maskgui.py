from __future__ import annotations

import argparse
import sys
from pathlib import Path

from PyQt5 import QtWidgets

from hapy.maskgui.mask_window import MaskWindow
from hapy.scripts.common_args import add_masking_args, add_object_args, add_weight_args


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_maskgui",
        description="Launch the interactive masking GUI.",
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
        "--title",
        default="Mask GUI",
        help="Window title",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print basic launch info",
    )

    add_masking_args(parser, include_force_mask=False, include_gaia_dir=True)
    add_weight_args(parser)
    add_object_args(parser)

    return parser


def parse_objparams(args: argparse.Namespace):
    vals = (args.objra, args.objdec, args.objsma, args.objba, args.objpa)

    if all(v is None for v in vals):
        return None

    if not all(v is not None for v in vals):
        raise ValueError(
            "If any object parameter is provided, all of "
            "--objra --objdec --objsma --objba --objpa must be provided."
        )

    # convert args.objpa to CCW from +x axis
    theta_deg = float(args.objpa) - 90.
    return [
        float(args.objra),
        float(args.objdec),
        float(args.objsma),
        float(args.objba),
        theta_deg,
    ]


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    image_path = Path(args.image)
    if not image_path.exists():
        parser.error(f"Image not found: {image_path}")

    if args.haimage is not None and not Path(args.haimage).exists():
        parser.error(f"Second image not found: {args.haimage}")

    try:
        objparams = parse_objparams(args)
    except ValueError as e:
        parser.error(str(e))

    if args.verbose:
        print(f"Launching mask GUI for {image_path}")
        if args.haimage is not None:
            print(f"Secondary image: {args.haimage}")
        if objparams is not None:
            print(f"Object parameters: {objparams}")

    app = QtWidgets.QApplication.instance()
    owns_app = app is None
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    mainwin = QtWidgets.QMainWindow()
    form = QtWidgets.QWidget(mainwin)
    mainwin.setCentralWidget(form)

    mw = MaskWindow(
        form,
        logger=None,
        image=str(image_path),
        haimage=args.haimage,
        sepath=args.sepath,
        config=args.seconfig,
        threshold=args.sethreshold,
        snr=args.sesnr,
        minarea=args.seminarea,
        objparams=objparams,
        auto=False,
        ngrow=args.grow_iterations,
        weightim=args.weightim,
        weight_threshold=args.weight_threshold,
        gaia_catalog=args.gaia_catalog,
        gaia_min_radius=args.gaia_min_radius,
        addgaia=not args.no_gaia,
    )

    title = args.title
    if title == "Mask GUI":
        title = f"Mask GUI - {image_path.name}"

    mainwin.setWindowTitle(title)
    mainwin.resize(1200, 800)
    mainwin.show()

    if owns_app:
        return app.exec_()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
