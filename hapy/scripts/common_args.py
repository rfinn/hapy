from __future__ import annotations


def add_masking_args(parser, *, include_force_mask: bool = False, include_gaia_dir: bool = True):
    """
    Add a consistent set of masking-related CLI args.

    Intended to be shared by:
      - run_analysis
      - build_mask
      - run_maskgui
    """
    g_mask = parser.add_argument_group("Masking Options")

    g_mask.add_argument(
        "--sepath",
        default="sex",
        help="Path to SExtractor executable",
    )
    g_mask.add_argument(
        "--seconfig",
        default="default.sex.HDI.mask",
        help="SExtractor config file path",
    )
    g_mask.add_argument(
        "--sethreshold",
        type=float,
        default=0.005,
        help="SExtractor detection/deblend threshold. Default is 0.005.",
    )
    g_mask.add_argument(
        "--sesnr",
        type=float,
        default=5.0,
        help="SExtractor SNR threshold. Default is 5.",
    )
    g_mask.add_argument(
        "--seminarea",
        type=int,
        default=5,
        help="SExtractor minimum object area. Default is 5.",
    )
    g_mask.add_argument(
        "--grow-size",
        type=int,
        default=7,
        help="Grow size in mask expansion. Default is 7.",
    )
    g_mask.add_argument(
        "--grow-iterations",
        type=int,
        default=4,
        help="Number of mask-growth iterations. Default is 4.",
    )

    if include_gaia_dir:
        g_mask.add_argument(
            "--gaia-dir",
            default="gaia_catalogs",
            help="Directory containing precomputed Gaia catalogs (default: gaia_catalogs)",
        )

    g_mask.add_argument(
        "--no-gaia",
        action="store_true",
        help="Disable Gaia star masking",
    )

    if include_force_mask:
        g_mask.add_argument(
            "--force-mask",
            action="store_true",
            help="Rebuild mask even if an existing mask file is present",
        )

    return g_mask


def add_object_args(parser):
    g_obj = parser.add_argument_group("Target Galaxy Options")

    g_obj.add_argument("--objra", type=float, default=None, help="Target RA in deg")
    g_obj.add_argument("--objdec", type=float, default=None, help="Target Dec in deg")
    g_obj.add_argument("--objsma", type=float, default=None, help="Target semi-major axis in arcsec")
    g_obj.add_argument("--objba", type=float, default=None, help="Target axis ratio b/a")
    g_obj.add_argument("--objpa", type=float, default=None, help="Target PA in deg (astronomy convention)")

    return g_obj


def add_weight_args(parser):
    g_wt = parser.add_argument_group("Weight Image Options")

    g_wt.add_argument(
        "--weightim",
        default=None,
        help="Optional weight image",
    )
    g_wt.add_argument(
        "--weight-threshold",
        type=float,
        default=1.0,
        help="Weight threshold passed to SExtractor / masking",
    )

    return g_wt
