#!/usr/bin/env python

"""
python ~/github/hapy/scripts/run_photometry.py --r VFID2507-UGC08693-HDI-20170523-p007-R.fits --mask VFID2507-UGC08693-HDI-20170523-p007-R-mask.fits --prefix test1
"""

import argparse
from hapy.ellipse.photometry import run_ellipse_photometry

def main():
    p = argparse.ArgumentParser(description="Run elliptical photometry (headless)")
    p.add_argument("--r", dest="r_fits", required=True, help="R-band FITS")
    p.add_argument("--cs", dest="cs_fits", default=None, help="Continuum-subtracted FITS (optional)")
    p.add_argument("--mask", dest="mask_fits", default=None, help="Mask FITS (optional)")
    p.add_argument("--image2-filter", dest="image2_filter", default=None, help="e.g. 4, 8, 12, 16, inthalpha, intha6657")
    p.add_argument("--filter-ratio", dest="filter_ratio", type=float, default=None, help="flux ratio image2/image1 (optional)")
    p.add_argument("--objra", type=float, default=None)
    p.add_argument("--objdec", type=float, default=None)
    p.add_argument("--fixcenter", action="store_true")
    p.add_argument("--statmorph", action="store_true")
    p.add_argument("--prefix", default=None, help="Prefix tag for output tables")

    args = p.parse_args()

    print("input statmorph = ",args.statmorph)
    e = run_ellipse_photometry(
        r_fits=args.r_fits,
        cs_fits=args.cs_fits,
        mask_fits=args.mask_fits,
        image2_filter=args.image2_filter,
        filter_ratio=args.filter_ratio,
        objra=args.objra,
        objdec=args.objdec,
        fixcenter=args.fixcenter,
        write_prefix=args.prefix,
        run_statmorph=args.statmorph
    )

    # make plot of profiles
    e.plot_fancy_profiles()
    e.draw_phot_results_mpl()

    # for building tables
    print("EllipsePhotometry ATTRIBUTES")
    print(type(e), e.__dict__.keys())    

    if args.statmorph:
        print()
        print('STATMORPH IMAGE 1 ATTRIBUTES')
        print(type(e.morph),e.morph.__dict__.keys())

        print()
        print('STATMORPH IMAGE 2 ATTRIBUTES')        
        print(type(e.morph2),e.morph2.__dict__.keys())        
if __name__ == "__main__":
    main()
