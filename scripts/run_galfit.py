"""
How to run:

python test_rungalfit.py \
  --image test_images/VFID2507-UGC08693-HDI-20170523-p007-R.fits \
  --mask_image test_images/VFID2507-UGC08693-HDI-20170523-p007-R-mask.fits \
  --psf_image test_images/MKW8_R.coadd-psf.fits \
  --sigma_image test_images/VFID2507-UGC08693-HDI-20170523-p007-R.mask.fits \
  --ncomp 1
"""

import argparse
from astropy.io import fits

from hapy.galfittools.rungalfit import RunGalfit, parse_galfit_results_dc
from hapy.imagetools.imutils import get_pixel_scale_from_filename


def main():
    ap = argparse.ArgumentParser(description="Run GALFIT on a test image")
    ap.add_argument("--image", required=True, help="input FITS image")
    ap.add_argument("--sigma_image", default=None, help="sigma/RMS image (optional)")
    ap.add_argument("--psf_image", default=None, help="PSF image (optional)")
    ap.add_argument("--psf_oversampling", type=int, default=2, help="PSF oversampling (often 2)")
    ap.add_argument("--mask_image", default=None, help="mask FITS (optional)")
    ap.add_argument("--convflag", type=int, default=1, help="1 to convolve with PSF, 0 otherwise")
    ap.add_argument("--ncomp", type=int, default=1, choices=[1, 2], help="number of sersic components")
    ap.add_argument("--magzp", type=float, default=None, help="mag zeropoint (default: PHOTZP header or 25.0)")
    ap.add_argument("--sky", type=float, default=0.0, help="initial sky guess (ADU)")
    args = ap.parse_args()

    galname = args.image.replace(".fits", "")
    pscale = get_pixel_scale_from_filename(args.image)

    # image shape -> default fit region = full image
    data, hdr = fits.getdata(args.image, header=True)
    ny, nx = data.shape
    xminfit, xmaxfit = 1, nx
    yminfit, ymaxfit = 1, ny

    # magzp default
    if args.magzp is not None:
        magzp = args.magzp
    else:
        magzp = float(hdr.get("PHOTZP", 25.0))

    convflag = bool(args.convflag)

    rg = RunGalfit(
        galname=galname,
        image=args.image,
        sigma_image=args.sigma_image,
        psf_image=args.psf_image,
        psf_oversampling=args.psf_oversampling,
        mask_image=args.mask_image,
        xminfit=xminfit,
        yminfit=yminfit,
        xmaxfit=xmaxfit,
        ymaxfit=ymaxfit,
        convolution_size=min(nx, ny),
        magzp=magzp,
        pscale=pscale,
        convflag=convflag,
        fitallflag=False,
        ncomp=args.ncomp,
        asym=False,
    )

    # crude initial guess at center
    xc = nx / 2
    yc = ny / 2

    # very basic initial guesses (fine for smoke testing)
    rg.set_sersic_params(
        xobj=xc, yobj=yc,
        mag=15.0,
        rad=10.0,
        nsersic=2.0,
        BA=0.7,
        PA=0.0,
        fitmag=1, fitcenter=1, fitrad=1, fitBA=1, fitPA=1, fitn=1,
        first_time=1,
    )
    rg.set_sky(args.sky)

    rg.run_galfit(displayflag=False)

    if getattr(rg, "galfit_flag", 0) == 0:
        raise RuntimeError("GALFIT did not complete (galfit_flag=0). Check stdout and fit.log")

    # parse & print dataclass results
    res = parse_galfit_results_dc(rg.output_image, ncomp=args.ncomp, asymflag=False)
    print("\n=== Parsed results (dataclass) ===")
    print("CHI2NU:", res.chi2nu)
    print("SKY:", res.sky, "+/-", res.sky_err)
    print("COMP1:", res.comp1)
    if res.comp2 is not None:
        print("COMP2:", res.comp2)

    print("\nOutput FITS:", rg.output_image)
    print("Log file:", rg.galfit_log)
    print("Galfit .01 renamed to:", rg.galfit_out)


if __name__ == "__main__":
    main()
