#!/usr/bin/env python
"""



USAGE:
python ~/github/hapy/scripts/get_cutouts.py --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo 
"""
from astropy.io import fits
#from astropy.table import Table
import numpy as np
import os
from hapy.hatools import GalaxyCatalog, CoaddImage, HalphaImageSet, FilterTrace
from hapy.hatools.utils import parse_coadd_name, build_cutout_name, get_survey_vectors
from pathlib import Path
from astropy.table import Table

import getpass
from datetime import datetime
import json

def main(args=None):

    import argparse

    parser = argparse.ArgumentParser(description ='create psf image from image that contains stars')

    #parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    parser.add_argument('--rimage', help='r-band image name')
    parser.add_argument("--outdir",type=str,default=None,
                            help="Directory where cutouts/ will be created (default: current working directory).")
    parser.add_argument('--psfdir', help='set to coadd directory')
    parser.add_argument('--cutout_scale',type=float, default=2, help='multiplicative scale factor for increasing the size of cutout images')
    parser.add_argument('--overwrite_metadata',default=False, action='store_true', help='Set this to overwrite metadata.json.  Will store a *.bak file.')
    parser.add_argument("--no-skysub", action="store_true",help="Disable local sky subtraction in cutouts (default: sky is subtracted).")
    parser.add_argument('--catalog',
                            help='full path to galaxy catalog to use for cutouts.  ')
    #parser.add_argument('--outdir',  default='cutouts',
    #                        help='base output directory for cutouts (default: cutouts)')
    parser.add_argument(
    "--scheme",
    choices=["generic", "virgo", "agc"],
    default="generic",
    help="Filename parsing scheme for coadd images.")
   
    parser.add_argument('--maxcorrection',dest='maxcorrection', default=3, help='maximum filter correction for galaxies in FOV.  default is 3, so galaxies whose redshift falls where filter transmission < 33 percent will be skipped.')        
    #parser.add_argument('--oneimage',dest = 'oneimage',default=None, help='give full path to the r-band image name to run on just one image')
    
    args = parser.parse_args()
    
    if args.outdir is None:
        outdir = os.getcwd()
    else:
        outdir = args.outdir
    base_outdir = Path(outdir).resolve() if args.outdir else Path.cwd()
    cutouts_dir = base_outdir / "cutouts"
    cutouts_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing cutouts to: {cutouts_dir}")
    
    try:
        rheader = fits.getheader(args.rimage)
        himage = rheader['HAIMAGE']
    except KeyError:
        print(f"WARNING: could not get HAIMAGE in rimage header {args.rimage}")
        himage = None
    try:
        rheader = fits.getheader(args.rimage)
        filter_ratio = float(rheader['FLTRATIO'])
    except KeyError:
        print(f"WARNING: could not get FLTRATIO in rimage header {args.rimage}.  Make sure you ran filter ratio program!")
        filter_ratio = None
    image_set = HalphaImageSet(args.rimage, himage, psfdir=args.psfdir)
    image_set.load_coadds()


    subtract_sky = not args.no_skysub
    
    ###################################################    
    # get galaxy catalog
    ###################################################
    agcflag = args.scheme == 'agc'
    virgoflag = args.scheme == 'virgo'
    gcat = GalaxyCatalog(args.catalog,nsa=False,agc=agcflag,virgo=virgoflag,sizecat=None, verbose=False)

    
    gcat.galaxies_in_fov(image_set.h.wcs,zmin=None,zmax=None,image_name=himage, agcflag=agcflag,virgoflag=virgoflag)

    ###################################################
    # get redshift from filter trace module
    ###################################################    


    #redshift_full, galid_full = get_survey_vectors(gcat, args.scheme)
    #redshift = None if redshift_full is None else redshift_full[gcat.keepflag]
    #galid = np.asarray(galid_full)[gcat.keepflag]


    redshift_full, galid_full = get_survey_vectors(gcat, args.scheme)
    redshift = None if redshift_full is None else redshift_full[gcat.keepflag]

    if redshift is not None:
        myfilter = FilterTrace(image_set.h.filter, instrument=image_set.h.instrument)
        corrections = myfilter.get_trans_correction(redshift,outfile=None)
        filter_keepflag = corrections < args.maxcorrection # this is a crazy big cut, but we can adjust with halphagui
        gcat.keepflag[gcat.keepflag] = filter_keepflag
        filter_corrections = corrections[filter_keepflag]

    print(f"number of galaxies in FOV = {np.sum(gcat.keepflag)}")
    
    ###################################################
    # get galaxy cutouts
    ###################################################    
    gra = gcat.RA[gcat.keepflag]
    gdec = gcat.DEC[gcat.keepflag]
    gradius = gcat.radius_arcsec[gcat.keepflag]
    gBA = gcat.BA[gcat.keepflag]
    gPA = gcat.PA_DEG[gcat.keepflag]
    galid = np.asarray(galid_full)[gcat.keepflag]

    rows = []
    tokens = parse_coadd_name(args.rimage, scheme=args.scheme)
    outbase = Path(outdir)
    outbase.mkdir(parents=True, exist_ok=True)
    #tokens = parse_coadd_name(args.rimage, scheme=args.scheme)
    for i in range(len(gra)):
        

        #rootname = build_cutout_name(tokens, galid[i], args.outdir)
        rootname = build_cutout_name(tokens, galid[i], cutouts_dir)

        # Write mask ellipse params for downstream masking
        cutdir = Path(rootname).parent
        params_path = cutdir / "metadata.json"
        params = dict(
                objid=str(galid[i]),
                tag=Path(rootname).name,          # e.g., VFID3084-NGC3512-HDI-20200226-p012
                root=str(rootname),              # full cutout root path
                ra=float(gra[i]),
                dec=float(gdec[i]),
                sma_arcsec=float(gradius[i]),   # semi-major = radius
                ba=float(gBA[i]),
                pa_deg=float(gPA[i]),
                telescope=tokens.get("telescope"),
                dateobs=tokens.get("dateobs"),
                pointing=tokens.get("pointing"),
                scheme=args.scheme,
                parent_rimage=Path(args.rimage).name,
                parent_haimage=Path(himage).name if himage else None,
                hafilter=image_set.h.filter if himage else None,
                filter_correction = float(filter_corrections[i]),
                rimage_psf= image_set.r.psf_image_name,
                himage_psf= image_set.h.psf_image_name if himage else None,                
                rimage_fwhm_arcsec= float(image_set.r.fwhm_arcsec) if image_set.r.fwhm_arcsec is not None else None, 
                rimage_fwhm_pixels= float(image_set.h.fwhm_pixels) if image_set.r.fwhm_pixels is not None else None,
                himage_fwhm_arcsec= float(image_set.r.fwhm_arcsec) if image_set.h.fwhm_arcsec is not None else None, 
                himage_fwhm_pixels= float(image_set.h.fwhm_pixels) if image_set.h.fwhm_pixels is not None else None,
                cutout_scale = float(args.cutout_scale),
                filter_ratio = filter_ratio,
                cutout_sky_subtracted = subtract_sky,

        )
        if params_path.exists() and args.overwrite_metadata:
            backup = params_path.with_suffix(".json.bak")
            params_path.replace(backup)
        params_path.write_text(json.dumps(params, indent=2))

        # commenting the next line for testing
        image_set.get_cutout_all_filters(gra[i], gdec[i], args.cutout_scale*2*gradius[i], rootname, subtract_sky=subtract_sky)
        # parent pixel position
        x, y = image_set.h.wcs.world_to_pixel_values(gra[i], gdec[i])
        #print(f"{rootname}: ra={gra[i]:.6f}, dec={gdec[i]:.6f}, radius_arcsec={gradius[i]:6.2f}, BA={gBA[i]:.2f}, PA={gPA[i]:5.1f}, x={x:.1f}, y={y:.1f}")
        rows.append(dict(
        objid=str(galid[i]),
        parent_rimage=Path(args.rimage).name,
        parent_haimage=Path(himage).name if himage is not None else None,
        telescope=tokens.get("telescope"),
        dateobs=tokens.get("dateobs"),
        pointing=tokens.get("pointing"),
        ra=float(gra[i]),
        dec=float(gdec[i]),
        size_arcsec=float(2 * gradius[i]),
        cutout_root=str(rootname),
        x_parent=float(x),
        y_parent=float(y),
            ))
    tab = Table(rows=rows)

    user = getpass.getuser()
    #ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    ts = datetime.now().strftime("%Y%m%d")

    stem = Path(rootname).stem
    summary_path = Path(outdir) / f"cutouts_summary-{stem}-{user}-{ts}.fits"
    tab.write(summary_path, overwrite=True)
    print(f"Wrote summary table: {summary_path}")
if __name__ == '__main__':


    main()
