#!/usr/bin/env python
"""
USAGE:
python ~/github/hapy/scripts/get_cutouts.py --rimage VF-165.869+28.044-HDI-20200226-p012-r.fits --catalog /Users/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo 
"""
import numpy as np
import os
from pathlib import Path
from astropy.table import Table
from astropy.io import fits
#from astropy.table import Table

import getpass
from datetime import datetime
import json

from hapy.hatools import GalaxyCatalog, CoaddImage, HalphaImageSet, FilterTrace
from hapy.hatools.utils import parse_coadd_name, build_cutout_name, get_survey_vectors
from hapy.utils.logging_utils import setup_logging


def resolve_sibling_path(base_image, sibling_name):
    """
    Resolve sibling_name relative to base_image directory unless already absolute.
    """
    sibling = Path(sibling_name)
    if sibling.is_absolute():
        return str(sibling)
    return str(Path(base_image).resolve().parent / sibling)

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
    
    parser.add_argument(
    "--overwrite",
    action="store_true",
    default = False,
    help="Overwrite existing cutouts if they already exist"
    )
   
    parser.add_argument('--maxcorrection', dest='maxcorrection', default=5., help='maximum filter correction for galaxies in FOV.  default is 5, so galaxies whose redshift falls where filter transmission < 20 percent will be skipped.')        
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

    # set up logging


    image_id = Path(args.rimage).stem

    outdir = args.outdir or "."

    log = setup_logging(
       outdir=outdir,
       tag=image_id,
       script_name="cutouts"
       )

    log.info("Starting get_cutouts with rimage=%s scheme=%s", args.rimage, args.scheme)


        

    # get images
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

    himage_full_path = resolve_sibling_path(args.rimage, himage)
    image_set = HalphaImageSet(args.rimage, himage_full_path, psfdir=args.psfdir)
    image_set.load_coadds()


    # default is to subtract sky locally in the r and halpha cutouts
    # then make CS cutout from those
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
        filter_keepflag = corrections < float(args.maxcorrection) # this is a crazy big cut, but we can adjust downstream
        gcat.keepflag[gcat.keepflag] = filter_keepflag
        filter_corrections = corrections[filter_keepflag]
        redshift = redshift[filter_keepflag]
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

    for i in range(len(gra)):

        size_arcsec = args.cutout_scale * 2 * gradius[i]

        # --------------------------------------------------
        # Pre-check validity on R and Ha weight images
        # --------------------------------------------------
        ok, status = image_set.cutout_location_is_valid(
            gra[i],
            gdec[i],
            size_arcsec
        )

        if not ok:
            #if "MOS" in args.rimage:
            #    print(f"Invalid regions in weight file for {args.rimage} - making cutout anyway")
            #else:
            print(f"Skipping {galid[i]}: invalid cutout region ({status}); ra={gra[i]:.6f},dec={gdec[i]:.6f}")
            continue

        # --------------------------------------------------
        # Build rootname only after validity check
        # --------------------------------------------------
        rootname = build_cutout_name(tokens, galid[i], cutouts_dir)
        cutdir = Path(rootname).parent
        params_path = cutdir / "metadata.json"

        # --------------------------------------------------
        # Metadata
        # --------------------------------------------------
        meta_redshift = None if redshift[i] is None else round(float(redshift[i]),6)
        if meta_redshift is not None:
            vr = round(float(meta_redshift*3.e5),1)
        else:
            vr = None
        params = dict(
            objid=str(galid[i]),
            tag=Path(rootname).name,
            root=str(rootname),
            ra=float(gra[i]),
            dec=float(gdec[i]),
            sma_arcsec=float(gradius[i]),
            ba=float(gBA[i]),
            pa_deg=float(gPA[i]),
            telescope=tokens.get("telescope"),
            dateobs=tokens.get("dateobs"),
            pointing=tokens.get("pointing"),
            scheme=args.scheme,
            parent_rimage=Path(args.rimage).name,
            parent_haimage=Path(himage).name if himage else None,
            hafilter=image_set.h.filter if himage else None,
            filter_correction=float(filter_corrections[i]),
            rimage_psf=image_set.r.psf_image_name,
            himage_psf=image_set.h.psf_image_name if himage else None,
            rimage_fwhm_se_arcsec=round(float(image_set.r.fwhm_se_arcsec),2) if image_set.r.fwhm_se_arcsec is not None else None,
            rimage_fwhm_psf_arcsec=round(float(image_set.r.fwhm_psf_arcsec),2) if image_set.r.fwhm_psf_arcsec is not None else None,
            himage_fwhm_se_arcsec=round(float(image_set.h.fwhm_se_arcsec),2) if (himage and image_set.h.fwhm_se_arcsec is not None) else None,
            himage_fwhm_psf_arcsec=round(float(image_set.h.fwhm_psf_arcsec),2) if (himage and image_set.h.fwhm_psf_arcsec is not None) else None,
            cutout_scale=float(args.cutout_scale),
            filter_ratio=float(filter_ratio) if filter_ratio is not None else None,
            cutout_sky_subtracted=bool(subtract_sky),
            valid_region=True,
            valid_status=str(status),
            rfilter_name = image_set.r.filter_file,
            rfilter_center_A = image_set.r.filter_center,
            rfilter_width_A = image_set.r.filter_width,
            hafilter_name = image_set.h.filter_file,
            hafilter_center_A = image_set.h.filter_center,
            hafilter_width_A = image_set.h.filter_width,
            redshift=meta_redshift,
            vr=vr,
        )

        if params_path.exists() and args.overwrite_metadata:
            backup = params_path.with_suffix(".json.bak")
            params_path.replace(backup)

        params_path.write_text(json.dumps(params, indent=2))

        # --------------------------------------------------
        # Make cutouts
        # --------------------------------------------------
        result = image_set.get_cutout_all_filters(
            gra[i],
            gdec[i],
            size_arcsec,
            rootname,
            subtract_sky=subtract_sky,
            overwrite=args.overwrite
        )

        # If get_cutout_all_filters later returns an invalid/failure flag,
        # you can handle it here. For now, assume pre-check was sufficient.

        # --------------------------------------------------
        # Parent pixel position
        # --------------------------------------------------
        x, y = image_set.h.wcs.world_to_pixel_values(gra[i], gdec[i])

        filter_warning = float(filter_corrections[i]) > 2
        rows.append(dict(
            objid=str(galid[i]),
            tag=Path(rootname).name,
            parent_rimage=Path(args.rimage).name,
            parent_haimage=Path(himage).name if himage is not None else None,
            telescope=tokens.get("telescope"),
            dateobs=tokens.get("dateobs"),
            pointing=tokens.get("pointing"),
            scheme=args.scheme,
            ra=float(gra[i]),
            dec=float(gdec[i]),
            hafilter=image_set.h.filter if himage else None,
            filter_correction=float(filter_corrections[i]),
            filter_warning=filter_warning,
            size_arcsec=float(size_arcsec),
            cutout_root=str(rootname),
            x_parent=float(x),
            y_parent=float(y),
            valid_region=True,
            valid_status=str(status),
        ))


    tab = Table(rows=rows)

    user = getpass.getuser()
    #ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    ts = datetime.now().strftime("%Y%m%d")

    #stem = Path(rootname).stem
    tag= Path(args.rimage).name.replace("-R.fits","").replace("-r.fits","")
    #summary_path = Path(outdir) / f"cutouts_summary-{tag}-{user}-{ts}.fits"
    #tab.write(summary_path, overwrite=True)
    summary_dir = Path(outdir) / "cutouts_summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    summary_file = f"cutouts_summary-{tag}-{user}-{ts}.ecsv"
    summary_path = summary_dir / summary_file
    #summary_path = Path(outdir) / f"cutouts_summary-{tag}-{user}-{ts}.ecsv"
    tab.write(summary_path, overwrite=True, format="ascii.ecsv")
    print(f"Wrote summary table: {summary_path}")
if __name__ == '__main__':


    main()
