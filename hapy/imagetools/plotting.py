import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from astropy.visualization import simple_norm
from astropy.stats import sigma_clip
from astropy.io import fits

from matplotlib.patches import Ellipse
def ellipse_patch(xc, yc, sma_pix, ba, pa_deg, **kwargs):
    # matplotlib Ellipse needs angle relative to +x axis
    return Ellipse(
        (xc, yc),
        width=2 * sma_pix,
        height=2 * sma_pix * ba,
        angle=pa_deg,
        fill=False,
        **kwargs,
    )


def display_image(image, percent=99.9, lowrange=False, sigclip=True,mask=None,cmap='gray_r',zoom=None):
    masked=False
    if mask is not None:
        #image = np.ma.array(image, mask=mask)
        masked=True
        statarray = image.copy()
        statarray[mask] = np.nan
    if sigclip and (mask is None):
        clipped_data = sigma_clip(image, sigma_lower=5, sigma_upper=5)
    else:
        clipped_data = image

    stretch = "linear" if lowrange else "asinh"
    if mask is not None:
        norm = simple_norm(statarray, stretch=stretch, percent=percent)
        masked_image = np.ma.array(image, mask=mask)
        plt.imshow(masked_image, norm=norm,origin="lower", cmap=cmap)#, origin="lower")
    else:
        norm = simple_norm(clipped_data, stretch=stretch, percent=percent)
        plt.imshow(image, norm=norm,origin="lower", cmap=cmap)#, origin="lower")


def display_unwise(ra,dec,galname,imsize_arcsec=60):
    # get unwise images

    imsize_pixels_legacy = round(imsize_arcsec/LEGACY_PIXSCALE)
    imsize_pixels_unwise = round(imsize_arcsec/UNWISE_PIXSCALE)
    
    t = get_unwise_image(ra,dec,galid=galname,makeplots=False,imsize=str(imsize_pixels_unwise),verbose=verbose)
    imagefiles = t[0]
    noisefiles = t[1]
    imagefiles.sort()
    noisefiles.sort()
    plt.figure(figsize=(12,4))

    # concatinate lists

    # plot WISE images
    imnames = [galname+' W1','W2','W3','W4']
    for i,im in enumerate(imagefiles):
        plt.subplot(2,4,i+1)
        data = fits.getdata(im)
        display_image(data,percent=92)
        plt.title(imnames[i],fontsize=14)    
    
def display_legacy_unwise(ra,dec,galname,imsize_arcsec=60,plotdir=None,verbose=False):
    """

    download and display legacy and unwise images

    INPUT:
    ra in deg
    dec in deg
    galid = galaxy name, this will be prefix of output images
    imsize = length/width of image in arcsec

    RETURN:
    list of legacy images, include jpg as the first element of this list

    list of wise image names

    """

    if plotdir is None:
        plotdir = 'plots/'
        
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)
    
    imsize_pixels_legacy = round(imsize_arcsec/LEGACY_PIXSCALE)
    imsize_pixels_unwise = round(imsize_arcsec/UNWISE_PIXSCALE)

    # get unwise images
    t = get_unwise_image(ra,dec,galid=galname,makeplots=False,imsize=str(imsize_pixels_unwise))
    imagefiles = t[0]
    noisefiles = t[1]
    imagefiles.sort()
    noisefiles.sort()

    # get legacy images

    bands = ['g','r','z']
    legimfiles = []
    plot_legacy = True
    for i,b in enumerate(bands):
        print(i,b)
        if i == 0:
            # only need to download this once
            getjpg = True
        else:
            getjpg = False
        try:
            t = get_legacy_images(ra,dec,galid=galname,band=b,makeplots=False,imsize=str(imsize_pixels_legacy),verbose=verbose)
        except:
            print(f"WARNING: got a http error when downloading {b} image for {galname}")
            continue

        if t is None:
            print(f"WARNING: GALAXY {galname} IS OUTSIDE LEGACY FOOTPRINT")
            plot_legacy = False
            break
        print('return from get_legacy_images = ',t,b,galname)
        if i == 0:
            legimfiles.append(t[0])
            legjpgfile = t[1]
        else:
            legimfiles.append(t[0])

    if plot_legacy:
        # make a plot
        plt.figure(figsize=(12,6.5))
        plt.gca().set_facecolor("white")
        # concatinate lists
        legacy_images = [legjpgfile]+legimfiles
        imnames = [galname+' grz','g','r','z']
        # plot legacy images in top row
        for i,im in enumerate(legacy_images):
            plt.subplot(2,4,i+1)
            if not os.path.exists(im):
                continue
            if i == 0:
                # display jpg
                t = Image.open(im)
                plt.imshow(t,origin='upper')
            else:
                data = fits.getdata(im)
                display_image(data,lowrange=False,percent=95)
        
            plt.title(imnames[i],fontsize=14)
        Noffset = 4
        Nrow=2
        
    else:
        plt.figure(figsize=(12,3))
        plt.gca().set_facecolor("white")        
        Noffset = 0
        Nrow=1    
    # plot WISE images
    imnames = [galname+' W1','W2','W3','W4']
    for i,im in enumerate(imagefiles):
        plt.subplot(Nrow,4,Noffset+i+1)
        data = fits.getdata(im)
        display_image(data,percent=92)
        plt.title(imnames[i],fontsize=14)
    plt.savefig(plotdir+galname+"cutouts.png",transparent=False)
    plt.close()
    return legacy_images, imagefiles



def plot_mask_ellipse_diagnostic(
    r_fits,
    mask_fits,
    e0,
    eph,
    outfile,
    row
):
    Xsize=10
    r_data, r_hdr = fits.getdata(r_fits, header=True)
    m_data = fits.getdata(mask_fits)
    mmask = m_data > 0
    rmasked = np.ma.array(r_data, mask=mmask)
    #plt.figure()
    #plt.imshow(rmasked, origin="lower")
    #plt.savefig("debug-plot.png")
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    plt.sca(ax[0])
    display_image(r_data, mask=mmask)
    #display_image(rmasked)
    #plt.colorbar()
    ax[0].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
                                  edgecolor="cyan", linewidth=2))
    ax[0].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
                                  edgecolor="magenta", linewidth=2))

    ax[0].plot(e0.xc, e0.yc,'cX',markersize=Xsize)
    ax[0].plot(eph.xc, eph.yc,'mX',markersize=Xsize)    
    objid = row.get("OBJID")
    #print("DEBUG: in plotting.plot_mask_ellipse_diagnostic, objid=",objid)
     
    ax[0].set_title(f"R Image: {objid}")

    ax[1].imshow(m_data, origin="lower")
    #plt.sca(ax[1])
    #display_image(m_data)
    ax[1].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
                                  edgecolor="cyan", linewidth=2))
    ax[1].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
                                  edgecolor="magenta", linewidth=2))

    # mark center
    ax[1].plot(e0.xc, e0.yc,'cX',markersize=Xsize)
    ax[1].plot(eph.xc, eph.yc,'mX',markersize=Xsize)        
    #ax[1].add_patch(ellipse_patch(xc, yc, sma_pix, params["ba"], theta_photutils, edgecolor="cyan", linewidth=2))
    ax[1].set_title("Mask")

    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close(fig)
    print("Wrote:", outfile)


def plot_segmentation_diagnostic(
    r_fits,
    se_seg_fits,
    mask_fits,
    phot_seg_fits,
    e0,
    eph,
    outfile,
    row,
):
    Xsize = 10

    r_data, r_hdr = fits.getdata(r_fits, header=True)
    se_seg = fits.getdata(se_seg_fits)
    m_data = fits.getdata(mask_fits)
    phot_seg = fits.getdata(phot_seg_fits)

    mmask = m_data > 0
    objid = row.get("OBJID", "")

    fig, ax = plt.subplots(1, 4, figsize=(20, 5))

    # Panel 1: R image
    plt.sca(ax[0])
    display_image(r_data, mask=mmask)
    ax[0].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
                                  edgecolor="cyan", linewidth=2))
    ax[0].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
                                  edgecolor="magenta", linewidth=2))
    ax[0].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
    ax[0].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
    ax[0].set_title(f"R image: {objid}")

    # Panel 2: SE segmentation
    ax[1].imshow(se_seg, origin="lower", interpolation="nearest")
    ax[1].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
                                  edgecolor="cyan", linewidth=2))
    ax[1].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
                                  edgecolor="magenta", linewidth=2))
    ax[1].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
    ax[1].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
    ax[1].set_title("SE segmentation")

    # Panel 3: mask
    ax[2].imshow(m_data, origin="lower", interpolation="nearest")
    ax[2].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
                                  edgecolor="cyan", linewidth=2))
    ax[2].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, e0.ba, eph.theta_deg,
                                  edgecolor="magenta", linewidth=2))
    ax[2].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
    ax[2].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
    ax[2].set_title("Mask from SE")

    # Panel 4: photutils segmentation
    ax[3].imshow(phot_seg, origin="lower", interpolation="nearest")
    ax[3].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
                                  edgecolor="cyan", linewidth=2))
    ax[3].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
                                  edgecolor="magenta", linewidth=2))
    ax[3].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
    ax[3].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
    ax[3].set_title("Photutils segmentation")

    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close(fig)
    print("Wrote:", outfile)
