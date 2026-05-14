import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from astropy.visualization import simple_norm
from astropy.stats import sigma_clip
from astropy.io import fits

from matplotlib.patches import Ellipse
from pathlib import Path

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


# def plot_segmentation_diagnostic(
#     r_fits,
#     se_seg_fits,
#     mask_fits,
#     phot_seg_fits,
#     e0,
#     eph,
#     outfile,
#     row,
# ):
#     Xsize = 10

#     r_data, r_hdr = fits.getdata(r_fits, header=True)
#     se_seg = fits.getdata(se_seg_fits)
#     m_data = fits.getdata(mask_fits)
#     phot_seg = fits.getdata(phot_seg_fits)

#     mmask = m_data > 0
#     objid = row.get("OBJID", "")

#     fig, ax = plt.subplots(1, 4, figsize=(20, 5))

#     # Panel 1: R image
#     plt.sca(ax[0])
#     display_image(r_data, mask=mmask)
#     ax[0].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
#                                   edgecolor="cyan", linewidth=2))
#     ax[0].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
#                                   edgecolor="magenta", linewidth=2))
#     ax[0].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
#     ax[0].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
#     ax[0].set_title(f"R image: {objid}")

#     # Panel 2: SE segmentation
#     from matplotlib.colors import LogNorm
#     seg_plot = np.array(se_seg, dtype=float)
#     seg_plot[seg_plot <= 0] = np.nan
#     ax[1].imshow(seg_plot, origin="lower", interpolation="nearest", norm=LogNorm(vmin=1, vmax=np.nanmax(seg_plot)),)
#     ax[1].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
#                                   edgecolor="cyan", linewidth=2))
#     ax[1].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
#                                   edgecolor="magenta", linewidth=2))
#     ax[1].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
#     ax[1].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
#     ax[1].set_title("SE segmentation")

#     # Panel 3: mask
#     ax[2].imshow(m_data, origin="lower", interpolation="nearest")
#     ax[2].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
#                                   edgecolor="cyan", linewidth=2))
#     ax[2].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, e0.ba, eph.theta_deg,
#                                   edgecolor="magenta", linewidth=2))
#     ax[2].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
#     ax[2].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
#     ax[2].set_title("Mask from SE")

#     # Panel 4: photutils segmentation
#     ax[3].imshow(phot_seg, origin="lower", interpolation="nearest")
#     ax[3].add_patch(ellipse_patch(e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
#                                   edgecolor="cyan", linewidth=2))
#     ax[3].add_patch(ellipse_patch(eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
#                                   edgecolor="magenta", linewidth=2))
#     ax[3].plot(e0.xc, e0.yc, 'cX', markersize=Xsize)
#     ax[3].plot(eph.xc, eph.yc, 'mX', markersize=Xsize)
#     ax[3].set_title("Photutils segmentation")

#     plt.tight_layout()
#     plt.savefig(outfile, dpi=150)
#     plt.close(fig)
#     print("Wrote:", outfile)


def _read_fits_if_exists(path, header=False):
    if path is None:
        return (None, None) if header else None
    path = Path(path)
    if not path.exists():
        return (None, None) if header else None
    return fits.getdata(path, header=header)

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
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from astropy.io import fits

    Xsize = 10

    r_data, r_hdr = _read_fits_if_exists(r_fits, header=True)
    se_seg = _read_fits_if_exists(se_seg_fits)
    m_data = _read_fits_if_exists(mask_fits)
    phot_seg = _read_fits_if_exists(phot_seg_fits)

    if r_data is None or m_data is None:
        print(f"Skipping diagnostic plot; missing required files for {outfile}")
        return

    mmask = m_data > 0
    objid = row.get("OBJID", "")

    fig, ax = plt.subplots(1, 5, figsize=(25, 5))

    def add_overlays(this_ax):
        this_ax.add_patch(
            ellipse_patch(
                e0.xc, e0.yc, e0.sma_pix, e0.ba, e0.theta_deg,
                edgecolor="cyan", linewidth=2
            )
        )
        this_ax.add_patch(
            ellipse_patch(
                eph.xc, eph.yc, eph.sma_pix, eph.ba, eph.theta_deg,
                edgecolor="magenta", linewidth=2
            )
        )
        this_ax.plot(e0.xc, e0.yc, "cX", markersize=Xsize)
        this_ax.plot(eph.xc, eph.yc, "mX", markersize=Xsize)

    # Panel 1: R image (unmasked)
    plt.sca(ax[0])
    display_image(r_data)
    add_overlays(ax[0])
    ax[0].set_title(f"R image: {objid}")

    # Panel 2: R image with mask applied
    plt.sca(ax[1])
    display_image(r_data, mask=mmask)
    add_overlays(ax[1])
    ax[1].set_title("R image + mask")

    # Colormap for segmentation images
    cmap = plt.cm.viridis.copy()
    cmap.set_bad(color="black")

    finite = np.isfinite(se_plot) & (se_plot > 0)

    if np.any(finite):
        vmax = np.nanmax(se_plot[finite])
        if vmax > 1:
            seg_norm = LogNorm(vmin=1, vmax=vmax)
        else:
            seg_norm = None
    else:
        seg_norm = None
    
    if se_seg is not None:
    # Panel 3: SE segmentation
        se_plot = np.array(se_seg, dtype=float)
        se_plot[se_plot <= 0] = np.nan
        ax[2].imshow(
            se_plot,
            origin="lower",
            interpolation="nearest",
            norm=seg_norm,
            cmap=cmap,
            )
        add_overlays(ax[2])
    else:
        ax[2].text(0.5, 0.5, "SE segmentation\nnot found",
               ha="center", va="center", transform=ax[2].transAxes)
    ax[2].set_title("SE segmentation")

    # Panel 4: mask
    ax[3].imshow(m_data, origin="lower", interpolation="nearest")
    add_overlays(ax[3])
    ax[3].set_title("Mask from SE")

    # Panel 5: photutils segmentation
    if phot_seg is not None:
        phot_plot = np.array(phot_seg, dtype=float)
        phot_plot[phot_plot <= 0] = np.nan
        ax[4].imshow(
            phot_plot,
            origin="lower",
            interpolation="nearest",
            norm=LogNorm(vmin=1, vmax=np.nanmax(phot_plot)),
            cmap=cmap,
            )
        add_overlays(ax[4])
    else:
        ax[4].text(0.5, 0.5, "Photutils segmentation\nnot found",
               ha="center", va="center", transform=ax[2].transAxes)
    ax[4].set_title("Photutils segmentation")

    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close(fig)
    print("Wrote:", outfile)
