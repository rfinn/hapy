from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Any

import numpy as np

# import statmorph ONLY here
import statmorph
from statmorph.utils.image_diagnostics import make_figure

from astropy.utils import lazyproperty
import scipy.ndimage as ndi


class MyStatmorph(statmorph.SourceMorphology):
    """Statmorph subclass that forces gini segmap behavior."""

    @lazyproperty
    def _segmap_gini(self):
        segmap = np.array(self._segmap.data == 1, "i")
        return segmap[self._slice_stamp]

    def print(self):
        ''' adding a print method to print out the instance variables '''
        for k in self.__dict__.keys():
            if k.startswith('_'):
                continue
            print(f"{k}: {self.__dict__[k]}")

    #def plot_segmap_weightmap(self):
    #    plt.figure()
        
    import warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.utils.exceptions import AstropyUserWarning


class MyStatmorph(statmorph.SourceMorphology):
    """Statmorph subclass that forces gini segmap behavior."""

    @lazyproperty
    def _segmap_gini(self):
        segmap = np.array(self._segmap.data == 1, dtype=int)
        return segmap[self._slice_stamp]

    def print(self):
        """Print public instance variables."""
        for k in self.__dict__.keys():
            if k.startswith('_'):
                continue
            print(f"{k}: {self.__dict__[k]}")

    def _get_stamp_image(self):
        for name in ["_cutout_stamp_maskzeroed", "_cutout_stamp", "_cutout"]:
            if hasattr(self, name):
                arr = getattr(self, name)
                if arr is not None:
                    return np.asarray(arr)
        return None

    def _get_segmap_stamp(self):
        try:
            return np.asarray(self._segmap.data[self._slice_stamp])
        except Exception:
            return None

    def _get_sigma_map(self):
        """
        statmorph 'weightmap' = sigma image.
        Try likely internal names.
        """
        for name in ["_weightmap_stamp", "_weightmap", "weightmap"]:
            if hasattr(self, name):
                arr = getattr(self, name)
                if arr is None:
                    continue
                arr = np.asarray(arr)

                # If full-size image, crop to stamp when possible.
                if hasattr(self, "_slice_stamp") and arr.ndim == 2:
                    img = self._get_stamp_image()
                    if img is not None and arr.shape != img.shape:
                        try:
                            arr = arr[self._slice_stamp]
                        except Exception:
                            pass
                return arr
        return None

    def _get_mask_stamp(self):
        for name in ["_mask_stamp", "_mask"]:
            if hasattr(self, name):
                arr = getattr(self, name)
                if arr is None:
                    continue
                arr = np.asarray(arr)
                img = self._get_stamp_image()
                if img is not None and arr.shape != img.shape and hasattr(self, "_slice_stamp"):
                    try:
                        arr = arr[self._slice_stamp]
                    except Exception:
                        pass
                return arr.astype(bool)
        return None

    def _get_center_xy_stamp(self):
        for xname, yname in [
            ("_xc_stamp", "_yc_stamp"),
            ("xc_centroid", "yc_centroid"),
            ("xcentroid", "ycentroid"),
        ]:
            if hasattr(self, xname) and hasattr(self, yname):
                try:
                    return float(getattr(self, xname)), float(getattr(self, yname))
                except Exception:
                    pass
        return None

    def plot_diagnostic(self, outfile=None, title=None, show=False, dpi=150):
        """
        6-panel statmorph QC plot.

        Panels:
          1. science stamp
          2. segmentation map
          3. Gini segmap
          4. sigma image used by statmorph
          5. valid pixels used for sn_per_pixel diagnostic
          6. science image with overlays
        """
        img = self._get_stamp_image()
        seg = self._get_segmap_stamp()
        seg_gini = np.asarray(self._segmap_gini) if hasattr(self, "_segmap_gini") else None
        sigma = self._get_sigma_map()
        mask = self._get_mask_stamp()

        if img is None:
            raise RuntimeError("Could not find statmorph stamp image.")

        seg_gini_bool = np.zeros_like(img, dtype=bool) if seg_gini is None else seg_gini.astype(bool)

        finite_img = np.isfinite(img)

        if sigma is not None:
            finite_sigma = np.isfinite(sigma)
            pos_sigma = sigma > 0
            locs = seg_gini_bool & finite_img & finite_sigma & pos_sigma
            sn_map = np.full_like(img, np.nan, dtype=float)
            if np.any(locs):
                sn_map[locs] = img[locs] / sigma[locs]
                snp = np.nanmean(sn_map[locs])
            else:
                warnings.warn("Invalid sn_per_pixel: no valid pixels.", AstropyUserWarning)
                try:
                    self.flag = max(getattr(self, "flag", 0), 2)
                except Exception:
                    pass
                snp = np.nan
        else:
            locs = seg_gini_bool & finite_img
            sn_map = np.full_like(img, np.nan, dtype=float)
            snp = np.nan

        fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
        ax = axes.ravel()

        finite_vals = img[np.isfinite(img)]
        if finite_vals.size > 0:
            vmin = np.nanpercentile(finite_vals, 5)
            vmax = np.nanpercentile(finite_vals, 99)
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
                vmin, vmax = np.nanmin(finite_vals), np.nanmax(finite_vals)
        else:
            vmin, vmax = 0, 1

        # 1. Science image
        im0 = ax[0].imshow(img, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        ax[0].set_title("Science stamp")
        plt.colorbar(im0, ax=ax[0], fraction=0.046)

        # 2. Segmentation map
        if seg is not None:
            im1 = ax[1].imshow(seg, origin="lower", cmap="viridis")
            ax[1].set_title("Segmentation map")
            plt.colorbar(im1, ax=ax[1], fraction=0.046)
        else:
            ax[1].text(0.5, 0.5, "No segmap", ha="center", va="center", transform=ax[1].transAxes)
            ax[1].set_title("Segmentation map")

        # 3. Gini segmap
        if seg_gini is not None:
            im2 = ax[2].imshow(seg_gini, origin="lower", cmap="viridis")
            ax[2].set_title("Gini segmap")
            plt.colorbar(im2, ax=ax[2], fraction=0.046)
        else:
            ax[2].text(0.5, 0.5, "No Gini segmap", ha="center", va="center", transform=ax[2].transAxes)
            ax[2].set_title("Gini segmap")

        # 4. Sigma image
        if sigma is not None:
            sigma_plot = np.array(sigma, dtype=float)
            sigma_plot[~np.isfinite(sigma_plot)] = np.nan
            im3 = ax[3].imshow(sigma_plot, origin="lower", cmap="magma")
            ax[3].set_title("Sigma image")
            plt.colorbar(im3, ax=ax[3], fraction=0.046)
        else:
            ax[3].text(0.5, 0.5, "No sigma image found", ha="center", va="center", transform=ax[3].transAxes)
            ax[3].set_title("Sigma image")

        # 5. Valid S/N pixels
        im4 = ax[4].imshow(locs.astype(int), origin="lower", cmap="gray_r", vmin=0, vmax=1)
        ax[4].set_title(f"Valid S/N pixels\nmean S/N per pix = {snp:.3g}" if np.isfinite(snp)
                        else "Valid S/N pixels\nmean S/N per pix = nan")
        plt.colorbar(im4, ax=ax[4], fraction=0.046)

        # 6. Overlay
        ax[5].imshow(img, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        ax[5].set_title("Science + overlays")

        if seg_gini is not None and np.any(seg_gini_bool):
            try:
                ax[5].contour(seg_gini_bool.astype(float), levels=[0.5], colors="cyan", linewidths=1.0)
            except Exception:
                pass

        if mask is not None and np.any(mask):
            try:
                ax[5].contour(mask.astype(float), levels=[0.5], colors="yellow", linewidths=1.0)
            except Exception:
                pass

        if np.any(locs):
            try:
                ax[5].contour(locs.astype(float), levels=[0.5], colors="red", linewidths=1.0)
            except Exception:
                pass

        center = self._get_center_xy_stamp()
        if center is not None:
            xc, yc = center
            ax[5].plot(xc, yc, marker="+", ms=12, mew=1.5)

        text_lines = []
        if title:
            text_lines.append(title)

        for name in ["flag", "sn_per_pixel", "gini", "m20", "concentration",
                     "asymmetry", "smoothness", "rhalf_ellip"]:
            if hasattr(self, name):
                try:
                    val = getattr(self, name)
                    if isinstance(val, (int, float, np.floating)):
                        text_lines.append(f"{name}={val:.3g}")
                    else:
                        text_lines.append(f"{name}={val}")
                except Exception:
                    pass

        if sigma is not None:
            good = np.isfinite(sigma) & (sigma > 0)
            if np.any(good):
                text_lines.append(f"sigma_med={np.nanmedian(sigma[good]):.3g}")
                text_lines.append(f"sigma_min={np.nanmin(sigma[good]):.3g}")
                text_lines.append(f"sigma_max={np.nanmax(sigma[good]):.3g}")

        text_lines.append(f"Nseg_gini={np.sum(seg_gini_bool)}")
        text_lines.append(f"Nvalid={np.sum(locs)}")

        ax[5].text(
            0.02, 0.98,
            "\n".join(text_lines),
            transform=ax[5].transAxes,
            va="top", ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.75),
        )

        for a in ax:
            a.set_xlabel("x [pix]")
            a.set_ylabel("y [pix]")

        if outfile is not None:
            fig.savefig(outfile, dpi=dpi, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig, axes, snp
            

@dataclass
class MorphologyResult:
    morph_r: Any
    morph_img2: Optional[Any] = None
    fig_r: Optional[Any] = None
    fig_img2: Optional[Any] = None
    diag_fig_r: Optional[Any] = None
    diag_fig_img2: Optional[Any] = None
    snp_r: Optional[float] = None
    snp_img2: Optional[float] = None
    
def make_object_segmap(segmentation_data: np.ndarray, object_label: int, smooth_size: int = 10) -> np.ndarray:
    """Return a 0/1 segmap for the central object only."""
    segmap = segmentation_data == object_label
    segmap_float = ndi.uniform_filter(np.float64(segmap), size=smooth_size)
    return np.array(segmap_float > 0.5, "i")


def run_statmorph_single(
    image: np.ndarray,
    segmap: np.ndarray,
    gain: float = 1.0,
    mask: Optional[np.ndarray] = None,
    weightmap=None,
    psf: Optional[np.ndarray] = None,
    cutout_extent: float = 1.5,
    make_fig: bool = True,
    make_diag: bool = False,
    diag_outfile: Optional[str] = None,
    diag_label: Optional[str] = None,
) -> tuple[Any, Optional[Any], Optional[Any]]:
    """
    Run statmorph on a single image and segmap (0/1).

    Returns
    -------
    morph, fig, diag_fig
    """
    label = 1

    morph = MyStatmorph(
        image,
        segmap,
        label,
        gain=gain,
        mask=mask,
        weightmap=weightmap,
        psf=psf,
        cutout_extent=cutout_extent,
    )

    fig = make_figure(morph) if make_fig else None
    diag_fig = None
    diag_snp = None
    if make_diag:
        diag_fig, _, diag_snp = morph.plot_diagnostic(
            outfile=diag_outfile,
            title=diag_label,
            show=False,
            )

        return morph, fig, diag_fig, diag_snp



def run_statmorph_for_photometry(
    *,
    image: np.ndarray,
    segmentation_data: np.ndarray,
    object_label: int,
    gain: float = 1.0,
    mask: Optional[np.ndarray] = None,
    weightmap=None,
    psf: Optional[np.ndarray] = None,
    image2: Optional[np.ndarray] = None,
    psf2: Optional[np.ndarray] = None,
    weightmap2=None,
    make_fig: bool = True,
    make_diag: bool = False,
    diag_prefix: Optional[str] = None,
) -> MorphologyResult:
    """
    High-level helper: build single-object segmap and run statmorph
    on image and optional image2.
    """
    segmap = make_object_segmap(segmentation_data, object_label)

    label1 = f"{diag_prefix or 'obj'} image1"
    morph_r, fig_r, diag_r, diag_snp = run_statmorph_single(
        image=image,
        segmap=segmap,
        gain=gain,
        mask=mask,
        weightmap=weightmap,
        psf=psf,
        make_fig=make_fig,
        make_diag=make_diag,
        diag_outfile=None,
        diag_label=label1,
    )

    print("STATMORPH RESULTS IMAGE1")
    morph_r.print()
    morph_r.print_diagnostic_summary()

    morph2 = fig2 = diag2 = None
    if image2 is not None:
        label2 = f"{diag_prefix or 'obj'} image2"
        morph2, fig2, diag2, diag_snp2 = run_statmorph_single(
            image=image2,
            segmap=segmap,
            gain=gain,
            mask=mask,
            weightmap=weightmap2,
            psf=psf2,
            make_fig=make_fig,
            make_diag=make_diag,
            diag_outfile=None,
            diag_label=label2,
        )
        print("STATMORPH RESULTS IMAGE2")
        morph2.print()
        morph2.print_diagnostic_summary()

    return MorphologyResult(
        morph_r=morph_r,
        morph_img2=morph2,
        fig_r=fig_r,
        fig_img2=fig2,
        diag_fig_r=diag_r,
        diag_fig_img2=diag2,
        diag_snr_r=diag_snp,
        diag_snr_img2=diag_snp2,
    )




def m20_from_sourcecatalog(cat, idx):
    """
    Compute M20 for a single object in a photutils SourceCatalog.

    Uses only pixels in the object's segment, ignoring masked pixels.
    """
    label = cat.label[idx]
    seg = (cat.segment[idx] == label)

    data_ma = cat.data_ma[idx]          # MaskedArray
    data = np.asarray(data_ma.filled(0.0), dtype=float)
    # exclude masked pixels
    good = seg & ~data_ma.mask

    if not np.any(good):
        return np.nan

    # centroid in cutout coordinates
    xc, yc = cat.cutout_centroid[idx]

    yy, xx = np.indices(data.shape)
    r2 = (xx - xc) ** 2 + (yy - yc) ** 2

    flux = data[good]
    r2_good = r2[good]

    # If flux can be negative (sky-subtracted), M20 becomes ill-defined.
    # Common approach: restrict to positive-flux pixels.
    pos = flux > 0
    if not np.any(pos):
        return np.nan

    flux = flux[pos]
    r2_good = r2_good[pos]

    mtot = np.sum(flux * r2_good)
    if mtot <= 0:
        return np.nan

    # Find the subset of pixels containing the brightest 20% of the total flux.
    order = np.argsort(flux)[::-1]  # descending
    flux_sorted = flux[order]
    r2_sorted = r2_good[order]

    cumsum = np.cumsum(flux_sorted)
    frac = cumsum / cumsum[-1]

    top20 = frac <= 0.20
    # ensure at least one pixel
    if not np.any(top20):
        top20[np.argmax(flux_sorted)] = True

    m20 = np.log10(np.sum(flux_sorted[top20] * r2_sorted[top20]) / mtot)
    return float(m20)
