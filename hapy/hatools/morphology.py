from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Any


import numpy as np

# import statmorph ONLY here
import statmorph
from statmorph.utils.image_diagnostics import make_figure

from astropy.utils import lazyproperty
import scipy.ndimage as ndi
from scipy.ndimage import shift

import warnings
import matplotlib.pyplot as plt
from astropy.utils.exceptions import AstropyUserWarning

#     @lazyproperty
#     def _segmap_gini(self):
#         segmap = np.array(self._segmap.data == 1, "i")
#         return segmap[self._slice_stamp]

from hapy.imagetools.imutils import get_bbox_from_mask, make_limits_square

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

@dataclass
class HapyMorphology:
    """
    Container for custom HAPY morphology measurements.

    Notes
    -----
    This dataclass is intended to hold both the scalar morphology results
    and the diagnostic arrays/products needed for QC and plotting.
    """

    # -----------------------------
    # Status / QC flags
    # -----------------------------
    ok: bool = False
    flag: int = 0

    # -----------------------------
    # Core morphology metrics
    # -----------------------------
    r_gini: float = np.nan
    h_gini: float = np.nan

    r_m20: float = np.nan
    h_m20: float = np.nan

    r_asym: float = np.nan
    h_asym: float = np.nan

    # -----------------------------
    # Optional associated errors
    # -----------------------------
    r_asym_err: float = np.nan
    h_asym_err: float = np.nan

    # -----------------------------
    # Second-moment bookkeeping, m20
    # -----------------------------
    r_mtot: float = np.nan
    h_mtot: float = np.nan

    r_m20sum: float = np.nan
    h_m20sum: float = np.nan

    # -----------------------------
    # Second-moment bookkeeping
    # -----------------------------
    r_mtot2: float = np.nan
    h_mtot2: float = np.nan

    r_flux_seg: float = np.nan
    h_flux_seg: float = np.nan

    r_rmom: float = np.nan
    h_rmom: float = np.nan
    

    # -----------------------------
    # Pixel bookkeeping
    # -----------------------------
    r_npix: int = 0
    h_npix: int = 0
    h_fillfrac: float = np.nan

    # -----------------------------
    # Signal-to-noise
    # -----------------------------
    r_snp_all: float = np.nan
    h_snp_det: float = np.nan
    h_snp_all: float = np.nan

    # -----------------------------
    # Threshold / sky values
    # -----------------------------
    r_sigma_sky: float = np.nan
    h_sigma_sky: float = np.nan

    r_threshold: float = np.nan
    h_threshold: float = np.nan

    # -----------------------------
    # Centers
    # -----------------------------
    xc: float = np.nan
    yc: float = np.nan

    r_asym_center: Optional[np.ndarray] = None
    h_asym_center: Optional[np.ndarray] = None

    # -----------------------------
    # Diagnostic arrays / masks
    # -----------------------------
    r_mask: Optional[np.ndarray] = None
    r_seg: Optional[Any] = None
    h_detect_mask: Optional[np.ndarray] = None
    h_gini_image: Optional[np.ndarray] = None

    # Optional full asymmetry grids for debugging
    r_asym_grid: Optional[np.ndarray] = None
    h_asym_grid: Optional[np.ndarray] = None

    # -----------------------------
    # Optional figure handles / paths
    # -----------------------------
    diag_fig: Optional[Any] = None
    diag_outfile: Optional[str] = None

    # -----------------------------
    # Optional free-form notes
    # -----------------------------
    note: str = ""

    def as_dict(self) -> dict:
        """
        Return a dictionary of scalar values for table writing.

        Arrays, masks, figures, and other large diagnostic products are excluded.
        """
        return {
            "HAPY_MORPH_OK": self.ok,
            "HAPY_MORPH_FLAG": self.flag,
            "R_HAPY_GINI": self.r_gini,
            "H_HAPY_GINI": self.h_gini,
            "R_HAPY_M20": self.r_m20,
            "H_HAPY_M20": self.h_m20,
            "R_HAPY_ASYM": self.r_asym,
            "H_HAPY_ASYM": self.h_asym,
            "R_HAPY_ASYM_ERR": self.r_asym_err,
            "H_HAPY_ASYM_ERR": self.h_asym_err,
            "R_HAPY_MTOT": self.r_mtot,
            "H_HAPY_MTOT": self.h_mtot,
            "R_HAPY_M20SUM": self.r_m20sum,
            "H_HAPY_M20SUM": self.h_m20sum,
            "R_HAPY_NPIX": self.r_npix,
            "H_HAPY_NPIX": self.h_npix,
            "H_HAPY_FILLFRAC": self.h_fillfrac,
            "R_HAPY_SNP_ALL": self.r_snp_all,
            "H_HAPY_SNP_DET": self.h_snp_det,
            "H_HAPY_SNP_ALL": self.h_snp_all,
            "R_HAPY_SIGMA_SKY": self.r_sigma_sky,
            "H_HAPY_SIGMA_SKY": self.h_sigma_sky,
            "R_HAPY_THRESHOLD": self.r_threshold,
            "H_HAPY_THRESHOLD": self.h_threshold,
            "R_HAPY_XC": self.xc,
            "R_HAPY_YC": self.yc,
        }
    
class MyStatmorph(statmorph.SourceMorphology):
    """Overiding some of statmorph's behavior """

    @lazyproperty
    def sn_per_pixel(self):
        sigma = np.asarray(self._weightmap_stamp)
        img = np.asarray(self._cutout_stamp_maskzeroed)
        seg = np.asarray(self._segmap_gini).astype(bool)

        #print("DEBUG sn_per_pixel")
        #print("  img shape:", img.shape)
        #print("  sigma shape:", sigma.shape)
        #print("  seg npix:", np.sum(seg))

        locs_raw = seg & (img >= 0) & (sigma > 0)
        #print("  raw locs:", np.sum(locs_raw))

        if np.any(locs_raw):
            raw_ratio = img[locs_raw] / sigma[locs_raw]
            #print("  raw finite img:", np.all(np.isfinite(img[locs_raw])))
            #print("  raw finite sigma:", np.all(np.isfinite(sigma[locs_raw])))
            #print("  raw finite ratio:", np.all(np.isfinite(raw_ratio)))
            #print("  raw n bad ratio:", np.sum(~np.isfinite(raw_ratio)))

        locs = seg & np.isfinite(img) & np.isfinite(sigma) & (img >= 0) & (sigma > 0)
        #print("  safe locs:", np.sum(locs))

        if np.sum(locs) == 0:
            self.flag = 2
            return -99.0

        ratio = img[locs] / sigma[locs]
        #print("  safe finite ratio:", np.all(np.isfinite(ratio)))

        return float(np.mean(ratio))
    
    def print_diagnostic_summary(self):
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

            locs2 = (
                seg_gini.astype(bool)
                & (self._cutout_stamp_maskzeroed >= 0)
                & (sigma > 0)
            )

            locs = (
                seg_gini.astype(bool)
                & np.isfinite(self._cutout_stamp_maskzeroed)
                & np.isfinite(sigma)
                & (sigma > 0)
            )

            #print("Npix locs2 (statmorph):", np.sum(locs2))
            #print("Npix locs  (safe)     :", np.sum(locs))

            if np.any(locs2):
                vals_img2 = self._cutout_stamp_maskzeroed[locs2]
                vals_sigma2 = sigma[locs2]
                ratio2 = vals_img2 / vals_sigma2

                #print("locs2: finite img   :", np.all(np.isfinite(vals_img2)))
                #print("locs2: finite sigma :", np.all(np.isfinite(vals_sigma2)))
                #print("locs2: finite ratio :", np.all(np.isfinite(ratio2)))
                #print("locs2: n bad ratio  :", np.sum(~np.isfinite(ratio2)))

                if not np.all(np.isfinite(vals_img2)):
                    print("WARNING: non-finite values in img[locs2]")
                if not np.all(np.isfinite(vals_sigma2)):
                    print("WARNING: non-finite values in sigma[locs2]")

                snp2 = np.mean(ratio2)
            else:
                snp2 = np.nan

            sn_map = np.full_like(img, np.nan, dtype=float)
            if np.any(locs):
                vals_img = img[locs]
                vals_sigma = sigma[locs]
                ratio = vals_img / vals_sigma
                sn_map[locs] = ratio
                snp = np.nanmean(ratio)
            else:
                warnings.warn("Invalid sn_per_pixel: no valid pixels.", AstropyUserWarning)
                try:
                    self.flag = max(getattr(self, "flag", 0), 2)
                except Exception:
                    pass
                snp = np.nan

            print("snp statmorph-style:", snp2)
            print("snp safe version   :", snp)

    
        # if sigma is not None:
        #     finite_sigma = np.isfinite(sigma)
        #     pos_sigma = sigma > 0
        #     # statmorph criteria
        #     locs2 = (seg_gini & (self._cutout_stamp_maskzeroed >= 0) & (sigma > 0))
        #     # chatgpt version
        #     locs = seg_gini_bool & finite_img & finite_sigma & pos_sigma
        #     sn_map = np.full_like(img, np.nan, dtype=float)
        #     if np.any(locs):
        #         # check img[locs] for nans
        #         if not np.all(np.isfinite(img[locs2])):
        #             print(f"WARNING: found nans in img[locs]")
        #         if not np.all(np.isfinite(sigma[locs2])):
        #             print(f"WARNING: found nans in sigma[locs]")
        #         sn_map[locs] = img[locs] / sigma[locs]
        #         snp = np.nanmean(sn_map[locs])
        #     else:
        #         warnings.warn("Invalid sn_per_pixel: no valid pixels.", AstropyUserWarning)
        #         try:
        #             self.flag = max(getattr(self, "flag", 0), 2)
        #         except Exception:
        #             pass
        #         snp = np.nan
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


        locs2 = (seg_gini.astype(bool) &
                 (self._cutout_stamp_maskzeroed >= 0) &
                 (sigma > 0))

        print("npix locs2 =", np.sum(locs2))

        if np.any(locs2):
            bad_img = locs2 & ~np.isfinite(self._cutout_stamp_maskzeroed)
            bad_sigma = locs2 & ~np.isfinite(sigma)
            zero_sigma = locs2 & (sigma == 0)

            #print("bad img in locs2:", np.sum(bad_img))
            #print("bad sigma in locs2:", np.sum(bad_sigma))
            #print("zero sigma in locs2:", np.sum(zero_sigma))

            vals_img = self._cutout_stamp_maskzeroed[locs2]
            vals_sigma = sigma[locs2]

            #print("img finite?", np.all(np.isfinite(vals_img)))
            #print("sigma finite?", np.all(np.isfinite(vals_sigma)))
            #print("sigma min/max:", np.nanmin(vals_sigma), np.nanmax(vals_sigma))

            ratio = vals_img / vals_sigma
            #print("ratio finite?", np.all(np.isfinite(ratio)))
            #print("n bad ratio:", np.sum(~np.isfinite(ratio)))

            only_in_locs2 = locs2 & ~locs
            only_in_locs = locs & ~locs2

            #print("npix locs2:", np.sum(locs2))
            #print("npix locs :", np.sum(locs))
            #print("npix only_in_locs2:", np.sum(only_in_locs2))
            #print("npix only_in_locs :", np.sum(only_in_locs))

            if np.any(only_in_locs2):
                print("nonfinite img among only_in_locs2:", np.sum(~np.isfinite(img[only_in_locs2])))
                print("nonfinite sigma among only_in_locs2:",np.sum(~np.isfinite(sigma[only_in_locs2])))

            print(type(self._cutout_stamp_maskzeroed))
            print(type(sigma))

            #print("img is masked array:", np.ma.isMaskedArray(self._cutout_stamp_maskzeroed))
            #print("sigma is masked array:", np.ma.isMaskedArray(sigma))
            #print()
            #print()
            return fig, axes, snp
            




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
    diag_outfile: Optional[str] = None,    
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
        diag_outfile=diag_outfile,
        diag_label=label1,
    )

    #print("STATMORPH RESULTS IMAGE1")
    #morph_r.print()
    #morph_r.print_diagnostic_summary()

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
            diag_outfile=diag_outfile.replace('.pdf','2.pdf'),
            diag_label=label2,
        )
        #print("STATMORPH RESULTS IMAGE2")
        #morph2.print()
        #morph2.print_diagnostic_summary()


    return MorphologyResult(
        morph_r=morph_r,
        morph_img2=morph2,
        fig_r=fig_r,
        fig_img2=fig2,
        diag_fig_r=diag_r,
        diag_fig_img2=diag2,
        snp_r=diag_snp,
        snp_img2=diag_snp2,)




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




def compute_gini(values, allow_negative=False):
    """
    Compute the Gini coefficient for a 1D array of pixel values.

    Parameters
    ----------
    values : array-like
        Input pixel values.
    allow_negative : bool, optional
        If False, negative values are clipped to zero before computing Gini.
        If True, negative values are kept. In most astronomical imaging
        applications for flux-like quantities, False is the safer choice.

    Returns
    -------
    gini : float
        Gini coefficient in the range [0, 1] for non-negative inputs.
        Returns np.nan if no finite values remain or if the sum is zero.

    Notes
    -----
    For non-negative values, the Gini coefficient is:

        G = sum_i (2i - n - 1) x_i / [n * sum_i x_i]

    where x_i are the sorted values in ascending order.
    """
    x = np.asarray(values, dtype=float).ravel()
    x = x[np.isfinite(x)]

    if x.size == 0:
        return np.nan

    if not allow_negative:
        x = np.clip(x, 0.0, None)

    if x.size == 0:
        return np.nan

    total = np.sum(x)
    if total <= 0:
        return np.nan

    x = np.sort(x)
    n = x.size
    index = np.arange(1, n + 1, dtype=float)

    gini = np.sum((2.0 * index - n - 1.0) * x) / (n * total)
    return float(gini)

def compute_hapy_gini(
    image,
    r_segmap,
    sigma_sky=None,
    nsigma=3.0,
    detection_segmap=None,
    zero_below_threshold=False,
):
    """
    Compute HAPY Gini over the r-band segmentation region.

    Parameters
    ----------
    image : 2D array
        Image on which to compute the Gini coefficient.
        This can be the r-band or Halpha image.
    r_segmap : 2D array-like
        Boolean or integer segmentation map defining the stellar disk region.
        Nonzero values are treated as inside the region.
    sigma_sky : float, optional
        Sky RMS for thresholding. Required if zero_below_threshold=True.
    nsigma : float, optional
        Threshold multiplier for sigma_sky. Default is 3.
    detection_segmap : 2D array-like, optional
        Optional second segmentation map describing detected emission.
        If provided, pixels inside r_segmap but outside detection_segmap
        are set to zero before computing Gini.
    zero_below_threshold : bool, optional
        If True, pixels inside r_segmap with image < nsigma*sigma_sky
        are set to zero before computing Gini.

    Returns
    -------
    gini : float
        Gini coefficient computed over the r-band segmentation region.

    processed_values : 1D ndarray
        The values actually used in the Gini calculation.
    """
    img = np.asarray(image, dtype=float)
    rmask = np.asarray(r_segmap).astype(bool)

    if img.shape != rmask.shape:
        raise ValueError("image and r_segmap must have the same shape")

    vals = np.array(img[rmask], dtype=float)

    if detection_segmap is not None:
        dmask = np.asarray(detection_segmap).astype(bool)
        if dmask.shape != img.shape:
            raise ValueError("detection_segmap must have the same shape as image")
        detvals = dmask[rmask]
        vals[~detvals] = 0.0

    if zero_below_threshold:
        if sigma_sky is None:
            raise ValueError("sigma_sky is required when zero_below_threshold=True")
        threshold = nsigma * sigma_sky
        vals[vals < threshold] = 0.0

    gini = compute_gini(vals, allow_negative=False)
    return gini, vals


def compute_r_hapy_gini(r_image, r_segmap):
    """
    R_HAPY_GINI:
    Gini of r-band image over the r-band segmentation region.
    """
    gini, _ = compute_hapy_gini(
        image=r_image,
        r_segmap=r_segmap,
        zero_below_threshold=False,
    )
    return gini


def compute_halpha_hapy_gini(ha_image, r_segmap, sigma_sky, nsigma=3.0):
    """
    H_HAPY_GINI:
    Gini of Halpha image over the r-band segmentation region,
    after setting pixels below nsigma*sigma_sky to zero.
    """
    gini, _ = compute_hapy_gini(
        image=ha_image,
        r_segmap=r_segmap,
        sigma_sky=sigma_sky,
        nsigma=nsigma,
        zero_below_threshold=True,
    )
    return gini

def compute_halpha_hapy_gini_with_segmaps(ha_image, r_segmap, ha_det_segmap):
    """
    H_HAPY_GINI using an explicit Halpha detection map.

    Pixels inside r_segmap but outside ha_det_segmap are set to zero.
    """
    gini, _ = compute_hapy_gini(
        image=ha_image,
        r_segmap=r_segmap,
        detection_segmap=ha_det_segmap,
        zero_below_threshold=False,
    )
    return gini



def compute_flux_centroid(image, segmap):
    img = np.asarray(image, dtype=float)
    mask = np.asarray(segmap).astype(bool)

    vals = np.array(img, copy=True)
    vals[~np.isfinite(vals)] = 0.0
    vals[vals < 0] = 0.0
    vals[~mask] = 0.0

    total = np.sum(vals[mask])
    if total <= 0:
        return np.nan, np.nan

    y, x = np.indices(vals.shape)
    xc = np.sum(vals[mask] * x[mask]) / total
    yc = np.sum(vals[mask] * y[mask]) / total
    return float(xc), float(yc)


def decode_hapy_flag(flag: int) -> list[str]:
    meanings = {
        1: "invalid or empty r-band morphology mask",
        2: "no Halpha image available",
        4: "no Halpha pixels above threshold",
        8: "invalid/non-finite core morphology metric",
        16: "exception during morphology calculation",
    }

    out = [text for bit, text in meanings.items() if flag & bit]

    return out if out else ["OK"]




def compute_m20(image, segmap, xc=None, yc=None):
    """
    Compute M20 on a specified image and segmentation mask.

    Parameters
    ----------
    image : 2D ndarray
        Image values used for the calculation.
    segmap : 2D ndarray or bool array
        Boolean or integer mask defining the pixels to include.
        Nonzero values are treated as inside the object.
    xc, yc : float, optional
        Center coordinates in full-image pixel coordinates.
        If not provided, the flux-weighted centroid is computed from the
        positive flux inside the segmap.

    Returns
    -------
    m20 : float
        M20 value.
    second_moment_tot : float
        Total second-order moment.
    second_moment_20 : float
        Second-order moment of brightest 20% of the flux.

    Notes
    -----
    M20 is defined as:

        M20 = log10(sum(M_i brightest 20%) / M_tot)

    where M_i = f_i * [(x_i - x_c)^2 + (y_i - y_c)^2]

    This implementation uses only pixels inside `segmap`.
    Negative values are clipped to zero.
    """
    img = np.asarray(image, dtype=float)
    mask = np.asarray(segmap).astype(bool)

    if img.shape != mask.shape:
        raise ValueError("image and segmap must have the same shape")

    vals = np.array(img, copy=True)
    vals[~np.isfinite(vals)] = 0.0
    vals[vals < 0] = 0.0
    vals[~mask] = 0.0

    if np.sum(mask) == 0:
        return np.nan, np.nan, np.nan

    total_flux = np.sum(vals[mask])
    if total_flux <= 0:
        return np.nan, np.nan, np.nan

    y, x = np.indices(vals.shape)

    # Compute centroid from masked positive flux if not supplied
    if xc is None or yc is None:
        flux = vals[mask]
        xmask = x[mask]
        ymask = y[mask]
        xc = np.sum(flux * xmask) / np.sum(flux)
        yc = np.sum(flux * ymask) / np.sum(flux)

    distsq = (x - xc) ** 2 + (y - yc) ** 2
    mi = vals * distsq

    second_moment_tot = np.sum(mi[mask])
    if second_moment_tot <= 0:
        return np.nan, np.nan, np.nan

    # Find pixels containing brightest 20% of total flux
    pixvals = vals[mask]
    pixmom = mi[mask]

    order = np.argsort(pixvals)[::-1]   # descending
    pixvals_sorted = pixvals[order]
    pixmom_sorted = pixmom[order]

    cumsum_flux = np.cumsum(pixvals_sorted)
    flux_limit = 0.2 * np.sum(pixvals_sorted)

    use = cumsum_flux <= flux_limit

    # Ensure at least one pixel is included
    if not np.any(use):
        use[0] = True
    else:
        # Include the first pixel that crosses the threshold
        first_over = np.argmax(cumsum_flux >= flux_limit)
        use[first_over] = True

    second_moment_20 = np.sum(pixmom_sorted[use])

    m20 = np.log10(second_moment_20 / second_moment_tot)
    return float(m20), float(second_moment_tot), float(second_moment_20)

def compute_second_moment_metrics(image, segmap, xc=None, yc=None):
    """
    Compute flux, total second moment, and second-moment radius
    on a specified image and segmentation mask.

    Uses the same conventions as compute_m20:
      - nonzero segmap pixels define the object
      - non-finite values are set to zero
      - negative values are clipped to zero
      - if xc/yc are not supplied, use the flux-weighted centroid

    Returns
    -------
    total_flux : float
        Total positive flux inside segmap.
    second_moment_tot : float
        Total second-order moment.
    rmom : float
        Flux-normalized second-moment radius:
            sqrt(second_moment_tot / total_flux)
    xc : float
        Center used for calculation.
    yc : float
        Center used for calculation.
    """
    img = np.asarray(image, dtype=float)
    mask = np.asarray(segmap).astype(bool)

    if img.shape != mask.shape:
        raise ValueError("image and segmap must have the same shape")

    vals = np.array(img, copy=True)
    vals[~np.isfinite(vals)] = 0.0
    vals[vals < 0] = 0.0
    vals[~mask] = 0.0

    if np.sum(mask) == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    total_flux = np.sum(vals[mask])
    if total_flux <= 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    y, x = np.indices(vals.shape)

    if xc is None or yc is None:
        flux = vals[mask]
        xmask = x[mask]
        ymask = y[mask]
        xc = np.sum(flux * xmask) / total_flux
        yc = np.sum(flux * ymask) / total_flux

    distsq = (x - xc) ** 2 + (y - yc) ** 2
    mi = vals * distsq

    second_moment_tot = np.sum(mi[mask])
    if second_moment_tot <= 0:
        return float(total_flux), np.nan, np.nan, float(xc), float(yc)

    rmom = np.sqrt(second_moment_tot / total_flux)

    return (
        float(total_flux),
        float(second_moment_tot),
        float(rmom),
        float(xc),
        float(yc),
    )


def quick_asym_test(image, mask, xc, yc):
    ypix, xpix = np.where(mask)
    x_ref = np.rint(2 * xc - xpix).astype(int)
    y_ref = np.rint(2 * yc - ypix).astype(int)

    ny, nx = image.shape
    good = (
        (y_ref >= 0) & (y_ref < ny) &
        (x_ref >= 0) & (x_ref < nx)
    )

    orig = image[ypix[good], xpix[good]]
    refl = image[y_ref[good], x_ref[good]]

    num = np.sum(np.abs(orig - refl))
    den = np.sum(np.abs(image[mask]))
    return num / den if den > 0 else np.nan

def compute_asymmetry(image, segmap, xc=None, yc=None, search_radius=1, step=1.0):
    """
    Compute 180-degree rotational asymmetry on a specified image and mask.

    Parameters
    ----------
    image : 2D ndarray
        Image used for the asymmetry calculation.
    segmap : 2D ndarray or bool array
        Boolean or integer mask defining the pixels to include.
        Nonzero values are treated as inside the object.
    xc, yc : float, optional
        Center coordinates in full-image pixel coordinates.
        If not provided, the flux-weighted centroid is computed from the
        positive flux inside the mask.
    search_radius : int, optional
        Radius in pixels around the input center over which to search
        for the minimum asymmetry. Default is 1, giving a 3x3 grid.
    step : float, optional
        Step size in pixels for the center search. Default is 1.0.

    Returns
    -------
    asym : float
        Minimum asymmetry value.
    asym_err : float
        Standard deviation of asymmetry values over the search grid.
    best_center : ndarray, shape (2,)
        Best-fit center as [yc, xc].
    asym_grid : 2D ndarray
        Grid of asymmetry values evaluated around the center.

    Notes
    -----
    Asymmetry is calculated as:

        A = sum(|I - I_180|) / sum(|I|)

    where the sums are performed only over the supplied mask.

    This implementation rotates the image by 180 degrees about each trial
    center using interpolation via scipy.ndimage.shift.
    """
    from scipy.ndimage import map_coordinates

    img = np.asarray(image, dtype=float)
    mask = np.asarray(segmap).astype(bool)


    if img.shape != mask.shape:
        raise ValueError("image and segmap must have the same shape")


    vals = np.array(img, copy=True)
    vals[~np.isfinite(vals)] = 0.0
    #vals[~mask] = 0.0

    if np.sum(mask) == 0:
        return np.nan, np.nan, np.array([np.nan, np.nan]), np.full((1, 1), np.nan)

    # Compute default center from positive flux if not supplied
    if xc is None or yc is None:
        pos = np.array(vals, copy=True)
        pos[pos < 0] = 0.0
        total = np.sum(pos[mask])

        if total <= 0:
            # fallback to geometric center of the mask
            yy, xx = np.where(mask)
            yc = np.mean(yy)
            xc = np.mean(xx)
        else:
            y, x = np.indices(vals.shape)
            yc = np.sum(pos[mask] * y[mask]) / total
            xc = np.sum(pos[mask] * x[mask]) / total

    print("DEBUG: quick_asym = ",quick_asym_test(image, mask, xc, yc))
    
    offsets = np.arange(-search_radius, search_radius + step, step, dtype=float)
    asym_grid = np.full((len(offsets), len(offsets)), np.nan, dtype=float)

    vals_pos = np.clip(vals, 0, None)
    norm = np.sum(np.abs(vals[mask]))
    if not np.isfinite(norm) or norm <= 0:
        return np.nan, np.nan, np.array([yc, xc]), asym_grid

    for iy, dy in enumerate(offsets):
        for ix, dx in enumerate(offsets):
            yc_trial = yc + dy
            xc_trial = xc + dx

            # Shift so that trial center is at origin, rotate by 180 via flip,
            # then shift back.
            #shifted = shift(vals, shift=(-yc_trial, -xc_trial), order=1, mode="constant", cval=0.0)
            #rotated = np.flipud(np.fliplr(shifted))
            #rotated_back = shift(rotated, shift=(yc_trial, xc_trial), order=1, mode="constant", cval=0.0)
            y, x = np.indices(vals.shape)
            x_rot = 2.0 * xc_trial - x
            y_rot = 2.0 * yc_trial - y

            rotated_back = map_coordinates(
                vals,
                [y_rot, x_rot],
                order=1,
                mode="constant",
                cval=0.0
                )
            #diff = np.abs(vals - rotated_back)
            #print(f"DEBUG: xc_trial={xc_trial}, yc_trial={yc_trial},diff[0]={diff[0]}, diff[1]={diff[1]}")
            #asym_grid[iy, ix] = np.sum(diff[mask]) / norm
            rotated_mask = map_coordinates(
                mask.astype(float),
                [y_rot, x_rot],
                order=0,
                mode="constant",
                cval=0.0
                ) > 0.5

            valid = mask & rotated_mask & np.isfinite(rotated_back)

            if np.any(valid):
                asym_grid[iy, ix] = np.sum(np.abs(vals[valid] - rotated_back[valid])) / np.sum(np.abs(vals[valid]))
            else:
                asym_grid[iy, ix] = np.nan

    if not np.any(np.isfinite(asym_grid)):
        return np.nan, np.nan, np.array([yc, xc]), asym_grid

    iy_min, ix_min = np.unravel_index(np.nanargmin(asym_grid), asym_grid.shape)
    best_yc = yc + offsets[iy_min]
    best_xc = xc + offsets[ix_min]

    asym = float(asym_grid[iy_min, ix_min])
    asym_err = float(np.nanstd(asym_grid))
    best_center = np.array([best_yc, best_xc], dtype=float)
    print(f"DEBUG in compute_asymmetry, asym={asym}")
    return asym, asym_err, best_center, asym_grid




def plot_hapy_morphology_diagnostic(
    r_image,
    ha_image,
    morph,
    outfile=None,
    show=False,
    dpi=150,
):
    """
    Make a 2x4 diagnostic plot for HAPY morphology measurements.

    Parameters
    ----------
    r_image : 2D ndarray
        R-band image.
    ha_image : 2D ndarray
        Halpha image.
    morph : HapyMorphology
        Dataclass containing masks, thresholded Halpha image, and scalar metrics.
    outfile : str, optional
        Output filename.
    show : bool, optional
        If True, call plt.show(); otherwise close the figure after saving.
    dpi : int, optional
        Save resolution.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import numpy as np
    import matplotlib.pyplot as plt

    r_image = np.asarray(r_image, dtype=float)
    ha_image = np.asarray(ha_image, dtype=float)

    r_gini_mask = np.asarray(morph.r_mask).astype(bool) if morph.r_mask is not None else np.zeros_like(r_image, dtype=bool)
    ha_detect_mask = np.asarray(morph.h_detect_mask).astype(bool) if morph.h_detect_mask is not None else np.zeros_like(r_image, dtype=bool)
    ha_gini_image = np.asarray(morph.h_gini_image, dtype=float) if morph.h_gini_image is not None else np.full_like(ha_image, np.nan, dtype=float)

    if r_image.shape != ha_image.shape:
        raise ValueError("r_image and ha_image must have the same shape")
    if r_gini_mask.shape != r_image.shape:
        raise ValueError("morph.r_mask must have same shape as images")
    if ha_detect_mask.shape != r_image.shape:
        raise ValueError("morph.h_detect_mask must have same shape as images")
    if ha_gini_image.shape != r_image.shape:
        raise ValueError("morph.h_gini_image must have same shape as images")

    r_vals = r_image[r_gini_mask]
    r_vals = r_vals[np.isfinite(r_vals)]

    ha_vals = ha_gini_image[r_gini_mask]
    ha_vals = ha_vals[np.isfinite(ha_vals)]

    def _get_vrange(arr, positive_only=False, symmetric=False):
        vals = arr[np.isfinite(arr)]
        if vals.size == 0:
            return 0.0, 1.0

        if positive_only:
            vals = vals[vals > 0]
            if vals.size == 0:
                vals = arr[np.isfinite(arr)]

        if symmetric:
            vmax = np.nanpercentile(np.abs(vals), 99.5)
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = np.nanmax(np.abs(vals))
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = 1.0
            return -vmax, vmax

        vmin = np.nanpercentile(vals, 5)
        vmax = np.nanpercentile(vals, 99.5)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin = np.nanmin(vals)
            vmax = np.nanmax(vals)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin, vmax = 0.0, 1.0
        return vmin, vmax

    r_vmin, r_vmax = _get_vrange(r_image)
    ha_vmin, ha_vmax = _get_vrange(ha_image, symmetric=True)

    hagi_vals = ha_gini_image[np.isfinite(ha_gini_image)]
    hagi_pos = hagi_vals[hagi_vals > 0]
    if hagi_pos.size > 0:
        hagi_vmin = 0.0
        hagi_vmax = np.nanpercentile(hagi_pos, 99)
        if not np.isfinite(hagi_vmax) or hagi_vmax <= 0:
            hagi_vmax = np.nanmax(hagi_pos)
        if not np.isfinite(hagi_vmax) or hagi_vmax <= 0:
            hagi_vmax = 1.0
    else:
        hagi_vmin, hagi_vmax = 0.0, 1.0

    xlim, ylim = get_bbox_from_mask(r_gini_mask, pad=75, min_size=75)
    xlim, ylim = make_limits_square(xlim, ylim, image_shape=r_gini_mask.shape)

    fig, axes = plt.subplots(2, 4, figsize=(18, 9), constrained_layout=True)
    ax = axes.ravel()

    # 1. r-band image
    im0 = ax[0].imshow(r_image, origin="lower", cmap="gray_r", vmin=r_vmin, vmax=r_vmax)
    ax[0].set_title("r-band image")
    plt.colorbar(im0, ax=ax[0], fraction=0.046)

    # 2. r-band morphology mask
    im1 = ax[1].imshow(r_gini_mask.astype(int), origin="lower", cmap="gray_r", vmin=0, vmax=1)
    ax[1].set_title(f"r morphology mask\nNpix = {morph.r_npix}")
    plt.colorbar(im1, ax=ax[1], fraction=0.046)

    # 3. r-band + mask overlay
    ax[2].imshow(r_image, origin="lower", cmap="gray_r", vmin=r_vmin, vmax=r_vmax)
    try:
        ax[2].contour(r_gini_mask.astype(float), levels=[0.5], colors="cyan", linewidths=1.0)
    except Exception:
        pass
    ax[2].set_title(
        f"r + mask\n"
        f"G={morph.r_gini:.3f}, M20={morph.r_m20:.3f}, A={morph.r_asym:.3f}"
    )

    # 4. histogram of r pixels used in morphology
    if r_vals.size > 0:
        ax[3].hist(r_vals, bins=50)
    ax[3].set_title("r-band morphology pixels")
    ax[3].set_xlabel("r-band pixel value")
    ax[3].set_ylabel("N")
    ax[3].set_yscale("log")    
    # 5. Halpha image
    im4 = ax[4].imshow(ha_image, origin="lower", cmap="gray_r", vmin=ha_vmin, vmax=ha_vmax)
    ax[4].set_title("Hα image")
    plt.colorbar(im4, ax=ax[4], fraction=0.046)

    # 6. thresholded Halpha image used for morphology
    hagi_plot = np.full_like(ha_gini_image, np.nan, dtype=float)
    hagi_plot[r_gini_mask] = ha_gini_image[r_gini_mask]
    im5 = ax[5].imshow(hagi_plot, origin="lower", cmap="viridis", vmin=hagi_vmin, vmax=hagi_vmax)
    if np.isfinite(morph.h_threshold) and np.isfinite(morph.h_sigma_sky):
        ax[5].set_title(
            f"Hα used for morphology\n"
            f"threshold={morph.h_threshold:.3g} ({morph.h_threshold/morph.h_sigma_sky:.1f}σ)"
        )
    elif np.isfinite(morph.h_threshold):
        ax[5].set_title(f"Hα used for morphology\nthreshold={morph.h_threshold:.3g}")
    else:
        ax[5].set_title("Hα used for morphology")
    try:
        ax[5].contour(ha_detect_mask.astype(float), levels=[0.5], colors="red", linewidths=1.0)
    except Exception:
        pass
    plt.colorbar(im5, ax=ax[5], fraction=0.046)

    # 7. Halpha + overlays
    ax[6].imshow(ha_image, origin="lower", cmap="gray_r", vmin=ha_vmin, vmax=ha_vmax)
    try:
        ax[6].contour(r_gini_mask.astype(float), levels=[0.5], colors="cyan", linewidths=1.0, alpha=0.6)
    except Exception:
        pass
    try:
        ax[6].contour(ha_detect_mask.astype(float), levels=[0.5], colors="red", linewidths=1.0)
    except Exception:
        pass

    title7 = (
        f"Hα + overlays\n"
        f"G={morph.h_gini:.3f}, M20={morph.h_m20:.3f}, A={morph.h_asym:.3f}\n"
        f"fill={morph.h_fillfrac:.3f}"
    )
    if np.isfinite(morph.h_snp_det):
        title7 += f", S/N(det)={morph.h_snp_det:.2f}"
    if np.isfinite(morph.h_snp_all):
        title7 += f"\nS/N(all)={morph.h_snp_all:.2f}"
    ax[6].set_title(title7)

    # 8. histogram of Halpha pixels used in morphology
    if ha_vals.size > 0:
        ax[7].hist(ha_vals, bins=50)
    ax[7].set_title("Hα morphology pixels")
    ax[7].set_xlabel("Hα pixel value (thresholded)")
    ax[7].set_ylabel("N")
    ax[7].set_yscale("log")    

    rtext_lines = []
    htext_lines = []    
    if morph.note:
        rtext_lines.append(morph.note)


    # Zoom image panels to the r-band segmentation bounding box.
    # Leave histogram panels unzoomed.
    image_axes = [ax[0], ax[1], ax[2], ax[4], ax[5], ax[6]]

    for a in image_axes:
        a.set_xlim(*xlim)
        a.set_ylim(*ylim)
        a.set_xlabel("x [pix]")
        a.set_ylabel("y [pix]")

    for a in [ax[3], ax[7]]:
        a.set_xlabel(a.get_xlabel())
        a.set_ylabel(a.get_ylabel())
    
    rtext_lines.extend([
        f"HAPY_MORPH_OK = {morph.ok}",
        f"HAPY_MORPH_FLAG = {morph.flag}",
        f"R_GINI = {morph.r_gini:.3f}",
        f"R_M20 = {morph.r_m20:.3f}",
        f"R_ASYM = {morph.r_asym:.3f}",
        f"R_NPIX = {morph.r_npix}",
    ])
    htext_lines.extend([
        f"H_GINI = {morph.h_gini:.3f}",
        f"H_M20 = {morph.h_m20:.3f}",
        f"H_ASYM = {morph.h_asym:.3f}",
        f"H_NPIX = {morph.h_npix}",
        f"H_FILLFRAC = {morph.h_fillfrac:.3f}",
    ])

    if np.isfinite(morph.r_snp_all):
        rtext_lines.append(f"R_SNP_ALL = {morph.r_snp_all:.3f}")
    if np.isfinite(morph.h_snp_det):
        htext_lines.append(f"H_SNP_DET = {morph.h_snp_det:.3f}")
    if np.isfinite(morph.h_snp_all):
        htext_lines.append(f"H_SNP_ALL = {morph.h_snp_all:.3f}")
    if np.isfinite(morph.h_sigma_sky):
        htext_lines.append(f"H sigma_sky = {morph.h_sigma_sky:.3g}")
    if np.isfinite(morph.h_threshold):
        htext_lines.append(f"H threshold = {morph.h_threshold:.3g}")



    ax[3].text(
        0.98, 0.98,
        "\n".join(rtext_lines),
        transform=ax[3].transAxes,
        ha="right", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    ax[7].text(
        0.98, 0.98,
        "\n".join(htext_lines),
        transform=ax[7].transAxes,
        ha="right", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )    

    for a in ax[:7]:
        a.set_xlabel("x [pix]")
        a.set_ylabel("y [pix]")

    if outfile is not None:
        fig.savefig(outfile, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig



