# hapy/io/schemas.py

import astropy.units as u

"""
Schema for per-aperture photometry profile tables written by HAPY.

Each row corresponds to one elliptical aperture, ordered by increasing
semi-major axis. Flux-like quantities are cumulative within the aperture.
Surface-brightness quantities are averaged over the unmasked aperture area.
"""
import astropy.units as u

PHOT_TABLE_SCHEMA = [
    ("ap_index", u.dimensionless_unscaled, "Integer aperture index, increasing with semi-major axis"),
    ("sma_arcsec", u.arcsec, "Elliptical aperture semi-major axis"),
    ("sma_pix", u.pixel, "Elliptical aperture semi-major axis in pixels"),
    ("area_total_pix", u.pixel**2, "Total number of pixels in the aperture"),
    ("area_unmasked_pix", u.pixel**2, "Number of unmasked pixels in the aperture"),
    ("masked_fraction", u.dimensionless_unscaled, "Fraction of aperture area masked"),
    ("flux_cum", u.adu / u.s, "Cumulative flux within the aperture"),
    ("flux_cum_err", u.adu / u.s, "Uncertainty in cumulative aperture flux"),
    ("sb_avg", u.adu / u.s / u.pixel**2, "Mean surface brightness in the annulus over unmasked area"),
    ("sb_avg_err", u.adu / u.s / u.pixel**2, "Uncertainty in mean annulus surface brightness"),
    ("sb_avg_snr", u.dimensionless_unscaled, "Signal-to-noise ratio of mean annulus surface brightness"),
    ("flux_cgs", u.erg / u.s / u.cm**2, "Cumulative calibrated flux within the aperture"),
    ("flux_cgs_err", u.erg / u.s / u.cm**2, "Uncertainty in cumulative calibrated flux"),
    ("mag_cum", u.mag, "Cumulative magnitude within the aperture"),
    ("mag_cum_err", u.mag, "Uncertainty in cumulative aperture magnitude"),
    ("sb_cgs_arcsec2", u.erg / u.s / u.cm**2 / u.arcsec**2, "Mean calibrated surface brightness per square arcsecond"),
    ("sb_cgs_arcsec2_err", u.erg / u.s / u.cm**2 / u.arcsec**2, "Uncertainty in calibrated surface brightness per square arcsecond"),
    ("sb_mag_arcsec2", u.mag / u.arcsec**2, "Surface brightness in magnitudes per square arcsecond"),
    ("sb_mag_arcsec2_err", u.mag / u.arcsec**2, "Uncertainty in surface brightness in magnitudes per square arcsecond"),
    ("snr_total", u.dimensionless_unscaled, "Signal-to-noise ratio of cumulative aperture flux"),
    ("snr_per_pixel", u.dimensionless_unscaled, "Signal-to-noise ratio per unmasked pixel"),
    ("snr_image_units", u.dimensionless_unscaled, "Surface-brightness-like signal-to-noise ratio in image units"),
]
