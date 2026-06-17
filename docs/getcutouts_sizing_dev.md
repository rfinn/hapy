# Notes

Working on an algorithm for setting the cutout size that includes a minimum size for the smallest galaxies, and a size + buffer for the larger galaxies.


# Testing

Testing directory:
```bash
/data-pool/Halpha/tests/getcutouts_size_buffer_test
```


### Test Images

```bash
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-v20260330/ --rimage /data-pool/Halpha/coadds-v20260330/VF-215.042+03.959-INT-20190208-p131-r.fits 
```

```bash
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-v20260330/ --rimage /data-pool/Halpha/coadds-v20260330/VF-208.876+05.130-HDI-20180313-p056-R.fits  --overwrite
```

```bash
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-v20260330/ --rimage /data-pool/Halpha/coadds-v20260330/VF-185.450+04.633-INT-20220505-VFID5966-r.fits --overwrite
```

```bash
get_cutouts --catalog ~/research/Virgo/tables-north/v2/vf_v2_main.fits --scheme virgo --maxcorrection 5 --psfdir /data-pool/Halpha/psf-images-v20260330/ --rimage /data-pool/Halpha/coadds-v20260330/VF-217.936+03.183-INT-20190210-p140-r.fits 
```


# Notes


I need to implement a minimum cutout size, and it probably should be around 40-50 arcsec on a side. I would also like to implement a max size, or a buffer that is applied to big galaxy b/c I don't need to double the size. for bigger galaxies, it would be better to have a buffer - maybe an arcmin?


Tried something like this:
```python
MIN_CUTOUT_SIZE = 50.0      # arcsec
EDGE_BUFFER = 60.0          # arcsec

diameter = 2 * gradius[i]

size_arcsec = max(
    MIN_CUTOUT_SIZE,
    diameter + EDGE_BUFFER
)
```

I realized that if the galaxy is inclined (smaller B/A), then a buffer size of 120 arcsec is fine.  If the galaxy is more face on, then the buffer does not seem big enough.  so it's almost like we need the area of sky to be way more than the area of the galaxy.

Tried something like:
```python
size_arcsec = max(
    75,
    2 * radius + 75 + 75 * np.sqrt(BA)
)
```

Also considered using a fixed ratio of sky_pixels/area_of_galaxy

```python
def get_cutout_size_arcsec(radius_arcsec, ba,
                           min_cutout_size=75.0,
                           min_ba=0.25,
                           sky_to_gal_area=5.0,
                           max_cutout_size=None):
    """
    Compute square cutout size from desired sky/galaxy area ratio.

    Parameters
    ----------
    radius_arcsec : float
        Galaxy semi-major axis radius in arcsec.
    ba : float
        Axis ratio b/a.
    min_cutout_size : float
        Minimum cutout side length in arcsec.
    min_ba : float
        Floor on b/a to avoid bogus tiny BA values.
    sky_to_gal_area : float
        Desired sky area divided by galaxy ellipse area.
        Example: 5 means sky area is 5 times galaxy area.
    max_cutout_size : float or None
        Optional maximum cutout size in arcsec.

    Returns
    -------
    size_arcsec : float
        Square cutout side length in arcsec.
    """

    radius_arcsec = float(radius_arcsec)

    try:
        ba = float(ba)
    except Exception:
        ba = np.nan

    if not np.isfinite(ba):
        ba = 1.0

    ba_eff = np.clip(ba, min_ba, 1.0)

    gal_area = np.pi * radius_arcsec**2 * ba_eff

    required_total_area = (1.0 + sky_to_gal_area) * gal_area
    area_based_size = np.sqrt(required_total_area)

    size_arcsec = max(min_cutout_size, area_based_size)

    if max_cutout_size is not None:
        size_arcsec = min(size_arcsec, max_cutout_size)

    return float(size_arcsec)
```
