# Transferring Legacy Images After Implementing a Min Cutout Size
```
rsync -av hapy-output-20260612/cutouts/ \
hapy-output-20260620/cutouts/ --include '*/' --include 'legacy/***' \
--exclude '*' --exclude '*logs*' --ignore-existing --prune-empty-dirs --dry-run
```


```
python ~/github/hapy/scripts/find_legacy_with_wrong_size.py
```

```bash
while read d; do
    rm -f "$d"/legacy/*legacy*.fits "$d"/legacy/*legacy*.jpg
done < legacy_redownload_needed.txt
```


ok, so just do the regular rsync first, then run this script, and then download the new legacy images. the only issue is that having two versions of the legacy images may throw a wrench in the scripts that create the cs-gr images.


It's saying that 849/849 need new legacy images.


```
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260620$ gethead NAXIS1 NAXIS2 CD1_1 CD2_2 PIXSCALE cutouts/VFID10*/legacy/*r.fits
VFID1008-NGC3990-INT-20220502-VFID1010-legacy-735-r.fits                  735 735 -9.25E-05 9.25E-05 
VFID1009-NGC3998-INT-20220502-VFID1010-legacy-1491-r.fits                 1491 1491 -9.25E-05 9.25E-05 
VFID1014-WISEAJ115701.80+552510.7-INT-20220502-VFID1010-legacy-297-r.fits 297 297 -9.25E-05 9.25E-05 
VFID1016-WISEAJ115813.70+552316.8-INT-20220502-VFID1010-legacy-226-r.fits 226 226 -9.25E-05 9.25E-05 
VFID1018-NGC3972-INT-20220502-VFID1010-legacy-2810-r.fits                 2810 2810 -9.25E-05 9.25E-05 
VFID1035-NGC3982-INT-20220502-VFID1010-legacy-820-r.fits                  820 820 -9.25E-05 9.25E-05 

```

ok, looks like I had legacy images with the pixel scale from each individual telescope/instrument (INT, BOK, HDI). so I am going to wipe these out and download a fresh set at the native pixel scale.



```
find cutouts -mindepth 2 -maxdepth 2 -type d -name legacy > legacy_dirs.txt
```

```
while read d; do echo rm -f "$d"/*legacy*.fits "$d"/*legacy*.jpg; done < legacy_dirs.txt
```

```
while read d; do rm -f "$d"/*legacy*.fits "$d"/*legacy*.jpg; done < legacy_dirs.txt
```

regenerate cutout list
```
find cutouts -mindepth 1 -maxdepth 1 -type d | sort > cutout_list.txt
```


```
parallel --bar -j 2 --joblog fetch_legacy_native.joblog python ~/github/hapy/scripts/fetch_legacy_cutouts.py --cutout-dir {} :::: cutout_list.txt
```


# Transferring Manual Masks After Implementing a Min Cutout Size

- I have a script to check the size of cutouts in the new directory `hapy-output-20260620` and the previous one `hapy-output-20260612`. 
- If the cutout size has changed, and the previous directory had a manual mask, then add file to a list: `needs_new_manual_mask.txt`.
- Copy the manual mask if the manual mask exists and the cutouts are the same size.


To run the mask copy script:
```
python ~/github/hapy/scripts/copy_manual_masks.py 
```


I got this output:
```
Copied manual masks: 58
Skipped existing masks: 0
Need new manual masks: 41
Wrote needs_new_manual_mask.txt
```

So.... need to remake a bunch of masks - woo hoo!

```bash
while read d; do
    tag=$(basename "$d")
    if [ ! -f "$d/${tag}-mask-manual.fits" ]; then
        echo "$d"
    fi
done < needs_new_manual_mask.txt > still_need_manual_mask.txt
```

Then count how many are left:
```bash
(hapy) rfinn@draco:/data-pool/Halpha/hapy-output-20260620$ wc -l still_need_manual_mask.txt 
40 still_need_manual_mask.txt
```


To get the next images to mask:
```
head still_need_manual_mask.txt
```

And run the mask
```

```
