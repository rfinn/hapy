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
