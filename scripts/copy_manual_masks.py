#!/usr/bin/env python

from pathlib import Path
from astropy.io import fits
import shutil

OLD_ROOT = Path("/data-pool/Halpha/hapy-output-20260612/cutouts")
NEW_ROOT = Path("/data-pool/Halpha/hapy-output-20260620/cutouts")

OUTFILE = Path("needs_new_manual_mask.txt")


def find_r_cutout(cutdir):
    tag = cutdir.name
    rfile = cutdir / f"{tag}-R.fits"
    if rfile.exists():
        return rfile

    matches = sorted(cutdir.glob("*-R.fits"))
    return matches[0] if matches else None


def get_shape(fitsfile):
    hdr = fits.getheader(fitsfile)
    return int(hdr["NAXIS1"]), int(hdr["NAXIS2"])


def main():
    needs = []
    copied = 0
    skipped_existing = 0

    for newdir in sorted(NEW_ROOT.iterdir()):
        if not newdir.is_dir():
            continue

        tag = newdir.name
        olddir = OLD_ROOT / tag

        if not olddir.exists():
            continue

        old_manual = olddir / f"{tag}-mask-manual.fits"
        new_manual = newdir / f"{tag}-mask-manual.fits"

        if not old_manual.exists():
            continue

        old_r = find_r_cutout(olddir)
        new_r = find_r_cutout(newdir)

        if old_r is None or new_r is None:
            print(f"WARNING: missing R cutout for {tag}")
            continue

        old_shape = get_shape(old_r)
        new_shape = get_shape(new_r)

        if old_shape == new_shape:
            if new_manual.exists():
                skipped_existing += 1
                continue

            print(f"copying manual mask: {tag}")
            shutil.copy2(old_manual, new_manual)
            copied += 1

        else:
            print(f"{tag}: old={old_shape}, new={new_shape}, manual mask needs remake")
            needs.append(str(newdir))

    with open(OUTFILE, "w") as out:
        for d in needs:
            out.write(d + "\n")

    print()
    print(f"Copied manual masks: {copied}")
    print(f"Skipped existing masks: {skipped_existing}")
    print(f"Need new manual masks: {len(needs)}")
    print(f"Wrote {OUTFILE}")


if __name__ == "__main__":
    main()
