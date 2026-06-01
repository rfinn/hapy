#!/usr/bin/env python

from pathlib import Path

CUTOUT_ROOT = Path("/data-pool/Halpha/hapy-output-20260519/cutouts")

missing_csgr = []
missing_csgr_phot = []
missing_any = []

for cutdir in sorted(CUTOUT_ROOT.iterdir()):

    if not cutdir.is_dir():
        continue

    if not cutdir.name.startswith("VF"):
        continue

    tag = cutdir.name

    csgr = cutdir / f"{tag}-CS-gr.fits"
    csgr_phot = cutdir / f"{tag}-CS-gr-phot.fits"

    have_csgr = csgr.exists()
    #have_phot = csgr_phot.exists()
    have_phot = csgr_phot.exists() and csgr_phot.stat().st_size > 0

    if not have_csgr:
        missing_csgr.append(cutdir)

    if not have_phot:
        missing_csgr_phot.append(cutdir)

    if (not have_csgr) or (not have_phot):
        missing_any.append(cutdir)

print()
print(f"Missing CS-gr image      : {len(missing_csgr)}")
print(f"Missing CS-gr photometry : {len(missing_csgr_phot)}")
print(f"Missing either           : {len(missing_any)}")
print()

with open("missing_csgr.txt", "w") as out:
    for d in missing_csgr:
        out.write(f"{d}\n")

with open("missing_csgr_phot.txt", "w") as out:
    for d in missing_csgr_phot:
        out.write(f"{d}\n")

with open("missing_csgr_or_phot.txt", "w") as out:
    for d in missing_any:
        out.write(f"{d}\n")

print("Wrote:")
print("  missing_csgr.txt")
print("  missing_csgr_phot.txt")
print("  missing_csgr_or_phot.txt")
