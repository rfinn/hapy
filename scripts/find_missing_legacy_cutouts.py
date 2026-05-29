#!/usr/bin/env python

from pathlib import Path

CUTOUT_ROOT = Path("/data-pool/Halpha/hapy-output-20260519-hybrid/cutouts")

missing = []

for cutdir in sorted(CUTOUT_ROOT.iterdir()):

    if not cutdir.is_dir():
        continue

    if not cutdir.name.startswith("VF"):
        continue

    legacy_dir = cutdir / "legacy"

    if not legacy_dir.exists():
        print(f"MISSING legacy dir: {cutdir.name}")
        missing.append(cutdir)
        continue

    gfile = legacy_dir / "legacy-g.fits"
    rfile = legacy_dir / "legacy-r.fits"
    zfile = legacy_dir / "legacy-z.fits"

    have_g = gfile.exists()
    have_r = rfile.exists()
    have_z = zfile.exists()

    if not (have_g and have_r and have_z):
        print(
            f"MISSING legacy files: {cutdir.name} "
            f"(g={have_g}, r={have_r}, z={have_z})"
        )
        missing.append(cutdir)

print()
print(f"Missing legacy images: {len(missing)}")

outfile = "missing_legacy_cutouts.txt"

with open(outfile, "w") as out:
    for d in missing:
        out.write(str(d) + "\n")

print(f"Wrote {outfile}")
