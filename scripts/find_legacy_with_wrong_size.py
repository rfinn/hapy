#!/usr/bin/env python

from pathlib import Path
from astropy.io import fits
import re

#OLD = Path("/data-pool/Halpha/hapy-output-20260612/cutouts")

# using the same path for both b/c I already rsynced the legacy images from 20260612
OLD = Path("/data-pool/Halpha/hapy-output-20260620/cutouts")
NEW = Path("/data-pool/Halpha/hapy-output-20260620/cutouts")

PIX_SCALE = 0.262  # Legacy pixscale used by fetch_legacy_cutouts.py, update if different

needs = []

for newdir in sorted(NEW.iterdir()):
    if not newdir.is_dir():
        continue

    tag = newdir.name
    rfile = newdir / f"{tag}-R.fits"
    legdir = newdir / "legacy"

    if not rfile.exists():
        continue

    data = fits.getdata(rfile)
    new_size_pix = max(data.shape)
    new_legacy_size = int(round(new_size_pix * fits.getheader(rfile).get("PIXSCALE", 0.4) / PIX_SCALE))

    old_leg = sorted((OLD / tag / "legacy").glob("*-legacy-*-r.fits"))
    if not old_leg:
        needs.append(newdir)
        continue

    m = re.search(r"-legacy-(\d+)-r\.fits$", old_leg[0].name)
    if not m:
        needs.append(newdir)
        continue

    old_legacy_size = int(m.group(1))

    if new_legacy_size > old_legacy_size:
        print(f"{tag}: old legacy={old_legacy_size}, new needed={new_legacy_size}")
        needs.append(newdir)

with open("legacy_redownload_needed.txt", "w") as out:
    for d in needs:
        out.write(str(d) + "\n")

print(f"\nNeed redownload: {len(needs)}")
print("Wrote legacy_redownload_needed.txt")
