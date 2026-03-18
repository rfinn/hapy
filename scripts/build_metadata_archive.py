#!/usr/bin/env python3
"""
build_metadata_archive.py

Prepare archival galaxy cutouts for hapy by:

1. Reading an authoritative image list (e.g. virgo.list)
2. Grouping files by object prefix (e.g. ic3392, n4064)
3. Selecting the correct CS, R, and optional mask image for each target
4. Matching each target to a Virgo / VFS catalog row
5. Building a standardized tag:
       {VFID}-{NEDname}-{instrument}-{YYYYMMDD}-{target}
6. Renaming the directory and selected files to:
       <tag>-R.fits
       <tag>-CS.fits
       <tag>-mask.fits
7. Writing metadata.json with the minimal fields needed by run_analysis.py

Notes
-----
- This script assumes the image list is authoritative.
- Extra FITS files in the original directories are ignored.
- The raw CS filenames may contain 'ha'; output is standardized to '-CS.fits'.
- WCS is assumed to be poor for these archive images, so metadata sets:
      scheme = "archive"
      has_bad_wcs = True

Example
-------
python build_metadata_archive.py \
    --archive-root /data/archive_cutouts \
    --catalog /data/vf_v2_main.fits \
    --image-list /data/virgo.list \
    --instrument HDI \
    --target-mode foldername \
    --dry-run
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from astropy.io import fits
from astropy.table import Table


FWHM_KEYS = ("FWHM",)
FILTER_KEYS = ("FILTER",)
# -----------------------------------------------------------------------------
# Catalog column candidates: adjust if needed for your actual Virgo/VFS tables
# -----------------------------------------------------------------------------
VFID_COL_CANDIDATES = ("VFID", "vfid")
RA_COL_CANDIDATES = ("RA", "ra", "RA_DEG", "RAdeg", "RAJ2000")
DEC_COL_CANDIDATES = ("DEC", "dec", "DEC_DEG", "DECdeg", "DEJ2000")
SMA_COL_CANDIDATES = ("SMA_ARCSEC", "sma_arcsec", "SMA", "radius")
BA_COL_CANDIDATES = ("BA", "ba", "B_A", "q")
PA_COL_CANDIDATES = ("PA", "pa", "PA_DEG", "pa_deg")

# For name matching / output GALNAME
NAME_COL_CANDIDATES = (
    "NEDname",
    "NEDNAME",
    "GALNAME",
    "NAME",
    "Object",
    "OBJECT",
)


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
def pick_col(tbl: Table, candidates: tuple[str, ...], required: bool = True) -> Optional[str]:
    """Return first matching column name from candidates."""
    for c in candidates:
        if c in tbl.colnames:
            return c
    if required:
        raise KeyError(f"Could not find any of {candidates} in catalog columns: {tbl.colnames}")
    return None


def sanitize_name(s: str) -> str:
    """
    Make a string safe for file/directory names.
    Keeps letters, numbers, underscore, hyphen.
    Removes spaces and most punctuation.
    """
    s = str(s).strip()
    s = s.replace(" ", "")
    s = s.replace("/", "-")
    s = re.sub(r"[^A-Za-z0-9_\-]+", "", s)
    return s or "UNKNOWN"


def normalize_match_key(s: str) -> str:
    """
    Normalize object names for matching.
    Example:
        'NGC 4354' -> 'ngc4354'
        'n4354'    -> 'n4354'
        'IC3392'   -> 'ic3392'
    """
    s = str(s).strip().lower()
    s = s.replace(" ", "")
    s = s.replace("_", "")
    s = s.replace("-", "")
    return s


def read_header_value(path: Path, keys: tuple[str, ...]) -> Optional[str]:
    """Read first non-empty FITS header key from path."""
    try:
        hdr = fits.getheader(path)
    except Exception:
        return None
    for k in keys:
        if k in hdr and hdr[k] not in ("", None):
            return str(hdr[k]).strip()
    return None

def read_header_float(path: Path, keys: tuple[str, ...]) -> Optional[float]:
    """Read first non-empty FITS header key from path and cast to float."""
    try:
        hdr = fits.getheader(path)
    except Exception:
        return None
    for k in keys:
        if k in hdr and hdr[k] not in ("", None):
            try:
                return float(hdr[k])
            except Exception:
                return None
    return None

def extract_yyyymmdd(date_str: Optional[str]) -> Optional[str]:
    """
    Convert DATE-OBS-ish strings to YYYYMMDD.
    Handles:
      YYYY-MM-DD
      YYYY-MM-DDThh:mm:ss
      YYYY/MM/DD
      YYYYMMDD
    """
    if not date_str:
        return None
    m = re.match(r"^\s*(\d{4})[-/]?(\d{2})[-/]?(\d{2})", date_str)
    if not m:
        return None
    return f"{m.group(1)}{m.group(2)}{m.group(3)}"


def ensure_unique_path(path: Path) -> Path:
    """Append suffix if needed to avoid collisions."""
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    parent = path.parent
    i = 2
    while True:
        trial = parent / f"{stem}_{i}{suffix}"
        if not trial.exists():
            return trial
        i += 1


def classify_list_entry(filename: str) -> str:
    """
    Classify an entry from the image list.

    Priority:
      1. mask  -> contains 'mask'
      2. cs    -> contains 'ha'
      3. r     -> contains 'r' (and not already classified)
      4. other
    """
    name = filename.lower()

    if "mask" in name:
        return "mask"
    if "ha" in name:
        return "cs"
    if "r" in name:
        return "r"
    return "other"


def infer_object_prefix(filename: str) -> str:
    """
    Infer object prefix from the beginning of the filename.

    Examples:
      ic3392ha.fits           -> ic3392
      ic3392r.mask.fits       -> ic3392
      n4064hacljccsasfc.fits  -> n4064
      n4192rcljcbc.fits       -> n4192

    Strategy:
      take basename without extension
      then capture leading [a-z]+[0-9]+
    """
    base = Path(filename).name.lower()
    base = re.sub(r"\.fits?$", "", base)
    m = re.match(r"^([a-z]+[0-9]+)", base)
    if m:
        return m.group(1)
    return base


def format_vfid(vfid_value) -> str:
    """
    Normalize VFID string.

    If already like VFID3084 -> keep uppercase.
    If numeric-ish, zero-pad to VFIDxxxxx.
    """
    s = str(vfid_value).strip()
    if re.fullmatch(r"VFID\d+", s, re.IGNORECASE):
        return s.upper()

    digits = re.sub(r"\D+", "", s)
    if digits:
        return f"VFID{int(digits):05d}"

    return sanitize_name(s).upper()


# -----------------------------------------------------------------------------
# Image-list parsing
# -----------------------------------------------------------------------------
def parse_image_list(image_list_path: Path) -> dict[str, dict[str, str]]:
    """
    Parse authoritative image list into:
        {
          "ic3392": {"cs": "...", "r": "...", "mask": "..."},
          "n4064":  {"cs": "...", "r": "..."},
          ...
        }

    Assumes at most one CS, one R, one mask per object in the list.
    """
    manifest: dict[str, dict[str, str]] = {}

    with open(image_list_path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            kind = classify_list_entry(line)
            if kind == "other":
                continue

            obj = infer_object_prefix(line)
            manifest.setdefault(obj, {})

            if kind in manifest[obj]:
                raise ValueError(
                    f"Duplicate {kind} entry for object '{obj}' in {image_list_path}: "
                    f"{manifest[obj][kind]} and {line}"
                )

            manifest[obj][kind] = line

    return manifest


# -----------------------------------------------------------------------------
# Catalog matching
# -----------------------------------------------------------------------------
def build_catalog_name_index(catalog: Table, name_cols: list[str]) -> dict[str, int]:
    """
    Build dict: normalized object name -> row index
    Uses first occurrence if duplicates exist.
    """
    idx: dict[str, int] = {}

    for i in range(len(catalog)):
        for col in name_cols:
            val = catalog[col][i]
            if val is None:
                continue
            key = normalize_match_key(val)
            if key and key not in idx:
                idx[key] = i

    return idx


def match_object_to_catalog(obj_prefix: str, catalog: Table, name_index: dict[str, int]) -> Optional[int]:
    """
    Match object prefix like 'ic3392' or 'n4064' to catalog row.

    Tries:
      exact normalized match
      ngc<n> expansion for n####
      ic<n> expansion for ic####
    """
    q = normalize_match_key(obj_prefix)

    if q in name_index:
        return name_index[q]

    # Expand "n4064" -> "ngc4064"
    m = re.fullmatch(r"n(\d+)", q)
    if m:
        alt = f"ngc{m.group(1)}"
        if alt in name_index:
            return name_index[alt]

    # Expand "i3392" not used here, but keep future flexibility if needed
    m = re.fullmatch(r"ic(\d+)", q)
    if m:
        alt = f"ic{m.group(1)}"
        if alt in name_index:
            return name_index[alt]

    return None


# -----------------------------------------------------------------------------
# Planning / execution
# -----------------------------------------------------------------------------
@dataclass
class Plan:
    obj_prefix: str
    old_dir: Path
    new_dir: Path

    cs_old: Path
    r_old: Path
    mask_old: Optional[Path]

    cs_new: Path
    r_new: Path
    mask_new: Optional[Path]

    metadata_path: Path
    metadata: dict


def find_exact_file_in_dir(dirpath: Path, filename: str) -> Path:
    """
    Find exact file in dirpath matching filename.
    First tries exact path. Then case-insensitive fallback by basename.
    """
    exact = dirpath / filename
    if exact.exists():
        return exact

    lower_target = filename.lower()
    matches = [p for p in dirpath.iterdir() if p.is_file() and p.name.lower() == lower_target]
    if len(matches) == 1:
        return matches[0]

    raise FileNotFoundError(f"Could not find listed file '{filename}' in directory '{dirpath}'")


def build_tag(
    vfid: str,
    galname: str,
    instrument: str,
    yyyymmdd: str,
    target: str,
) -> str:
    """Construct standardized tag."""
    return f"{vfid}-{sanitize_name(galname)}-{sanitize_name(instrument)}-{sanitize_name(yyyymmdd)}-{sanitize_name(target)}"


def build_plans(
    archive_root: Path,
    catalog: Table,
    manifest: dict[str, dict[str, str]],
    instrument_default: str,
    date_default: str,
    target_mode: str,
    scheme: str = "archive",
) -> tuple[list[Plan], list[str]]:
    """
    Build execution plans.

    Returns
    -------
    plans : list[Plan]
        Objects that matched catalog + had required listed files.
    skipped : list[str]
        Human-readable reasons for skipped objects.
    """
    vfid_col = pick_col(catalog, VFID_COL_CANDIDATES)
    ra_col = pick_col(catalog, RA_COL_CANDIDATES)
    dec_col = pick_col(catalog, DEC_COL_CANDIDATES)
    sma_col = pick_col(catalog, SMA_COL_CANDIDATES)
    ba_col = pick_col(catalog, BA_COL_CANDIDATES)
    pa_col = pick_col(catalog, PA_COL_CANDIDATES)

    available_name_cols = [c for c in NAME_COL_CANDIDATES if c in catalog.colnames]
    if not available_name_cols:
        raise KeyError(f"No usable name columns found. Need one of {NAME_COL_CANDIDATES}")

    galname_col = available_name_cols[0]
    name_index = build_catalog_name_index(catalog, available_name_cols)

    plans: list[Plan] = []
    skipped: list[str] = []

    # Stable ordering
    object_keys = sorted(manifest.keys())

    seq_counter = 1
    for obj_prefix in object_keys:
        info = manifest[obj_prefix]

        # Require both cs and r
        if "cs" not in info or "r" not in info:
            skipped.append(f"{obj_prefix}: missing required CS or R entry in image list")
            continue

        # Require matching directory
        old_dir = archive_root / obj_prefix
        if not old_dir.exists() or not old_dir.is_dir():
            skipped.append(f"{obj_prefix}: directory not found at {old_dir}")
            continue

        row_idx = match_object_to_catalog(obj_prefix, catalog, name_index)
        if row_idx is None:
            skipped.append(f"{obj_prefix}: no catalog match; skipping")
            continue

        row = catalog[row_idx]

        try:
            cs_old = find_exact_file_in_dir(old_dir, info["cs"])
            r_old = find_exact_file_in_dir(old_dir, info["r"])
            mask_old = find_exact_file_in_dir(old_dir, info["mask"]) if "mask" in info else None
        except FileNotFoundError as e:
            skipped.append(f"{obj_prefix}: {e}")
            continue

        # read info from header
        instrument = read_header_value(r_old, ("INSTRUME", "INSTRUMENT", "CAMERA", "DETECTOR")) or instrument_default
        yyyymmdd = extract_yyyymmdd(read_header_value(r_old, ("DATE-OBS", "DATEOBS", "DATE_OBS"))) or date_default
        
        r_fwhm_arcsec = read_header_float(r_old, FWHM_KEYS)        
        cs_fwhm_arcsec = read_header_float(cs_old, FWHM_KEYS)

        r_filter = read_header_value(r_old, FILTER_KEYS)        
        cs_filter = read_header_value(cs_old,FILTER_KEYS)

        
        # Gather tag components
        vfid = format_vfid(row[vfid_col])
        galname = str(row[galname_col]).strip()

        
        if target_mode == "foldername":
            target = obj_prefix
        else:
            target = f"p{seq_counter:03d}"
            seq_counter += 1

        tag = build_tag(vfid, galname, instrument, yyyymmdd, target)

        new_dir = archive_root / tag
        new_dir = ensure_unique_path(new_dir)

        tag_final = new_dir.name  # in case ensure_unique_path altered it

        cs_new = new_dir / f"{tag_final}-CS.fits"
        r_new = new_dir / f"{tag_final}-R.fits"
        mask_new = new_dir / f"{tag_final}-mask.fits" if mask_old else None

        metadata = {
            "VFID": vfid,
            "objid": vfid,
            "GALNAME": galname,
            "tag": tag_final,
            "ra": float(row[ra_col]),
            "dec": float(row[dec_col]),
            "sma_arcsec": float(row[sma_col]),
            "ba": float(row[ba_col]),
            "pa_deg": float(row[pa_col]),
            "scheme": scheme,
            "has_bad_wcs": True,
            "r_fits": r_new.name,
            "cs_fits": cs_new.name,
            "mask_fits": mask_new.name if mask_new else "",
            "rimage_fwhm_arcsec": r_fwhm_arcsec,
            "csimage_fwhm_arcsec": cs_fwhm_arcsec,
            "rimage_filter": r_filter or "",
            "csimage_filter": cs_filter or "",
            # TODO - GET FILTER CENTER WAVELENGTH AND WIDTH FOR ARCHIVE GALAXIES
            # Table 5 in KKY
            rfilter_name = image_set.r.filter_file,
            rfilter_center_A = image_set.r.filter_center,
            rfilter_width_A = image_set.r.filter_width,
            hafilter_name = image_set.h.filter_file,
            hafilter_center_A = image_set.h.filter_center,
            hafilter_width_A = image_set.h.filter_width,
        }

        plans.append(
            Plan(
                obj_prefix=obj_prefix,
                old_dir=old_dir,
                new_dir=new_dir,
                cs_old=cs_old,
                r_old=r_old,
                mask_old=mask_old,
                cs_new=cs_new,
                r_new=r_new,
                mask_new=mask_new,
                metadata_path=new_dir / "metadata.json",
                metadata=metadata,
            )
        )

    return plans, skipped


def move_or_link(src: Path, dst: Path, link: bool = False) -> None:
    """Move or symlink a file."""
    if link:
        dst.symlink_to(src.resolve())
    else:
        shutil.move(str(src), str(dst))


def apply_plan(plan: Plan, link: bool = False, keep_old_dir: bool = False) -> None:
    """
    Execute one rename/write plan.

    If link=True:
      - create new standardized dir
      - symlink selected files
      - leave old dir untouched

    If link=False:
      - create new standardized dir
      - move selected files into it
      - optionally remove old dir if empty and keep_old_dir is False
    """
    plan.new_dir.mkdir(parents=False, exist_ok=False)

    move_or_link(plan.cs_old, plan.cs_new, link=link)
    move_or_link(plan.r_old, plan.r_new, link=link)

    if plan.mask_old and plan.mask_new:
        move_or_link(plan.mask_old, plan.mask_new, link=link)

    with open(plan.metadata_path, "w") as f:
        json.dump(plan.metadata, f, indent=2)

    if not link and not keep_old_dir:
        try:
            next(plan.old_dir.iterdir())
        except StopIteration:
            plan.old_dir.rmdir()


def print_plan_summary(plans: list[Plan], skipped: list[str]) -> None:
    """Pretty-print planned actions."""
    print("\nPlanned objects:")
    print("-" * 80)
    for p in plans:
        print(f"{p.obj_prefix}  ->  {p.new_dir.name}")
        print(f"   CS:   {p.cs_old.name}   ->   {p.cs_new.name}")
        print(f"   R:    {p.r_old.name}    ->   {p.r_new.name}")
        if p.mask_old:
            print(f"   MASK: {p.mask_old.name} ->   {p.mask_new.name}")
        else:
            print("   MASK: none")
        print(f"   META: {p.metadata_path.name}")
        print()

    if skipped:
        print("\nSkipped objects:")
        print("-" * 80)
        for s in skipped:
            print(f" - {s}")
    else:
        print("\nSkipped objects: none")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Build archive metadata + standardized names for hapy.")
    parser.add_argument("--archive-root", type=Path, required=True,
                        help="Root directory containing per-object archive folders (e.g. ic3392/, n4064/)")
    parser.add_argument("--catalog", type=Path, required=True,
                        help="Virgo/VFS catalog readable by astropy Table.read")
    parser.add_argument("--image-list", type=Path, required=True,
                        help="Authoritative list of image filenames to use (e.g. virgo.list)")
    parser.add_argument("--instrument", default="ARCH",
                        help="Default instrument string if not found in FITS header")
    parser.add_argument("--date-default", default="19990101",
                        help="Fallback YYYYMMDD if DATE-OBS missing or unparsable")
    parser.add_argument("--target-mode", choices=("foldername", "sequential"), default="foldername",
                        help="How to populate the {target} field in the tag")
    parser.add_argument("--scheme", default="archive",
                        help="Metadata scheme name; default='archive'")
    parser.add_argument("--link", action="store_true",
                        help="Create new standardized directories with symlinks instead of moving files")
    parser.add_argument("--keep-old-dir", action="store_true",
                        help="When moving files, do not remove old directory even if it becomes empty")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print planned actions but do not modify filesystem")
    args = parser.parse_args()

    catalog = Table.read(args.catalog)
    manifest = parse_image_list(args.image_list)

    plans, skipped = build_plans(
        archive_root=args.archive_root,
        catalog=catalog,
        manifest=manifest,
        instrument_default=args.instrument,
        date_default=args.date_default,
        target_mode=args.target_mode,
        scheme=args.scheme,
    )

    print_plan_summary(plans, skipped)

    if args.dry_run:
        print("\nDry run only: no changes made.")
        return

    for p in plans:
        apply_plan(p, link=args.link, keep_old_dir=args.keep_old_dir)

    print(f"\nDone. Created {len(plans)} standardized archive directories.")


if __name__ == "__main__":
    main()
