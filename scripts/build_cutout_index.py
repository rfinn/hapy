#!/usr/bin/env python

"""
build_cutout_index.py

Build an index page for per-cutout HAPY webpages.

Expected layout:
    runroot/
        cutouts/
            <tag>/
                *-results.ecsv
        html/
            cutouts/
                <tag>/
                    <tag>.html
                    ...png/jpg assets...
"""

import argparse
from pathlib import Path

from build_web_common import (
    find_results_file,
    read_results_row,
    get_result,
    fmt_result,
    status_cell,
    combined_flag,
)


def find_local_legacy_jpg(html_cutout_dir):
    """
    Find a local legacy jpg already copied into the html/cutouts/<tag>/ directory.
    """
    patterns = ["*legacy*.jpg", "*legacy*.jpeg", "*.jpg", "*.jpeg"]
    for pattern in patterns:
        matches = sorted(html_cutout_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def collect_entries(runroot):
    runroot = Path(runroot).resolve()
    html_root = runroot / "html" / "cutouts"
    cutout_root = runroot / "cutouts"

    entries = []

    for subdir in sorted(html_root.iterdir()):
        if not subdir.is_dir():
            continue

        tag = subdir.name
        galname = tag.split('-')[0]
        html_file = subdir / f"{tag}.html"
        if not html_file.exists():
            continue

        cutout_dir = cutout_root / tag
        results_file = find_results_file(cutout_dir)
        results_row = read_results_row(results_file) if results_file is not None else None

        legacy_jpg = find_local_legacy_jpg(subdir)

        mask_ok = get_result(results_row, "MASK_OK", None)
        bright_star = get_result(results_row, "BRIGHT_STAR_FLAG", None)
        phot_ok = get_result(results_row, "PHOT_OK", None)
        psf_ok = get_result(results_row, "PSF_OK", None)

        r_prof_ok = get_result(results_row, "R_PROFILE_OK", None)
        h_prof_ok = get_result(results_row, "H_PROFILE_OK", None)
        profile_ok = combined_flag(r_prof_ok, h_prof_ok)

        r_sm_ok = get_result(results_row, "R_SM_FLAG", None)
        h_sm_ok = get_result(results_row, "H_SM_FLAG", None)
        statmorph_ok = combined_flag(r_sm_ok, h_sm_ok)

        gal_nc_ok = get_result(results_row, "GAL_NC_OK", None)
        gal_cv_ok = get_result(results_row, "GAL_CV_OK", None)
        galfit_ok = combined_flag(gal_nc_ok, gal_cv_ok)

        entry = dict(
            tag=tag,
            galname=galname,
            html_file=html_file,
            legacy_jpg=legacy_jpg,
            results_row=results_row,
            mask_ok=mask_ok,
            bright_star=bright_star,
            phot_ok=phot_ok,
            psf_ok=psf_ok,
            profile_ok=profile_ok,
            r_prof_ok=r_prof_ok,
            h_prof_ok=h_prof_ok,            
            statmorph_ok=statmorph_ok,
            r_sm_ok=r_sm_ok,
            h_sm_ok=h_sm_ok,            
            galfit_ok=galfit_ok,
            gal_nc_ok=gal_nc_ok,
            gal_cv_ok=gal_cv_ok,            
            status=get_result(results_row, "STATUS", ""),
            stage=get_result(results_row, "STAGE", ""),
            r_fwhm=fmt_result(results_row, "R_FWHM_PSF", "{:.2f}"),
            h_fwhm=fmt_result(results_row, "H_FWHM_PSF", "{:.2f}"),
        )
        entries.append(entry)

    return entries


def write_index(entries, outfile):
    outfile = Path(outfile)

    lines = []
    lines.append("<html><body>")
    lines.append("<style type='text/css'>")
    lines.append("body { font-family: Arial, sans-serif; margin: 20px; }")
    lines.append("table, td, th { padding: 6px; text-align: center; border: 1px solid black; border-collapse: collapse; }")
    lines.append("th { background-color: #f0f0f0; }")
    lines.append("img.thumb { max-width: 180px; height: auto; }")
    lines.append("a { text-decoration: none; color: #1565c0; }")
    lines.append("a:hover { text-decoration: underline; }")
    lines.append("</style>")

    lines.append("<h1>HAPY Cutout Index</h1>")
    lines.append(f"<p>Found {len(entries)} cutout pages.</p>")

    lines.append("<table width='95%'>")
    lines.append("<tr>")
    headers = [
        "Index",
        "GALID",
        "Legacy",
        "Cutout Page",
        "PSF",
        "R FWHM",
        "H FWHM",        
        "Mask",
        "Star Flag",
        "Phot",  
        "R Prof",
        "H Prof",        
        "R SM",
        "H SM",        
        "GAL NC",
        "GAL CV",        
        "Status",
        "Stage",
    ]
    for h in headers:
        lines.append(f"<th>{h}</th>")
    lines.append("</tr>")

    for i, e in enumerate(entries, start=1):
        lines.append("<tr>")
        lines.append(f"<td>{i}</td>")
        lines.append(f"<td>{e['galname']}</td>")

        if e["legacy_jpg"] is not None:
            rel_jpg = f"{e['tag']}/{e['legacy_jpg'].name}"
            lines.append(
                f"<td><a href='{rel_jpg}' target='_blank'>"
                f"<img class='thumb' src='{rel_jpg}' alt='legacy image for {e['tag']}'></a></td>"
            )
        else:
            lines.append("<td>Missing</td>")

        rel_html = f"{e['tag']}/{e['tag']}.html"
        lines.append(f"<td><a href='{rel_html}'>{e['tag']}</a></td>")
        
        lines.append(f"<td>{status_cell(e['psf_ok'])}</td>")
        lines.append(f"<td>{e['r_fwhm']}</td>")
        lines.append(f"<td>{e['h_fwhm']}</td>")        
        lines.append(f"<td>{status_cell(e['mask_ok'])}</td>")
        lines.append(f"<td>{status_cell((e['bright_star']==0))}</td>")
        lines.append(f"<td>{status_cell(e['phot_ok'])}</td>")

        lines.append(f"<td>{status_cell(e['r_prof_ok'])}</td>")
        lines.append(f"<td>{status_cell(e['h_prof_ok'])}</td>")
        
        lines.append(f"<td>{status_cell(e['r_sm_ok'])}</td>")
        lines.append(f"<td>{status_cell(e['h_sm_ok'])}</td>")        
        lines.append(f"<td>{status_cell(e['gal_nc_ok'])}</td>")
        lines.append(f"<td>{status_cell(e['gal_cv_ok'])}</td>")        

        lines.append(f"<td>{e['status']}</td>")
        lines.append(f"<td>{e['stage']}</td>")
        lines.append("</tr>")

    lines.append("</table>")
    lines.append("</body></html>")

    outfile.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {outfile}")


def main():
    parser = argparse.ArgumentParser(description="Build index.html for HAPY cutout webpages.")
    parser.add_argument(
        "--runroot",
        required=True,
        help="Run directory containing cutouts/ and html/cutouts/"
    )
    args = parser.parse_args()

    runroot = Path(args.runroot).resolve()
    entries = collect_entries(runroot)

    outfile = runroot / "html" / "cutouts" / "index.html"
    write_index(entries, outfile)


if __name__ == "__main__":
    main()
