# hapy/hatools/results.py
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable

from astropy.table import Table, vstack


def write_result_row_ecsv(path: str | Path, row: Dict[str, Any], overwrite: bool = True) -> Path:
    """
    Write a single-row ECSV results file.
    Intended for per-galaxy outputs: cutouts/<tag>/<tag>-results.ecsv
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Normalize to 1-row table
    t = Table({k: [v] for k, v in row.items()})
    t.write(path, format="ascii.ecsv", overwrite=overwrite)
    return path


def merge_result_rows_ecsv(paths: Iterable[str | Path], outpath: str | Path, overwrite: bool = True) -> Path:
    """
    Merge many single-row ECSV files into one ECSV table.
    """
    tables = []
    for p in paths:
        p = Path(p)
        tables.append(Table.read(p, format="ascii.ecsv"))

    if len(tables) == 0:
        raise ValueError("No result files provided to merge.")

    merged = vstack(tables, metadata_conflicts="silent")
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    merged.write(outpath, format="ascii.ecsv", overwrite=overwrite)
    return outpath
