# src/qs/metadata/augment.py
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from qs.metadata.read import load_metadata_table


def _has_types_row(lines: List[List[str]]) -> bool:
    return bool(lines) and lines[0] and lines[0][0].strip().lower() == "#q2:types"


def augment_metadata_with_runs_and_groups(
    metadata_in: Path,
    metadata_out: Path,
    *,
    run_by_sample_id: Dict[str, str],
    f_col: str = "__f_primer",
    r_col: str = "__r_primer",
    group_col_name: str = "primer_group",
    run_col_name: str = "run",
) -> Tuple[List[str], Dict[str, List[str]]]:
    """
    Add 'run' and 'primer_group' columns. Primer group is f"{F}|{R}" (blank if either missing).
    Returns (header, groups mapping: group_key -> list of sample-ids).
    """
    header, rows = load_metadata_table(metadata_in)

    # Build new header and possibly types row
    has_types = False
    raw_lines = [ln.rstrip("\n") for ln in metadata_in.read_text(encoding="utf-8").splitlines() if ln.strip()]
    matrix: List[List[str]] = [[c.strip() for c in raw_lines[0].split("\t")]]
    if len(raw_lines) > 1 and raw_lines[1].split("\t", 1)[0].strip().lower() == "#q2:types":
        matrix.append([c.strip() for c in raw_lines[1].split("\t")])
        has_types = True

    # Ensure columns exist in header
    new_header = header[:]
    if run_col_name not in new_header:
        new_header.append(run_col_name)
    if group_col_name not in new_header:
        new_header.append(group_col_name)

    # Write header (and types if present)
    out_rows: List[List[str]] = [new_header]
    if has_types:
        types = matrix[1][:]
        # extend types with 'categorical' for any new columns
        while len(types) < len(new_header):
            types.append("categorical")
        out_rows.append(types)

    # Fill rows
    groups: Dict[str, List[str]] = {}
    for r in rows:
        sid = (r.get("#SampleID") or r.get("sample-id") or r.get("SampleID") or r.get("#Sample ID") or "").lstrip("#")
        run = run_by_sample_id.get(sid, "")
        f = (r.get(f_col) or "").strip()
        rev = (r.get(r_col) or "").strip()
        group_key = f"{f}|{rev}" if f and rev else ""
        groups.setdefault(group_key, []).append(sid)

        line = [r.get(h, "") for h in header]
        if run_col_name not in header:
            line.append(run)
        else:
            line[header.index(run_col_name)] = run
        if group_col_name not in header:
            line.append(group_key)
        else:
            line[header.index(group_col_name)] = group_key
        out_rows.append(line)

    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    with metadata_out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for row in out_rows:
            w.writerow(row)
    return new_header, groups
