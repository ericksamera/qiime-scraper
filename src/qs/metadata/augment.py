# src/qs/metadata/augment.py
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple

from qs.metadata.read import load_metadata_table


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
    Add 'run' and 'primer_group' columns, normalize ID header to '#SampleID',
    and write a clean TSV (no quotes).
    Returns (new_header, groups mapping: group_key -> list of sample-ids).
    """
    header, rows = load_metadata_table(metadata_in)

    # Ensure first header is canonical '#SampleID'
    if header and header[0] != "#SampleID":
        header[0] = "#SampleID"

    # Ensure new columns exist (append if missing)
    new_header = list(header)
    if run_col_name not in new_header:
        new_header.append(run_col_name)
    if group_col_name not in new_header:
        new_header.append(group_col_name)

    groups: Dict[str, List[str]] = {}
    out_lines: List[List[str]] = [new_header]

    for r in rows:
        sid = (r.get("#SampleID") or "").lstrip("#").strip()

        # Compute primer group key
        f = (r.get(f_col) or "").strip()
        rev = (r.get(r_col) or "").strip()
        gkey = f"{f}|{rev}" if f and rev else ""

        groups.setdefault(gkey, []).append(sid)

        # Build row respecting new_header order
        line_map = dict(r)  # copy original cells
        line_map[run_col_name] = run_by_sample_id.get(sid, "")
        line_map[group_col_name] = gkey

        out_lines.append([line_map.get(col, "") for col in new_header])

    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    with metadata_out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        for row in out_lines:
            w.writerow(row)

    return new_header, groups
