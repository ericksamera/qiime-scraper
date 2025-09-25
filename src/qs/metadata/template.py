# src/qs/metadata/template.py
from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, List

from qs.metadata import PROTECTED_PREFIX, DEFAULT_PROTECTED_COLS


def normalize_protected(cols: Iterable[str]) -> List[str]:
    out: List[str] = []
    for c in cols:
        c = c.strip()
        if not c:
            continue
        if not c.startswith(PROTECTED_PREFIX):
            c = PROTECTED_PREFIX + c
        out.append(c)
    return out


def write_metadata_template(
    sample_ids: Iterable[str],
    output_file: Path,
    protected_cols: Iterable[str] = DEFAULT_PROTECTED_COLS,
    *,
    include_types_row: bool = True,
) -> None:
    protected = normalize_protected(protected_cols)
    header = ["#SampleID", *protected]
    rows: List[List[str]] = []

    if include_types_row:
        types = ["#q2:types", *["categorical"] * len(protected)]
        rows.append(types)

    for sid in sample_ids:
        rows.append([sid, *[""] * len(protected)])

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
