# src/qs/metadata/read.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple


def load_metadata_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """
    Load a QIIME-style metadata TSV.

    Returns (header, rows) where rows are dicts keyed by header names.
    Skips an optional '#q2:types' second row.
    """
    lines = [ln.rstrip("\n") for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"Empty metadata file: {path}")
    header = [h.strip() for h in lines[0].split("\t")]
    start = 1
    if len(lines) > 1 and lines[1].split("\t", 1)[0].strip().lower() == "#q2:types":
        start = 2
    rows: List[Dict[str, str]] = []
    for ln in lines[start:]:
        cols = ln.split("\t")
        row = {h: (cols[i] if i < len(cols) else "") for i, h in enumerate(header)}
        rows.append(row)
    return header, rows


def collect_unique_primers(
    rows: List[Dict[str, str]],
    f_col: str = "__f_primer",
    r_col: str = "__r_primer",
) -> Tuple[List[str], List[str]]:
    """Return sorted unique forward and reverse primers (non-empty)."""
    fset, rset = set(), set()
    for r in rows:
        f = (r.get(f_col) or "").strip()
        g = (r.get(r_col) or "").strip()
        if f:
            fset.add(f)
        if g:
            rset.add(g)
    return sorted(fset), sorted(rset)
