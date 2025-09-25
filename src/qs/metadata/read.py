# src/qs/metadata/read.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple


def _unquote(s: str) -> str:
    """
    Strip BOM, surrounding quotes, and outer whitespace from a single cell.
    """
    s = s.replace("\ufeff", "")
    s = s.strip()
    if len(s) >= 2 and ((s[0] == s[-1] == '"') or (s[0] == s[-1] == "'")):
        s = s[1:-1].strip()
    return s


# Header aliases that should normalize to '#SampleID'
_SAMPLE_ID_ALIASES = {
    "#sampleid", "#sample id", "sample id", "sample-id", "sampleid", "sample_name",
}


def _normalize_header(raw: List[str]) -> List[str]:
    """
    Normalize header cells (unquote + fix first column to '#SampleID' if it is an alias).
    """
    out: List[str] = []
    for i, cell in enumerate(raw):
        c = _unquote(cell)
        if i == 0:
            if c.lower() in _SAMPLE_ID_ALIASES or c in {"#SampleID", "#Sample ID"}:
                c = "#SampleID"
        out.append(c)
    return out


def load_metadata_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """
    Load a QIIME-style sample metadata TSV.

    Returns:
      header: List[str] (first element guaranteed to be '#SampleID')
      rows:   List[Dict[str, str]] (keys are header names)

    Skips an optional '#q2:types' second row if present.
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"Empty metadata file: {path}")

    raw_header = lines[0].split("\t")
    header = _normalize_header(raw_header)

    start = 1
    if len(lines) > 1:
        first_cell = _unquote(lines[1].split("\t", 1)[0])
        if first_cell.lower() == "#q2:types":
            start = 2

    rows: List[Dict[str, str]] = []
    for ln in lines[start:]:
        cols = [_unquote(c) for c in ln.split("\t")]
        row = {header[i]: (cols[i] if i < len(cols) else "") for i in range(len(header))}
        rows.append(row)
    return header, rows


def collect_unique_primers(
    rows: List[Dict[str, str]],
    f_col: str = "__f_primer",
    r_col: str = "__r_primer",
) -> Tuple[List[str], List[str]]:
    """
    Collect sorted unique forward and reverse primer sequences from metadata rows.
    Blank cells are ignored.
    """
    fset, rset = set(), set()
    for r in rows:
        f = (r.get(f_col) or "").strip()
        g = (r.get(r_col) or "").strip()
        if f:
            fset.add(f)
        if g:
            rset.add(g)
    return sorted(fset), sorted(rset)
