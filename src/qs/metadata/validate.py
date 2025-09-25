# src/qs/metadata/validate.py
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

from qs.metadata import PROTECTED_PREFIX, DEFAULT_PROTECTED_COLS


class MetadataError(Exception):
    pass


def _read_tsv(path: Path) -> Tuple[List[str], List[List[str]]]:
    lines = [ln.rstrip("\n") for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise MetadataError("Empty metadata file.")
    header = [h.strip() for h in lines[0].split("\t")]
    rows = [[c for c in ln.split("\t")] for ln in lines[1:]]
    return header, rows


def _is_types_row(row: List[str]) -> bool:
    return bool(row) and row[0].strip().lower() == "#q2:types"


def validate_metadata_file(
    metadata_path: Path,
    *,
    protected_prefix: str = PROTECTED_PREFIX,
    required_protected: Optional[List[str]] = None,
    require_types_row: bool = True,
    against_sample_ids: Optional[List[str]] = None,
) -> List[str]:
    """
    Return a list of warnings. Raise MetadataError for hard failures.
    """
    if required_protected is None:
        required_protected = DEFAULT_PROTECTED_COLS

    header, rows = _read_tsv(metadata_path)
    if header[0] not in {"#SampleID", "#Sample ID", "sample-id", "SampleID"}:
        raise MetadataError("First column must be '#SampleID' (or 'sample-id').")

    warnings: List[str] = []
    # Detect types row
    types_row_idx = 0
    has_types = False
    if rows and _is_types_row(rows[0]):
        has_types = True
        types_row_idx = 1
    elif require_types_row:
        warnings.append("No '#q2:types' row found; recommended for QIIME2.")

    # Protected columns
    protected_cols = [c for c in header if c.startswith(protected_prefix)]
    if required_protected:
        missing = [c for c in required_protected if c not in protected_cols]
        if missing:
            raise MetadataError(f"Missing required protected columns: {', '.join(missing)}")

    # Duplicate SampleIDs
    sample_ids = []
    for r in rows[types_row_idx:]:
        if not r:
            continue
        sample_ids.append(r[0].lstrip("#").strip())
    if len(sample_ids) != len(set(sample_ids)):
        raise MetadataError("Duplicate SampleIDs detected.")

    # Optional cross-check against expected IDs
    if against_sample_ids is not None:
        s_meta = set(sample_ids)
        s_ref = set(against_sample_ids)
        only_in_meta = sorted(s_meta - s_ref)
        only_in_ref = sorted(s_ref - s_meta)
        if only_in_meta or only_in_ref:
            msg = []
            if only_in_meta:
                msg.append(f"IDs only in metadata: {only_in_meta[:5]}{'...' if len(only_in_meta) > 5 else ''}")
            if only_in_ref:
                msg.append(f"IDs missing from metadata: {only_in_ref[:5]}{'...' if len(only_in_ref) > 5 else ''}")
            raise MetadataError("; ".join(msg))

    return warnings
