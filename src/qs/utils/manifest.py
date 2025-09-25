# src/qs/utils/manifest.py
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from qs.utils.logger import get_logger
from qs.utils.samples import discover_fastqs, collect_pairs, normalize_id

LOG = get_logger("manifest")


def _rows_from_pairs(
    pairs: Dict[str, Dict[str, str]],
    *,
    strip_illumina_suffix: bool,
    id_regex: Optional[str],
    allowed_sample_ids: Optional[Set[str]] = None,
) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for root, d in sorted(pairs.items()):
        if "1" not in d or "2" not in d:
            continue
        sid = normalize_id(root, strip_illumina_suffix=strip_illumina_suffix, id_regex=id_regex)
        if allowed_sample_ids is not None and sid not in allowed_sample_ids:
            continue
        rows.append({
            "sample-id": sid,
            "forward-absolute-filepath": d["1"],
            "reverse-absolute-filepath": d["2"],
        })
    return rows


def generate_manifest(
    fastq_dir: Path,
    manifest_path: Path,
    *,
    strip_illumina_suffix: bool = True,
    id_regex: Optional[str] = None,
    allowed_sample_ids: Optional[Iterable[str]] = None,
) -> Path:
    paths = discover_fastqs(fastq_dir)
    pairs = collect_pairs(paths)
    allowed_set = set(allowed_sample_ids) if allowed_sample_ids is not None else None
    rows = _rows_from_pairs(
        pairs,
        strip_illumina_suffix=strip_illumina_suffix,
        id_regex=id_regex,
        allowed_sample_ids=allowed_set,
    )
    if not rows:
        raise RuntimeError(f"No paired R1/R2 FASTQs found under {fastq_dir}")

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sample-id",
                "forward-absolute-filepath",
                "reverse-absolute-filepath",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    LOG.info("Manifest written → %s (%d rows)", manifest_path, len(rows))
    return manifest_path


def read_manifest_sample_ids(manifest_path: Path) -> List[str]:
    with manifest_path.open("r", encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    header = [h.strip() for h in lines[0].split("\t")]
    try:
        idx = header.index("sample-id")
    except ValueError as e:
        raise RuntimeError("Manifest missing 'sample-id' header.") from e
    sample_ids = []
    for ln in lines[1:]:
        cols = ln.split("\t")
        sample_ids.append(cols[idx].strip())
    return sample_ids
