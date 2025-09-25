# src/qs/utils/samples.py
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional

from qs.utils.logger import get_logger

LOG = get_logger("samples")

# Common Illumina-style paired-end filenames:
#   <sample>_S1_L001_R1_001.fastq.gz
#   <sample>_L001_R2_001.fastq.gz
#   <sample>_R1.fastq.gz
ILLUMINA_RE = re.compile(
    r"^(?P<root>.+?)_R(?P<read>[12])(?:_[0-9]{3})?\.fastq(?:\.gz)?$",
    re.IGNORECASE,
)
S_TOKEN_RE = re.compile(r"^S\d+$", re.IGNORECASE)       # e.g., S1, S18
L_TOKEN_RE = re.compile(r"^L\d{3}$", re.IGNORECASE)     # e.g., L001, L002


def discover_fastqs(fastq_dir: Path) -> List[Path]:
    if not fastq_dir.is_dir():
        raise NotADirectoryError(fastq_dir)
    paths = list(fastq_dir.glob("*.fastq*"))
    if not paths:
        LOG.info("No FASTQs at top level; searching recursively…")
        paths = list(fastq_dir.rglob("*.fastq*"))
    return paths


def collect_pairs(paths: List[Path]) -> Dict[str, Dict[str, str]]:
    """
    Map <root> (portion before _R1/_R2) -> {"1": path_to_R1, "2": path_to_R2}
    Multiple lanes appear as multiple distinct roots.
    """
    pairs: Dict[str, Dict[str, str]] = {}
    for fq in paths:
        m = ILLUMINA_RE.match(fq.name)
        if not m:
            continue
        root = m.group("root")
        read = m.group("read")
        pairs.setdefault(root, {})[read] = str(fq.resolve())
    return pairs


def normalize_id(
    root: str,
    *,
    strip_illumina_suffix: bool = True,
    id_regex: Optional[str] = None,
) -> str:
    """
    Turn an Illumina 'root' into a clean SampleID.

    Priority:
      1) If id_regex is given, use the first capturing group or the named 'id' group.
      2) Else, if strip_illumina_suffix is True (default), drop trailing _S<d+>/_L<ddd> tokens.
      3) Else, return root as-is.
    """
    if id_regex:
        rx = re.compile(id_regex)
        m = rx.search(root)
        if m:
            if "id" in m.groupdict():
                return m.group("id")
            if m.groups():
                return m.group(1)
        # if regex provided but no match, fall through to other logic

    if not strip_illumina_suffix:
        return root

    # remove S#/L### tokens anywhere in the root (typically at the end)
    tokens = root.split("_")
    filtered = [t for t in tokens if not (S_TOKEN_RE.fullmatch(t) or L_TOKEN_RE.fullmatch(t))]
    # preserve original underscores between remaining tokens
    return "_".join(filtered)


def paired_sample_ids(
    paths: List[Path],
    *,
    strip_illumina_suffix: bool = True,
    id_regex: Optional[str] = None,
) -> List[str]:
    """
    Return unique normalized sample IDs for which we observed both R1 and R2.
    """
    pairs = collect_pairs(paths)
    sids: set[str] = set()
    for root, reads in pairs.items():
        if {"1", "2"} <= set(reads):
            sid = normalize_id(root, strip_illumina_suffix=strip_illumina_suffix, id_regex=id_regex)
            sids.add(sid)
    return sorted(sids)
