# src/qs/utils/runs.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List


def discover_runs(fastq_dir: Path) -> Dict[str, Path]:
    """
    Treat immediate subdirectories with FASTQs as distinct runs.
    If none found, use fastq_dir itself as a single 'run1'.
    """
    fastq_dir = fastq_dir.resolve()
    runs: Dict[str, Path] = {}
    for d in sorted(p for p in fastq_dir.iterdir() if p.is_dir()):
        if any(d.rglob("*.fastq*")):
            runs[d.name] = d
    if not runs:
        runs["run1"] = fastq_dir
    return runs
