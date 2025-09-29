# src/qs/plan/types.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


@dataclass(frozen=True)
class GroupPlan:
    key: str              # original group key "F|R" or "all"
    slug: str             # filesystem-safe slug
    dir: Path             # project/groups/<slug> (or project/ when all)
    table_qza: Path
    rep_seqs_qza: Path
    taxonomy_qza: Optional[Path]
    classifier_source: Optional[Path]  # .qza file OR a directory for sweep (resolved)


@dataclass(frozen=True)
class Plan:
    project_dir: Path
    groups: List[GroupPlan]
    metadata_aug: Path
