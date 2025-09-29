# src/qs/plan/build.py
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from qs.plan.types import GroupPlan, Plan
from qs.utils.text import slugify


def _load_merged_index(project_dir: Path) -> Dict[str, dict]:
    idx = project_dir / "MERGED.json"
    if not idx.exists():
        raise FileNotFoundError(f"Missing {idx}")
    return json.loads(idx.read_text())


def _select_groups(merged_index: Dict[str, dict], wanted: Optional[Iterable[str]]) -> Set[str]:
    if not wanted:
        return set(merged_index.keys())
    want = {str(g).strip() for g in wanted if str(g).strip()}
    if not want:
        return set(merged_index.keys())
    slug_to_key = {slugify(k if k else "all"): k for k in merged_index}
    out: Set[str] = set()
    for g in want:
        if g in merged_index:
            out.add(g); continue
        if g in slug_to_key:
            out.add(slug_to_key[g])
    return out


def _resolve_classifier_source(
    *,
    group_key: str,
    group_slug: str,
    classifiers_map: Optional[Dict[str, str]],
    default_dir: Optional[Path],
) -> Optional[Path]:
    if classifiers_map:
        for k in (group_key, group_slug, "default"):
            v = classifiers_map.get(k)
            if v:
                return Path(v)
    return default_dir


def build_after_denoise(
    *,
    project_dir: Path,
    groups_requested: Optional[Iterable[str]],
    classifiers_map: Optional[Dict[str, str]],
    default_classifiers_dir: Optional[Path],
) -> Plan:
    mi = _load_merged_index(project_dir)
    selected = _select_groups(mi, groups_requested)
    if not selected:
        raise ValueError("no matching groups found to operate on")

    gps: List[GroupPlan] = []
    for key in sorted(selected):
        rec = mi[key]
        gdir = Path(rec["dir"])
        slug = gdir.name
        table = gdir / "table.qza"
        rep = gdir / "rep-seqs.qza"
        tax = gdir / "taxonomy.qza"
        tax = tax if tax.exists() else None
        src = _resolve_classifier_source(
            group_key=key,
            group_slug=slug,
            classifiers_map=classifiers_map,
            default_dir=default_classifiers_dir,
        )
        gps.append(GroupPlan(
            key=key,
            slug=slug,
            dir=gdir,
            table_qza=table,
            rep_seqs_qza=rep,
            taxonomy_qza=tax,
            classifier_source=src,
        ))
    return Plan(project_dir=project_dir, groups=gps, metadata_aug=project_dir / "metadata.augmented.tsv")
