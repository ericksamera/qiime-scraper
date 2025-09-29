# src/qs/commands/common.py
from __future__ import annotations

import json
import sys
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Optional, Set

from qs.utils.text import slugify


# -------- Namespace & arg helpers --------------------------------------------

def to_ns(args=None, **kwargs) -> SimpleNamespace:
    """
    Accept argparse-style positional or kwargs and normalize to a namespace.
    IMPORTANT: if a namespace (argparse.Namespace or SimpleNamespace) was
    already provided, just return it unchanged.
    """
    if args is not None:
        return args  # keep whatever namespace we were given
    return SimpleNamespace(**kwargs)


def derive_default_project_dir(fastq_dir: Path) -> Path:
    fd = fastq_dir.resolve()
    return fd.parent / f"{fd.name}-qs"


def ensure_fastq_dir(ns: Namespace) -> Path:
    """
    Guarantee we have a Path for fastq_dir.
    - If present on ns, use it.
    - Else, if ./fastq exists, use that (common case).
    - Else, print a clear error and exit (rather than raising AttributeError).
    """
    v = getattr(ns, "fastq_dir", None)
    if isinstance(v, Path):
        return v
    if isinstance(v, str) and v:
        return Path(v)
    default = Path("fastq")
    if default.is_dir():
        return default
    print("error: --fastq-dir is required. Pass it explicitly (e.g., --fastq-dir ./fastq).", file=sys.stderr)
    sys.exit(2)


# -------- Project index & selection ------------------------------------------

def load_merged_index(project_dir: Path) -> Dict[str, dict]:
    path = project_dir / "MERGED.json"
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def select_groups(merged_index: Dict[str, dict], groups_arg: Optional[str]) -> Set[str]:
    if not groups_arg:
        return set(merged_index.keys())
    want = {g.strip() for g in groups_arg.split(",") if g.strip()}
    slug_to_key = {slugify(k if k else "all"): k for k in merged_index}
    out: Set[str] = set()
    for g in want:
        if g in merged_index:
            out.add(g)
        elif g in slug_to_key:
            out.add(slug_to_key[g])
    return out


def winner_if_present(group_dir: Path) -> Optional[str]:
    rpt = group_dir / "classifiers" / "optimal_classifiers.json"
    if not rpt.exists():
        return None
    try:
        data = json.loads(rpt.read_text())
        for k in ("pct_depth≥7", "median_conf", "mean_conf"):
            picks = data.get("winners", {}).get(k, [])
            if picks:
                return str(picks[0]["classifier"])
    except Exception:
        return None
    return None


def resolve_classifier_source(
    *,
    group_key: str,
    group_slug: str,
    mapping: Optional[Dict[str, str]],
    fallback_dir: Optional[Path],
) -> Optional[Path]:
    if mapping:
        for k in (group_key, group_slug, "default"):
            v = mapping.get(k)
            if v:
                return Path(v)
    return fallback_dir
