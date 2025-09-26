# src/qs/commands/core_metrics.py
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, Optional, Set

from qs.utils.logger import get_logger
from qs.analysis.metrics import stage_and_run_metrics_for_group, write_alpha_metrics_tsv
from qs.utils.text import slugify

LOG = get_logger("coremetrics")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "core-metrics", parents=[parent],
        help="Run phylogeny + core metrics for merged groups (no denoise/classify).",
        description=(
            "Reads MERGED.json, then for each selected group ensures rep-seqs.qza and table.qza "
            "exist (auto-repair from per-run outputs if needed), builds a phylogeny, "
            "chooses a depth from table.qzv (~retain_fraction samples), and runs core metrics."
        ),
    )
    p.add_argument("--project-dir", type=Path, required=True, help="Project directory (must contain MERGED.json).")
    p.add_argument("--metadata-file", type=Path, default=None,
                   help="Metadata TSV for core metrics (default: project/metadata.augmented.tsv).")
    p.add_argument("--groups", type=str, default=None,
                   help="Comma-separated group keys or slugs to operate on (default: all).")
    p.add_argument("--retain-fraction", type=float, default=0.90,
                   help="Fraction of samples to retain when choosing depth (default 0.90).")
    p.add_argument("--min-depth", type=int, default=1000,
                   help="Minimum sampling depth (default 1000).")
    p.add_argument("--if-exists", choices=("skip", "overwrite", "error", "new"), default="skip",
                   help="How to handle existing core-metrics-phylo directory (default: skip).")
    p.add_argument("--alpha-tsv", action="store_true",
                   help="Also write <group>/core-metrics-phylo/alpha-metrics.tsv.")
    p.set_defaults(func=run)


def _load_merged_index(project_dir: Path) -> Dict[str, dict]:
    idx = project_dir / "MERGED.json"
    if not idx.exists():
        raise FileNotFoundError(f"Missing {idx}")
    return json.loads(idx.read_text())


def _select_groups(merged_index: Dict[str, dict], groups_arg: Optional[str]) -> Set[str]:
    if not groups_arg:
        return set(merged_index.keys())
    want = {g.strip() for g in groups_arg.split(",") if g.strip()}
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


def run(args) -> None:
    project_dir: Path = args.project_dir
    merged_index = _load_merged_index(project_dir)
    groups = _select_groups(merged_index, args.groups)
    if not groups:
        print("error: no matching groups found to operate on.", file=sys.stderr)
        sys.exit(2)

    meta_aug = args.metadata_file or (project_dir / "metadata.augmented.tsv")
    if not meta_aug.exists():
        print(f"error: metadata file not found: {meta_aug}", file=sys.stderr)
        sys.exit(3)

    for key in sorted(groups):
        rec = merged_index[key]
        group_dir = Path(rec["dir"])
        core_dir = stage_and_run_metrics_for_group(
            group_dir=group_dir,
            metadata_augmented=meta_aug,
            retain_fraction=args.retain_fraction,
            min_depth=args.min_depth,
            if_exists=args.if_exists,
            dry_run=getattr(args, "dry_run", False),
            show_qiime=getattr(args, "show_qiime", True),
        )
        if args.alpha_tsv:
            try:
                path = write_alpha_metrics_tsv(core_dir)
                LOG.info("Alpha metrics → %s", path)
            except Exception as e:
                LOG.warning("Alpha TSV export skipped for %s: %s", key, e)
        LOG.info("Group %s → %s", key, core_dir)
        print(f"[ok] core metrics: {group_dir}")
