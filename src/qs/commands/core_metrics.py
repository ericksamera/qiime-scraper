# src/qs/commands/core_metrics.py
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Optional, Set

from qs.utils.logger import get_logger
from qs.analysis.metrics import stage_and_run_metrics_for_group, write_alpha_metrics_tsv
from qs.analysis.filtering import compute_min_freq_from_qzv, filter_table_and_seqs
from qs.qiime import commands as qiime
from qs.commands.common import load_merged_index, select_groups
from qs.analysis.stats import run_diversity_stats_for_group

LOG = get_logger("coremetrics")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "core-metrics", parents=[parent],
        help="Run filtering (optional) + phylogeny + core metrics for merged groups.",
        description=(
            "Reads MERGED.json; for each selected group can filter the table (rare features, contaminants, low-depth), "
            "build a phylogeny, choose a sampling depth (~retain_fraction), then run core metrics."
        ),
    )
    p.add_argument("--project-dir", type=Path, required=True, help="Project directory (must contain MERGED.json).")
    p.add_argument("--metadata-file", type=Path, default=None,
                   help="Metadata TSV for core metrics (default: project/metadata.augmented.tsv).")
    p.add_argument("--groups", type=str, default=None, help="Comma-separated group keys or slugs (default: all).")
    p.add_argument("--retain-fraction", type=float, default=0.90, help="Fraction of samples to retain when choosing depth.")
    p.add_argument("--min-depth", type=int, default=1000, help="Minimum sampling depth.")
    p.add_argument("--if-exists", choices=("skip", "overwrite", "error", "new"), default="skip",
                   help="How to handle existing core-metrics-phylo directory.")
    p.add_argument("--alpha-tsv", action="store_true", help="Also write <group>/core-metrics-phylo/alpha-metrics.tsv.")

    # Filtering toggles
    p.add_argument("--filter", dest="do_filter", action="store_true", help="Enable table filtering (default: on).")
    p.add_argument("--no-filter", dest="do_filter", action="store_false", help="Disable filtering.")
    p.set_defaults(do_filter=True)
    p.add_argument("--rare-freq-frac", type=float, default=0.001)
    p.add_argument("--rare-min-samples", type=int, default=1)
    p.add_argument("--contam-exclude", type=str, default="mitochondria,chloroplast")
    p.add_argument("--filter-unclassified-phylum", action="store_true")
    p.add_argument("--min-sample-depth", type=int, default=0)

    # Stats toggles
    p.add_argument("--run-stats", action="store_true", help="Run diversity stats after core metrics.")
    p.add_argument("--beta-group-cols", type=str, default=None, help="Comma-separated metadata columns for beta-group-significance.")
    p.add_argument("--beta-method", choices=("permanova", "anosim", "permdisp"), default="permanova")
    p.add_argument("--pairwise", action="store_true")
    p.add_argument("--permutations", type=int, default=999)
    p.add_argument("--adonis-formula", type=str, default=None)

    # Taxa barplot toggle
    p.add_argument("--taxa-barplot", dest="taxa_barplot", action="store_true", help="Render taxa-barplot.qzv if taxonomy exists.")
    p.add_argument("--no-taxa-barplot", dest="taxa_barplot", action="store_false")
    p.set_defaults(taxa_barplot=True)

    p.set_defaults(func=run)


def run(args) -> None:
    project_dir: Path = args.project_dir
    merged_index = load_merged_index(project_dir)
    groups = select_groups(merged_index, args.groups)
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
        table = group_dir / "table.qza"
        rep_seqs = group_dir / "rep-seqs.qza"
        taxonomy = group_dir / "taxonomy.qza" if (group_dir / "taxonomy.qza").exists() else None

        # ensure table.qzv for thresholding
        table_qzv = group_dir / "table.qzv"
        if not table_qzv.exists():
            qiime.feature_table_summarize(
                input_table=table,
                output=table_qzv,
                sample_metadata_file=meta_aug,
                dry_run=getattr(args, "dry_run", False),
                show_stdout=getattr(args, "show_qiime", True),
            )

        if args.do_filter:
            min_freq = compute_min_freq_from_qzv(table_qzv, args.rare_freq_frac, fallback_min=1)
            f = filter_table_and_seqs(
                group_dir=group_dir,
                table_qza=table,
                rep_seqs_qza=rep_seqs,
                table_qzv=table_qzv,
                taxonomy_qza=taxonomy,
                min_freq=min_freq,
                min_samples=args.rare_min_samples,
                contam_exclude=args.contam_exclude,
                keep_only_phylum_classified=args.filter_unclassified_phylum,
                min_sample_depth=args.min_sample_depth,
                metadata_augmented=meta_aug,
                dry_run=getattr(args, "dry_run", False),
                show_qiime=getattr(args, "show_qiime", True),
            )
            table = f["table_final"]
            rep_seqs = f["rep_seqs_final"]

        core_dir = stage_and_run_metrics_for_group(
            group_dir=group_dir,
            metadata_augmented=meta_aug,
            retain_fraction=args.retain_fraction,
            min_depth=args.min_depth,
            if_exists=args.if_exists,
            dry_run=getattr(args, "dry_run", False),
            show_qiime=getattr(args, "show_qiime", True),
            make_taxa_barplot=getattr(args, "taxa_barplot", True),
        )

        if args.run_stats:
            cols = [c.strip() for c in (args.beta_group_cols.split(",") if args.beta_group_cols else []) if c.strip()]
            run_diversity_stats_for_group(
                group_dir=group_dir,
                metadata_file=meta_aug,
                beta_group_cols=cols,
                beta_method=args.beta_method,
                pairwise=args.pairwise,
                permutations=args.permutations,
                adonis_formula=args.adonis_formula,
                run_alpha_correlation=True,
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
