# src/qs/commands/auto_run.py
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import List, Optional

from qs.utils.logger import get_logger
from qs.utils.samples import discover_fastqs, paired_sample_ids
from qs.metadata.template import write_metadata_template
from qs.utils.manifest import generate_manifest
from qs.commands import denoise_runs as cmd_denoise
from qs.commands import classify_sweep as cmd_cls
from qs.qiime import commands as qiime
from qs.analysis.depth import choose_sampling_depth

LOG = get_logger("auto")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "auto-run", parents=[parent],
        help="End-to-end pipeline: (optional init) → import+cutadapt+DADA2 per run/group → merge → classify sweep → core metrics.",
    )
    # Init / inputs
    p.add_argument("--fastq-dir", type=Path, required=True, help="Top FASTQ directory; subfolders are runs.")
    p.add_argument("--project-dir", type=Path, required=True, help="Project directory (outputs).")
    p.add_argument("--metadata-file", type=Path, default=None,
                   help="Metadata TSV (with primers). If missing and --init-metadata is set, will be generated.")
    p.add_argument("--init-metadata", action="store_true",
                   help="Generate metadata and a top-level manifest, then stop so you can fill primers.")
    p.add_argument("--protected-cols", type=str, default="__f_primer,__r_primer",
                   help="Protected columns for metadata template (used with --init-metadata).")
    p.add_argument("--keep-illumina-suffix", action="store_true", help="Keep _S#/_L### tokens in SampleIDs (default: strip).")
    p.add_argument("--id-regex", type=str, default=None, help="Optional regex to derive SampleID (named 'id' or group 1).")

    # Denoise pipeline knobs (forwarded to denoise-runs)
    p.add_argument("--split-by-primer-group", action="store_true", help="Per primer group within each run (recommended for mixed loci).")
    p.add_argument("--discard-untrimmed", action="store_true", help="cutadapt: discard reads without primer.")
    p.add_argument("--no-indels", action="store_true", help="cutadapt: exact matches only.")
    p.add_argument("--cores", type=int, default=0, help="Cores for cutadapt & DADA2.")
    # Auto truncation (forwarded)
    p.add_argument("--auto-trunc", action="store_true", help="Auto optimize trunc-len-f/r per run/group.")
    p.add_argument("--trunc-len-f", type=int, default=0, help="Override trunc-len-f (disables auto if both F/R set).")
    p.add_argument("--trunc-len-r", type=int, default=0, help="Override trunc-len-r.")
    p.add_argument("--trunc-lower-frac", type=float, default=0.80)
    p.add_argument("--trunc-min-len", type=int, default=50)
    p.add_argument("--trunc-step", type=int, default=10)
    p.add_argument("--trunc-refine-step", type=int, default=5)
    p.add_argument("--trunc-quick-learn", type=int, default=250000)
    # DADA2 misc (forwarded)
    p.add_argument("--trim-left-f", type=int, default=0)
    p.add_argument("--trim-left-r", type=int, default=0)
    p.add_argument("--max-ee-f", type=int, default=2)
    p.add_argument("--max-ee-r", type=int, default=2)
    p.add_argument("--trunc-q", type=int, default=2)
    p.add_argument("--min-overlap", type=int, default=12)
    p.add_argument("--pooling-method", type=str, default="independent")
    p.add_argument("--chimera-method", type=str, default="consensus")

    # Classifier sweep
    p.add_argument("--classifiers-dir", type=Path, default=None, help="Directory of sklearn classifier .qza files (if omitted, skip sweep).")

    # Core metrics (phylogeny & diversity)
    p.add_argument("--retain-fraction", type=float, default=0.90, help="Fraction of samples to retain when choosing depth.")
    p.add_argument("--min-depth", type=int, default=1000, help="Minimum sampling depth.")
    # (Optional) beta/alpha tests could be added later

    p.set_defaults(func=run)


def _init_and_exit(
    *,
    fastq_dir: Path,
    project_dir: Path,
    metadata_file: Optional[Path],
    protected_cols: str,
    keep_illumina_suffix: bool,
    id_regex: Optional[str],
) -> None:
    meta_path = metadata_file or (project_dir / "metadata.tsv")
    sids = paired_sample_ids(
        discover_fastqs(fastq_dir),
        strip_illumina_suffix=not keep_illumina_suffix,
        id_regex=id_regex,
    )
    if not sids:
        print("error: no paired FASTQs discovered; cannot init metadata.", file=sys.stderr)
        sys.exit(2)
    write_metadata_template(sids, meta_path, [c.strip() for c in protected_cols.split(",") if c.strip()])
    top_manifest = project_dir / "fastq.manifest"
    generate_manifest(fastq_dir, top_manifest, strip_illumina_suffix=not keep_illumina_suffix, id_regex=id_regex)
    LOG.info("Initialized metadata → %s", meta_path)
    LOG.info("Top-level manifest → %s", top_manifest)
    print(f"[ok] Wrote metadata + manifest. Fill primers, then re-run without --init-metadata.")
    sys.exit(0)


def _stage_three(rep_seqs: Path, table: Path, taxonomy: Optional[Path], meta_aug: Path, out_dir: Path, *, dry_run: bool, show_qiime: bool, retain_fraction: float, min_depth: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Symlink or copy inputs
    for src, name in ((rep_seqs, "rep-seqs.qza"), (table, "table.qza")):
        dst = out_dir / name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        try:
            dst.symlink_to(src)
        except Exception:
            import shutil
            shutil.copy2(src, dst)
    if taxonomy:
        dst = out_dir / "taxonomy.qza"
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        try:
            dst.symlink_to(taxonomy)
        except Exception:
            import shutil
            shutil.copy2(taxonomy, dst)

    # Build phylogeny
    rooted = out_dir / "rooted-tree.qza"
    if not rooted.exists():
        qiime.phylogeny_align_to_tree_mafft_fasttree(
            input_sequences=out_dir / "rep-seqs.qza",
            output_alignment=out_dir / "aligned-rep-seqs.qza",
            output_masked_alignment=out_dir / "masked-aligned-rep-seqs.qza",
            output_tree=out_dir / "unrooted-tree.qza",
            output_rooted_tree=rooted,
            dry_run=dry_run,
            show_stdout=show_qiime,
        )

    # choose depth (expects table.qzv nearby — created during merge; if not there, create it)
    table_qzv = out_dir / "table.qzv"
    if not table_qzv.exists():
        qiime.feature_table_summarize(input_table=out_dir / "table.qza", output=table_qzv, sample_metadata_file=meta_aug,
                                      dry_run=dry_run, show_stdout=show_qiime)
    depth = choose_sampling_depth(table_qzv, retain_fraction=retain_fraction, min_depth=min_depth)

    core_dir = out_dir / "core-metrics-phylo"
    qiime.diversity_core_metrics_phylogenetic(
        input_phylogeny=rooted,
        input_table=out_dir / "table.qza",
        sampling_depth=depth,
        metadata_file=meta_aug,
        output_dir=core_dir,
        dry_run=dry_run,
        show_stdout=show_qiime,
    )
    LOG.info("Core metrics done (depth=%d) → %s", depth, core_dir)


def run(args) -> None:
    fastq_dir: Path = args.fastq_dir
    project_dir: Path = args.project_dir
    meta_path: Optional[Path] = args.metadata_file

    project_dir.mkdir(parents=True, exist_ok=True)

    # Step 0: optional init and stop
    if args.init_metadata:
        _init_and_exit(
            fastq_dir=fastq_dir,
            project_dir=project_dir,
            metadata_file=meta_path,
            protected_cols=args.protected_cols,
            keep_illumina_suffix=args.keep_illumina_suffix,
            id_regex=args.id_regex,
        )

    # Require metadata to proceed
    meta_path = meta_path or (project_dir / "metadata.tsv")
    if not meta_path.exists():
        print(f"error: metadata file not found: {meta_path} (use --init-metadata first)", file=sys.stderr)
        sys.exit(2)

    # Step 1-2: run the per-run(/group) pipeline
    ns = SimpleNamespace(
        fastq_dir=fastq_dir,
        project_dir=project_dir,
        metadata_file=meta_path,
        keep_illumina_suffix=args.keep_illumina_suffix,
        id_regex=args.id_regex,
        front_f_col="__f_primer",
        front_r_col="__r_primer",
        split_by_primer_group=args.split_by_primer_group,
        cores=args.cores,
        discard_untrimmed=args.discard_untrimmed,
        no_indels=args.no_indels,
        auto_trunc=args.auto_trunc,
        trunc_len_f=args.trunc_len_f,
        trunc_len_r=args.trunc_len_r,
        trunc_lower_frac=args.trunc_lower_frac,
        trunc_min_len=args.trunc_min_len,
        trunc_step=args.trunc_step,
        trunc_refine_step=args.trunc_refine_step,
        trunc_quick_learn=args.trunc_quick_learn,
        trim_left_f=args.trim_left_f,
        trim_left_r=args.trim_left_r,
        max_ee_f=args.max_ee_f,
        max_ee_r=args.max_ee_r,
        trunc_q=args.trunc_q,
        min_overlap=args.min_overlap,
        pooling_method=args.pooling_method,
        chimera_method=args.chimera_method,
        dry_run=getattr(args, "dry_run", False),
        show_qiime=getattr(args, "show_qiime", True),
    )
    cmd_denoise.run(ns)

    # Where are the merged outputs?
    merged_index_path = project_dir / "MERGED.json"
    if not merged_index_path.exists():
        print("error: expected MERGED.json not found; denoise phase failed?", file=sys.stderr)
        sys.exit(3)
    merged_index = json.loads(merged_index_path.read_text())

    # Step 3: classifier sweep (optional)
    winners: dict[str, dict] = {}
    if args.classifiers_dir and Path(args.classifiers_dir).is_dir():
        for group_key, rec in merged_index.items():
            group_dir = Path(rec["dir"])
            rep_seqs = group_dir / "rep-seqs.qza"
            cls_out = group_dir / "classifiers"
            res = cmd_cls.sweep_and_pick(
                input_reads=rep_seqs,
                classifiers_dir=Path(args.classifiers_dir),
                out_dir=cls_out,
                dry_run=getattr(args, "dry_run", False),
                show_qiime=getattr(args, "show_qiime", True),
            )
            winners[group_key] = res
            print(f"[ok] group={group_key} classifier winner={res['classifier_tag']} by {res['priority']}")
    else:
        LOG.info("No classifiers-dir provided; skipping sweep.")

    # Step 4: stage + phylogeny + core metrics per group
    meta_aug = project_dir / "metadata.augmented.tsv"
    for group_key, rec in merged_index.items():
        group_dir = Path(rec["dir"])
        table = group_dir / "table.qza"
        rep_seqs = group_dir / "rep-seqs.qza"
        taxonomy = None
        if winners.get(group_key):
            tag = winners[group_key]["classifier_tag"]
            taxonomy = group_dir / "classifiers" / f"{tag}_classification.qza"
        _stage_three(
            rep_seqs=rep_seqs,
            table=table,
            taxonomy=taxonomy,
            meta_aug=meta_aug,
            out_dir=group_dir,  # keep analysis beside the merged outputs
            dry_run=getattr(args, "dry_run", False),
            show_qiime=getattr(args, "show_qiime", True),
            retain_fraction=args.retain_fraction,
            min_depth=args.min_depth,
        )

    print(f"[ok] auto-run complete → {project_dir}")
