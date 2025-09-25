# src/qs/commands/auto_run.py
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Set

from qs.utils.logger import get_logger
from qs.utils.samples import discover_fastqs, paired_sample_ids
from qs.metadata.template import write_metadata_template
from qs.utils.manifest import generate_manifest
from qs.commands import denoise_runs as cmd_denoise
from qs.commands import classify_sweep as cmd_cls
from qs.qiime import commands as qiime
from qs.analysis.depth import choose_sampling_depth
from qs.utils.text import slugify

LOG = get_logger("auto")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "auto-run", parents=[parent],
        help=(
            "End-to-end pipeline with resume: "
            "(optional init) → per-run/group import+cutadapt+DADA2(+auto trunc) → merge "
            "→ classifier sweep → core metrics."
        ),
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

    # Stage skipping / resume / selection
    p.add_argument("--skip-denoise", action="store_true",
                   help="Skip import/cutadapt/DADA2/merge and reuse existing MERGED.json.")
    p.add_argument("--skip-classify", action="store_true",
                   help="Skip classifier sweep.")
    p.add_argument("--skip-metrics", action="store_true",
                   help="Skip phylogeny + core metrics.")
    p.add_argument("--resume", action="store_true",
                   help="Auto-skip stages with completed outputs (idempotent resume).")
    p.add_argument("--groups", type=str, default=None,
                   help="Comma-separated group keys or slugs to operate on (default: all).")

    # Denoise pipeline knobs (forwarded to denoise-runs)
    p.add_argument("--split-by-primer-group", action="store_true",
                   help="Per primer group within each run (recommended for mixed loci).")
    p.add_argument("--discard-untrimmed", action="store_true",
                   help="cutadapt: discard reads without primer.")
    p.add_argument("--no-indels", action="store_true",
                   help="cutadapt: exact matches only.")
    p.add_argument("--cores", type=int, default=0, help="Cores for cutadapt & DADA2.")

    # Auto truncation (forwarded)
    p.add_argument("--auto-trunc", action="store_true",
                   help="Auto optimize trunc-len-f/r per run/group.")
    p.add_argument("--trunc-len-f", type=int, default=0,
                   help="Override trunc-len-f (disables auto if both F/R set).")
    p.add_argument("--trunc-len-r", type=int, default=0,
                   help="Override trunc-len-r.")
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
    p.add_argument("--classifiers-dir", type=Path, default=None,
                   help="Directory of sklearn classifier .qza files (if omitted, skip sweep).")

    # Core metrics (phylogeny & diversity)
    p.add_argument("--retain-fraction", type=float, default=0.90,
                   help="Fraction of samples to retain when choosing depth.")
    p.add_argument("--min-depth", type=int, default=1000,
                   help="Minimum sampling depth.")

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


def _stage_three(
    rep_seqs: Path,
    table: Path,
    taxonomy: Optional[Path],
    meta_aug: Path,
    out_dir: Path,
    *,
    dry_run: bool,
    show_qiime: bool,
    retain_fraction: float,
    min_depth: int
) -> None:
    """
    Prepare analysis in out_dir and run phylogeny + core metrics.

    - Avoids self-referential links
    - Auto-repairs missing/broken rep-seqs.qza and table.qza from per-run outputs
    - Resume-safe: keeps existing good outputs
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    def _is_broken_symlink(p: Path) -> bool:
        return p.is_symlink() and not p.exists()

    def _stage_artifact(src: Path, dest: Path) -> None:
        if str(src) == str(dest):
            # already in place
            return
        if dest.exists() or dest.is_symlink():
            try:
                if dest.resolve() == src.resolve():
                    return
            except Exception:
                pass
            try:
                dest.unlink()
            except FileNotFoundError:
                pass
        # prefer hardlink → symlink → copy
        try:
            os.link(src, dest)
            return
        except Exception:
            try:
                dest.symlink_to(src)
                return
            except Exception:
                import shutil
                shutil.copy2(src, dest)

    # locate project root and runs root
    group_slug = out_dir.name
    project_dir = out_dir.parent.parent  # .../project/groups/<slug> → project/
    runs_root = project_dir / "runs"

    # --- repair rep-seqs if missing/broken ---
    if (not rep_seqs.exists()) or _is_broken_symlink(rep_seqs):
        per_run_rep = sorted(runs_root.glob(f"*/{group_slug}/rep-seqs.qza"))
        if not per_run_rep:
            raise FileNotFoundError(
                f"rep-seqs.qza missing and no per-run rep-seqs found under {runs_root}/ */{group_slug}/"
            )
        if len(per_run_rep) == 1:
            import shutil
            rep_seqs.parent.mkdir(parents=True, exist_ok=True)
            try:
                rep_seqs.unlink()
            except Exception:
                pass
            shutil.copy2(per_run_rep[0], rep_seqs)
            LOG.warning("Repaired rep-seqs.qza by copying from %s", per_run_rep[0])
        else:
            tmp_merged = out_dir / "_rep-seqs.repaired.qza"
            qiime.merge_seqs(per_run_rep, tmp_merged, dry_run=dry_run, show_stdout=show_qiime)
            try:
                rep_seqs.unlink()
            except Exception:
                pass
            os.replace(tmp_merged, rep_seqs)
            LOG.warning("Repaired rep-seqs.qza by merging %d per-run files", len(per_run_rep))

    # --- repair table if missing/broken ---
    if (not table.exists()) or _is_broken_symlink(table):
        # most common (split-by-primer-group): tables are in runs/*/<group>/table.qza
        per_run_tbl = sorted(runs_root.glob(f"*/{group_slug}/table.qza"))
        # fallback: if user ran without split and group is 'all'
        if not per_run_tbl and group_slug == "all":
            per_run_tbl = sorted(runs_root.glob("*/table.qza"))
        if not per_run_tbl:
            raise FileNotFoundError(
                f"table.qza missing and no per-run tables found under {runs_root}/ */{group_slug}/"
            )
        if len(per_run_tbl) == 1:
            import shutil
            table.parent.mkdir(parents=True, exist_ok=True)
            try:
                table.unlink()
            except Exception:
                pass
            shutil.copy2(per_run_tbl[0], table)
            LOG.warning("Repaired table.qza by copying from %s", per_run_tbl[0])
        else:
            tmp_tbl = out_dir / "_table.repaired.qza"
            qiime.merge_tables(per_run_tbl, tmp_tbl, dry_run=dry_run, show_stdout=show_qiime)
            try:
                table.unlink()
            except Exception:
                pass
            os.replace(tmp_tbl, table)
            LOG.warning("Repaired table.qza by merging %d per-run files", len(per_run_tbl))

    # stage artifacts into out_dir if their source is elsewhere
    rep_dest = out_dir / "rep-seqs.qza"
    tbl_dest = out_dir / "table.qza"
    if rep_seqs.exists() and not _is_broken_symlink(rep_seqs):
        _stage_artifact(rep_seqs, rep_dest)
    if table.exists() and not _is_broken_symlink(table):
        _stage_artifact(table, tbl_dest)
    if taxonomy and taxonomy.exists():
        _stage_artifact(taxonomy, out_dir / "taxonomy.qza")

    # build phylogeny (only if needed)
    rooted = out_dir / "rooted-tree.qza"
    if not rooted.exists():
        qiime.phylogeny_align_to_tree_mafft_fasttree(
            input_sequences=rep_dest,
            output_alignment=out_dir / "aligned-rep-seqs.qza",
            output_masked_alignment=out_dir / "masked-aligned-rep-seqs.qza",
            output_tree=out_dir / "unrooted-tree.qza",
            output_rooted_tree=rooted,
            dry_run=dry_run,
            show_stdout=show_qiime,
        )

    # ensure we have a table summary to decide depth
    table_qzv = out_dir / "table.qzv"
    if not table_qzv.exists():
        qiime.feature_table_summarize(
            input_table=tbl_dest,
            output=table_qzv,
            sample_metadata_file=meta_aug,
            dry_run=dry_run,
            show_stdout=show_qiime,
        )
    depth = choose_sampling_depth(table_qzv, retain_fraction=retain_fraction, min_depth=min_depth)

    core_dir = out_dir / "core-metrics-phylo"
    if core_dir.exists():
        LOG.info("Core metrics already present → %s (resume)", core_dir)
        return
    qiime.diversity_core_metrics_phylogenetic(
        input_phylogeny=rooted,
        input_table=tbl_dest,
        sampling_depth=depth,
        metadata_file=meta_aug,
        output_dir=core_dir,
        dry_run=dry_run,
        show_stdout=show_qiime,
    )
    LOG.info("Core metrics done (depth=%d) → %s", depth, core_dir)


def _load_merged_index(project_dir: Path) -> Dict[str, dict]:
    idx_path = project_dir / "MERGED.json"
    if not idx_path.exists():
        raise FileNotFoundError(f"Missing {idx_path}")
    return json.loads(idx_path.read_text())


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
    if not out:
        LOG.warning("None of the requested groups matched: %s ; available: %s", sorted(want), sorted(merged_index.keys()))
        return set()
    return out


def _winner_if_present(group_dir: Path) -> Optional[str]:
    report = group_dir / "classifiers" / "optimal_classifiers.json"
    if not report.exists():
        return None
    try:
        data = json.loads(report.read_text())
        for key in ("pct_depth≥7", "median_conf", "mean_conf"):
            winners = data.get("winners", {}).get(key, [])
            if winners:
                return str(winners[0]["classifier"])
    except Exception:
        return None
    return None


def run(args) -> None:
    fastq_dir: Path = args.fastq_dir
    project_dir: Path = args.project_dir
    meta_path: Optional[Path] = args.metadata_file

    project_dir.mkdir(parents=True, exist_ok=True)

    if args.init_metadata:
        _init_and_exit(
            fastq_dir=fastq_dir,
            project_dir=project_dir,
            metadata_file=meta_path,
            protected_cols=args.protected_cols,
            keep_illumina_suffix=args.keep_illumina_suffix,
            id_regex=args.id_regex,
        )

    meta_path = meta_path or (project_dir / "metadata.tsv")
    if not meta_path.exists():
        print(f"error: metadata file not found: {meta_path} (use --init-metadata first)", file=sys.stderr)
        sys.exit(2)

    # Stage 1-2: denoise + merge unless skipped / already present
    merged_index: Dict[str, dict]
    merged_index_path = project_dir / "MERGED.json"

    if args.resume and merged_index_path.exists():
        LOG.info("Resume: MERGED.json found, skipping denoise/merge.")
        args.skip_denoise = True

    if not args.skip_denoise:
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

    try:
        merged_index = _load_merged_index(project_dir)
    except FileNotFoundError:
        print("error: expected MERGED.json not found; denoise phase required or use --skip-denoise only when MERGED.json exists.",
              file=sys.stderr)
        sys.exit(3)

    selected = _select_groups(merged_index, args.groups)
    if not selected:
        print("error: no matching groups found to operate on.", file=sys.stderr)
        sys.exit(4)

    # Stage 3: classifier sweep
    winners: dict[str, dict] = {}
    classifiers_ok = (args.classifiers_dir and Path(args.classifiers_dir).is_dir())

    for group_key in sorted(selected):
        rec = merged_index[group_key]
        group_dir = Path(rec["dir"])
        rep_seqs = group_dir / "rep-seqs.qza"

        if args.skip_classify or not classifiers_ok:
            LOG.info("Skipping classification for group=%s", group_key)
            continue

        if args.resume:
            prior = _winner_if_present(group_dir)
            if prior:
                LOG.info("Resume: winner already present for group=%s (%s); skipping sweep.", group_key, prior)
                winners[group_key] = {"classifier_tag": prior, "priority": "resume"}
                continue

        cls_out = group_dir / "classifiers"
        res = cmd_cls.sweep_and_pick(
            input_reads=rep_seqs,
            classifiers_dir=Path(args.classifiers_dir),  # type: ignore[arg-type]
            out_dir=cls_out,
            dry_run=getattr(args, "dry_run", False),
            show_qiime=getattr(args, "show_qiime", True),
        )
        winners[group_key] = res
        print(f"[ok] group={group_key} classifier winner={res['classifier_tag']} by {res['priority']}")

    if args.skip_metrics:
        LOG.info("Skipping core metrics per user request.")
        print(f"[ok] auto-run (skipped metrics) → {project_dir}")
        return

    # Stage 4: phylogeny + core metrics (resume-aware, and robust staging)
    meta_aug = project_dir / "metadata.augmented.tsv"
    for group_key in sorted(selected):
        rec = merged_index[group_key]
        group_dir = Path(rec["dir"])
        table = group_dir / "table.qza"
        rep_seqs = group_dir / "rep-seqs.qza"

        taxonomy = None
        if (group_key in winners) and winners[group_key].get("classifier_tag"):
            tag = winners[group_key]["classifier_tag"]
            taxonomy = group_dir / "classifiers" / f"{tag}_classification.qza"

        if args.resume and (group_dir / "core-metrics-phylo").exists():
            LOG.info("Resume: core-metrics present for group=%s; skipping.", group_key)
            continue

        _stage_three(
            rep_seqs=rep_seqs,
            table=table,
            taxonomy=taxonomy,
            meta_aug=meta_aug,
            out_dir=group_dir,
            dry_run=getattr(args, "dry_run", False),
            show_qiime=getattr(args, "show_qiime", True),
            retain_fraction=args.retain_fraction,
            min_depth=args.min_depth,
        )

    print(f"[ok] auto-run complete → {project_dir}")
