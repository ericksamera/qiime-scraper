# src/qs/commands/denoise_runs.py
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from qs.utils.logger import get_logger
from qs.utils.runs import discover_runs
from qs.utils.text import slugify
from qs.utils.manifest import generate_manifest, read_manifest_sample_ids
from qs.metadata.read import load_metadata_table
from qs.metadata.augment import augment_metadata_with_runs_and_groups
from qs.qiime import commands as qiime
from qs.optimize.truncation import find_optimal_truncation  # <— NEW

LOG = get_logger("denoise")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "denoise-runs", parents=[parent],
        help="Per-run (and per primer-group) import → cutadapt → DADA2 → merge.",
        description=(
            "Detect runs from subfolders of --fastq-dir. Add a 'run' column to metadata. "
            "Optionally split each run by primer pair (from metadata columns). "
            "Run cutadapt and DADA2 per run(/group), then merge across runs per primer-group."
        ),
    )
    p.add_argument("--fastq-dir", type=Path, required=True, help="Top directory of FASTQs. Subfolders are treated as runs.")
    p.add_argument("--project-dir", type=Path, required=True, help="Project output directory.")
    p.add_argument("--metadata-file", type=Path, required=True, help="Metadata TSV (with primer columns).")

    # ID normalization
    p.add_argument("--keep-illumina-suffix", action="store_true", help="Keep _S#/_L### tokens (default: strip).")
    p.add_argument("--id-regex", type=str, default=None, help="Optional regex to extract SampleID (named 'id' or group 1).")

    # Primer columns & grouping
    p.add_argument("--front-f-col", type=str, default="__f_primer", help="Metadata column with forward primers.")
    p.add_argument("--front-r-col", type=str, default="__r_primer", help="Metadata column with reverse primers.")
    p.add_argument("--split-by-primer-group", action="store_true",
                   help="Process each primer pair as its own sub-pipeline per run (recommended for mixed loci).")

    # cutadapt knobs
    p.add_argument("--cores", type=int, default=0, help="Cores for cutadapt and DADA2.")
    p.add_argument("--discard-untrimmed", action="store_true", help="cutadapt: --p-discard-untrimmed")
    p.add_argument("--no-indels", action="store_true", help="cutadapt: --p-no-indels")

    # DADA2 fixed or auto truncation
    p.add_argument("--trunc-len-f", type=int, default=0, help="DADA2 --p-trunc-len-f (override auto if set).")
    p.add_argument("--trunc-len-r", type=int, default=0, help="DADA2 --p-trunc-len-r (override auto if set).")
    p.add_argument("--auto-trunc", action="store_true",
                   help="Automatically choose trunc-len-f/r via coarse-to-fine optimization.")
    p.add_argument("--trunc-lower-frac", type=float, default=0.80,
                   help="Lower bound as fraction of read length (default 0.80).")
    p.add_argument("--trunc-min-len", type=int, default=50, help="Absolute minimum trunc length bound.")
    p.add_argument("--trunc-step", type=int, default=10, help="Coarse step size in bp (default 10).")
    p.add_argument("--trunc-refine-step", type=int, default=5, help="Refine step size in bp (default 5).")
    p.add_argument("--trunc-quick-learn", type=int, default=250000,
                   help="n-reads-learn for optimization trials (speeds up search).")

    # DADA2 other knobs
    p.add_argument("--trim-left-f", type=int, default=0, help="DADA2 --p-trim-left-f")
    p.add_argument("--trim-left-r", type=int, default=0, help="DADA2 --p-trim-left-r")
    p.add_argument("--max-ee-f", type=int, default=2, help="DADA2 --p-max-ee-f")
    p.add_argument("--max-ee-r", type=int, default=2, help="DADA2 --p-max-ee-r")
    p.add_argument("--trunc-q", type=int, default=2, help="DADA2 --p-trunc-q")
    p.add_argument("--min-overlap", type=int, default=12, help="DADA2 --p-min-overlap")
    p.add_argument("--pooling-method", type=str, default="independent", help="DADA2 --p-pooling-method")
    p.add_argument("--chimera-method", type=str, default="consensus", help="DADA2 --p-chimera-method")

    p.set_defaults(func=run)


def _group_key(f: str, r: str) -> str:
    return f"{f}|{r}" if f and r else ""


def _collect_group_primers(rows: List[dict], f_col: str, r_col: str) -> Dict[str, Tuple[List[str], List[str]]]:
    by_group: Dict[str, Tuple[set, set]] = {}
    for r in rows:
        f = (r.get(f_col) or "").strip()
        rv = (r.get(r_col) or "").strip()
        key = _group_key(f, rv)
        by_group.setdefault(key, (set(), set()))
        if f:
            by_group[key][0].add(f)
        if rv:
            by_group[key][1].add(rv)
    return {k: (sorted(v[0]), sorted(v[1])) for k, v in by_group.items() if k}


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _auto_or_fixed_trunc(
    *,
    run_path: Path,
    sub_root: Path,
    input_for_dada2: Path,
    args,
) -> Tuple[int, int, Dict[str, object]]:
    """
    Decide trunc-len-f/r. If --auto-trunc, run optimizer; else use provided values.
    Returns (f, r, opt_result_dict_or_empty).
    """
    if not args.auto_trunc and (args.trunc_len_f > 0 and args.trunc_len_r > 0):
        return args.trunc_len_f, args.trunc_len_r, {}

    # Auto mode
    opt_dir = sub_root / "optimize_trunc"
    result = find_optimal_truncation(
        run_path=run_path,
        input_seqs=input_for_dada2,
        outdir=opt_dir,
        lower_frac=args.trunc_lower_frac,
        abs_min_len=args.trunc_min_len,
        coarse_step=args.trunc_step,
        refine_step=args.trunc_refine_step,
        cores=args.cores,
        quick_learn=args.trunc_quick_learn,
        dry_run=getattr(args, "dry_run", False),
        show_qiime=getattr(args, "show_qiime", True),
        trim_left_f=args.trim_left_f,
        trim_left_r=args.trim_left_r,
        max_ee_f=args.max_ee_f,
        max_ee_r=args.max_ee_r,
        trunc_q=args.trunc_q,
        min_overlap=args.min_overlap,
        pooling_method=args.pooling_method,
        chimera_method=args.chimera_method,
    )
    best = result["best"]  # type: ignore[index]
    return int(best["trunc_len_f"]), int(best["trunc_len_r"]), result  # type: ignore[index]


def run(args) -> None:
    fastq_dir: Path = args.fastq_dir
    project_dir: Path = args.project_dir
    meta_in: Path = args.metadata_file
    project_dir.mkdir(parents=True, exist_ok=True)

    runs = discover_runs(fastq_dir)
    if not runs:
        print(f"error: no runs detected under {fastq_dir}", file=sys.stderr)
        sys.exit(2)
    LOG.info("Detected %d run(s): %s", len(runs), ", ".join(runs.keys()))

    # Per-run manifests to get normalized IDs
    run_sid_map: Dict[str, str] = {}
    for run_id, run_path in runs.items():
        run_root = project_dir / "runs" / run_id
        _ensure_dir(run_root)
        manifest_path = run_root / "fastq.manifest"
        generate_manifest(
            run_path,
            manifest_path,
            strip_illumina_suffix=not args.keep_illumina_suffix,
            id_regex=args.id_regex,
        )
        for sid in read_manifest_sample_ids(manifest_path):
            run_sid_map[sid] = run_id

    meta_aug = project_dir / "metadata.augmented.tsv"
    header, groups_map = augment_metadata_with_runs_and_groups(
        meta_in,
        meta_aug,
        run_by_sample_id=run_sid_map,
        f_col=args.front_f_col,
        r_col=args.front_r_col,
    )
    LOG.info("Wrote augmented metadata → %s", meta_aug)

    per_group_tables: Dict[str, List[Path]] = {}
    per_group_seqs: Dict[str, List[Path]] = {}
    _, all_rows = load_metadata_table(meta_aug)

    for run_id, run_path in runs.items():
        run_root = project_dir / "runs" / run_id
        manifest_path = run_root / "fastq.manifest"
        run_sids = set(read_manifest_sample_ids(manifest_path))
        run_rows = [r for r in all_rows if (r.get("#SampleID") or r.get("sample-id") or "").lstrip("#") in run_sids]

        if args.split_by_primer_group:
            group_primers = _collect_group_primers(run_rows, args.front_f_col, args.front_r_col)
            if not group_primers:
                LOG.warning("Run %s: no primer groups; processing all together without trimming.", run_id)
                group_primers = {}

            for gkey, (fwd_list, rev_list) in (group_primers.items() or {"all": ([], [])}.items()):
                gslug = slugify(gkey if gkey != "all" else "all")
                sub_root = run_root / gslug
                _ensure_dir(sub_root)

                # manifest limited to this group's SIDs
                group_sids = [r.get("#SampleID") or r.get("sample-id") for r in run_rows
                              if f"{(r.get(args.front_f_col) or '').strip()}|{(r.get(args.front_r_col) or '').strip()}" == gkey]
                group_sids = [s.lstrip("#") for s in group_sids if s]

                sub_manifest = sub_root / "fastq.manifest"
                generate_manifest(
                    run_path,
                    sub_manifest,
                    strip_illumina_suffix=not args.keep_illumina_suffix,
                    id_regex=args.id_regex,
                    allowed_sample_ids=group_sids if gkey != "all" else None,
                )

                demux_qza = sub_root / "demux.qza"
                qiime.import_data(
                    input_path=sub_manifest,
                    output_path=demux_qza,
                    dry_run=getattr(args, "dry_run", False),
                    show_stdout=getattr(args, "show_qiime", True),
                )

                trimmed_qza = sub_root / "trimmed.qza"
                input_for_dada2 = demux_qza
                if fwd_list and rev_list:
                    qiime.cutadapt_trim_paired(
                        input_seqs=demux_qza,
                        forward_primers=fwd_list,
                        reverse_primers=rev_list,
                        output_path=trimmed_qza,
                        cores=args.cores,
                        discard_untrimmed=args.discard_untrimmed,
                        no_indels=args.no_indels,
                        dry_run=getattr(args, "dry_run", False),
                        show_stdout=getattr(args, "show_qiime", True),
                    )
                    input_for_dada2 = trimmed_qza

                # >>> AUTO/FIXED TRUNCATION HERE <<<
                trunc_f, trunc_r, opt_meta = _auto_or_fixed_trunc(
                    run_path=run_path,
                    sub_root=sub_root,
                    input_for_dada2=input_for_dada2,
                    args=args,
                )

                table_qza = sub_root / "table.qza"
                seqs_qza = sub_root / "rep-seqs.qza"
                stats_qza = sub_root / "dada2-stats.qza"
                qiime.dada2_denoise_paired(
                    input_seqs=input_for_dada2,
                    trunc_len_f=trunc_f,
                    trunc_len_r=trunc_r,
                    output_table=table_qza,
                    output_rep_seqs=seqs_qza,
                    output_stats=stats_qza,
                    trim_left_f=args.trim_left_f,
                    trim_left_r=args.trim_left_r,
                    max_ee_f=args.max_ee_f,
                    max_ee_r=args.max_ee_r,
                    trunc_q=args.trunc_q,
                    min_overlap=args.min_overlap,
                    pooling_method=args.pooling_method,
                    chimera_method=args.chimera_method,
                    n_threads=args.cores,
                    dry_run=getattr(args, "dry_run", False),
                    show_stdout=getattr(args, "show_qiime", True),
                )

                per_group_tables.setdefault(gkey or "all", []).append(table_qza)
                per_group_seqs.setdefault(gkey or "all", []).append(seqs_qza)

        else:
            # No split
            fset, rset = set(), set()
            for r in run_rows:
                f = (r.get(args.front_f_col) or "").strip()
                g = (r.get(args.front_r_col) or "").strip()
                if f: fset.add(f)
                if g: rset.add(g)
            fwd_list, rev_list = sorted(fset), sorted(rset)

            demux_qza = run_root / "demux.qza"
            qiime.import_data(
                input_path=manifest_path,
                output_path=demux_qza,
                dry_run=getattr(args, "dry_run", False),
                show_stdout=getattr(args, "show_qiime", True),
            )

            input_for_dada2 = demux_qza
            if fwd_list and rev_list:
                trimmed_qza = run_root / "trimmed.qza"
                qiime.cutadapt_trim_paired(
                    input_seqs=demux_qza,
                    forward_primers=fwd_list,
                    reverse_primers=rev_list,
                    output_path=trimmed_qza,
                    cores=args.cores,
                    discard_untrimmed=args.discard_untrimmed,
                    no_indels=args.no_indels,
                    dry_run=getattr(args, "dry_run", False),
                    show_stdout=getattr(args, "show_qiime", True),
                )
                input_for_dada2 = trimmed_qza

            # >>> AUTO/FIXED TRUNCATION HERE <<<
            trunc_f, trunc_r, opt_meta = _auto_or_fixed_trunc(
                run_path=run_path,
                sub_root=run_root,
                input_for_dada2=input_for_dada2,
                args=args,
            )

            table_qza = run_root / "table.qza"
            seqs_qza = run_root / "rep-seqs.qza"
            stats_qza = run_root / "dada2-stats.qza"
            qiime.dada2_denoise_paired(
                input_seqs=input_for_dada2,
                trunc_len_f=trunc_f,
                trunc_len_r=trunc_r,
                output_table=table_qza,
                output_rep_seqs=seqs_qza,
                output_stats=stats_qza,
                trim_left_f=args.trim_left_f,
                trim_left_r=args.trim_left_r,
                max_ee_f=args.max_ee_f,
                max_ee_r=args.max_ee_r,
                trunc_q=args.trunc_q,
                min_overlap=args.min_overlap,
                pooling_method=args.pooling_method,
                chimera_method=args.chimera_method,
                n_threads=args.cores,
                dry_run=getattr(args, "dry_run", False),
                show_stdout=getattr(args, "show_qiime", True),
            )

            per_group_tables.setdefault("all", []).append(table_qza)
            per_group_seqs.setdefault("all", []).append(seqs_qza)

    # Merge and summarize (unchanged)
    merged_index = {}
    for gkey, tables in per_group_tables.items():
        seqs = per_group_seqs.get(gkey, [])
        gslug = slugify(gkey if gkey else "all")
        out_root = project_dir / ("groups" if args.split_by_primer_group else ".")
        if args.split_by_primer_group:
            out_root = project_dir / "groups" / gslug
        out_root.mkdir(parents=True, exist_ok=True)

        merged_table = out_root / "table.qza"
        merged_seqs = out_root / "rep-seqs.qza"
        qiime.merge_tables(tables, merged_table, dry_run=getattr(args, "dry_run", False), show_stdout=getattr(args, "show_qiime", True))
        qiime.merge_seqs(seqs, merged_seqs, dry_run=getattr(args, "dry_run", False), show_stdout=getattr(args, "show_qiime", True))

        table_qzv = out_root / "table.qzv"
        rep_qzv = out_root / "rep-seqs.qzv"
        qiime.feature_table_summarize(input_table=merged_table, output=table_qzv, sample_metadata_file=meta_aug,
                                      dry_run=getattr(args, "dry_run", False), show_stdout=getattr(args, "show_qiime", True))
        qiime.feature_table_tabulate_seqs(input_data=merged_seqs, output=rep_qzv,
                                          dry_run=getattr(args, "dry_run", False), show_stdout=getattr(args, "show_qiime", True))

        merged_index[gkey or "all"] = {
            "dir": str(out_root),
            "table": str(merged_table),
            "rep_seqs": str(merged_seqs),
            "table_qzv": str(table_qzv),
            "rep_seqs_qzv": str(rep_qzv),
            "runs_contributed": len(tables),
        }

    (project_dir / "MERGED.json").write_text(json.dumps(merged_index, indent=2) + "\n")
    LOG.info("Merging complete. Index → %s", project_dir / "MERGED.json")
    print(f"[ok] denoise pipeline complete → {project_dir}")
