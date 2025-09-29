# src/qs/commands/init.py
from __future__ import annotations
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from qs.utils.logger import get_logger
from qs.utils.samples import discover_fastqs, paired_sample_ids
from qs.metadata import DEFAULT_PROTECTED_COLS
from qs.metadata.template import write_metadata_with_primers, write_metadata_template
from qs.primers.detect import detect_primers_for_samples
from qs.primers.grouping import canonicalize_primer_pairs
from qs.primers import degen
from qs.utils.text import slugify
from qs.config.io import write_params_and_plan

LOG = get_logger("init")

def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "init", parents=[parent],
        help="Scan FASTQs, write metadata.tsv, and generate params.yaml with primer pairs & per-group classifier map.",
    )
    p.add_argument("--fastq-dir", type=Path, required=True)
    p.add_argument("--output-file", type=Path, required=True)
    p.add_argument("--protected-cols", type=str, default="__f_primer,__r_primer")
    p.add_argument("--keep-illumina-suffix", action="store_true")
    p.add_argument("--id-regex", type=str, default=None)
    p.add_argument("--force", action="store_true")

    # Scan controls
    p.add_argument("--auto-detect-primers", dest="auto_detect_primers", action="store_true")
    p.add_argument("--no-auto-detect-primers", dest="auto_detect_primers", action="store_false")
    p.set_defaults(auto_detect_primers=True)
    p.add_argument("--scan-max-reads", type=int, default=4000)
    p.add_argument("--scan-kmin", type=int, default=16)
    p.add_argument("--scan-kmax", type=int, default=24)
    p.add_argument("--scan-min-frac", type=float, default=0.30)

    # Indel/mismatch tolerant grouping
    p.add_argument("--group-edit-max", type=int, default=1,
                   help="Collapse primer pairs that differ by ≤ this many edits (IUPAC-aware). 0 = exact only.")

    # YAML expansion limit
    p.add_argument("--expand-limit", type=int, default=64)

    # Params output
    p.add_argument("--emit-params", dest="emit_params", action="store_true")
    p.add_argument("--no-emit-params", dest="emit_params", action="store_false")
    p.set_defaults(emit_params=True)
    p.add_argument("--params-file", type=Path, default=None)

    p.set_defaults(func=run)

def _parse_protected(s: str) -> List[str]:
    cols = [c.strip() for c in s.split(",") if c.strip()]
    return cols or list(DEFAULT_PROTECTED_COLS)

def _summaries_from_clusters(
    clusters, expand_limit: int
) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for c in clusters:
        fsum = degen.summarize_variants(c.fwd_variants, expand_limit=expand_limit)
        rsum = degen.summarize_variants(c.rev_variants, expand_limit=expand_limit)
        out.append({
            "group_key": c.key,
            "group_slug": c.slug,
            "n_samples": len(c.samples),
            "forward": fsum,
            "reverse": rsum,
        })
    return out

def run(args) -> None:
    fastq_dir: Path = args.fastq_dir
    meta_out: Path = args.output_file
    strip = not args.keep_illumina_suffix

    if meta_out.exists() and not args.force:
        LOG.error("Refusing to overwrite: %s (use --force)", meta_out)
        print(f"error: {meta_out} exists (use --force)", file=sys.stderr)
        sys.exit(1)

    paths = discover_fastqs(fastq_dir)
    if not paths:
        print(f"error: no FASTQ files under {fastq_dir}", file=sys.stderr)
        sys.exit(3)

    sample_ids = paired_sample_ids(paths, strip_illumina_suffix=strip, id_regex=args.id_regex)
    if not sample_ids:
        print("error: no paired R1/R2 FASTQs detected; check filenames.", file=sys.stderr)
        sys.exit(4)

    # 1) Detect primers per sample
    groups: Dict[str, List[str]] = {}
    per_sample: Dict[str, Tuple[str, str]] = {}
    if args.auto_detect_primers:
        per_sample = detect_primers_for_samples(
            fastq_dir, sample_ids,
            strip_illumina_suffix=strip, id_regex=args.id_regex,
            kmin=args.scan_kmin, kmax=args.scan_kmax, max_reads=args.scan_max_reads, min_fraction=args.scan_min_frac,
        )
        LOG.info("Primer detection: confident calls for %d/%d samples.", len(per_sample), len(sample_ids))

    # 2) Write metadata (filled if possible)
    if per_sample:
        write_metadata_with_primers(sample_ids, meta_out, per_sample, _parse_protected(args.protected_cols))
    else:
        write_metadata_template(sample_ids, meta_out, _parse_protected(args.protected_cols))
    print(f"[ok] metadata → {meta_out}")

    # 3) Build INDDEL-tolerant groups + params plan
    if args.emit_params:
        params_path = args.params_file or (meta_out.parent / "params.yaml")

        if per_sample:
            groups, sample_to_key, per_group_variants, clusters = canonicalize_primer_pairs(
                per_sample, max_edits=max(0, int(args.group_edit_max))
            )
        else:
            groups, clusters = {}, []

        params = {
            "resume": True,
            "split_by_primer_group": (len([k for k in groups if k]) > 1),
            "cores": 0,
            "discard_untrimmed": False,
            "no_indels": False,
            "auto_trunc": False,
            "trunc_len_f": 0, "trunc_len_r": 0,
            "trunc_lower_frac": 0.80, "trunc_min_len": 50, "trunc_step": 10, "trunc_refine_step": 5, "trunc_quick_learn": 250000,
            "trim_left_f": 0, "trim_left_r": 0, "max_ee_f": 2, "max_ee_r": 2, "trunc_q": 2, "min_overlap": 12,
            "pooling_method": "independent", "chimera_method": "consensus",
            "do_filter": True, "rare_freq_frac": 0.001, "rare_min_samples": 1,
            "contam_exclude": "mitochondria,chloroplast", "filter_unclassified_phylum": False, "min_sample_depth": 0,
            "retain_fraction": 0.90, "min_depth": 1000, "alpha_tsv": False, "taxa_barplot": True,
            "run_stats": False, "beta_group_cols": ["primer_group"], "beta_method": "permanova", "pairwise": False,
            "permutations": 999, "adonis_formula": None,
            "show_qiime": True,
            # Per-group classifier map scaffold
            "classifiers_map": { (slugify(k) if k else "all"): None for k in groups } | {"default": None},
            # Canonical IUPAC consensus per group (what trimming expects at the start of reads)
            "primer_pairs": {
                (slugify(k) if k else "all"): {
                    "forward_iupac": k.split("|", 1)[0] if k else "",
                    "reverse_iupac": k.split("|", 1)[1] if (k and "|" in k) else "",
                } for k in groups
            },
        }

        plan = {
            "generated": datetime.utcnow().isoformat(timespec="seconds") + "Z",
            "detected_groups_n": len([k for k in groups if k]),
            "primer_pairs": _summaries_from_clusters(clusters, args.expand_limit),
            "detected_primers_per_sample": per_sample,
            "suggested_command": "qs auto-run --fastq-dir <FASTQ_DIR> --params " + params_path.name,
            "notes": f"Primer groups were clustered with IUPAC-aware edit distance ≤ {max(0, int(args.group_edit_max))}.",
        }

        write_params_and_plan(params_path, params, plan)
        print(f"[ok] params → {params_path}")
