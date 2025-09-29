# src/qs/commands/auto_run.py
from __future__ import annotations
import json, sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Optional, Set

from qs.utils.logger import get_logger
from qs.utils.samples import discover_fastqs, paired_sample_ids
from qs.metadata.template import write_metadata_template
from qs.utils.manifest import generate_manifest
from qs.commands import denoise_runs as cmd_denoise
from qs.commands import classify_sweep as cmd_cls
from qs.analysis.metrics import stage_and_run_metrics_for_group
from qs.analysis.filtering import compute_min_freq_from_qzv, filter_table_and_seqs
from qs.analysis.stats import run_diversity_stats_for_group
from qs.utils.text import slugify
from qs.config.io import load_params_typed, apply_params_defaults
from qs.config.schema import Params

LOG = get_logger("auto")

_DEFAULTS: Dict[str, Any] = {
    "groups": None, "skip_denoise": False, "skip_classify": False, "skip_metrics": False, "resume": False,
    "split_by_primer_group": True, "discard_untrimmed": False, "no_indels": False, "cores": 0,
    "auto_trunc": False, "trunc_len_f": 0, "trunc_len_r": 0, "trunc_lower_frac": 0.80,
    "trunc_min_len": 50, "trunc_step": 10, "trunc_refine_step": 5, "trunc_quick_learn": 250000,
    "trim_left_f": 0, "trim_left_r": 0, "max_ee_f": 2, "max_ee_r": 2, "trunc_q": 2, "min_overlap": 12,
    "pooling_method": "independent", "chimera_method": "consensus",
    "classifiers_dir": None,
    "do_filter": True, "rare_freq_frac": 0.001, "rare_min_samples": 1,
    "contam_exclude": "mitochondria,chloroplast", "filter_unclassified_phylum": False, "min_sample_depth": 0,
    "retain_fraction": 0.90, "min_depth": 1000, "alpha_tsv": False, "taxa_barplot": True,
    "run_stats": False, "beta_group_cols": None, "beta_method": "permanova", "pairwise": False, "permutations": 999,
    "adonis_formula": None,
    "dry_run": False, "show_qiime": True,
    # NEW default used when forwarding to denoise
    "group_edit_max": 1,
}

def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "auto-run", parents=[parent],
        help="End-to-end: init (optional) → per-run/group import+cutadapt+DADA2 → merge → classify (per-group map or sweep) → filter → metrics (+stats).",
    )
    p.add_argument("--fastq-dir", type=Path, required=True)
    p.add_argument("--project-dir", type=Path, default=None)
    p.add_argument("--metadata-file", type=Path, default=None)
    p.add_argument("--params", type=Path, default=None, help="YAML/JSON with 'params' and optional 'plan'")
    p.add_argument("--auto-init", dest="auto_init", action="store_true")
    p.add_argument("--no-auto-init", dest="auto_init", action="store_false")
    p.set_defaults(auto_init=True)

    p.add_argument("--init-metadata", action="store_true")
    p.add_argument("--protected-cols", type=str, default="__f_primer,__r_primer")
    p.add_argument("--keep-illumina-suffix", action="store_true")
    p.add_argument("--id-regex", type=str, default=None)

    p.add_argument("--skip-denoise", action="store_true")
    p.add_argument("--skip-classify", action="store_true")
    p.add_argument("--skip-metrics", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--groups", type=str, default=None)

    p.add_argument("--split-by-primer-group", dest="split_by_primer_group", action="store_true")
    p.add_argument("--no-split-by-primer-group", dest="split_by_primer_group", action="store_false")
    p.set_defaults(split_by_primer_group=True)

    # Forwarded to cutadapt/DADA2
    p.add_argument("--discard-untrimmed", action="store_true")
    p.add_argument("--no-indels", action="store_true")
    p.add_argument("--cores", type=int, default=0)

    p.add_argument("--auto-trunc", action="store_true")
    p.add_argument("--trunc-len-f", type=int, default=0)
    p.add_argument("--trunc-len-r", type=int, default=0)
    p.add_argument("--trunc-lower-frac", type=float, default=0.80)
    p.add_argument("--trunc-min-len", type=int, default=50)
    p.add_argument("--trunc-step", type=int, default=10)
    p.add_argument("--trunc-refine-step", type=int, default=5)
    p.add_argument("--trunc-quick-learn", type=int, default=250000)

    p.add_argument("--trim-left-f", type=int, default=0)
    p.add_argument("--trim-left-r", type=int, default=0)
    p.add_argument("--max-ee-f", type=int, default=2)
    p.add_argument("--max-ee-r", type=int, default=2)
    p.add_argument("--trunc-q", type=int, default=2)
    p.add_argument("--min-overlap", type=int, default=12)
    p.add_argument("--pooling-method", type=str, default="independent")
    p.add_argument("--chimera-method", type=str, default="consensus")

    # Classifiers
    p.add_argument("--classifiers-dir", type=Path, default=None)
    p.add_argument("--classifiers-map", type=Path, default=None, help="YAML/JSON mapping slug|key|default -> dir|.qza")

    # Filtering / metrics / stats
    p.add_argument("--filter", dest="do_filter", action="store_true")
    p.add_argument("--no-filter", dest="do_filter", action="store_false")
    p.set_defaults(do_filter=True)
    p.add_argument("--rare-freq-frac", type=float, default=0.001)
    p.add_argument("--rare-min-samples", type=int, default=1)
    p.add_argument("--contam-exclude", type=str, default="mitochondria,chloroplast")
    p.add_argument("--filter-unclassified-phylum", action="store_true")
    p.add_argument("--min-sample-depth", type=int, default=0)

    p.add_argument("--retain-fraction", type=float, default=0.90)
    p.add_argument("--min-depth", type=int, default=1000)
    p.add_argument("--alpha-tsv", action="store_true")

    p.add_argument("--run-stats", action="store_true")
    p.add_argument("--beta-group-cols", type=str, default=None)
    p.add_argument("--beta-method", choices=("permanova", "anosim", "permdisp"), default="permanova")
    p.add_argument("--pairwise", action="store_true")
    p.add_argument("--permutations", type=int, default=999)
    p.add_argument("--adonis-formula", type=str, default=None)

    p.add_argument("--taxa-barplot", dest="taxa_barplot", action="store_true")
    p.add_argument("--no-taxa-barplot", dest="taxa_barplot", action="store_false")
    p.set_defaults(taxa_barplot=True)

    # NEW: expose grouping tolerance here too (forwarded to denoise)
    p.add_argument("--group-edit-max", type=int, default=_DEFAULTS["group_edit_max"],
                   help="Collapse primer pairs that differ by ≤ this many edits (IUPAC-aware). 0 = exact only.")

    p.set_defaults(func=run)

def _derive_default_project_dir(fastq_dir: Path) -> Path:
    fd = fastq_dir.resolve()
    return fd.parent / f"{fd.name}-qs"

def _generate_metadata_and_manifest(*, fastq_dir: Path, project_dir: Path, metadata_file: Path,
                                    protected_cols: str, keep_illumina_suffix: bool, id_regex: Optional[str]) -> None:
    sids = paired_sample_ids(
        discover_fastqs(fastq_dir), strip_illumina_suffix=not keep_illumina_suffix, id_regex=id_regex,
    )
    if not sids:
        print("error: no paired FASTQs discovered; cannot init metadata.", file=sys.stderr)
        sys.exit(2)
    write_metadata_template(sids, metadata_file, [c.strip() for c in protected_cols.split(",") if c.strip()])
    top_manifest = project_dir / "fastq.manifest"
    generate_manifest(fastq_dir, top_manifest, strip_illumina_suffix=not keep_illumina_suffix, id_regex=id_regex)

def _load_merged_index(project_dir: Path) -> Dict[str, dict]:
    p = project_dir / "MERGED.json"
    if not p.exists():
        raise FileNotFoundError(p)
    return json.loads(p.read_text())

def _select_groups(idx: Dict[str, dict], groups_arg: Optional[str]) -> Set[str]:
    if not groups_arg:
        return set(idx.keys())
    want = {g.strip() for g in groups_arg.split(",") if g.strip()}
    slugs = {slugify(k if k else "all"): k for k in idx}
    out: Set[str] = set()
    for g in want:
        if g in idx: out.add(g)
        elif g in slugs: out.add(slugs[g])
    return out

def _winner_if_present(group_dir: Path) -> Optional[str]:
    rpt = group_dir / "classifiers" / "optimal_classifiers.json"
    if not rpt.exists(): return None
    try:
        data = json.loads(rpt.read_text())
        for k in ("pct_depth≥7", "median_conf", "mean_conf"):
            picks = data.get("winners", {}).get(k, [])
            if picks: return str(picks[0]["classifier"])
    except Exception:
        return None
    return None

def _resolve_classifier_source(group_key: str, slug: str, *, mapping: Optional[Dict[str, str]],
                               fallback_dir: Optional[Path]) -> Optional[Path]:
    if mapping:
        for k in (group_key, slug, "default"):
            v = mapping.get(k)
            if v: return Path(v)
    return fallback_dir

def _load_mapping_file(path: Path) -> Dict[str, str]:
    from qs.config.io import load_params_raw
    data = load_params_raw(path)
    if "params" in data and isinstance(data["params"], dict) and isinstance(data["params"].get("classifiers_map"), dict):
        return {str(k): str(v) for k, v in data["params"]["classifiers_map"].items() if v is not None}
    return {str(k): str(v) for k, v in data.items() if v is not None}

def _to_ns(args=None, **kwargs) -> SimpleNamespace:
    """Accept argparse-style positional or kwargs and normalize to a SimpleNamespace."""
    if args is not None and not isinstance(args, SimpleNamespace):
        # Usually argparse.Namespace; that's fine too
        return args  # type: ignore[return-value]
    return SimpleNamespace(**kwargs)

def run(args=None, **kwargs) -> None:
    # Accept both call styles: func(args) or func(**vars(args))
    args = _to_ns(args, **kwargs)

    # params file → typed + apply to CLI only where at defaults
    params_typed: Optional[Params] = None
    params_map: Optional[Dict[str, Any]] = None
    if getattr(args, "params", None):
        try:
            params_typed, plan = load_params_typed(args.params)
            params_map = (params_typed.model_dump() if params_typed else None)
            apply_params_defaults(args, params_map or {}, _DEFAULTS)
        except Exception as e:
            print(f"error: reading params {args.params}: {e}", file=sys.stderr)
            sys.exit(2)

    fastq_dir: Path = args.fastq_dir
    project_dir: Path = args.project_dir or _derive_default_project_dir(fastq_dir)
    meta_path: Path = args.metadata_file or (project_dir / "metadata.tsv")
    project_dir.mkdir(parents=True, exist_ok=True)

    if getattr(args, "init_metadata", False):
        _generate_metadata_and_manifest(
            fastq_dir=fastq_dir, project_dir=project_dir, metadata_file=meta_path,
            protected_cols=args.protected_cols, keep_illumina_suffix=args.keep_illumina_suffix, id_regex=args.id_regex,
        )
        print("[ok] Wrote metadata + manifest. Fill primers, then re-run without --init-metadata.")
        return

    if not meta_path.exists():
        if getattr(args, "auto_init", True):
            _generate_metadata_and_manifest(
                fastq_dir=fastq_dir, project_dir=project_dir, metadata_file=meta_path,
                protected_cols=args.protected_cols, keep_illumina_suffix=args.keep_illumina_suffix, id_regex=args.id_regex,
            )
        else:
            print(f"error: metadata not found: {meta_path}", file=sys.stderr); sys.exit(2)

    merged_idx = project_dir / "MERGED.json"
    if getattr(args, "resume", False) and merged_idx.exists():
        args.skip_denoise = True

    if not getattr(args, "skip_denoise", False):
        ns = SimpleNamespace(
            fastq_dir=fastq_dir, project_dir=project_dir, metadata_file=meta_path,
            keep_illumina_suffix=args.keep_illumina_suffix, id_regex=args.id_regex,
            front_f_col="__f_primer", front_r_col="__r_primer",
            split_by_primer_group=args.split_by_primer_group,
            cores=args.cores, discard_untrimmed=args.discard_untrimmed, no_indels=args.no_indels,
            auto_trunc=args.auto_trunc, trunc_len_f=args.trunc_len_f, trunc_len_r=args.trunc_len_r,
            trunc_lower_frac=args.trunc_lower_frac, trunc_min_len=args.trunc_min_len,
            trunc_step=args.trunc_step, trunc_refine_step=args.trunc_refine_step, trunc_quick_learn=args.trunc_quick_learn,
            trim_left_f=args.trim_left_f, trim_left_r=args.trim_left_r, max_ee_f=args.max_ee_f, max_ee_r=args.max_ee_r,
            trunc_q=args.trunc_q, min_overlap=args.min_overlap, pooling_method=args.pooling_method, chimera_method=args.chimera_method,
            dry_run=getattr(args, "dry_run", False), show_qiime=getattr(args, "show_qiime", True),
            auto_detect_primers=True, scan_max_reads=4000, scan_kmin=16, scan_kmax=24, scan_min_frac=0.30,
            # NEW: forward grouping tolerance (default 1 if not supplied)
            group_edit_max=getattr(args, "group_edit_max", _DEFAULTS["group_edit_max"]),
        )
        cmd_denoise.run(ns)

    try:
        merged_index = _load_merged_index(project_dir)
    except FileNotFoundError:
        print("error: MERGED.json not found; run denoise or use --resume only with existing project.", file=sys.stderr)
        sys.exit(3)

    selected = _select_groups(merged_index, getattr(args, "groups", None))
    if not selected:
        print("error: no matching groups found.", file=sys.stderr); sys.exit(4)

    # Resolve per-group classifier sources
    classifiers_map: Optional[Dict[str, str]] = None
    if getattr(args, "classifiers_map", None):
        try:
            classifiers_map = _load_mapping_file(args.classifiers_map)
        except Exception as e:
            print(f"error: reading classifiers map: {e}", file=sys.stderr); sys.exit(5)
    elif params_typed and params_typed.classifiers_map:
        classifiers_map = {k: v for k, v in (params_typed.classifiers_map or {}).items() if v}

    winners: Dict[str, Dict[str, str]] = {}
    if not getattr(args, "skip_classify", False):
        for g in sorted(selected):
            gdir = Path(merged_index[g]["dir"])
            rep = gdir / "rep-seqs.qza"
            slug = slugify(g if g else "all")

            if getattr(args, "resume", False):
                prior = _winner_if_present(gdir)
                if prior:
                    winners[g] = {"classifier_tag": prior, "priority": "resume"}
                    continue

            src = _resolve_classifier_source(g, slug, mapping=classifiers_map, fallback_dir=args.classifiers_dir)
            if not src:
                LOG.info("No classifier for group=%s; skipping classification.", g); continue
            out = gdir / "classifiers"
            if src.is_file() and src.suffix.lower() == ".qza":
                res = cmd_cls.classify_with_single(
                    input_reads=rep, classifier_qza=src, out_dir=out,
                    dry_run=getattr(args, "dry_run", False), show_qiime=getattr(args, "show_qiime", True),
                )
            elif src.is_dir():
                res = cmd_cls.sweep_and_pick(
                    input_reads=rep, classifiers_dir=src, out_dir=out,
                    dry_run=getattr(args, "dry_run", False), show_qiime=getattr(args, "show_qiime", True),
                )
            else:
                LOG.warning("Invalid classifier source for %s: %s", g, src); continue
            winners[g] = res  # type: ignore[assignment]
            print(f"[ok] group={g} classifier winner={res['classifier_tag']} by {res['priority']}")

    if getattr(args, "skip_metrics", False):
        print(f"[ok] auto-run (skipped metrics) → {project_dir}"); return

    meta_aug = project_dir / "metadata.augmented.tsv"
    for g in sorted(selected):
        rec = merged_index[g]; gdir = Path(rec["dir"])
        table = gdir / "table.qza"
        rep = gdir / "rep-seqs.qza"

        taxonomy = None
        if (g in winners) and winners[g].get("classifier_tag"):
            tag = winners[g]["classifier_tag"]
            taxonomy = gdir / "classifiers" / f"{tag}_classification.qza"
        elif (gdir / "taxonomy.qza").exists():
            taxonomy = gdir / "taxonomy.qza"

        table_qzv = gdir / "table.qzv"
        if not table_qzv.exists():
            from qs.qiime import commands as qiime
            qiime.feature_table_summarize(table, table_qzv, sample_metadata_file=meta_aug,
                                          dry_run=getattr(args, "dry_run", False), show_stdout=getattr(args, "show_qiime", True))

        if getattr(args, "do_filter", True):
            min_freq = compute_min_freq_from_qzv(table_qzv, args.rare_freq_frac, fallback_min=1)
            f = filter_table_and_seqs(
                group_dir=gdir, table_qza=table, rep_seqs_qza=rep, table_qzv=table_qzv, taxonomy_qza=taxonomy,
                min_freq=min_freq, min_samples=args.rare_min_samples, contam_exclude=args.contam_exclude,
                keep_only_phylum_classified=args.filter_unclassified_phylum, min_sample_depth=args.min_sample_depth,
                metadata_augmented=meta_aug, dry_run=getattr(args, "dry_run", False), show_qiime=getattr(args, "show_qiime", True),
            )
            table = f["table_final"]; rep = f["rep_seqs_final"]

        stage_and_run_metrics_for_group(
            group_dir=gdir, metadata_augmented=meta_aug, retain_fraction=args.retain_fraction, min_depth=args.min_depth,
            if_exists="skip", dry_run=getattr(args, "dry_run", False), show_qiime=getattr(args, "show_qiime", True),
            make_taxa_barplot=getattr(args, "taxa_barplot", True),
        )

        if getattr(args, "run_stats", False):
            cols = [c.strip() for c in (args.beta_group_cols.split(",") if args.beta_group_cols else []) if c.strip()]
            run_diversity_stats_for_group(
                group_dir=gdir, metadata_file=meta_aug,
                beta_group_cols=cols, beta_method=args.beta_method, pairwise=args.pairwise, permutations=args.permutations,
                adonis_formula=args.adonis_formula, run_alpha_correlation=True,
                dry_run=getattr(args, "dry_run", False), show_qiime=getattr(args, "show_qiime", True),
            )

    print(f"[ok] auto-run complete → {project_dir}")
