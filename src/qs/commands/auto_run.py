# src/qs/commands/auto_run.py
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Optional, Set, Any

from qs.utils.logger import get_logger
from qs.utils.samples import discover_fastqs, paired_sample_ids
from qs.metadata.template import write_metadata_template
from qs.utils.manifest import generate_manifest
from qs.commands import denoise_runs as cmd_denoise
from qs.commands import classify_sweep as cmd_cls
from qs.analysis.metrics import stage_and_run_metrics_for_group
from qs.analysis.staging import _stage_artifact
from qs.analysis.filtering import compute_min_freq_from_qzv, filter_table_and_seqs
from qs.analysis.stats import run_diversity_stats_for_group
from qs.utils.text import slugify

LOG = get_logger("auto")

# Defaults snapshot used to decide whether a CLI value was left at default
_DEFAULTS: Dict[str, Any] = {
    # selection / resume
    "groups": None, "skip_denoise": False, "skip_classify": False, "skip_metrics": False, "resume": False,
    # denoise knobs
    "split_by_primer_group": True, "discard_untrimmed": False, "no_indels": False, "cores": 0,
    # auto-truncation
    "auto_trunc": False, "trunc_len_f": 0, "trunc_len_r": 0, "trunc_lower_frac": 0.80,
    "trunc_min_len": 50, "trunc_step": 10, "trunc_refine_step": 5, "trunc_quick_learn": 250000,
    # dada2 misc
    "trim_left_f": 0, "trim_left_r": 0, "max_ee_f": 2, "max_ee_r": 2, "trunc_q": 2, "min_overlap": 12,
    "pooling_method": "independent", "chimera_method": "consensus",
    # classify
    "classifiers_dir": None,
    # filtering
    "do_filter": True, "rare_freq_frac": 0.001, "rare_min_samples": 1,
    "contam_exclude": "mitochondria,chloroplast", "filter_unclassified_phylum": False, "min_sample_depth": 0,
    # metrics
    "retain_fraction": 0.90, "min_depth": 1000, "alpha_tsv": False, "taxa_barplot": True,
    # stats
    "run_stats": False, "beta_group_cols": None, "beta_method": "permanova", "pairwise": False, "permutations": 999,
    "adonis_formula": None,
    # display
    "dry_run": False, "show_qiime": True,
}

def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "auto-run", parents=[parent],
        help=(
            "End-to-end pipeline with resume: "
            "(optional init) → per-run/group import+cutadapt+DADA2(+auto trunc) → merge "
            "→ classifier sweep → filtering → core metrics (+optional stats)."
        ),
    )

    # Inputs (project-dir now optional; default derived from fastq-dir)
    p.add_argument("--fastq-dir", type=Path, required=True, help="Top FASTQ directory; subfolders are runs.")
    p.add_argument("--project-dir", type=Path, default=None,
                   help="Project directory (default: <fastq-parent>/<fastq-name>-qs).")
    p.add_argument("--metadata-file", type=Path, default=None,
                   help="Metadata TSV (default: <project-dir>/metadata.tsv).")
    p.add_argument("--params", type=Path, default=None,
                   help="YAML/JSON file with parameter defaults to avoid many CLI flags.")
    p.add_argument("--auto-init", dest="auto_init", action="store_true",
                   help="If metadata.tsv is missing, create it and continue (default).")
    p.add_argument("--no-auto-init", dest="auto_init", action="store_false",
                   help="Do not auto-create metadata; require it to exist.")
    p.set_defaults(auto_init=True)

    # Optional one-shot metadata initialization then exit (unchanged)
    p.add_argument("--init-metadata", action="store_true",
                   help="Generate metadata and a top-level manifest, then stop so you can fill primers.")

    p.add_argument("--protected-cols", type=str, default="__f_primer,__r_primer",
                   help="Protected columns for metadata template (used with --init-metadata/--auto-init).")
    p.add_argument("--keep-illumina-suffix", action="store_true",
                   help="Keep _S#/_L### tokens in SampleIDs (default: strip).")
    p.add_argument("--id-regex", type=str, default=None, help="Optional regex to derive SampleID (named 'id' or group 1).")

    # Stage skipping / resume / selection
    p.add_argument("--skip-denoise", action="store_true", help="Skip import/cutadapt/DADA2/merge and reuse existing MERGED.json.")
    p.add_argument("--skip-classify", action="store_true", help="Skip classifier sweep.")
    p.add_argument("--skip-metrics", action="store_true", help="Skip phylogeny + core metrics.")
    p.add_argument("--resume", action="store_true", help="Auto-skip stages with completed outputs (idempotent resume).")
    p.add_argument("--groups", type=str, default=None, help="Comma-separated group keys or slugs to operate on (default: all).")

    # Denoise knobs (forwarded)
    p.add_argument("--split-by-primer-group", dest="split_by_primer_group", action="store_true",
                   help="Per primer group within each run (default: on).")
    p.add_argument("--no-split-by-primer-group", dest="split_by_primer_group", action="store_false",
                   help="Process all primer pairs together as a single group ('all').")
    p.set_defaults(split_by_primer_group=True)

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
    p.add_argument("--classifiers-dir", type=Path, default=None, help="Directory of sklearn classifier .qza files.")

    # Filtering (default ON). NOTE: escape % in help string to avoid argparse crash.
    p.add_argument("--filter", dest="do_filter", action="store_true",
                   help="Enable table filtering: rare ASVs, contaminants, low-depth samples (default: on).")
    p.add_argument("--no-filter", dest="do_filter", action="store_false", help="Disable filtering stage.")
    p.set_defaults(do_filter=True)

    p.add_argument("--rare-freq-frac", type=float, default=0.001,
                   help="Min ASV frequency = fraction of mean sample depth (e.g., 0.001 = 0.1%%).")
    p.add_argument("--rare-min-samples", type=int, default=1, help="Keep features present in ≥ this many samples.")
    p.add_argument("--contam-exclude", type=str, default="mitochondria,chloroplast",
                   help="Comma list to exclude as contaminants (taxonomy contains any of these).")
    p.add_argument("--filter-unclassified-phylum", action="store_true",
                   help="Keep only features with 'p__' in taxonomy (phylum-level classified). Off by default.")
    p.add_argument("--min-sample-depth", type=int, default=0,
                   help="Drop samples with total counts < this value after filtering (0 = keep all).")

    # Core metrics (phylogeny & diversity)
    p.add_argument("--retain-fraction", type=float, default=0.90, help="Fraction of samples to retain when choosing depth.")
    p.add_argument("--min-depth", type=int, default=1000, help="Minimum sampling depth.")
    p.add_argument("--alpha-tsv", action="store_true", help="Also write <group>/core-metrics-phylo/alpha-metrics.tsv.")

    # Stats (optional)
    p.add_argument("--run-stats", action="store_true",
                   help="Run diversity stats after core metrics (alpha-group-significance, alpha-correlation, beta-group-significance, ADONIS).")
    p.add_argument("--beta-group-cols", type=str, default=None,
                   help="Comma-separated categorical metadata columns for beta-group-significance (e.g., 'primer_group,treatment').")
    p.add_argument("--beta-method", choices=("permanova", "anosim", "permdisp"), default="permanova")
    p.add_argument("--pairwise", action="store_true", help="Pairwise tests for beta-group-significance.")
    p.add_argument("--permutations", type=int, default=999)
    p.add_argument("--adonis-formula", type=str, default=None,
                   help='Optional ADONIS formula, e.g., "treatment+block" or "treatment*block".')

    # Taxa barplot toggle (default ON)
    p.add_argument("--taxa-barplot", dest="taxa_barplot", action="store_true",
                   help="If taxonomy.qza is available, render taxa-barplot.qzv (default: on).")
    p.add_argument("--no-taxa-barplot", dest="taxa_barplot", action="store_false",
                   help="Disable taxa barplot generation.")
    p.set_defaults(taxa_barplot=True)

    p.set_defaults(func=run)


def _derive_default_project_dir(fastq_dir: Path) -> Path:
    fd = fastq_dir.resolve()
    return fd.parent / f"{fd.name}-qs"


def _generate_metadata_and_manifest(
    *,
    fastq_dir: Path,
    project_dir: Path,
    metadata_file: Path,
    protected_cols: str,
    keep_illumina_suffix: bool,
    id_regex: Optional[str],
) -> None:
    meta_path = metadata_file
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
    LOG.info("[auto-init] Metadata → %s", meta_path)
    LOG.info("[auto-init] Top-level manifest → %s", top_manifest)


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
    _generate_metadata_and_manifest(
        fastq_dir=fastq_dir,
        project_dir=project_dir,
        metadata_file=meta_path,
        protected_cols=protected_cols,
        keep_illumina_suffix=keep_illumina_suffix,
        id_regex=id_regex,
    )
    print(f"[ok] Wrote metadata + manifest. Fill primers, then re-run without --init-metadata.")
    sys.exit(0)


def _load_params_file(path: Path) -> Dict[str, Any]:
    text = path.read_text()
    try:
        import yaml  # type: ignore
        data = yaml.safe_load(text)
    except Exception:
        data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError("Params file must contain a mapping/object at the top level.")
    return data


def _flatten(d: Dict[str, Any]) -> Dict[str, Any]:
    flat: Dict[str, Any] = {}
    for k, v in d.items():
        if isinstance(v, dict):
            for k2, v2 in v.items():
                flat[k2] = v2
        else:
            flat[k] = v
    return flat


def _apply_params_defaults(args, params: Dict[str, Any]) -> None:
    flat = _flatten(params)
    def _norm_cols(x: Any) -> Optional[str]:
        if x is None:
            return None
        if isinstance(x, (list, tuple)):
            return ",".join(map(str, x))
        return str(x)

    if "beta_group_cols" in flat:
        flat["beta_group_cols"] = _norm_cols(flat["beta_group_cols"])
    if "groups" in flat:
        flat["groups"] = _norm_cols(flat["groups"])

    for name, value in flat.items():
        if not hasattr(args, name):
            continue
        cur = getattr(args, name)
        default = _DEFAULTS.get(name, None)
        if cur == default:
            setattr(args, name, value)


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
    if getattr(args, "params", None):
        try:
            params = _load_params_file(args.params)
            _apply_params_defaults(args, params)
            LOG.info("Loaded defaults from params file: %s", args.params)
        except Exception as e:
            print(f"error: failed to read params file {args.params}: {e}", file=sys.stderr)
            sys.exit(2)

    fastq_dir: Path = args.fastq_dir
    project_dir: Path = args.project_dir or _derive_default_project_dir(fastq_dir)
    meta_path: Optional[Path] = args.metadata_file or (project_dir / "metadata.tsv")

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

    if not meta_path.exists():
        if args.auto_init:
            _generate_metadata_and_manifest(
                fastq_dir=fastq_dir,
                project_dir=project_dir,
                metadata_file=meta_path,
                protected_cols=args.protected_cols,
                keep_illumina_suffix=args.keep_illumina_suffix,
                id_regex=args.id_regex,
            )
        else:
            print(f"error: metadata file not found: {meta_path} (use --init-metadata or enable --auto-init)", file=sys.stderr)
            sys.exit(2)

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

    # Stage 3: classifier sweep (optional)
    winners: dict[str, dict] = {}
    classifiers_ok = (args.classifiers_dir and Path(args.classifiers_dir).is_dir())

    for group_key in sorted(selected):
        rec = merged_index[group_key]
        group_dir = Path(rec["dir"])
        rep_seqs = group_dir / "rep-seqs.qza"

        if args.skip_classify or not classifiers_ok:
            LOG.info("Skipping classification for group=%s", group_key)
        else:
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

    # Stage 4: filtering → phylogeny + core metrics (+ optional stats)
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
        elif (group_dir / "taxonomy.qza").exists():
            taxonomy = group_dir / "taxonomy.qza"

        table_qzv = group_dir / "table.qzv"
        if not table_qzv.exists():
            from qs.qiime import commands as qiime
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
            # metrics stage will now prefer *_final.qza if present (see metrics.py)
        core_dir = stage_and_run_metrics_for_group(
            group_dir=group_dir,
            metadata_augmented=meta_aug,
            retain_fraction=args.retain_fraction,
            min_depth=args.min_depth,
            if_exists="skip",
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

        if getattr(args, "alpha_tsv", False):
            from qs.analysis.metrics import write_alpha_metrics_tsv
            try:
                path = write_alpha_metrics_tsv(core_dir)
                LOG.info("Alpha metrics → %s", path)
            except Exception as e:
                LOG.warning("Alpha TSV export skipped: %s", e)

        LOG.info("Group %s → %s", group_key, core_dir)

    print(f"[ok] auto-run complete → {project_dir}")
