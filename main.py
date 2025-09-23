# ./main.py
import argparse
import json
from pathlib import Path

from modules import io_utils, iterators, logger, qiime_wrapper, downstream

__VERSION__="1.0.0"
__AUTHOR__="Erick Samera (erick.samera@kpu.ca)"

def parse_args():
    p = argparse.ArgumentParser(prog="pipeline", epilog=f"{__VERSION__} | {__AUTHOR__}")

    # ---------------- Global options (available on all subcommands) ---------------
    common = argparse.ArgumentParser(add_help=False)
    # Stream QIIME logs live by default; allow turning off.
    common.add_argument("--show-qiime", dest="show_qiime", action="store_true",
                        help="Stream QIIME logs live (default).")
    common.add_argument("--no-show-qiime", dest="show_qiime", action="store_false",
                        help="Do not stream QIIME logs live; capture and print on error.")
    common.set_defaults(show_qiime=True)
    # Global dry-run parity
    common.add_argument("--dry-run", action="store_true", help="Print commands without executing")

    sp = p.add_subparsers(dest="cmd", required=True)

    # shared opts for steps that still touch FASTQs
    def add_fastq_ops(sub):
        sub.add_argument("--fastq-dir", type=Path, required=True, help="Project FASTQ dir")
        sub.add_argument(
            "--project-dir",
            type=Path,
            required=True,
            help="Named output directory to create/use (pipeline writes here)",
        )

    # full run
    run = sp.add_parser("run", help="Run full pipeline", parents=[common])
    add_fastq_ops(run)
    run.add_argument("--classifiers-dir", type=Path, default=Path("/home/erick/qiime"))

    # import-only
    imp = sp.add_parser("import", help="Only import FASTQs to a QIIME artifact", parents=[common])
    add_fastq_ops(imp)

    # trimming-only
    trm = sp.add_parser("trim", help="Only search optimal trimming", parents=[common])
    add_fastq_ops(trm)

    # classify-only (NO fastq-dir required)
    cls = sp.add_parser("classify", help="Run classifier sweep on rep_seqs in <project-dir>", parents=[common])
    cls.add_argument(
        "--project-dir",
        type=Path,
        required=True,
        help="Directory containing work/optimal_trimming and rep_seqs artifact",
    )
    cls.add_argument(
        "--input-reads",
        type=Path,
        default=None,
        help="Optional FeatureData[Sequence] .qza to classify (overrides best rep_seqs)",
    )
    cls.add_argument("--classifiers-dir", type=Path, default=Path("/home/erick/qiime"))

    # ---------------- Downstream (Moving Pictures style) ----------------

    stg = sp.add_parser("stage", help="Stage winners into analysis/ and build summaries", parents=[common])
    stg.add_argument("--project-dir", type=Path, required=True)
    stg.add_argument("--metadata-file", type=Path, required=False)
    stg.add_argument(
        "--preferred-metric",
        type=str,
        default="pct_depth≥7",
        help="Which classifier metric to prefer (pct_depth≥7 | pct_depth>=7 | median_conf | mean_conf)",
    )

    phy = sp.add_parser("phylogeny", help="Build phylogeny (MAFFT→FastTree→root) in analysis/", parents=[common])
    phy.add_argument("--project-dir", type=Path, required=True)

    div = sp.add_parser("diversity", help="Run core phylogenetic diversity + tests + Emperor", parents=[common])
    div.add_argument("--project-dir", type=Path, required=True)
    div.add_argument("--metadata-file", type=Path, required=True)
    div.add_argument("--sampling-depth", type=int, default=None)
    div.add_argument(
        "--beta-cols",
        type=str,
        default="body-site,subject",
        help="Comma-separated metadata columns for beta-group-significance",
    )
    div.add_argument("--time-column", type=str, default=None)
    div.add_argument("--skip-alpha-tests", action="store_true")
    div.add_argument("--skip-beta-tests", action="store_true")
    div.add_argument("--skip-emperor", action="store_true")
    div.add_argument(
        "--if-exists",
        choices=["skip", "overwrite", "new", "error"],
        default="skip",
        help="What to do if analysis/core-metrics-phylo already exists (default: skip).",
    )

    arf = sp.add_parser("alpha-rarefaction", help="Alpha rarefaction visualization", parents=[common])
    arf.add_argument("--project-dir", type=Path, required=True)
    arf.add_argument("--metadata-file", type=Path, required=True)
    arf.add_argument("--max-depth", type=int, default=None)

    tx = sp.add_parser("taxonomy", help="Taxonomy bar plots", parents=[common])
    tx.add_argument("--project-dir", type=Path, required=True)
    tx.add_argument("--metadata-file", type=Path, required=True)

    dwn = sp.add_parser(
        "downstream",
        parents=[common],
        help="Stage winners → phylogeny → core diversity → alpha-rarefaction → taxa barplots",
    )
    dwn.add_argument("--project-dir", type=Path, required=True)
    dwn.add_argument("--metadata-file", type=Path, required=True)
    dwn.add_argument("--preferred-metric", type=str, default="pct_depth≥7")
    dwn.add_argument("--sampling-depth", type=int, default=None)
    dwn.add_argument("--beta-cols", type=str, default="body-site,subject")
    dwn.add_argument("--time-column", type=str, default=None)
    dwn.add_argument("--no-taxa-barplots", action="store_true")

    # ---------------- Diversity sweep (nested or one-level) ----------------
    sweep = sp.add_parser(
        "diversity-sweep",
        parents=[common],
        help="Iterate core diversity over metadata subsets (nested by/within, or one-level by-only).",
    )
    sweep.add_argument("--project-dir", type=Path, required=True)
    sweep.add_argument("--metadata-file", type=Path, required=True)
    sweep.add_argument("--by", type=str, default=None,
                       help="Comma-separated outer loop columns (or all categorical if omitted).")
    sweep.add_argument("--within", type=str, default=None,
                       help="Comma-separated inner loop columns (ignored if --by-only is set).")
    sweep.add_argument("--by-only", action="store_true",
                       help="Use a single-level sweep: only iterate values of --by (or all categorical if --by omitted).")
    sweep.add_argument("--min-samples", type=int, default=5,
                       help="Skip subsets with < N samples (default: 5).")
    sweep.add_argument("--retain-fraction", type=float, default=0.90,
                       help="Sampling depth retains ~this fraction of samples per subset (default: 0.90).")
    sweep.add_argument("--beta-cols", type=str, default="",
                       help="Additional comma-separated columns for beta-group-significance.")
    sweep.add_argument("--time-column", type=str, default=None,
                       help="Metadata column for Emperor custom axis, if present.")
    sweep.add_argument("--if-exists", choices=["skip", "overwrite", "new", "error"], default="skip",
                       help="Handle existing subset outputs (default: skip).")

    return p.parse_args()


def _paths(project_fastq_dir: Path, project_dir: Path):
    work = project_dir / "work"
    imported_qza = project_dir / "output.qza"
    trim_dir = work / "optimal_trimming"
    cls_dir = work / "optimal_classifier"
    return work, imported_qza, trim_dir, cls_dir


def cmd_import(project_fastq_dir: Path, project_dir: Path, dry_run: bool, show_qiime: bool):
    io_utils.generate_manifest(project_fastq_dir, project_dir)
    qiime_wrapper.import_data(
        input_path=project_dir / "fastq.manifest",
        output_path=project_dir / "output.qza",
        dry_run=dry_run,
        show_qiime=show_qiime,
    )


def cmd_trim(imported_qza: Path, dry_run: bool, show_qiime: bool):
    # Expected to return (trunc_len_f, trunc_len_r) from your iterator
    return iterators.get_optimal_trimming(imported_qza=imported_qza, dry_run=dry_run, show_qiime=show_qiime)


def _load_best_trim(trim_dir: Path) -> tuple[int, int]:
    meta = json.loads((trim_dir / "optimal_trimming.json").read_text())
    return int(meta["trunc_len_f"]), int(meta["trunc_len_r"])


def _best_rep_seqs(project_dir: Path) -> Path:
    trim_dir = project_dir / "work" / "optimal_trimming"
    f, r = _load_best_trim(trim_dir)
    rep = trim_dir / f"{f}-{r}_output_rep_seqs.qza"
    if not rep.exists():
        raise FileNotFoundError(f"rep_seqs not found: {rep}")
    return rep


def cmd_classify(project_dir: Path, classifiers_dir: Path, input_reads: Path | None, dry_run: bool, show_qiime: bool):
    rep_seqs = input_reads if input_reads else _best_rep_seqs(project_dir)
    return iterators.get_optimal_classifier(
        imported_qza=rep_seqs,
        classifiers_dir=classifiers_dir,
        dry_run=dry_run,
        show_qiime=show_qiime,
    )

def main():
    args = parse_args()
    logger.setup_logger()

    # Steps that still need FASTQs
    if args.cmd in {"run", "import", "trim"}:
        project_fastq_dir = args.fastq_dir.resolve()
        project_dir = args.project_dir.resolve()
        project_dir.mkdir(parents=True, exist_ok=True)
        (project_dir / "work").mkdir(parents=True, exist_ok=True)
        work, imported_qza, trim_dir, _ = _paths(project_fastq_dir, project_dir)

        if args.cmd == "import":
            cmd_import(project_fastq_dir, project_dir, args.dry_run, args.show_qiime)
            return

        if args.cmd == "trim":
            cmd_trim(imported_qza, args.dry_run, args.show_qiime)
            return

        if args.cmd == "run":
            cmd_import(project_fastq_dir, project_dir, args.dry_run, args.show_qiime)
            cmd_trim(imported_qza, args.dry_run, args.show_qiime)
            cmd_classify(project_dir, args.classifiers_dir, None, args.dry_run, args.show_qiime)
            return

    # Classify does NOT require --fastq-dir
    if args.cmd == "classify":
        project_dir = args.project_dir.resolve()
        cmd_classify(project_dir, args.classifiers_dir, args.input_reads, args.dry_run, args.show_qiime)
        return

    # ---------------- Downstream command handlers ----------------

    if args.cmd == "stage":
        proj = args.project_dir.resolve()
        downstream.stage_winners(
            proj,
            preferred_metric=args.preferred_metric,
            metadata_file=args.metadata_file,
            show_qiime=args.show_qiime,
            dry_run=args.dry_run,
        )
        return

    if args.cmd == "phylogeny":
        proj = args.project_dir.resolve()
        ana = downstream.analysis_dir(proj)
        downstream.build_phylogeny(ana, show_qiime=args.show_qiime, dry_run=args.dry_run)
        return

    if args.cmd == "diversity":
        proj = args.project_dir.resolve()
        ana = downstream.analysis_dir(proj)
        downstream.run_core_diversity(
            ana,
            args.metadata_file,
            sampling_depth=args.sampling_depth,
            beta_cols=[c.strip() for c in args.beta_cols.split(",") if c.strip()],
            time_column=args.time_column,
            include_alpha_tests=not args.skip_alpha_tests if hasattr(args, "skip_alpha_tests") else True,
            include_beta_tests=not args.skip_beta_tests if hasattr(args, "skip_beta_tests") else True,
            include_emperor=not args.skip_emperor if hasattr(args, "skip_emperor") else True,
            auto_build_tree=True,
            show_qiime=args.show_qiime,
            if_exists=args.if_exists,
            dry_run=args.dry_run,
        )
        return

    if args.cmd == "alpha-rarefaction":
        proj = args.project_dir.resolve()
        ana = downstream.analysis_dir(proj)
        downstream.run_alpha_rarefaction(ana, args.metadata_file, max_depth=args.max_depth,
                                         show_qiime=args.show_qiime, dry_run=args.dry_run)
        return

    if args.cmd == "taxonomy":
        proj = args.project_dir.resolve()
        ana = downstream.analysis_dir(proj)
        downstream.run_taxa_barplots(ana, args.metadata_file, show_qiime=args.show_qiime, dry_run=args.dry_run)
        return

    if args.cmd == "downstream":
        proj = args.project_dir.resolve()
        downstream.run_downstream(
            proj,
            args.metadata_file,
            preferred_metric=args.preferred_metric,
            sampling_depth=args.sampling_depth,
            beta_cols=[c.strip() for c in args.beta_cols.split(",") if c.strip()],
            time_column=args.time_column,
            include_taxa_barplots=not args.no_taxa_barplots,
            show_qiime=args.show_qiime,
            dry_run=args.dry_run,
        )
        return

    if args.cmd == "diversity-sweep":
        proj = args.project_dir.resolve()
        # Parse column lists
        by_cols = [c.strip() for c in (args.by or "").split(",") if c.strip()] or None
        within_cols = [c.strip() for c in (args.within or "").split(",") if c.strip()] or None
        extra_beta = [c.strip() for c in (args.beta_cols or "").split(",") if c.strip()]
        downstream.run_diversity_sweep(
            proj,
            args.metadata_file,
            by=by_cols,
            within=within_cols,
            min_samples=args.min_samples,
            retain_fraction=args.retain_fraction,
            beta_cols=extra_beta,
            time_column=args.time_column,
            by_only=args.by_only,  # <-- pass through
            if_exists=args.if_exists,
            show_qiime=args.show_qiime,
            dry_run=args.dry_run,
        )
        return


if __name__ == "__main__":
    main()