import argparse
from pathlib import Path

# Import refactored pipeline components
import sys
sys.path.append(str(Path(__file__).parent))
from qiime_pipeline.models import ProjectConfig
from qiime_pipeline.pipeline_runner import PipelineRunner
from qiime_pipeline.importing import Importer
from qiime_pipeline.denoising import Denoiser
from qiime_pipeline.classification import Classifier
from qiime_pipeline.analysis import AnalysisRunner
from qiime_pipeline.sweep import DiversitySweeper
from qiime_pipeline.utils import logger

__VERSION__ = "1.0.0"
__AUTHOR__ = "Erick Samera (erick.samera@kpu.ca)"

def parse_args():
    p = argparse.ArgumentParser(prog="pipeline", epilog=f"{__VERSION__} | {__AUTHOR__}")
    # Global options
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--show-qiime", dest="show_qiime", action="store_true", help="Stream QIIME logs live (default).")
    common.add_argument("--no-show-qiime", dest="show_qiime", action="store_false", help="Do not stream QIIME logs live; capture and print on error.")
    common.set_defaults(show_qiime=True)
    common.add_argument("--dry-run", action="store_true", help="Print commands without executing")

    # Shared options for analysis commands
    analysis_common = argparse.ArgumentParser(add_help=False)
    analysis_common.add_argument("--project-dir", type=Path, required=True, help="Project directory for analysis results")
    analysis_common.add_argument("--metadata-file", type=Path, required=True, help="Sample metadata TSV file for the project")

    sp = p.add_subparsers(dest="cmd", required=True)

    # Helper to add FASTQ and project dir options (for commands that need raw FASTQs)
    def add_fastq_ops(sub):
        sub.add_argument("--fastq-dir", type=Path, required=True, help="Project FASTQ directory")
        sub.add_argument("--project-dir", type=Path, required=True, help="Output project directory (pipeline writes here)")

    # Subcommands definitions
    run = sp.add_parser("run", help="Run full pipeline", parents=[common])
    add_fastq_ops(run)
    run.add_argument("--classifiers-dir", type=Path, default=Path("/home/erick/qiime"), help="Directory containing classifier .qza files")
    run.add_argument("--no-trim-primers", action="store_true", help="Skip primer trimming even if primers are present in metadata")

    imp = sp.add_parser("import", help="Only import FASTQs to a QIIME artifact", parents=[common])
    add_fastq_ops(imp)

    trm = sp.add_parser("trim", help="Only search for optimal trimming parameters", parents=[common])
    add_fastq_ops(trm)

    cls = sp.add_parser("classify", help="Run classifier sweep on rep_seqs in <project-dir>", parents=[common])
    cls.add_argument("--project-dir", type=Path, required=True, help="Project directory containing rep_seqs artifact (from trimming)")
    cls.add_argument("--input-reads", type=Path, default=None, help="Optional FeatureData[Sequence] .qza to classify (overrides best rep_seqs)")
    cls.add_argument("--classifiers-dir", type=Path, default=Path("/home/erick/qiime"))

    stg = sp.add_parser("stage", help="Stage winners into analysis/ and build summaries", parents=[common])
    stg.add_argument("--project-dir", type=Path, required=True)
    stg.add_argument("--metadata-file", type=Path, required=False)
    stg.add_argument("--preferred-metric", type=str, default="pct_depth≥7", help="Which classifier metric to prefer for choosing the winner")

    phy = sp.add_parser("phylogeny", help="Build phylogeny (MAFFT→FastTree→root) in analysis/", parents=[common])
    phy.add_argument("--project-dir", type=Path, required=True)

    div = sp.add_parser("diversity", help="Run core phylogenetic diversity + tests + Emperor", parents=[common, analysis_common])
    div.add_argument("--sampling-depth", type=int, default=None)
    div.add_argument("--beta-cols", type=str, default="body-site,subject", help="Comma-separated metadata columns for beta-group-significance")
    div.add_argument("--time-column", type=str, default=None)
    div.add_argument("--skip-alpha-tests", action="store_true", help="Skip alpha diversity significance tests")
    div.add_argument("--skip-beta-tests", action="store_true", help="Skip beta diversity significance tests")
    div.add_argument("--skip-emperor", action="store_true", help="Skip Emperor plot generation")
    div.add_argument("--if-exists", choices=["skip", "overwrite", "new", "error"], default="skip",
                     help="What to do if core-metrics output already exists (default: skip)")

    arf = sp.add_parser("alpha-rarefaction", help="Alpha rarefaction visualization", parents=[common, analysis_common])
    arf.add_argument("--max-depth", type=int, default=None)

    tx = sp.add_parser("taxonomy", help="Taxonomy bar plots", parents=[common, analysis_common])
    # (No additional arguments needed for taxonomy beyond project-dir and metadata-file)

    dwn = sp.add_parser("downstream", help="Stage winners → phylogeny → core diversity → rarefaction → taxa barplots", parents=[common, analysis_common])
    dwn.add_argument("--preferred-metric", type=str, default="pct_depth≥7")
    dwn.add_argument("--sampling-depth", type=int, default=None)
    dwn.add_argument("--beta-cols", type=str, default="body-site,subject")
    dwn.add_argument("--time-column", type=str, default=None)
    dwn.add_argument("--no-taxa-barplots", action="store_true", help="Skip generating taxa barplots")

    sweep = sp.add_parser("diversity-sweep", help="Iterate core diversity over metadata subsets (nested or one-level)", parents=[common, analysis_common])
    sweep.add_argument("--by", type=str, default=None, help="Comma-separated outer loop columns (if omitted, use all categorical columns)")
    sweep.add_argument("--within", type=str, default=None, help="Comma-separated inner loop columns (ignored if --by-only is set)")
    sweep.add_argument("--by-only", action="store_true", help="Use a single-level sweep (iterate only over --by columns)")
    sweep.add_argument("--min-samples", type=int, default=5, help="Skip subsets with < N samples (default: 5)")
    sweep.add_argument("--retain-fraction", type=float, default=0.90, help="Fraction of samples to retain when choosing rarefaction depth")
    sweep.add_argument("--beta-cols", type=str, default=None, help="Comma-separated metadata columns to use for beta diversity significance")
    sweep.add_argument("--time-column", type=str, default=None)
    sweep.add_argument("--if-exists", choices=["skip", "overwrite", "new", "error"], default="skip",
                       help="How to handle existing subset outputs (default: skip)")

    return p.parse_args()

def main():
    args = parse_args()
    logger.setup_logger()  # initialize logging to console and file

    # Commands that require FASTQ inputs
    if args.cmd in {"run", "import", "trim"}:
        # Prepare project configuration
        config = ProjectConfig(
            fastq_dir=args.fastq_dir.resolve(),
            project_dir=args.project_dir.resolve(),
            classifiers_dir=(args.classifiers_dir.resolve() if hasattr(args, "classifiers_dir") else None),
            dry_run=args.dry_run,
            show_qiime=args.show_qiime,
            no_trim_primers=(args.no_trim_primers if hasattr(args, "no_trim_primers") else False),
        )

        if args.cmd == "import":
            from qiime_pipeline.importing import Importer
            Importer(config).run()
            return
        if args.cmd == "trim":
            from qiime_pipeline.trimming import TrimmingOptimizer
            TrimmingOptimizer(config).find_optimal_trimming()
            return
        if args.cmd == "run":
            PipelineRunner(config).run()
            return

    # Classification only (no FASTQ required if rep-seqs provided)
    if args.cmd == "classify":
        config = ProjectConfig(
            project_dir=args.project_dir.resolve(),
            classifiers_dir=args.classifiers_dir.resolve(),
            input_reads=(args.input_reads.resolve() if args.input_reads else None),
            dry_run=args.dry_run,
            show_qiime=args.show_qiime
        )
        Classifier(config).find_optimal_classifier()
        return

    # Downstream analysis commands
    if args.cmd == "stage":
        config = ProjectConfig(project_dir=args.project_dir.resolve(), metadata_file=args.metadata_file, dry_run=args.dry_run, show_qiime=args.show_qiime)
        AnalysisRunner(config).stage_winners(preferred_metric=args.preferred_metric)
        return

    if args.cmd == "phylogeny":
        config = ProjectConfig(project_dir=args.project_dir.resolve(), dry_run=args.dry_run, show_qiime=args.show_qiime)
        AnalysisRunner(config).build_phylogeny()
        return

    if args.cmd == "diversity":
        config = ProjectConfig(project_dir=args.project_dir.resolve(), metadata_file=args.metadata_file, dry_run=args.dry_run, show_qiime=args.show_qiime)
        beta_list = [c.strip() for c in args.beta_cols.split(",") if c.strip()]
        AnalysisRunner(config).run_core_diversity(
            sampling_depth=args.sampling_depth,
            beta_cols=beta_list,
            time_column=args.time_column,
            include_alpha_tests=(not args.skip_alpha_tests),
            include_beta_tests=(not args.skip_beta_tests),
            include_emperor=(not args.skip_emperor),
            if_exists=args.if_exists
        )
        return

    if args.cmd == "alpha-rarefaction":
        config = ProjectConfig(project_dir=args.project_dir.resolve(), metadata_file=args.metadata_file, dry_run=args.dry_run, show_qiime=args.show_qiime)
        AnalysisRunner(config).run_alpha_rarefaction(max_depth=args.max_depth)
        return

    if args.cmd == "taxonomy":
        config = ProjectConfig(project_dir=args.project_dir.resolve(), metadata_file=args.metadata_file, dry_run=args.dry_run, show_qiime=args.show_qiime)
        AnalysisRunner(config).run_taxa_barplots()
        return

    if args.cmd == "downstream":
        config = ProjectConfig(project_dir=args.project_dir.resolve(), metadata_file=args.metadata_file, dry_run=args.dry_run, show_qiime=args.show_qiime)
        AnalysisRunner(config).run_full_analysis(
            preferred_metric=args.preferred_metric,
            sampling_depth=args.sampling_depth,
            beta_cols=args.beta_cols,
            time_column=args.time_column,
            include_taxa_barplots=(not args.no_taxa_barplots)
        )
        return

    if args.cmd == "diversity-sweep":
        config = ProjectConfig(project_dir=args.project_dir.resolve(), metadata_file=args.metadata_file, dry_run=args.dry_run, show_qiime=args.show_qiime)
        # Parse comma-separated lists for columns
        by_cols = [c.strip() for c in (args.by or "").split(",") if c.strip()] or None
        within_cols = [c.strip() for c in (args.within or "").split(",") if c.strip()] or None
        beta_list = [c.strip() for c in (args.beta_cols or "").split(",") if c.strip()]
        DiversitySweeper(config).run_diversity_sweep(
            by=by_cols,
            within=within_cols,
            min_samples=args.min_samples,
            retain_fraction=args.retain_fraction,
            beta_cols=(beta_list or None),
            time_column=args.time_column,
            by_only=args.by_only,
            if_exists=args.if_exists
        )
        return

if __name__ == "__main__":
    main()
