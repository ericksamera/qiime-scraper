import logging
from pathlib import Path
from qiime_pipeline.utils import io, qiime
from qiime_pipeline.importing import Importer
from qiime_pipeline.preprocessing import Preprocessor
from qiime_pipeline.denoising import Denoiser
from qiime_pipeline.classification import Classifier
from qiime_pipeline.analysis import AnalysisRunner
from qiime_pipeline.sweep import DiversitySweeper

logger = logging.getLogger("qiime_pipeline")

class PipelineRunner:
    """Coordinates the end-to-end execution of the QIIME2 pipeline."""
    def __init__(self, config,
                 importer=Importer, denoiser=Denoiser,
                 classifier=Classifier, analysis=AnalysisRunner, sweeper=DiversitySweeper):
        self.config = config
        self.importer = importer(self.config)
        self.preprocessor = Preprocessor(self.config) if (Preprocessor and not getattr(self.config, "no_trim_primers", False)) else None
        self.denoiser = denoiser(self.config)
        self.classifier = classifier(self.config) if classifier else None
        self.analysis = analysis(self.config) if analysis else None
        self.sweeper = sweeper(self.config) if sweeper else None

    def run(self):
        Path(self.config.project_dir).mkdir(parents=True, exist_ok=True)
        fastq_dir = Path(self.config.fastq_dir)

        # Detect multiple runs by subfolders
        run_dirs = [d for d in fastq_dir.iterdir() if d.is_dir()]
        multi_run = len(run_dirs) > 0
        if not multi_run:
            run_dirs = [fastq_dir]

        tables, rep_seqs, stats = [], [], []
        for run_dir in run_dirs:
            run_name = run_dir.name
            logger.info(f"Processing run: {run_name}")

            # Import FASTQs (per-run demux if multi-run)
            if multi_run:
                manifest_path = io.generate_manifest(run_dir, Path(self.config.project_dir))
                demux_path = Path(self.config.project_dir) / f"{run_name}_demux.qza"
                qiime.import_data(
                    input_path=manifest_path,
                    output_path=demux_path,
                    import_type="SampleData[PairedEndSequencesWithQuality]",
                    input_format="PairedEndFastqManifestPhred33V2",
                    dry_run=self.config.dry_run,
                    show_qiime=self.config.show_qiime
                )
            else:
                self.importer.run()
                demux_path = Path(self.config.project_dir) / "output.qza"

            # Optional primer trimming if per-sample primers available
            if self.preprocessor:
                demux_path = self.preprocessor.trim_if_needed(demux_path, run_dir)

            # Denoising per run
            table_qza = Path(self.config.project_dir) / (f"{run_name}_table.qza" if multi_run else "table.qza")
            rep_seqs_qza = Path(self.config.project_dir) / (f"{run_name}_rep-seqs.qza" if multi_run else "rep-seqs.qza")
            stats_qza = Path(self.config.project_dir) / (f"{run_name}_denoising-stats.qza" if multi_run else "denoising-stats.qza")
            self.denoiser.run(input_path=demux_path, output_table=table_qza,
                               output_rep_seqs=rep_seqs_qza, output_stats=stats_qza)
            tables.append(table_qza); rep_seqs.append(rep_seqs_qza); stats.append(stats_qza)

        # Merge results if multi-run
        if multi_run:
            merged_table = Path(self.config.project_dir) / "table.qza"
            merged_rep_seqs = Path(self.config.project_dir) / "rep-seqs.qza"
            self.denoiser.merge_results(tables, merged_table, rep_seqs, merged_rep_seqs)
            qiime.feature_table_summarize(
                input_table=merged_table,
                output=Path(self.config.project_dir) / "table.qzv",
                sample_metadata_file=getattr(self.config, 'metadata_file', None),
                dry_run=self.config.dry_run,
                show_qiime=self.config.show_qiime
            )
            qiime.feature_table_tabulate_seqs(
                input_data=merged_rep_seqs,
                output=Path(self.config.project_dir) / "rep-seqs.qzv",
                dry_run=self.config.dry_run,
                show_qiime=self.config.show_qiime
            )
        else:
            # Summaries for convenience
            final_table = tables[0]; final_rep_seqs = rep_seqs[0]
            if final_table.exists():
                qiime.feature_table_summarize(
                    input_table=final_table,
                    output=Path(self.config.project_dir) / "table.qzv",
                    sample_metadata_file=getattr(self.config, 'metadata_file', None),
                    dry_run=self.config.dry_run,
                    show_qiime=self.config.show_qiime
                )
            if final_rep_seqs.exists():
                qiime.feature_table_tabulate_seqs(
                    input_data=final_rep_seqs,
                    output=Path(self.config.project_dir) / "rep-seqs.qzv",
                    dry_run=self.config.dry_run,
                    show_qiime=self.config.show_qiime
                )

        # Continue with downstream analysis if configured (not explicitly asked in this runner)
        # For example, one might call AnalysisRunner or DiversitySweeper after pipeline run.
        # (No auto-trigger of downstream steps in PipelineRunner.run by default)

