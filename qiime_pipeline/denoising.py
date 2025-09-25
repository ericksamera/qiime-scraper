import logging
from pathlib import Path
from qiime_pipeline.utils import qiime

logger = logging.getLogger("qiime_pipeline")

class Denoiser:
    """Handles DADA2 denoising of demultiplexed sequences (single run)."""
    def __init__(self, config):
        self.config = config

    def run(self, input_path: Path = None,
            output_table: Path = None,
            output_rep_seqs: Path = None,
            output_stats: Path = None,
            **kwargs):
        """Run DADA2 denoise-paired to produce feature table, rep-seqs, and stats artifacts."""
        project_dir = Path(self.config.project_dir)
        if input_path is None:
            input_path = getattr(self.config, 'demux_artifact', project_dir / "output.qza")
        if output_table is None:
            output_table = project_dir / "table.qza"
        if output_rep_seqs is None:
            output_rep_seqs = project_dir / "rep-seqs.qza"
        if output_stats is None:
            output_stats = project_dir / "denoising-stats.qza"

        trunc_len_f = kwargs.get('trunc_len_f', getattr(self.config, 'trunc_len_f', 0))
        trunc_len_r = kwargs.get('trunc_len_r', getattr(self.config, 'trunc_len_r', 0))
        trim_left_f = kwargs.get('trim_left_f', getattr(self.config, 'trim_left_f', 0))
        trim_left_r = kwargs.get('trim_left_r', getattr(self.config, 'trim_left_r', 0))
        max_ee_f = kwargs.get('max_ee_f', getattr(self.config, 'max_ee_f', 2))
        max_ee_r = kwargs.get('max_ee_r', getattr(self.config, 'max_ee_r', 2))
        trunc_q   = kwargs.get('trunc_q',   getattr(self.config, 'trunc_q', 2))

        logger.info(f"Running DADA2 denoise-paired on {input_path} (trunc_len_f={trunc_len_f}, trunc_len_r={trunc_len_r})")
        qiime.dada2_denoise_paired(
            input_seqs=input_path,
            trunc_len_f=trunc_len_f,
            trunc_len_r=trunc_len_r,
            trim_left_f=trim_left_f,
            trim_left_r=trim_left_r,
            max_ee_f=max_ee_f,
            max_ee_r=max_ee_r,
            trunc_q=trunc_q,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            output_denoising_stats=output_stats,
            dry_run=getattr(self.config, 'dry_run', False),
            show_qiime=getattr(self.config, 'show_qiime', False)
        )
        logger.info(f"Denoising complete. Outputs: {output_table.name}, {output_rep_seqs.name}")
        self.config.table_artifact = output_table
        self.config.rep_seqs_artifact = output_rep_seqs
        self.config.denoising_stats_artifact = output_stats

    def merge_results(self, table_artifacts: list[Path], merged_table_path: Path,
                      rep_seqs_artifacts: list[Path], merged_rep_seqs_path: Path):
        """Merge feature tables and rep-seqs from multiple runs into single artifacts."""
        logger.info(f"Merging {len(table_artifacts)} tables into {merged_table_path}")
        from qiime_pipeline.utils import qiime as _q
        _q.merge_tables(table_artifacts, merged_table_path,
                        dry_run=getattr(self.config, 'dry_run', False),
                        show_qiime=getattr(self.config, 'show_qiime', False))
        logger.info(f"Merging {len(rep_seqs_artifacts)} representative sequence sets into {merged_rep_seqs_path}")
        _q.merge_seqs(rep_seqs_artifacts, merged_rep_seqs_path,
                      dry_run=getattr(self.config, 'dry_run', False),
                      show_qiime=getattr(self.config, 'show_qiime', False))
        logger.info(f"Merged outputs: table -> {merged_table_path}, rep-seqs -> {merged_rep_seqs_path}")
        self.config.table_artifact = merged_table_path
        self.config.rep_seqs_artifact = merged_rep_seqs_path
