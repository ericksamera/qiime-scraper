import logging
from pathlib import Path
from typing import Optional, Union, Sequence  # <- Sequence added
from qiime_pipeline.utils.io import run_command
from qiime_pipeline.utils.logger import log_success

logger = logging.getLogger("qiime_pipeline")

def import_data(input_path: Path, output_path: Path, import_type: str = "SampleData[PairedEndSequencesWithQuality]",
               input_format: str = "PairedEndFastqManifestPhred33V2", dry_run: bool = False, show_qiime: bool = False):
    if output_path.exists():
        logger.info("%s already exists.", output_path)
        return
    cmd = [
        "qiime", "tools", "import",
        "--type", import_type,
        "--input-format", input_format,
        "--input-path", str(input_path),
        "--output-path", str(output_path)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
    log_success("Import completed.")

def merge_tables(input_tables: Sequence[Path], output_table: Path,
                 dry_run: bool = False, show_qiime: bool = False):
    cmd = ["qiime", "feature-table", "merge"]
    for table in input_tables:
        cmd += ["--i-tables", str(table)]
    cmd += ["--o-merged-table", str(output_table)]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
    log_success("Feature tables merged.")

def merge_seqs(input_sequences: Sequence[Path], output_sequence: Path,
               dry_run: bool = False, show_qiime: bool = False):
    cmd = ["qiime", "feature-table", "merge-seqs"]
    for seq in input_sequences:
        cmd += ["--i-data", str(seq)]
    cmd += ["--o-merged-data", str(output_sequence)]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
    log_success("Representative sequences merged.")

def cutadapt_trim_paired(input_seqs: Path, forward_primers: Sequence[str],
                         reverse_primers: Sequence[str], output_path: Path,
                         cores: int = 0, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "cutadapt", "trim-paired",
        "--i-demultiplexed-sequences", str(input_seqs),
        "--o-trimmed-sequences", str(output_path),
        "--p-cores", str(cores)
    ]
    for primer in forward_primers:
        cmd += ["--p-front-f", primer]
    for primer in reverse_primers:
        cmd += ["--p-front-r", primer]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
    log_success("Primer trimming completed.")


def dada2_denoise_paired(input_seqs: Path, trunc_len_f: int, trunc_len_r: int, output_table: Path,
                         output_rep_seqs: Path, output_denoising_stats: Path, trim_left_f: int = 0, trim_left_r: int = 0,
                         max_ee_f: int = 2, max_ee_r: int = 2, trunc_q: int = 2, min_overlap: int = 12,
                         pooling_method: str = "independent", chimera_method: str = "consensus",
                         min_fold_parent_over_abundance: int = 1, allow_one_off: bool = False,
                         n_threads: int = 0, n_reads_learn: int = 1000000,
                         hashed_feature_ids: bool = True, retain_all_samples: bool = True,
                         dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", str(input_seqs),
        "--p-trunc-len-f", str(trunc_len_f),
        "--p-trunc-len-r", str(trunc_len_r),
        "--p-trim-left-f", str(trim_left_f),
        "--p-trim-left-r", str(trim_left_r),
        "--p-max-ee-f", str(max_ee_f),
        "--p-max-ee-r", str(max_ee_r),
        "--p-trunc-q", str(trunc_q),
        "--p-min-overlap", str(min_overlap),
        "--p-pooling-method", pooling_method,
        "--p-chimera-method", chimera_method,
        "--p-min-fold-parent-over-abundance", str(min_fold_parent_over_abundance),
        "--p-allow-one-off" if allow_one_off else "--p-no-allow-one-off",
        "--p-n-threads", str(n_threads),
        "--p-n-reads-learn", str(n_reads_learn),
        "--p-hashed-feature-ids" if hashed_feature_ids else "--p-no-hashed-feature-ids",
        "--p-retain-all-samples" if retain_all_samples else "--p-no-retain-all-samples",
        "--o-table", str(output_table),
        "--o-representative-sequences", str(output_rep_seqs),
        "--o-denoising-stats", str(output_denoising_stats)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
    log_success("DADA2 denoise completed.")

def feature_table_summarize(input_table: Path, output: Path, sample_metadata_file: Optional[Path] = None,
                            dry_run: bool = False, show_qiime: bool = False):
    cmd = ["qiime", "feature-table", "summarize", "--i-table", str(input_table), "--o-visualization", str(output)]
    if sample_metadata_file:
        cmd += ["--m-sample-metadata-file", str(sample_metadata_file)]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def feature_table_tabulate_seqs(input_data: Path, output: Path, dry_run: bool = False, show_qiime: bool = False):
    cmd = ["qiime", "feature-table", "tabulate-seqs", "--i-data", str(input_data), "--o-visualization", str(output)]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def feature_table_filter_samples(input_table: Path, metadata_file: Path, where: str,
                                 output_table: Path, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "feature-table", "filter-samples",
        "--i-table", str(input_table),
        "--m-metadata-file", str(metadata_file),
        "--p-where", where,
        "--o-filtered-table", str(output_table)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def classify_sklearn(input_reads: Path, input_classifier: Path, output_classification: Path,
                     reads_per_batch: Union[str, int] = "auto", n_jobs: int = 0, pre_dispatch: str = "2*n_jobs",
                     confidence: float = 0.7, read_orientation: str = "auto",
                     dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-reads", str(input_reads),
        "--i-classifier", str(input_classifier),
        "--o-classification", str(output_classification),
        "--p-reads-per-batch", str(reads_per_batch),
        "--p-n-jobs", str(n_jobs),
        "--p-pre-dispatch", str(pre_dispatch),
        "--p-confidence", str(confidence),
        "--p-read-orientation", str(read_orientation)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def metadata_tabulate(input_file: Path, output_visualization: Path, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "metadata", "tabulate",
        "--m-input-file", str(input_file),
        "--o-visualization", str(output_visualization)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def phylogeny_align_to_tree_mafft_fasttree(input_sequences: Path, output_alignment: Path,
                                          output_masked_alignment: Path, output_tree: Path,
                                          output_rooted_tree: Path, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
        "--i-sequences", str(input_sequences),
        "--o-alignment", str(output_alignment),
        "--o-masked-alignment", str(output_masked_alignment),
        "--o-tree", str(output_tree),
        "--o-rooted-tree", str(output_rooted_tree)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
    log_success("Phylogeny (MAFFT→FastTree) completed.")

def diversity_core_metrics_phylogenetic(input_phylogeny: Path, input_table: Path, sampling_depth: int,
                                       metadata_file: Path, output_dir: Path, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "diversity", "core-metrics-phylogenetic",
        "--i-phylogeny", str(input_phylogeny),
        "--i-table", str(input_table),
        "--p-sampling-depth", str(sampling_depth),
        "--m-metadata-file", str(metadata_file),
        "--output-dir", str(output_dir)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
    log_success("Core diversity metrics completed.")

def diversity_alpha_group_significance(input_alpha_vector: Path, metadata_file: Path, output_visualization: Path,
                                       dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "diversity", "alpha-group-significance",
        "--i-alpha-diversity", str(input_alpha_vector),
        "--m-metadata-file", str(metadata_file),
        "--o-visualization", str(output_visualization)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def diversity_beta_group_significance(input_distance: Path, metadata_file: Path, metadata_column: str,
                                      output_visualization: Path, pairwise: bool = True,
                                      dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "diversity", "beta-group-significance",
        "--i-distance-matrix", str(input_distance),
        "--m-metadata-file", str(metadata_file),
        "--m-metadata-column", metadata_column,
        "--o-visualization", str(output_visualization)
    ]
    if pairwise:
        cmd.append("--p-pairwise")
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def emperor_plot(input_pcoa: Path, metadata_file: Path, output_visualization: Path,
                 custom_axes: Optional[str] = None, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "emperor", "plot",
        "--i-pcoa", str(input_pcoa),
        "--m-metadata-file", str(metadata_file),
        "--o-visualization", str(output_visualization)
    ]
    if custom_axes:
        cmd += ["--p-custom-axes", str(custom_axes)]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def diversity_alpha_rarefaction(input_table: Path, input_phylogeny: Path, max_depth: int, metadata_file: Path,
                                output_visualization: Path, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "diversity", "alpha-rarefaction",
        "--i-table", str(input_table),
        "--i-phylogeny", str(input_phylogeny),
        "--p-max-depth", str(max_depth),
        "--m-metadata-file", str(metadata_file),
        "--o-visualization", str(output_visualization)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)

def taxa_barplot(input_table: Path, input_taxonomy: Path, metadata_file: Path,
                 output_visualization: Path, dry_run: bool = False, show_qiime: bool = False):
    cmd = [
        "qiime", "taxa", "barplot",
        "--i-table", str(input_table),
        "--i-taxonomy", str(input_taxonomy),
        "--m-metadata-file", str(metadata_file),
        "--o-visualization", str(output_visualization)
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_qiime)
