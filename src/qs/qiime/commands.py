# src/qs/qiime/commands.py
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Union, Mapping

from qs.utils.runner import run_command


# ---------------------------
# Imports & metadata helpers
# ---------------------------

def import_data(
    input_path: Path,
    output_path: Path,
    import_type: str = "SampleData[PairedEndSequencesWithQuality]",
    input_format: Optional[str] = "PairedEndFastqManifestPhred33V2",
    *,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "tools", "import",
        "--type", import_type,
        "--input-path", str(input_path),
        "--output-path", str(output_path),
    ]
    if input_format:
        cmd += ["--input-format", input_format]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def metadata_tabulate(
    input_file: Path,
    output_visualization: Path,
    *,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "metadata", "tabulate",
        "--m-input-file", str(input_file),
        "--o-visualization", str(output_visualization),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


# ---------------------------
# Cutadapt
# ---------------------------

def cutadapt_trim_paired(
    *,
    input_seqs: Path,
    forward_primers: Sequence[str],
    reverse_primers: Sequence[str],
    output_path: Path,
    cores: int = 0,
    discard_untrimmed: bool = False,
    no_indels: bool = False,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "cutadapt", "trim-paired",
        "--i-demultiplexed-sequences", str(input_seqs),
        "--o-trimmed-sequences", str(output_path),
        "--p-cores", str(cores),
    ]
    for p in forward_primers:
        cmd += ["--p-front-f", p]
    for p in reverse_primers:
        cmd += ["--p-front-r", p]
    if discard_untrimmed:
        cmd.append("--p-discard-untrimmed")
    if no_indels:
        cmd.append("--p-no-indels")
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


# ---------------------------
# DADA2
# ---------------------------

def dada2_denoise_paired(
    *,
    input_seqs: Path,
    trunc_len_f: int,
    trunc_len_r: int,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    trim_left_f: int = 0,
    trim_left_r: int = 0,
    max_ee_f: int = 2,
    max_ee_r: int = 2,
    trunc_q: int = 2,
    min_overlap: int = 12,
    pooling_method: str = "independent",
    chimera_method: str = "consensus",
    n_threads: int = 0,
    n_reads_learn: Optional[int] = None,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
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
        "--p-pooling-method", str(pooling_method),
        "--p-chimera-method", str(chimera_method),
        "--p-n-threads", str(n_threads),
        "--o-table", str(output_table),
        "--o-representative-sequences", str(output_rep_seqs),
        "--o-denoising-stats", str(output_stats),
    ]
    if n_reads_learn is not None:
        cmd += ["--p-n-reads-learn", str(n_reads_learn)]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


# ---------------------------
# Feature-table utilities
# ---------------------------

def feature_table_summarize(
    *,
    input_table: Path,
    output: Path,
    sample_metadata_file: Optional[Path] = None,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "feature-table", "summarize",
        "--i-table", str(input_table),
        "--o-visualization", str(output),
    ]
    if sample_metadata_file:
        cmd += ["--m-sample-metadata-file", str(sample_metadata_file)]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def feature_table_tabulate_seqs(
    *,
    input_data: Path,
    output: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "feature-table", "tabulate-seqs",
        "--i-data", str(input_data),
        "--o-visualization", str(output),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def merge_tables(
    input_tables: Sequence[Path],
    output_table: Path,
    *,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = ["qiime", "feature-table", "merge"]
    for t in input_tables:
        cmd += ["--i-tables", str(t)]
    cmd += ["--o-merged-table", str(output_table)]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def merge_seqs(
    input_sequences: Sequence[Path],
    output_sequence: Path,
    *,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = ["qiime", "feature-table", "merge-seqs"]
    for s in input_sequences:
        cmd += ["--i-data", str(s)]
    cmd += ["--o-merged-data", str(output_sequence)]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


# ---------------------------
# Filtering (NEW)
# ---------------------------

def feature_table_filter_features(
    *,
    input_table: Path,
    output_table: Path,
    min_frequency: int = 0,
    min_samples: int = 1,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "feature-table", "filter-features",
        "--i-table", str(input_table),
        "--o-filtered-table", str(output_table),
        "--p-min-frequency", str(int(min_frequency)),
        "--p-min-samples", str(int(min_samples)),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def taxa_filter_table(
    *,
    input_table: Path,
    input_taxonomy: Path,
    output_table: Path,
    include: Optional[str] = None,
    exclude: Optional[str] = None,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "taxa", "filter-table",
        "--i-table", str(input_table),
        "--i-taxonomy", str(input_taxonomy),
        "--o-filtered-table", str(output_table),
    ]
    if include:
        cmd += ["--p-include", include]
    if exclude:
        cmd += ["--p-exclude", exclude]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def feature_table_filter_samples(
    *,
    input_table: Path,
    output_table: Path,
    min_frequency: int,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "feature-table", "filter-samples",
        "--i-table", str(input_table),
        "--o-filtered-table", str(output_table),
        "--p-min-frequency", str(int(min_frequency)),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def feature_table_filter_seqs(
    *,
    input_data: Path,
    input_table: Path,
    output_data: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "feature-table", "filter-seqs",
        "--i-data", str(input_data),
        "--i-table", str(input_table),
        "--o-filtered-data", str(output_data),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


# ---------------------------
# Classification
# ---------------------------

def classify_sklearn(
    *,
    input_reads: Path,
    input_classifier: Path,
    output_classification: Path,
    reads_per_batch: Union[str, int] = 1000,
    n_jobs: int = 1,
    pre_dispatch: Optional[str] = None,
    confidence: float = 0.7,
    read_orientation: str = "auto",
    dry_run: bool = False,
    show_stdout: bool = False,
    extra_env: Optional[Mapping[str, str]] = None,
) -> None:
    if pre_dispatch is None:
        pre_dispatch = "1*n_jobs"
    cmd: list[str] = [
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-reads", str(input_reads),
        "--i-classifier", str(input_classifier),
        "--o-classification", str(output_classification),
        "--p-reads-per-batch", str(reads_per_batch),
        "--p-n-jobs", str(n_jobs),
        "--p-pre-dispatch", str(pre_dispatch),
        "--p-confidence", str(confidence),
        "--p-read-orientation", str(read_orientation),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout, env=dict(extra_env) if extra_env else None)


# ---------------------------
# Phylogeny & diversity
# ---------------------------

def phylogeny_align_to_tree_mafft_fasttree(
    *,
    input_sequences: Path,
    output_alignment: Path,
    output_masked_alignment: Path,
    output_tree: Path,
    output_rooted_tree: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
        "--i-sequences", str(input_sequences),
        "--o-alignment", str(output_alignment),
        "--o-masked-alignment", str(output_masked_alignment),
        "--o-tree", str(output_tree),
        "--o-rooted-tree", str(output_rooted_tree),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def diversity_core_metrics_phylogenetic(
    *,
    input_phylogeny: Path,
    input_table: Path,
    sampling_depth: int,
    metadata_file: Path,
    output_dir: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "diversity", "core-metrics-phylogenetic",
        "--i-phylogeny", str(input_phylogeny),
        "--i-table", str(input_table),
        "--p-sampling-depth", str(sampling_depth),
        "--m-metadata-file", str(metadata_file),
        "--output-dir", str(output_dir),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


# ---------------------------
# Taxonomy visuals & stats
# ---------------------------

def taxa_barplot(
    *,
    input_table: Path,
    input_taxonomy: Path,
    metadata_file: Path,
    output_visualization: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd: list[str] = [
        "qiime", "taxa", "barplot",
        "--i-table", str(input_table),
        "--i-taxonomy", str(input_taxonomy),
        "--m-metadata-file", str(metadata_file),
        "--o-visualization", str(output_visualization),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def feature_table_rarefy(
    *,
    input_table: Path,
    sampling_depth: int,
    output_table: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "feature-table", "rarefy",
        "--i-table", str(input_table),
        "--p-sampling-depth", str(int(sampling_depth)),
        "--o-rarefied-table", str(output_table),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


# ---------------------------
# Alpha diversity (NEW)
# ---------------------------

def diversity_alpha(
    *,
    input_table: Path,
    metric: str,
    output_vector: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "diversity", "alpha",
        "--i-table", str(input_table),
        "--p-metric", str(metric),
        "--o-alpha-diversity", str(output_vector),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def diversity_alpha_phylogenetic(
    *,
    input_table: Path,
    input_phylogeny: Path,
    metric: str,
    output_vector: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "diversity", "alpha-phylogenetic",
        "--i-table", str(input_table),
        "--i-phylogeny", str(input_phylogeny),
        "--p-metric", str(metric),
        "--o-alpha-diversity", str(output_vector),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def diversity_alpha_group_significance(
    *,
    alpha_diversity: Path,
    metadata_file: Path,
    output_visualization: Path,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "diversity", "alpha-group-significance",
        "--i-alpha-diversity", str(alpha_diversity),
        "--m-metadata-file", str(metadata_file),
        "--o-visualization", str(output_visualization),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def diversity_alpha_correlation(
    *,
    alpha_diversity: Path,
    metadata_file: Path,
    output_visualization: Path,
    method: str = "spearman",
    intersect_ids: bool = False,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "diversity", "alpha-correlation",
        "--i-alpha-diversity", str(alpha_diversity),
        "--m-metadata-file", str(metadata_file),
        "--p-method", method,
        "--o-visualization", str(output_visualization),
    ]
    if intersect_ids:
        cmd.append("--p-intersect-ids")
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def diversity_beta_group_significance(
    *,
    distance_matrix: Path,
    metadata_file: Path,
    metadata_column: str,
    output_visualization: Path,
    method: str = "permanova",
    permutations: int = 999,
    pairwise: bool = False,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "diversity", "beta-group-significance",
        "--i-distance-matrix", str(distance_matrix),
        "--m-metadata-file", str(metadata_file),
        "--m-metadata-column", metadata_column,
        "--p-method", method,
        "--p-permutations", str(permutations),
        "--o-visualization", str(output_visualization),
    ]
    if pairwise:
        cmd.append("--p-pairwise")
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


def diversity_adonis(
    *,
    distance_matrix: Path,
    metadata_file: Path,
    formula: str,
    output_visualization: Path,
    permutations: int = 999,
    n_jobs: int = 1,
    dry_run: bool = False,
    show_stdout: bool = False,
) -> None:
    cmd = [
        "qiime", "diversity", "adonis",
        "--i-distance-matrix", str(distance_matrix),
        "--m-metadata-file", str(metadata_file),
        "--p-formula", formula,
        "--p-permutations", str(permutations),
        "--p-n-jobs", str(n_jobs),
        "--o-visualization", str(output_visualization),
    ]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)
