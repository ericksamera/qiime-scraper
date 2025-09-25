# src/qs/qiime/commands.py
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Union, Mapping

from qs.utils.runner import run_command


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
        "--p-pooling-method", pooling_method,
        "--p-chimera-method", chimera_method,
        "--p-n-threads", str(n_threads),
        "--o-table", str(output_table),
        "--o-representative-sequences", str(output_rep_seqs),
        "--o-denoising-stats", str(output_stats),
    ]
    if n_reads_learn is not None:
        cmd += ["--p-n-reads-learn", str(n_reads_learn)]
    run_command(cmd, dry_run=dry_run, capture=not show_stdout)


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


def classify_sklearn(
    *,
    input_reads: Path,
    input_classifier: Path,
    output_classification: Path,
    reads_per_batch: Union[str, int] = 1000,     # safer default than 'auto'
    n_jobs: int = 1,                              # safer default than 0 (all cores)
    pre_dispatch: Optional[str] = None,           # default to '1*n_jobs'
    confidence: float = 0.7,
    read_orientation: str = "auto",
    dry_run: bool = False,
    show_stdout: bool = False,
    extra_env: Optional[Mapping[str, str]] = None,
) -> None:
    """
    Memory-safer wrapper for sklearn classifier.
    Set n_jobs small, reads_per_batch small, and pin BLAS threads via extra_env.
    """
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
