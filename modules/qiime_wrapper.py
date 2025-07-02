# modules/qiime_wrapper.py

import logging
from pathlib import Path

from typing import Union

from .io_utils import run_command
from .logger import log_success

logger = logging.getLogger("qiime_pipeline")

def import_data(
        input_path: Path,
        output_path: Path,
        import_type: str = "SampleData[PairedEndSequencesWithQuality]",
        input_format: str = "PairedEndFastqManifestPhred33V2",
        dry_run: bool = False,
        show_qiime: bool = False,
    ) -> None:

    if input_path is None or output_path is None:
        raise ValueError("Both input_path and output_path must be provided.")

    if output_path.exists() and output_path.is_dir():
        raise ValueError(f"Expected output_path to be a file, not a directory: {output_path}")

    if output_path.exists():
        logger.info(f"{output_path} already exists.")
        return

    run_command([
        "qiime", "tools", "import",
        "--type",         import_type,
        "--input-format", input_format,
        "--input-path",   str(input_path),
        "--output-path",  str(output_path)
    ], dry_run, capture=not show_qiime)

    log_success("Import completed.")

def dada2_denoise_paired(
        input_seqs: Path,
        trunc_len_f: int,
        trunc_len_r: int,

        output_table: Path,
        output_rep_seqs: Path,
        output_denoising_stats: Path,

        trim_left_f: int = 0,
        trim_left_r: int = 0,

        max_ee_f: int = 2,
        max_ee_r: int = 2,
        trunc_q: int = 2,

        min_overlap: int = 12,
        pooling_method: str = "independent",
        chimera_method: str = "consensus",

        min_fold_parent_over_abundance: int = 1,
        allow_one_off: bool = False,

        n_threads: int = 0,
        n_reads_learn: int = 1_000_000,

        hashed_feature_ids: bool = True,
        retain_all_samples: bool = True,

        dry_run: bool = False,
        show_qiime: bool = False
    ) -> None:
    """
    """

    if any(x is None for x in [input_seqs, trunc_len_f, trunc_len_r,]):
        raise ValueError("All of these [input_seqs, trunc_len_f, trunc_len_r] must be provided.")
    
    if any(x is None for x in [output_table, output_rep_seqs, output_denoising_stats,]):
        raise ValueError("All of these [output_table, output_rep_seqs, output_denoising_stats] must be provided.")

    if min_overlap < 4:
        raise ValueError("min_overlap must be at least 4.")

    if pooling_method not in ('independent', 'pseudo'):
        raise ValueError("pooling_method must be in ['independent', 'pseudo'].")

    if chimera_method not in ('consensus', 'none'):
        raise ValueError("chimera_method must be in ['consensus', 'none].")
    
    if min_fold_parent_over_abundance < 1:
        raise ValueError("min_fold_parent_over_abundance should be greater than or equal to 1")

    run_command([
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs",             str(input_seqs),
        "--p-trunc-len-f",                    str(trunc_len_f),
        "--p-trunc-len-r",                    str(trunc_len_r),

        "--p-trim-left-f",                    str(trim_left_f),
        "--p-trim-left-r",                    str(trim_left_r),
        "--p-max-ee-f",                       str(max_ee_f),
        "--p-max-ee-r",                       str(max_ee_r),
        "--p-trunc-q",                        str(trunc_q),
        "--p-min-overlap",                    str(min_overlap),
        "--p-pooling-method",                 pooling_method,
        "--p-chimera-method",                 chimera_method,
        "--p-min-fold-parent-over-abundance", str(min_fold_parent_over_abundance),
        "--p-allow-one-off" if allow_one_off else "--p-no-allow-one-off",
        "--p-n-threads",                      str(n_threads),
        "--p-n-reads-learn",                  str(n_reads_learn),
        "--p-hashed-feature-ids" if hashed_feature_ids else "--p-no-hashed-feature-ids",
        "--p-retain-all-samples" if retain_all_samples else "--p-no-retain-all-samples",

        "--o-table",                          str(output_table), 
        "--o-representative-sequences",       str(output_rep_seqs),
        "--o-denoising-stats",                str(output_denoising_stats)
    ], dry_run, capture=not show_qiime)

    log_success("Import completed.")

def feature_table_summarize(
        input_table: Path,
        output: Path,
        dry_run: bool = False,
        show_qiime: bool = False
    ) -> None:

    run_command([
        "qiime", "feature-table", "summarize",
        "--i-table",         str(input_table),
        "--o-visualization", str(output),
        ], dry_run, capture=not show_qiime)

def classify_sklearn(
        input_reads: Path,
        input_classifier: Path,
        output_classification: Path,

        reads_per_batch: str | int = 'auto',
        n_jobs: int = 0,
        pre_dispatch: str = '2*n_jobs',
        confidence: float = 0.7,
        read_orientation: str = 'auto',

        dry_run: bool = False,
        show_qiime: bool = False
    ) -> None:

    run_command([
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-reads",             str(input_reads),
        "--i-classifier",        str(input_classifier),
        "--o-classification",    str(output_classification),
        "--p-reads-per-batch",   str(reads_per_batch),
        "--p-n-jobs",            str(n_jobs),
        "--p-pre-dispatch",      str(pre_dispatch),
        "--p-confidence",        str(confidence),
        "--p-read-orientation",  str(read_orientation)
        ], dry_run, capture=not show_qiime)

# ---