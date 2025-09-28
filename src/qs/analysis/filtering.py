# src/qs/analysis/filtering.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

from qs.utils.logger import get_logger
from qs.analysis.depth import mean_sample_depth
from qs.analysis.staging import _stage_artifact
from qs.qiime import commands as qiime

LOG = get_logger("filtering")


def compute_min_freq_from_qzv(table_qzv: Path, frac_of_mean: float, fallback_min: int = 1) -> int:
    """
    Compute a min-feature frequency threshold as fraction of mean sample depth.
    E.g., frac_of_mean=0.001 ~ 0.1% of mean depth (Microbiome Helper heuristic).
    """
    md = mean_sample_depth(table_qzv)
    if md <= 0:
        return max(1, int(fallback_min))
    return max(1, int(round(md * max(0.0, frac_of_mean))))


def filter_table_and_seqs(
    *,
    group_dir: Path,
    table_qza: Path,
    rep_seqs_qza: Path,
    table_qzv: Path,
    taxonomy_qza: Optional[Path],
    # knobs
    min_freq: int = 0,
    min_samples: int = 1,
    contam_exclude: Optional[str] = "mitochondria,chloroplast",
    keep_only_phylum_classified: bool = False,   # keep entries containing 'p__'
    min_sample_depth: int = 0,
    metadata_augmented: Optional[Path] = None,
    dry_run: bool = False,
    show_qiime: bool = False,
) -> Dict[str, Path]:
    """
    Run a basic filtering cascade and return paths to final artifacts:
      - table_filt.qza                 (rare features removed)
      - table_filt_contam.qza          (contaminants removed, optional)
      - table_final.qza                (low-depth samples removed, optional)
      - rep_seqs_final.qza             (subset representative sequences)
      - table_final_summary.qzv
    """
    group_dir.mkdir(parents=True, exist_ok=True)

    # 1) Rare features
    tbl_rares = group_dir / "table_filt.qza"
    if min_freq > 0 or min_samples > 1:
        qiime.feature_table_filter_features(
            input_table=table_qza,
            output_table=tbl_rares,
            min_frequency=int(min_freq),
            min_samples=int(min_samples),
            dry_run=dry_run,
            show_stdout=show_qiime,
        )
    else:
        _stage_artifact(table_qza, tbl_rares)

    # 2) Contaminants / phylum-classified filter via taxonomy (if available)
    tbl_tax = group_dir / "table_filt_contam.qza"
    if taxonomy_qza and taxonomy_qza.exists():
        include = "p__" if keep_only_phylum_classified else None
        qiime.taxa_filter_table(
            input_table=tbl_rares,
            input_taxonomy=taxonomy_qza,
            output_table=tbl_tax,
            include=include,
            exclude=(contam_exclude or None),
            dry_run=dry_run,
            show_stdout=show_qiime,
        )
    else:
        _stage_artifact(tbl_rares, tbl_tax)

    # 3) Low-depth samples
    tbl_final = group_dir / "table_final.qza"
    if min_sample_depth and int(min_sample_depth) > 0:
        qiime.feature_table_filter_samples(
            input_table=tbl_tax,
            output_table=tbl_final,
            min_frequency=int(min_sample_depth),
            dry_run=dry_run,
            show_stdout=show_qiime,
        )
    else:
        _stage_artifact(tbl_tax, tbl_final)

    # 4) Subset representative sequences to retained features
    rep_final = group_dir / "rep_seqs_final.qza"
    qiime.feature_table_filter_seqs(
        input_data=rep_seqs_qza,
        input_table=tbl_final,
        output_data=rep_final,
        dry_run=dry_run,
        show_stdout=show_qiime,
    )

    # 5) Summarize final table
    qzv = group_dir / "table_final_summary.qzv"
    qiime.feature_table_summarize(
        input_table=tbl_final,
        output=qzv,
        sample_metadata_file=metadata_augmented,
        dry_run=dry_run,
        show_stdout=show_qiime,
    )

    LOG.info("Filtering complete (min_freq=%s, contam_exclude=%s, keep_p__=%s, min_sample_depth=%s)",
             min_freq, contam_exclude, keep_only_phylum_classified, min_sample_depth)

    return {
        "table_filt": tbl_rares,
        "table_filt_contam": tbl_tax,
        "table_final": tbl_final,
        "rep_seqs_final": rep_final,
        "table_final_summary": qzv,
    }
