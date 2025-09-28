# src/qs/analysis/stats.py
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from qs.utils.logger import get_logger
from qs.qiime import commands as qiime

LOG = get_logger("stats")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def run_diversity_stats_for_group(
    *,
    group_dir: Path,
    metadata_file: Path,
    beta_group_cols: Optional[Sequence[str]] = None,
    beta_method: str = "permanova",   # 'permanova' | 'anosim' | 'permdisp'
    pairwise: bool = False,
    permutations: int = 999,
    adonis_formula: Optional[str] = None,  # e.g. "treatment+block" or "treatment*block"
    run_alpha_correlation: bool = True,
    overwrite: bool = False,
    dry_run: bool = False,
    show_qiime: bool = False,
) -> Path:
    """
    Run:
      - alpha-group-significance (for each alpha vector present)
      - alpha-correlation       (optional, for each alpha vector)
      - beta-group-significance (for each DM present × each requested metadata column)
      - adonis                  (optional, on Bray-Curtis by default)

    Writes outputs to <group_dir>/stats/.
    """
    core = group_dir / "core-metrics-phylo"
    out_root = group_dir / "stats"
    _ensure_dir(out_root)

    # --- Alpha stats ---
    alpha_vectors = {
        "observed_features": core / "observed_features_vector.qza",
        "evenness":          core / "evenness_vector.qza",
        "faith_pd":          core / "faith_pd_vector.qza",
        "shannon":           core / "shannon_vector.qza",
    }

    for name, qza in alpha_vectors.items():
        if not qza.exists():
            continue
        ags_qzv = out_root / f"alpha-group-significance_{name}.qzv"
        if overwrite or not ags_qzv.exists():
            qiime.diversity_alpha_group_significance(
                alpha_diversity=qza,
                metadata_file=metadata_file,
                output_visualization=ags_qzv,
                dry_run=dry_run,
                show_stdout=show_qiime,
            )
        if run_alpha_correlation:
            ac_qzv = out_root / f"alpha-correlation_{name}.qzv"
            if overwrite or not ac_qzv.exists():
                qiime.diversity_alpha_correlation(
                    alpha_diversity=qza,
                    metadata_file=metadata_file,
                    method="spearman",
                    intersect_ids=False,
                    output_visualization=ac_qzv,
                    dry_run=dry_run,
                    show_stdout=show_qiime,
                )

    # --- Beta group significance ---
    beta_mats = {
        "bray_curtis":        core / "bray_curtis_distance_matrix.qza",
        "jaccard":            core / "jaccard_distance_matrix.qza",
        "unweighted_unifrac": core / "unweighted_unifrac_distance_matrix.qza",
        "weighted_unifrac":   core / "weighted_unifrac_distance_matrix.qza",
    }

    cols: List[str] = [c.strip() for c in (beta_group_cols or []) if c.strip()]
    for metric, dm in beta_mats.items():
        if not dm.exists() or not cols:
            continue
        for col in cols:
            bg_qzv = out_root / f"beta-group-significance_{metric}_{col}.qzv"
            if overwrite or not bg_qzv.exists():
                qiime.diversity_beta_group_significance(
                    distance_matrix=dm,
                    metadata_file=metadata_file,
                    metadata_column=col,
                    method=beta_method,
                    permutations=permutations,
                    pairwise=pairwise,
                    output_visualization=bg_qzv,
                    dry_run=dry_run,
                    show_stdout=show_qiime,
                )

    # --- ADONIS (PERMANOVA with formula) on Bray-Curtis by default ---
    if adonis_formula:
        dm = beta_mats["bray_curtis"]
        if dm.exists():
            ad_qzv = out_root / "adonis_bray_curtis.qzv"
            if overwrite or not ad_qzv.exists():
                qiime.diversity_adonis(
                    distance_matrix=dm,
                    metadata_file=metadata_file,
                    formula=adonis_formula,
                    permutations=permutations,
                    n_jobs=1,
                    output_visualization=ad_qzv,
                    dry_run=dry_run,
                    show_stdout=show_qiime,
                )

    LOG.info("Stats complete → %s", out_root)
    return out_root
