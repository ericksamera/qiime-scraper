# src/qs/analysis/metrics.py
from __future__ import annotations

import zipfile
from pathlib import Path
from typing import Dict, Optional, Sequence

from qs.analysis.staging import (
    _is_broken_symlink,
    _stage_artifact,
    _repair_rep_seqs,
    _repair_table,
)
from qs.qiime import commands as qiime
from qs.analysis.depth import choose_sampling_depth


def stage_and_run_metrics_for_group(
    *,
    group_dir: Path,
    metadata_augmented: Path,
    retain_fraction: float = 0.90,
    min_depth: int = 1000,
    if_exists: str = "skip",  # skip|overwrite|error|new
    dry_run: bool = False,
    show_qiime: bool = False,
    make_taxa_barplot: bool = True,
) -> Path:
    """
    Ensure rep-seqs/table exist in group_dir (prefer filtered *_final.qza if present),
    optionally emit taxa-barplot if taxonomy.qza is present,
    then build phylogeny and run core metrics.
    """
    group_dir.mkdir(parents=True, exist_ok=True)
    slug = group_dir.name

    # Prefer filtered artifacts when present
    rep_src = group_dir / "rep_seqs_final.qza"
    if not rep_src.exists():
        rep_src = group_dir / "rep-seqs.qza"
    tbl_src = group_dir / "table_final.qza"
    if not tbl_src.exists():
        tbl_src = group_dir / "table.qza"

    if (not rep_src.exists()) or _is_broken_symlink(rep_src):
        _repair_rep_seqs(slug, group_dir, rep_src, dry_run=dry_run, show_qiime=show_qiime)
    if (not tbl_src.exists()) or _is_broken_symlink(tbl_src):
        _repair_table(slug, group_dir, tbl_src, dry_run=dry_run, show_qiime=show_qiime)

    # Stage to canonical names used by downstream steps
    _stage_artifact(rep_src, group_dir / "rep-seqs.qza")
    _stage_artifact(tbl_src, group_dir / "table.qza")

    # taxa barplot (fast) if taxonomy present
    tax = group_dir / "taxonomy.qza"
    barplot_qzv = group_dir / "taxa-barplot.qzv"
    if make_taxa_barplot and tax.exists() and not barplot_qzv.exists():
        qiime.taxa_barplot(
            input_table=group_dir / "table.qza",
            input_taxonomy=tax,
            metadata_file=metadata_augmented,
            output_visualization=barplot_qzv,
            dry_run=dry_run,
            show_stdout=show_qiime,
        )

    rooted = group_dir / "rooted-tree.qza"
    if not rooted.exists():
        qiime.phylogeny_align_to_tree_mafft_fasttree(
            input_sequences=group_dir / "rep-seqs.qza",
            output_alignment=group_dir / "aligned-rep-seqs.qza",
            output_masked_alignment=group_dir / "masked-aligned-rep-seqs.qza",
            output_tree=group_dir / "unrooted-tree.qza",
            output_rooted_tree=rooted,
            dry_run=dry_run,
            show_stdout=show_qiime,
        )

    table_qzv = group_dir / "table.qzv"
    if not table_qzv.exists():
        qiime.feature_table_summarize(
            input_table=group_dir / "table.qza",
            output=table_qzv,
            sample_metadata_file=metadata_augmented,
            dry_run=dry_run,
            show_stdout=show_qiime,
        )
    depth = choose_sampling_depth(table_qzv, retain_fraction=retain_fraction, min_depth=min_depth)

    core_dir = group_dir / "core-metrics-phylo"
    if core_dir.exists():
        if if_exists == "skip":
            return core_dir
        if if_exists == "overwrite":
            import shutil
            shutil.rmtree(core_dir)
        elif if_exists == "new":
            from time import strftime
            core_dir = group_dir / f"core-metrics-phylo-{strftime('%Y%m%d-%H%M%S')}"
        elif if_exists == "error":
            raise FileExistsError(core_dir)

    qiime.diversity_core_metrics_phylogenetic(
        input_phylogeny=rooted,
        input_table=group_dir / "table.qza",
        sampling_depth=depth,
        metadata_file=metadata_augmented,
        output_dir=core_dir,
        dry_run=dry_run,
        show_stdout=show_qiime,
    )
    return core_dir


def write_alpha_metrics_tsv(core_metrics_dir: Path, out_tsv: Optional[Path] = None) -> Path:
    alpha_qzas = {
        "observed_features": core_metrics_dir / "observed_features_vector.qza",
        "evenness":          core_metrics_dir / "evenness_vector.qza",
        "faith_pd":          core_metrics_dir / "faith_pd_vector.qza",
        "shannon":           core_metrics_dir / "shannon_vector.qza",
    }

    def _read_alpha_vec(qza: Path) -> Dict[str, float]:
        if not qza.exists():
            return {}
        try:
            with zipfile.ZipFile(qza) as z:
                tsv_name = next((p for p in z.namelist() if p.endswith("alpha-diversity.tsv")), None)
                if not tsv_name:
                    tsv_name = next((p for p in z.namelist() if p.endswith(".tsv")), None)
                if not tsv_name:
                    return {}
                lines = z.read(tsv_name).decode("utf-8").splitlines()
        except Exception:
            return {}
        if not lines:
            return {}
        header = [h.strip().lower() for h in lines[0].split("\t")]
        try:
            sid_idx = next(i for i, h in enumerate(header) if h in {"#sampleid", "sample-id", "sample id", "sampleid"})
        except StopIteration:
            sid_idx = 0
        val_idx = next((i for i, h in enumerate(header) if h in {"value", "alpha-diversity", "alpha_diversity"}), 1 if len(header) > 1 else 0)
        out: Dict[str, float] = {}
        for ln in lines[1:]:
            cols = ln.split("\t")
            if len(cols) <= max(sid_idx, val_idx):
                continue
            sid = cols[sid_idx].lstrip("#").strip()
            try:
                out[sid] = float(cols[val_idx])
            except Exception:
                continue
        return out

    columns: Dict[str, Dict[str, float]] = {k: _read_alpha_vec(p) for k, p in alpha_qzas.items()}
    if not any(columns.values()):
        raise FileNotFoundError(f"No alpha diversity vectors found under {core_metrics_dir}")

    all_sids = sorted({sid for col in columns.values() for sid in col.keys()})
    out_tsv = out_tsv or (core_metrics_dir / "alpha-metrics.tsv")
    with out_tsv.open("w", encoding="utf-8") as fh:
        fh.write("sample-id\tobserved_features\tevenness\tfaith_pd\tshannon\n")
        for sid in all_sids:
            fh.write(
                f"{sid}\t"
                f"{columns.get('observed_features', {}).get(sid, '')}\t"
                f"{columns.get('evenness', {}).get(sid, '')}\t"
                f"{columns.get('faith_pd', {}).get(sid, '')}\t"
                f"{columns.get('shannon', {}).get(sid, '')}\n"
            )
    return out_tsv
