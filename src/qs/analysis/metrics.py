# src/qs/analysis/metrics.py
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Sequence

from qs.qiime import commands as qiime
from qs.analysis.depth import choose_sampling_depth


def _is_broken_symlink(p: Path) -> bool:
    return p.is_symlink() and not p.exists()


def _stage_artifact(src: Path, dest: Path) -> None:
    # no-ops if same path or already the same target
    if str(src) == str(dest):
        return
    if dest.exists() or dest.is_symlink():
        try:
            if dest.resolve() == src.resolve():
                return
        except Exception:
            pass
        try:
            dest.unlink()
        except FileNotFoundError:
            pass
    # prefer hardlink → symlink → copy
    try:
        os.link(src, dest)
        return
    except Exception:
        try:
            dest.symlink_to(src)
            return
        except Exception:
            import shutil
            shutil.copy2(src, dest)


def _repair_rep_seqs(group_slug: str, group_dir: Path, rep_seqs: Path,
                     *, dry_run: bool, show_qiime: bool) -> None:
    project_dir = group_dir.parent.parent  # .../project/groups/<slug> → project/
    runs_root = project_dir / "runs"
    per_run_rep = sorted(runs_root.glob(f"*/{group_slug}/rep-seqs.qza"))
    if not per_run_rep:
        raise FileNotFoundError(
            f"rep-seqs.qza missing and no per-run rep-seqs found under {runs_root}/ */{group_slug}/"
        )
    if len(per_run_rep) == 1:
        import shutil
        rep_seqs.parent.mkdir(parents=True, exist_ok=True)
        try:
            rep_seqs.unlink()
        except Exception:
            pass
        shutil.copy2(per_run_rep[0], rep_seqs)
    else:
        tmp_merged = group_dir / "_rep-seqs.repaired.qza"
        qiime.merge_seqs(per_run_rep, tmp_merged, dry_run=dry_run, show_stdout=show_qiime)
        try:
            rep_seqs.unlink()
        except Exception:
            pass
        os.replace(tmp_merged, rep_seqs)


def _repair_table(group_slug: str, group_dir: Path, table: Path,
                  *, dry_run: bool, show_qiime: bool) -> None:
    project_dir = group_dir.parent.parent
    runs_root = project_dir / "runs"
    per_run_tbl = sorted(runs_root.glob(f"*/{group_slug}/table.qza"))
    if not per_run_tbl and group_slug == "all":
        per_run_tbl = sorted(runs_root.glob("*/table.qza"))
    if not per_run_tbl:
        raise FileNotFoundError(
            f"table.qza missing and no per-run tables found under {runs_root}/ */{group_slug}/"
        )
    if len(per_run_tbl) == 1:
        import shutil
        table.parent.mkdir(parents=True, exist_ok=True)
        try:
            table.unlink()
        except Exception:
            pass
        shutil.copy2(per_run_tbl[0], table)
    else:
        tmp_tbl = group_dir / "_table.repaired.qza"
        qiime.merge_tables(per_run_tbl, tmp_tbl, dry_run=dry_run, show_stdout=show_qiime)
        try:
            table.unlink()
        except Exception:
            pass
        os.replace(tmp_tbl, table)


def stage_and_run_metrics_for_group(
    *,
    group_dir: Path,
    metadata_augmented: Path,
    retain_fraction: float = 0.90,
    min_depth: int = 1000,
    if_exists: str = "skip",  # skip|overwrite|error|new
    dry_run: bool = False,
    show_qiime: bool = False,
) -> Path:
    """
    Ensure rep-seqs.qza & table.qza exist in group_dir (repair if needed),
    then build phylogeny and run core metrics.

    Returns path to the core-metrics-phylo directory (or the new one if 'new').
    """
    group_dir.mkdir(parents=True, exist_ok=True)
    slug = group_dir.name

    rep_src = group_dir / "rep-seqs.qza"
    tbl_src = group_dir / "table.qza"

    if (not rep_src.exists()) or _is_broken_symlink(rep_src):
        _repair_rep_seqs(slug, group_dir, rep_src, dry_run=dry_run, show_qiime=show_qiime)
    if (not tbl_src.exists()) or _is_broken_symlink(tbl_src):
        _repair_table(slug, group_dir, tbl_src, dry_run=dry_run, show_qiime=show_qiime)

    # re-stage to local filenames (no-ops if already correct)
    _stage_artifact(rep_src, group_dir / "rep-seqs.qza")
    _stage_artifact(tbl_src, group_dir / "table.qza")

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
    # if_exists policy
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
