# src/qs/analysis/staging.py
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from qs.qiime import commands as qiime
from qs.utils.logger import get_logger

LOG = get_logger("staging")


def _is_broken_symlink(p: Path) -> bool:
    """Return True if *p* is a symlink whose target does not exist."""
    return p.is_symlink() and not p.exists()


def _stage_artifact(src: Path, dest: Path) -> None:
    """Place *src* at *dest* using hardlink → symlink → copy, idempotently."""
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
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)


def _repair_rep_seqs(group_slug: str, group_dir: Path, rep_seqs: Path, *,
                     dry_run: bool, show_qiime: bool) -> None:
    """Recreate group-level rep-seqs.qza by copying/merging per-run outputs."""
    project_dir = group_dir.parent.parent  # .../project/groups/<slug> → project/
    runs_root = project_dir / "runs"
    per_run_rep = sorted(runs_root.glob(f"*/{group_slug}/rep-seqs.qza"))
    if not per_run_rep and group_slug == "all":
        # fallback if user ran without split-by-primer-group
        per_run_rep = sorted(runs_root.glob("*/rep-seqs.qza"))
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
        LOG.warning("Repaired rep-seqs.qza by copying from %s", per_run_rep[0])
    else:
        tmp_merged = group_dir / "_rep-seqs.repaired.qza"
        qiime.merge_seqs(per_run_rep, tmp_merged, dry_run=dry_run, show_stdout=show_qiime)
        try:
            rep_seqs.unlink()
        except Exception:
            pass
        os.replace(tmp_merged, rep_seqs)
        LOG.warning("Repaired rep-seqs.qza by merging %d per-run files", len(per_run_rep))


def _repair_table(group_slug: str, group_dir: Path, table: Path, *,
                  dry_run: bool, show_qiime: bool) -> None:
    """Recreate group-level table.qza by copying/merging per-run outputs."""
    project_dir = group_dir.parent.parent  # .../project/groups/<slug> → project/
    runs_root = project_dir / "runs"
    per_run_tbl = sorted(runs_root.glob(f"*/{group_slug}/table.qza"))
    if not per_run_tbl and group_slug == "all":
        # fallback if user ran without split-by-primer-group
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
        LOG.warning("Repaired table.qza by copying from %s", per_run_tbl[0])
    else:
        tmp_tbl = group_dir / "_table.repaired.qza"
        qiime.merge_tables(per_run_tbl, tmp_tbl, dry_run=dry_run, show_stdout=show_qiime)
        try:
            table.unlink()
        except Exception:
            pass
        os.replace(tmp_tbl, table)
        LOG.warning("Repaired table.qza by merging %d per-run files", len(per_run_tbl))
