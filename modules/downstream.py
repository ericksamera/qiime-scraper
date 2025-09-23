# modules/downstream.py

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import zipfile
from pathlib import Path
from typing import Iterable, Optional

import shutil, time

from . import qiime_wrapper

logger = logging.getLogger("qiime_pipeline")

_SCRIPT_RE = re.compile(
    r'<script[^>]*id=["\']table-data["\'][^>]*>\s*(\{.*?\})\s*</script>',
    re.I | re.S,
)


def analysis_dir(project_dir: Path) -> Path:
    ana = project_dir / "analysis"
    ana.mkdir(parents=True, exist_ok=True)
    return ana


def _read_json(p: Path) -> dict:
    return json.loads(p.read_text())


def _stage(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    try:
        dest.symlink_to(src)
        how = "symlinked"
    except Exception:
        shutil.copy2(src, dest)
        how = "copied"
    logger.info("%s -> %s", how.capitalize(), dest)


def _pick_classifier(cls_json: Path, preferred: str) -> str:
    data = _read_json(cls_json)
    # try preferred, then fallbacks
    for k in (preferred, "median_conf", "mean_conf"):
        if k in data and data[k]:
            return data[k][0]["classifier"]  # e.g., "silva-138-99-515-806.qza"
    raise RuntimeError("No winners found in optimal_classifiers.json")


def stage_winners(
    project_dir: Path,
    *,
    preferred_metric: str = "pct_depth≥7",
    metadata_file: Optional[Path] = None,
) -> dict[str, Path]:
    """
    Link/copy winners into analysis/ with tutorial-style names:
    rep-seqs.qza, table.qza, taxonomy.qza (+ helpful .qzv summaries).
    """
    proj = project_dir.resolve()
    ana = analysis_dir(proj)

    trim_dir = proj / "work" / "optimal_trimming"
    cls_dir = proj / "work" / "optimal_classifier"

    # --- trimming winner (table + rep_seqs) ---
    t_json = trim_dir / "optimal_trimming.json"
    if not t_json.exists():
        raise FileNotFoundError(f"Missing: {t_json}")
    t = _read_json(t_json)
    f, r = int(t["trunc_len_f"]), int(t["trunc_len_r"])
    tag = f"{f}-{r}"
    rep = trim_dir / f"{tag}_output_rep_seqs.qza"
    tbl = trim_dir / f"{tag}_output_table.qza"
    tbl_qzv = trim_dir / f"{tag}_output_table.qzv"
    for p in (rep, tbl, tbl_qzv):
        if not p.exists():
            raise FileNotFoundError(f"Missing trimming output: {p}")

    # --- classifier winner (taxonomy) ---
    c_json = cls_dir / "optimal_classifiers.json"
    if not c_json.exists():
        raise FileNotFoundError(f"Missing: {c_json}")
    winner_name = _pick_classifier(c_json, preferred_metric)
    tax = cls_dir / f"{Path(winner_name).stem}_classification.qza"
    if not tax.exists():
        raise FileNotFoundError(f"Missing classification: {tax}")

    # --- stage into analysis/ ---
    rep_out = ana / "rep-seqs.qza"
    tbl_out = ana / "table.qza"
    tax_out = ana / "taxonomy.qza"
    _stage(rep, rep_out)
    _stage(tbl, tbl_out)
    _stage(tax, tax_out)
    _stage(tbl_qzv, ana / "table.qzv")  # keep the pipeline-made summary, too

    # Re-summarize table with metadata + tabulate sequences
    if metadata_file is not None:
        qiime_wrapper.feature_table_summarize(
            input_table=tbl_out,
            output=ana / "table-with-metadata.qzv",
            sample_metadata_file=metadata_file,
        )
    qiime_wrapper.feature_table_tabulate_seqs(
        input_data=rep_out,
        output=ana / "rep-seqs.qzv",
    )

    # Write a tiny winners record
    env = ana / "WINNERS.env"
    env.write_text(f"WIN_TRUNC_F={f}\nWIN_TRUNC_R={r}\nWIN_CLASSIFIER={winner_name}\n")

    return {
        "rep_seqs": rep_out,
        "table": tbl_out,
        "taxonomy": tax_out,
        "analysis_dir": ana,
    }


def _extract_freqs_from_table_qzv(qzv: Path) -> list[int]:
    if not qzv.exists():
        return []
    with zipfile.ZipFile(qzv) as z:
        try:
            html = z.read(
                next(p for p in z.namelist() if Path(p).name == "sample-frequency-detail.html")
            ).decode()
        except (KeyError, StopIteration):
            return []
    m = _SCRIPT_RE.search(html)
    if not m:
        return []
    try:
        data = json.loads(m.group(1))
        vals = list(map(int, data["Frequency"].values()))
        return sorted(vals)
    except Exception:
        return []


def choose_sampling_depth(
    ana: Path,
    *,
    retain_fraction: float = 0.90,
    min_depth: int = 1000,
) -> int:
    """
    Pick a depth that retains about *retain_fraction* of samples (default 90%).
    Falls back to *min_depth* if extraction fails.
    """
    qzv = ana / "table-with-metadata.qzv"
    if not qzv.exists():
        qzv = ana / "table.qzv"
    freqs = _extract_freqs_from_table_qzv(qzv)
    if not freqs:
        return min_depth
    k = int((1.0 - retain_fraction) * len(freqs))
    k = min(max(k, 0), len(freqs) - 1)
    return max(min_depth, freqs[k])


def build_phylogeny(ana: Path, *, show_qiime: bool = False) -> dict[str, Path]:
    rooted = ana / "rooted-tree.qza"
    if rooted.exists():
        logger.info("Rooted tree already exists: %s", rooted)
        return {"rooted_tree": rooted}

    req = ana / "rep-seqs.qza"
    if not req.exists():
        raise FileNotFoundError(
            f"Missing {req}. Run 'python main.py stage --project-dir <dir> --metadata-file <tsv>' first."
        )

    qiime_wrapper.phylogeny_align_to_tree_mafft_fasttree(
        input_sequences=req,
        output_alignment=ana / "aligned-rep-seqs.qza",
        output_masked_alignment=ana / "masked-aligned-rep-seqs.qza",
        output_tree=ana / "unrooted-tree.qza",
        output_rooted_tree=rooted,
        show_qiime=show_qiime,
    )
    return {"rooted_tree": rooted}


def _metadata_columns(metadata_file: Path) -> set[str]:
    try:
        header = metadata_file.read_text().splitlines()[0]
        return set([h.strip() for h in header.split("\t")])
    except Exception:
        return set()


def _core_outputs_ready(out_dir: Path) -> bool:
    need = [
        "unweighted_unifrac_distance_matrix.qza",
        "faith_pd_vector.qza",
        "evenness_vector.qza",
    ]
    return all((out_dir / n).exists() for n in need)

def _unique_dir(base: Path, name: str) -> Path:
    ts = time.strftime("%Y%m%d-%H%M%S")
    return base / f"{name}-{ts}"

def run_core_diversity(
    ana: Path,
    metadata_file: Path,
    *,
    sampling_depth: Optional[int] = None,
    beta_cols: Iterable[str] = ("body-site", "subject"),
    time_column: Optional[str] = None,
    include_alpha_tests: bool = True,
    include_beta_tests: bool = True,
    include_emperor: bool = True,
    auto_build_tree: bool = True,
    show_qiime: bool = False,
    if_exists: str = "skip",
) -> dict[str, Path]:

    table = ana / "table.qza"
    tree  = ana / "rooted-tree.qza"
    if not table.exists():
        raise FileNotFoundError(f"Missing {table}. Run 'python main.py stage ...' first.")
    if not tree.exists():
        if auto_build_tree:
            logger.info("No rooted-tree.qza found; building it now…")
            build_phylogeny(ana, show_qiime=show_qiime)
        else:
            raise FileNotFoundError(f"Missing {tree}. Run 'python main.py phylogeny ...' first.")

    if sampling_depth is None:
        sampling_depth = choose_sampling_depth(ana)
    logger.info("Using --p-sampling-depth=%d", sampling_depth)

    out_dir = ana / "core-metrics-phylo"

    # Handle existing output directory
    if out_dir.exists():
        if if_exists == "skip":
            if _core_outputs_ready(out_dir):
                logger.info("Reusing existing core metrics at: %s", out_dir)
            else:
                raise SystemExit(
                    "core-metrics-phylo exists but looks incomplete. "
                    "Rerun with --if-exists overwrite or --if-exists new."
                )
        elif if_exists == "overwrite":
            shutil.rmtree(out_dir)
            logger.info("Deleted existing directory: %s", out_dir)
        elif if_exists == "new":
            out_dir = _unique_dir(ana, "core-metrics-phylo")
            logger.info("Writing core metrics to new directory: %s", out_dir)
        elif if_exists == "error":
            raise SystemExit(
                f"Output directory exists: {out_dir}. Use --if-exists overwrite|skip|new."
            )
        else:
            raise ValueError("--if-exists must be one of: skip|overwrite|new|error")

    # Run core metrics only if we didn't choose 'skip' with a ready dir
    if not (out_dir.exists() and if_exists == "skip" and _core_outputs_ready(out_dir)):
        qiime_wrapper.diversity_core_metrics_phylogenetic(
            input_phylogeny=tree,
            input_table=table,
            sampling_depth=sampling_depth,
            metadata_file=metadata_file,
            output_dir=out_dir,
            show_qiime=show_qiime,
        )

    # Continue with group tests / Emperor as before … (your existing code here)
    cols_in_meta = _metadata_columns(metadata_file)
    if include_alpha_tests:
        qiime_wrapper.diversity_alpha_group_significance(
            input_alpha_vector=out_dir / "faith_pd_vector.qza",
            metadata_file=metadata_file,
            output_visualization=ana / "faith-pd-group-significance.qzv",
            show_qiime=show_qiime,
        )
        qiime_wrapper.diversity_alpha_group_significance(
            input_alpha_vector=out_dir / "evenness_vector.qza",
            metadata_file=metadata_file,
            output_visualization=ana / "evenness-group-significance.qzv",
            show_qiime=show_qiime,
        )

    if include_beta_tests:
        for col in (beta_cols or []):
            if col in cols_in_meta:
                qiime_wrapper.diversity_beta_group_significance(
                    input_distance=out_dir / "unweighted_unifrac_distance_matrix.qza",
                    metadata_file=metadata_file,
                    metadata_column=col,
                    output_visualization=ana / f"unweighted-unifrac-{col}-group-significance.qzv",
                    pairwise=True,
                    show_qiime=show_qiime,
                )
            else:
                logger.info("Skipping beta-group-significance: '%s' not in metadata.", col)

    if include_emperor and time_column and time_column in cols_in_meta:
        qiime_wrapper.emperor_plot(
            input_pcoa=out_dir / "unweighted_unifrac_pcoa_results.qza",
            metadata_file=metadata_file,
            output_visualization=ana / f"unweighted-unifrac-emperor-{time_column}.qzv",
            custom_axes=time_column,
            show_qiime=show_qiime,
        )
        qiime_wrapper.emperor_plot(
            input_pcoa=out_dir / "bray_curtis_pcoa_results.qza",
            metadata_file=metadata_file,
            output_visualization=ana / f"bray-curtis-emperor-{time_column}.qzv",
            custom_axes=time_column,
            show_qiime=show_qiime,
        )

    return {"core_dir": out_dir, "sampling_depth": Path(str(sampling_depth))}

def run_alpha_rarefaction(
    ana: Path,
    metadata_file: Path,
    *,
    max_depth: Optional[int] = None,
    sampling_depth_hint: Optional[int] = None,
) -> Path:
    if max_depth is None:
        if sampling_depth_hint is None:
            sampling_depth_hint = choose_sampling_depth(ana)
        max_depth = int(sampling_depth_hint * 2)
    out = ana / "alpha-rarefaction.qzv"
    qiime_wrapper.diversity_alpha_rarefaction(
        input_table=ana / "table.qza",
        input_phylogeny=ana / "rooted-tree.qza",
        max_depth=max_depth,
        metadata_file=metadata_file,
        output_visualization=out,
    )
    return out


def run_taxa_barplots(ana: Path, metadata_file: Path) -> Path:
    out = ana / "taxa-bar-plots.qzv"
    qiime_wrapper.taxa_barplot(
        input_table=ana / "table.qza",
        input_taxonomy=ana / "taxonomy.qza",
        metadata_file=metadata_file,
        output_visualization=out,
    )
    return out


def run_downstream(
    project_dir: Path,
    metadata_file: Path,
    *,
    preferred_metric: str = "pct_depth≥7",
    sampling_depth: Optional[int] = None,
    beta_cols: Iterable[str] = ("body-site", "subject"),
    time_column: Optional[str] = None,
    include_taxa_barplots: bool = True,
) -> dict[str, Path]:
    """Full flow: stage winners → tree → core metrics → rarefaction → (optional) taxa barplots."""
    st = stage_winners(project_dir, preferred_metric=preferred_metric, metadata_file=metadata_file)
    ana = st["analysis_dir"]
    build_phylogeny(ana)
    div = run_core_diversity(
        ana,
        metadata_file,
        sampling_depth=sampling_depth,
        beta_cols=beta_cols,
        time_column=time_column,
    )
    run_alpha_rarefaction(ana, metadata_file, sampling_depth_hint=int(str(div["sampling_depth"])))
    if include_taxa_barplots:
        run_taxa_barplots(ana, metadata_file)
    return {"analysis_dir": ana}
