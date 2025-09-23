# modules/downstream.py

from __future__ import annotations

import json
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, TypedDict

from . import qiime_wrapper

logger = logging.getLogger("qiime_pipeline")

_SCRIPT_RE = re.compile(
    r'<script[^>]*id=["\']table-data["\'][^>]*>\s*(\{.*?\})\s*</script>',
    re.I | re.S,
)

# ------------------- basics -------------------

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

# ------------------- provenance / RUN.json -------------------

def _now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")

def _try_qiime_info(dry_run: bool) -> str | None:
    if dry_run:
        return None
    try:
        out = subprocess.run(["qiime", "info"], check=True, capture_output=True, text=True)
        return out.stdout.strip()
    except Exception:
        return None

def write_run_provenance(ana: Path, command: str, payload: Mapping[str, object], *, dry_run: bool) -> None:
    p = ana / "RUN.json"
    rec = {
        "command": command,
        "started_at": _now_iso(),
        "args": payload,
        "env": {
            "python": sys.version,
            "platform": platform.platform(),
            "cwd": str(Path.cwd()),
            "qiime_info": _try_qiime_info(dry_run),
        },
    }
    try:
        winners = (ana / "WINNERS.env").read_text().splitlines()
        rec["winners_env"] = dict(kv.split("=", 1) for kv in winners if "=" in kv)
    except Exception:
        rec["winners_env"] = {}
    p.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")

# ------------------- helpers: metrics & depths -------------------

def _normalize_metric_key(k: str) -> str:
    """
    Normalize common variants to the JSON key used in reports.
    Accepts 'pct_depth>=7' as a synonym for 'pct_depth≥7'.
    """
    if not k:
        return k
    s = k.strip().lower()
    s = s.replace(" ", "")
    s = s.replace(">=", "≥")
    return s

def _pick_classifier(cls_json: Path, preferred: str) -> str:
    data = _read_json(cls_json)
    pref = _normalize_metric_key(preferred)
    # try preferred, then fallbacks
    for k in (pref, "pct_depth≥7", "median_conf", "mean_conf"):
        if k in data and data[k]:
            return data[k][0]["classifier"]  # e.g., "silva-138-99-515-806.qza"
    raise RuntimeError("No winners found in optimal_classifiers.json")

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

def choose_sampling_depth_from_qzv(qzv: Path, *, retain_fraction: float = 0.90, min_depth: int = 1000) -> int:
    freqs = _extract_freqs_from_table_qzv(qzv)
    if not freqs:
        return min_depth
    k = int((1.0 - retain_fraction) * len(freqs))
    k = min(max(k, 0), len(freqs) - 1)
    return max(min_depth, freqs[k])

# ------------------- core downstream steps -------------------

def stage_winners(
    project_dir: Path,
    *,
    preferred_metric: str = "pct_depth≥7",
    metadata_file: Optional[Path] = None,
    show_qiime: bool = False,
    dry_run: bool = False,
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
            dry_run=dry_run,
            show_qiime=show_qiime,
        )
    qiime_wrapper.feature_table_tabulate_seqs(
        input_data=rep_out,
        output=ana / "rep-seqs.qzv",
        dry_run=dry_run,
        show_qiime=show_qiime,
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

def build_phylogeny(ana: Path, *, show_qiime: bool = False, dry_run: bool = False) -> dict[str, Path]:
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
        dry_run=dry_run,
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

class CoreDiversityResult(TypedDict):
    core_dir: Path
    sampling_depth: int

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
    dry_run: bool = False,
) -> CoreDiversityResult:

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
            dry_run=dry_run,
        )

    # Continue with group tests / Emperor
    cols_in_meta = _metadata_columns(metadata_file)
    if include_alpha_tests:
        qiime_wrapper.diversity_alpha_group_significance(
            input_alpha_vector=out_dir / "faith_pd_vector.qza",
            metadata_file=metadata_file,
            output_visualization=ana / "faith-pd-group-significance.qzv",
            show_qiime=show_qiime,
            dry_run=dry_run,
        )
        qiime_wrapper.diversity_alpha_group_significance(
            input_alpha_vector=out_dir / "evenness_vector.qza",
            metadata_file=metadata_file,
            output_visualization=ana / "evenness-group-significance.qzv",
            show_qiime=show_qiime,
            dry_run=dry_run,
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
                    dry_run=dry_run,
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
            dry_run=dry_run,
        )
        qiime_wrapper.emperor_plot(
            input_pcoa=out_dir / "bray_curtis_pcoa_results.qza",
            metadata_file=metadata_file,
            output_visualization=ana / f"bray-curtis-emperor-{time_column}.qzv",
            custom_axes=time_column,
            show_qiime=show_qiime,
            dry_run=dry_run,
        )

    return {"core_dir": out_dir, "sampling_depth": int(sampling_depth)}

def run_alpha_rarefaction(
    ana: Path,
    metadata_file: Path,
    *,
    max_depth: Optional[int] = None,
    sampling_depth_hint: Optional[int] = None,
    show_qiime: bool = False,
    dry_run: bool = False,
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
        show_qiime=show_qiime,
        dry_run=dry_run
    )
    return out

def run_taxa_barplots(ana: Path, metadata_file: Path, *, show_qiime: bool = False, dry_run: bool = False) -> Path:
    out = ana / "taxa-bar-plots.qzv"
    qiime_wrapper.taxa_barplot(
        input_table=ana / "table.qza",
        input_taxonomy=ana / "taxonomy.qza",
        metadata_file=metadata_file,
        output_visualization=out,
        dry_run=dry_run,
        show_qiime=show_qiime,
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
    show_qiime: bool = False,
    dry_run: bool = False,
) -> dict[str, Path]:
    """Full flow: stage winners → tree → core metrics → rarefaction → (optional) taxa barplots."""
    st = stage_winners(project_dir, preferred_metric=preferred_metric, metadata_file=metadata_file,
                       show_qiime=show_qiime, dry_run=dry_run)
    ana = st["analysis_dir"]
    build_phylogeny(ana, show_qiime=show_qiime, dry_run=dry_run)
    div = run_core_diversity(
        ana,
        metadata_file,
        sampling_depth=sampling_depth,
        beta_cols=beta_cols,
        time_column=time_column,
        show_qiime=show_qiime,
        dry_run=dry_run,
    )
    run_alpha_rarefaction(ana, metadata_file, sampling_depth_hint=div["sampling_depth"],
                          show_qiime=show_qiime, dry_run=dry_run)
    if include_taxa_barplots:
        run_taxa_barplots(ana, metadata_file, show_qiime=show_qiime, dry_run=dry_run)
    return {"analysis_dir": ana}

# ------------------- NEW: Diversity sweep -------------------

_SLUG_RE = re.compile(r"[^A-Za-z0-9._-]+")

def _slug(s: str) -> str:
    return _SLUG_RE.sub("-", (s or "").strip()).strip("-._")

def _parse_metadata_table(metadata_file: Path) -> tuple[List[str], Dict[str, str], List[Dict[str, str]]]:
    """
    Returns (headers, types_map, rows)
    - headers: list of column names
    - types_map: column -> 'categorical' | 'numeric' | ''
    - rows: list of dicts for each sample row
    """
    lines = [ln.rstrip("\n") for ln in metadata_file.read_text(encoding="utf-8").splitlines() if ln.strip() != ""]
    if not lines:
        raise SystemExit(f"Empty metadata file: {metadata_file}")
    headers = [h.strip() for h in lines[0].split("\t")]
    types_map: Dict[str, str] = {h: "" for h in headers}
    start_row = 1
    if len(lines) > 1 and lines[1].startswith("#q2:types"):
        type_vals = [v.strip() for v in lines[1].split("\t")]
        for i, h in enumerate(headers):
            if i < len(type_vals):
                types_map[h] = type_vals[i].lower()
        start_row = 2
    rows: List[Dict[str, str]] = []
    for ln in lines[start_row:]:
        parts = ln.split("\t")
        row = {}
        for i, h in enumerate(headers):
            row[h] = parts[i] if i < len(parts) else ""
        rows.append(row)
    return headers, types_map, rows

def _categorical_columns(headers: Sequence[str], types_map: Mapping[str, str]) -> List[str]:
    cats = []
    for h in headers:
        if h.lower() in {"sample-id", "sampleid", "#sampleid"}:
            cats.append(h)  # keep for reference but we won't loop on it
            continue
        t = (types_map.get(h) or "").lower()
        if t == "categorical" or (t == "" and h.lower() not in {"date", "time"}):
            cats.append(h)
    return cats

def _value_counts(rows: Sequence[Mapping[str, str]], column: str) -> Counter:
    c = Counter()
    for r in rows:
        v = (r.get(column) or "").strip()
        if v != "":
            c[v] += 1
    return c

def _require_columns(cols_in_meta: Set[str], needed: Iterable[str], *, context: str) -> None:
    missing = [c for c in (needed or []) if c and c not in cols_in_meta]
    if not missing:
        return
    advice = f"Available columns: {', '.join(sorted(cols_in_meta))}"
    raise SystemExit(f"[metadata check] {context}: missing columns: {', '.join(missing)}.\n{advice}")

def _sql_quote(val: str) -> str:
    # quote single quotes for SQLite-style WHERE
    return "'" + val.replace("'", "''") + "'"

def run_diversity_sweep(
    project_dir: Path,
    metadata_file: Path,
    *,
    by: Iterable[str] | None,
    within: Iterable[str] | None,
    min_samples: int = 5,
    retain_fraction: float = 0.90,
    beta_cols: Iterable[str] | None = None,
    time_column: Optional[str] = None,
    by_only: bool = False,
    if_exists: str = "skip",
    show_qiime: bool = True,
    dry_run: bool = False,
) -> None:
    """
    Run diversity across metadata-defined subsets.

    Modes
    -----
    • Nested (default): outer loop over `by` columns then inner loop over `within` columns.
      For each (outer=value, inner=value2) subset, compute a per-subset sampling depth and run core metrics.
      Group significance tests use [inner_col] + beta_cols.

    • One level (`by_only=True`): loop *only* over `by` columns.
      For each (by=value) subset, compute depth and run core metrics.
      Group significance tests use only `beta_cols` that have ≥2 levels within the subset
      (the outer column itself is skipped because it is constant in the subset).

    Notes
    -----
    - If `by` (or `within`) is None, we use all categorical columns.
    - Subsets with < min_samples (by metadata row counts) are skipped.
    - Per-subset depth is chosen from that subset's table.qzv using `retain_fraction`.
    """
    proj = project_dir.resolve()
    ana = analysis_dir(proj)

    # Validate prerequisites
    table = ana / "table.qza"
    if not table.exists():
        raise SystemExit(
            f"Missing {table}. Run 'python main.py stage --project-dir <dir> --metadata-file <tsv>' first."
        )

    # Ensure rooted tree available
    build_phylogeny(ana, show_qiime=show_qiime, dry_run=dry_run)

    # Parse metadata and plan loops
    headers, types_map, rows = _parse_metadata_table(metadata_file)
    cols_in_meta = set(headers)
    cats = [
        c for c in _categorical_columns(headers, types_map)
        if c.lower() not in {"sample-id", "sampleid", "#sampleid"}
    ]

    by_cols = list(by) if by else list(cats)
    within_cols = list(within) if within else list(cats)

    # Friendly validation
    _require_columns(cols_in_meta, by_cols, context="--by")
    if not by_only:
        _require_columns(cols_in_meta, within_cols, context="--within")
    if time_column:
        _require_columns(cols_in_meta, [time_column], context="--time-column")
    if beta_cols:
        _require_columns(cols_in_meta, beta_cols, context="--beta-cols")

    # Run-level provenance
    write_run_provenance(
        ana,
        "diversity-sweep",
        {
            "project_dir": str(proj),
            "metadata_file": str(metadata_file),
            "by": by_cols,
            "within": None if by_only else within_cols,
            "by_only": by_only,  # <-- NEW
            "min_samples": min_samples,
            "retain_fraction": retain_fraction,
            "beta_cols": list(beta_cols or []),
            "time_column": time_column,
            "if_exists": if_exists,
            "show_qiime": show_qiime,
            "dry_run": dry_run,
        },
        dry_run=dry_run,
    )

    # Precompute counts per column/value for outer loop(s)
    value_counts: Dict[str, Counter] = {c: _value_counts(rows, c) for c in cats}

    subsets_root = ana / "subsets"
    subsets_root.mkdir(parents=True, exist_ok=True)
    index_records: List[Dict[str, object]] = []

    # Helper: compute distinct levels of a column inside an outer subset
    def _distinct_levels_in_subset(filter_col: str, filter_val: str, test_col: str) -> set[str]:
        levels: set[str] = set()
        for r in rows:
            if (r.get(filter_col) or "").strip() == filter_val:
                v = (r.get(test_col) or "").strip()
                if v != "":
                    levels.add(v)
        return levels

    # --------------------- MODE: BY-ONLY (one level) ---------------------
    if by_only:
        for outer_col in by_cols:
            for outer_val, n_outer in sorted(value_counts.get(outer_col, Counter()).items()):
                if n_outer < min_samples:
                    logger.info(
                        "[skip] %s=%s has only %d samples (< %d)",
                        outer_col, outer_val, n_outer, min_samples
                    )
                    continue

                outer_where = f'"{outer_col}" = {_sql_quote(outer_val)}'
                sub_dir = subsets_root / _slug(f"{outer_col}") / _slug(f"{outer_val}")
                sub_dir.mkdir(parents=True, exist_ok=True)

                # 1) Filter table
                sub_table = sub_dir / "table.qza"
                if sub_table.exists() and if_exists == "skip":
                    logger.info("Reusing filtered table: %s", sub_table)
                else:
                    if sub_table.exists() and if_exists == "overwrite":
                        sub_table.unlink()
                    elif sub_table.exists() and if_exists == "new":
                        sub_dir = _unique_dir(sub_dir.parent, "subset")
                        sub_dir.mkdir(parents=True, exist_ok=True)
                        sub_table = sub_dir / "table.qza"
                    qiime_wrapper.feature_table_filter_samples(
                        input_table=table,
                        metadata_file=metadata_file,
                        where=outer_where,
                        output_table=sub_table,
                        dry_run=dry_run,
                        show_qiime=show_qiime,
                    )

                # 2) Summarize subset table + choose depth
                sub_qzv = sub_dir / "table.qzv"
                qiime_wrapper.feature_table_summarize(
                    input_table=sub_table,
                    output=sub_qzv,
                    sample_metadata_file=metadata_file,
                    dry_run=dry_run,
                    show_qiime=show_qiime,
                )
                depth = choose_sampling_depth_from_qzv(sub_qzv, retain_fraction=retain_fraction)

                # 3) Core metrics
                sub_core = sub_dir / "core-metrics-phylo"
                root_tree = ana / "rooted-tree.qza"
                if not root_tree.exists():
                    build_phylogeny(ana, show_qiime=show_qiime, dry_run=dry_run)

                if sub_core.exists():
                    if if_exists == "skip":
                        logger.info("Reusing existing core metrics at: %s", sub_core)
                    elif if_exists == "overwrite":
                        shutil.rmtree(sub_core)
                    elif if_exists == "new":
                        sub_core = _unique_dir(sub_dir, "core-metrics-phylo")
                        sub_core.mkdir(parents=True, exist_ok=True)
                    elif if_exists == "error":
                        raise SystemExit(f"Subset output exists: {sub_core}")
                    else:
                        raise ValueError("--if-exists must be one of: skip|overwrite|new|error")

                if not (sub_core.exists() and if_exists == "skip"):
                    qiime_wrapper.diversity_core_metrics_phylogenetic(
                        input_phylogeny=root_tree,
                        input_table=sub_table,
                        sampling_depth=depth,
                        metadata_file=metadata_file,
                        output_dir=sub_core,
                        dry_run=dry_run,
                        show_qiime=show_qiime,
                    )

                # 4) Group tests: only beta_cols that have ≥2 levels inside this subset
                executed_tests: List[str] = []
                for col in (beta_cols or []):
                    if col not in cols_in_meta:
                        continue
                    if col == outer_col:
                        # Outer column is constant by construction; skip testing it here.
                        logger.info(
                            "Skipping beta-group-significance on '%s' (constant within %s=%s).",
                            col, outer_col, outer_val
                        )
                        continue
                    levels = _distinct_levels_in_subset(outer_col, outer_val, col)
                    if len(levels) < 2:
                        logger.info(
                            "Skipping beta-group-significance: '%s' has < 2 levels within %s=%s.",
                            col, outer_col, outer_val
                        )
                        continue
                    qiime_wrapper.diversity_beta_group_significance(
                        input_distance=sub_core / "unweighted_unifrac_distance_matrix.qza",
                        metadata_file=metadata_file,
                        metadata_column=col,
                        output_visualization=sub_dir / f"unweighted-unifrac-{_slug(col)}-group-significance.qzv",
                        pairwise=True,
                        dry_run=dry_run,
                        show_qiime=show_qiime,
                    )
                    executed_tests.append(col)

                # 5) Emperor (optional)
                if time_column and time_column in cols_in_meta:
                    qiime_wrapper.emperor_plot(
                        input_pcoa=sub_core / "unweighted_unifrac_pcoa_results.qza",
                        metadata_file=metadata_file,
                        output_visualization=sub_dir / f"unweighted-unifrac-emperor-{_slug(time_column)}.qzv",
                        custom_axes=time_column,
                        dry_run=dry_run,
                        show_qiime=show_qiime,
                    )

                # 6) Provenance
                subset_json = {
                    "mode": "by-only",
                    "outer": {"column": outer_col, "value": outer_val, "n_samples": n_outer},
                    "inner": None,
                    "where": outer_where,
                    "sampling_depth": depth,
                    "paths": {"dir": str(sub_dir), "table": str(sub_table), "core": str(sub_core)},
                    "retain_fraction": retain_fraction,
                    "min_samples": min_samples,
                    "time_column": time_column,
                    "beta_cols": executed_tests,
                    "timestamp": _now_iso(),
                }
                (sub_dir / "subset.json").write_text(json.dumps(subset_json, indent=2) + "\n", encoding="utf-8")
                index_records.append(subset_json)

        (subsets_root / "index.json").write_text(json.dumps(index_records, indent=2) + "\n", encoding="utf-8")
        logger.info("Diversity sweep (by-only) completed across %d subsets.", len(index_records))
        return

    # --------------------- MODE: NESTED (existing behavior) ---------------------
    for outer_col in by_cols:
        for outer_val, n_outer in sorted(value_counts.get(outer_col, Counter()).items()):
            if n_outer < min_samples:
                logger.info("[skip] %s=%s has only %d samples (< %d)", outer_col, outer_val, n_outer, min_samples)
                continue

            outer_where = f'"{outer_col}" = {_sql_quote(outer_val)}'
            outer_dir = subsets_root / _slug(f"{outer_col}") / _slug(f"{outer_val}")

            for inner_col in within_cols:
                if inner_col == outer_col:
                    # Avoid redundant (outer==inner) second loop; user can request explicitly if desired
                    continue

                # Compute inner counts restricted to the outer subset
                inner_counter: Counter = Counter()
                for r in rows:
                    if (r.get(outer_col) or "").strip() == outer_val:
                        iv = (r.get(inner_col) or "").strip()
                        if iv != "":
                            inner_counter[iv] += 1

                for inner_val, n_inner in sorted(inner_counter.items()):
                    if n_inner < min_samples:
                        logger.info(
                            "[skip] %s=%s AND %s=%s → %d samples (< %d)",
                            outer_col, outer_val, inner_col, inner_val, n_inner, min_samples
                        )
                        continue

                    where = f'{outer_where} AND "{inner_col}" = {_sql_quote(inner_val)}'
                    sub_dir = outer_dir / _slug(f"{inner_col}") / _slug(f"{inner_val}")
                    sub_dir.mkdir(parents=True, exist_ok=True)

                    # 1) Filter table
                    sub_table = sub_dir / "table.qza"
                    if sub_table.exists() and if_exists == "skip":
                        logger.info("Reusing filtered table: %s", sub_table)
                    else:
                        if sub_table.exists() and if_exists == "overwrite":
                            sub_table.unlink()
                        elif sub_table.exists() and if_exists == "new":
                            sub_dir = _unique_dir(outer_dir / _slug(f"{inner_col}") / _slug(f"{inner_val}"), "subset")
                            sub_dir.mkdir(parents=True, exist_ok=True)
                            sub_table = sub_dir / "table.qza"
                        qiime_wrapper.feature_table_filter_samples(
                            input_table=table,
                            metadata_file=metadata_file,
                            where=where,
                            output_table=sub_table,
                            dry_run=dry_run,
                            show_qiime=show_qiime,
                        )

                    # 2) Summarize + choose depth
                    sub_qzv = sub_dir / "table.qzv"
                    qiime_wrapper.feature_table_summarize(
                        input_table=sub_table,
                        output=sub_qzv,
                        sample_metadata_file=metadata_file,
                        dry_run=dry_run,
                        show_qiime=show_qiime,
                    )
                    depth = choose_sampling_depth_from_qzv(sub_qzv, retain_fraction=retain_fraction)

                    # 3) Core metrics
                    sub_core = sub_dir / "core-metrics-phylo"
                    root_tree = ana / "rooted-tree.qza"
                    if not root_tree.exists():
                        build_phylogeny(ana, show_qiime=show_qiime, dry_run=dry_run)

                    if sub_core.exists():
                        if if_exists == "skip":
                            logger.info("Reusing existing core metrics at: %s", sub_core)
                        elif if_exists == "overwrite":
                            shutil.rmtree(sub_core)
                        elif if_exists == "new":
                            sub_core = _unique_dir(sub_dir, "core-metrics-phylo")
                            sub_core.mkdir(parents=True, exist_ok=True)
                        elif if_exists == "error":
                            raise SystemExit(f"Subset output exists: {sub_core}")
                        else:
                            raise ValueError("--if-exists must be one of: skip|overwrite|new|error")

                    if not (sub_core.exists() and if_exists == "skip"):
                        qiime_wrapper.diversity_core_metrics_phylogenetic(
                            input_phylogeny=root_tree,
                            input_table=sub_table,
                            sampling_depth=depth,
                            metadata_file=metadata_file,
                            output_dir=sub_core,
                            dry_run=dry_run,
                            show_qiime=show_qiime,
                        )

                    # 4) Group tests: inner column + any extra beta columns
                    test_cols = [inner_col] + list(beta_cols or [])
                    for col in test_cols:
                        if col in cols_in_meta:
                            qiime_wrapper.diversity_beta_group_significance(
                                input_distance=sub_core / "unweighted_unifrac_distance_matrix.qza",
                                metadata_file=metadata_file,
                                metadata_column=col,
                                output_visualization=sub_dir / f"unweighted-unifrac-{_slug(col)}-group-significance.qzv",
                                pairwise=True,
                                dry_run=dry_run,
                                show_qiime=show_qiime,
                            )

                    # 5) Emperor (optional)
                    if time_column and time_column in cols_in_meta:
                        qiime_wrapper.emperor_plot(
                            input_pcoa=sub_core / "unweighted_unifrac_pcoa_results.qza",
                            metadata_file=metadata_file,
                            output_visualization=sub_dir / f"unweighted-unifrac-emperor-{_slug(time_column)}.qzv",
                            custom_axes=time_column,
                            dry_run=dry_run,
                            show_qiime=show_qiime,
                        )

                    # 6) Per-subset provenance
                    subset_json = {
                        "mode": "nested",
                        "outer": {"column": outer_col, "value": outer_val, "n_samples": n_outer},
                        "inner": {"column": inner_col, "value": inner_val, "n_samples": n_inner},
                        "where": where,
                        "sampling_depth": depth,
                        "paths": {"dir": str(sub_dir), "table": str(sub_table), "core": str(sub_core)},
                        "retain_fraction": retain_fraction,
                        "min_samples": min_samples,
                        "time_column": time_column,
                        "beta_cols": test_cols,
                        "timestamp": _now_iso(),
                    }
                    (sub_dir / "subset.json").write_text(json.dumps(subset_json, indent=2) + "\n", encoding="utf-8")
                    index_records.append(subset_json)

    # Write an index of all subset runs
    (subsets_root / "index.json").write_text(json.dumps(index_records, indent=2) + "\n", encoding="utf-8")
    logger.info("Diversity sweep (nested) completed across %d subsets.", len(index_records))
