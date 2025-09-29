# src/qs/analysis/stats.py
from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from qs.utils.logger import get_logger
from qs.qiime import commands as qiime

LOG = get_logger("stats")


# ------------------------------- helpers --------------------------------------

def _read_metadata_simple(path: Path) -> tuple[List[str], List[Dict[str, str]], Dict[str, str]]:
    """
    Minimal metadata reader that preserves '#q2:types' when present.
    Returns (header, rows, types_map).
    - rows: list of dicts for real samples (skips the '#q2:types' row)
    - types_map: column -> 'categorical'/'numeric'/other when types row is present
    """
    header: List[str] = []
    rows: List[Dict[str, str]] = []
    types_map: Dict[str, str] = {}

    with path.open("r", encoding="utf-8", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        for i, raw in enumerate(rdr):
            if i == 0:
                header = [c.strip() for c in raw]
                continue
            # if there's a types row, record it and continue
            if raw and str(raw[0]).strip().lower().startswith("#q2:types"):
                for col, val in zip(header, raw):
                    v = (val or "").strip().lower()
                    if v:
                        types_map[col] = v
                continue
            if not raw:
                continue
            d = {col: (val or "").strip() for col, val in zip(header, raw)}
            # skip empty rows or commented rows
            sid = (d.get("#SampleID") or d.get("sample-id") or "").strip()
            if not sid or sid.startswith("#"):
                continue
            rows.append(d)
    return header, rows, types_map


def _unique_non_empty(values: Sequence[str]) -> Set[str]:
    return {v for v in values if v not in {"", "NA", "NaN", "nan", "None"}}


def _col_values(rows: List[Dict[str, str]], col: str, sample_ids: Set[str]) -> List[str]:
    out: List[str] = []
    for r in rows:
        sid = (r.get("#SampleID") or r.get("sample-id") or "").lstrip("#")
        if sid in sample_ids:
            out.append(r.get(col, ""))
    return out


def _is_numeric(values: Sequence[str]) -> bool:
    any_num = False
    for v in values:
        if v == "" or v in {"NA", "NaN", "nan", "None"}:
            continue
        try:
            float(v)
            any_num = True
        except Exception:
            return False
    return any_num


def _eligible_columns(
    header: List[str],
    rows: List[Dict[str, str]],
    types_map: Dict[str, str],
    sample_ids: Set[str],
) -> tuple[List[str], List[str]]:
    """
    Return (categorical_cols, numeric_cols) that are valid **for the given samples**.
    - categorical: at least 2 unique values and either types row marks it categorical,
      or it's non-numeric by inspection.
    - numeric: at least 2 unique numeric values.
    Excludes '#SampleID'.
    """
    cats: List[str] = []
    nums: List[str] = []
    for col in header:
        if col.strip() in {"#SampleID", "sample-id", "#sampleid"}:
            continue
        vals = _col_values(rows, col, sample_ids)
        u = _unique_non_empty(vals)
        if len(u) < 2:
            # not enough variation -> skip for both cat/numeric uses
            continue
        ty = types_map.get(col, "")
        if ty.startswith("categorical"):
            cats.append(col); continue
        if ty in {"numeric", "number"}:
            if _is_numeric(vals):
                nums.append(col)
            continue
        # no explicit type: infer
        if _is_numeric(vals):
            nums.append(col)
        else:
            cats.append(col)
    return cats, nums


def _alpha_vec_sample_ids(qza: Path) -> Set[str]:
    """
    Read alpha-diversity vector.qza and return the set of sample IDs present.
    """
    sids: Set[str] = set()
    if not qza.exists():
        return sids
    try:
        with zipfile.ZipFile(qza) as z:
            # typical name: data/alpha-diversity.tsv (or any *.tsv)
            name = next((p for p in z.namelist() if p.endswith("alpha-diversity.tsv")), None)
            if not name:
                name = next((p for p in z.namelist() if p.endswith(".tsv")), None)
            if not name:
                return sids
            lines = z.read(name).decode("utf-8").splitlines()
            if not lines:
                return sids
            for ln in lines[1:]:
                sid = (ln.split("\t", 1)[0] if "\t" in ln else ln).lstrip("#").strip()
                if sid:
                    sids.add(sid)
    except Exception:
        return set()
    return sids


def _dm_sample_ids(qza: Path) -> Set[str]:
    """
    Read distance-matrix.qza to get sample IDs (header row/col labels).
    """
    sids: Set[str] = set()
    if not qza.exists():
        return sids
    try:
        with zipfile.ZipFile(qza) as z:
            name = next((p for p in z.namelist() if p.endswith("distance-matrix.tsv")), None)
            if not name:
                name = next((p for p in z.namelist() if p.endswith(".tsv")), None)
            if not name:
                return sids
            lines = z.read(name).decode("utf-8").splitlines()
            if not lines:
                return sids
            hdr = [c.strip() for c in lines[0].split("\t")]
            # first cell blank, then sample IDs
            for c in hdr[1:]:
                if c:
                    sids.add(c)
    except Exception:
        return set()
    return sids


def _safe_call(cmd, **kwargs) -> bool:
    """
    Call a qiime.* wrapper; if it throws, log and return False rather than crash.
    """
    try:
        cmd(**kwargs)
        return True
    except Exception as e:
        LOG.warning("Skip stat step (%s): %s", getattr(cmd, "__name__", "qiime"), e)
        return False


# ------------------------------- main -----------------------------------------

def run_diversity_stats_for_group(
    *,
    group_dir: Path,
    metadata_file: Path,
    beta_group_cols: Sequence[str] | None = None,
    beta_method: str = "permanova",
    pairwise: bool = False,
    permutations: int = 999,
    adonis_formula: Optional[str] = None,
    run_alpha_correlation: bool = True,
    dry_run: bool = False,
    show_qiime: bool = False,
) -> Dict[str, object]:
    """
    Run alpha-group-significance / alpha-correlation (when possible),
    beta-group-significance (when distance matrices exist and columns are valid),
    and optionally ADONIS. Skips gracefully when preconditions aren't met.
    Returns a summary dict of what was produced or skipped.
    """
    core = group_dir / "core-metrics-phylo"
    out = {"alpha_group": {}, "alpha_corr": {}, "beta_group": {}, "adonis": {}}

    # Load metadata once
    if not metadata_file.exists():
        LOG.warning("Metadata file missing; skipping all stats for %s", group_dir)
        return out
    header, rows, types_map = _read_metadata_simple(metadata_file)

    # Alpha vectors that may exist (alpha-only or full core-metrics)
    alpha_vecs = {
        "observed_features": core / "observed_features_vector.qza",
        "evenness":          core / "evenness_vector.qza",
        "faith_pd":          core / "faith_pd_vector.qza",
        "shannon":           core / "shannon_vector.qza",
    }
    available_alpha = {k: p for k, p in alpha_vecs.items() if p.exists()}

    # ---------- ALPHA GROUP SIGNIFICANCE ----------
    if available_alpha:
        # use the intersection of sample IDs across the available alpha vectors
        sample_sets = [ _alpha_vec_sample_ids(p) for p in available_alpha.values() ]
        alpha_sids = set.intersection(*sample_sets) if sample_sets else set()
        if len(alpha_sids) < 2:
            LOG.info("Not enough samples for alpha-group-significance in %s (n=%d). Skipping.", group_dir, len(alpha_sids))
            out["alpha_group"]["skipped"] = "too_few_samples"
        else:
            cat_cols, num_cols = _eligible_columns(header, rows, types_map, alpha_sids)
            if not cat_cols:
                LOG.info("No eligible categorical metadata columns for alpha-group-significance in %s. Skipping.", group_dir)
                out["alpha_group"]["skipped"] = "no_valid_categorical_cols"
            else:
                # run once per alpha vector
                ok = True
                for tag, vec in available_alpha.items():
                    viz = core / f"alpha-groups-{tag}.qzv"
                    ok = _safe_call(
                        qiime.diversity_alpha_group_significance,
                        input_alpha_diversity=vec,
                        metadata_file=metadata_file,
                        output_visualization=viz,
                        dry_run=dry_run,
                        show_qiime=show_qiime,
                    ) and ok
                out["alpha_group"]["status"] = "ok" if ok else "partial"
                out["alpha_group"]["columns_used"] = cat_cols
    else:
        out["alpha_group"]["skipped"] = "no_alpha_vectors"

    # ---------- ALPHA CORRELATION ----------
    if run_alpha_correlation and available_alpha:
        sample_sets = [ _alpha_vec_sample_ids(p) for p in available_alpha.values() ]
        alpha_sids = set.intersection(*sample_sets) if sample_sets else set()
        if len(alpha_sids) < 3:
            LOG.info("Not enough samples for alpha-correlation in %s (n=%d). Skipping.", group_dir, len(alpha_sids))
            out["alpha_corr"]["skipped"] = "too_few_samples"
        else:
            cat_cols, num_cols = _eligible_columns(header, rows, types_map, alpha_sids)
            if not num_cols:
                LOG.info("No eligible numeric metadata columns for alpha-correlation in %s. Skipping.", group_dir)
                out["alpha_corr"]["skipped"] = "no_valid_numeric_cols"
            else:
                ok = True
                for tag, vec in available_alpha.items():
                    viz = core / f"alpha-corr-{tag}.qzv"
                    ok = _safe_call(
                        qiime.diversity_alpha_correlation,
                        input_alpha_diversity=vec,
                        metadata_file=metadata_file,
                        output_visualization=viz,
                        dry_run=dry_run,
                        show_qiime=show_qiime,
                    ) and ok
                out["alpha_corr"]["status"] = "ok" if ok else "partial"
                out["alpha_corr"]["columns_used"] = num_cols
    else:
        out["alpha_corr"]["skipped"] = out["alpha_group"].get("skipped", "no_alpha_vectors")

    # ---------- BETA GROUP SIGNIFICANCE / ADONIS ----------
    # Distance matrices present only when full core-metrics succeeded
    dms = {
        "unweighted_unifrac": core / "unweighted_unifrac_distance_matrix.qza",
        "weighted_unifrac":   core / "weighted_unifrac_distance_matrix.qza",
        "bray_curtis":        core / "bray_curtis_distance_matrix.qza",
        "jaccard":            core / "jaccard_distance_matrix.qza",
    }
    dms = {k: p for k, p in dms.items() if p.exists()}

    if not dms:
        out["beta_group"]["skipped"] = "no_distance_matrices"
        out["adonis"]["skipped"] = "no_distance_matrices"
    else:
        # For each DM, compute eligible categorical cols among its samples
        for dm_tag, dm_path in dms.items():
            dm_sids = _dm_sample_ids(dm_path)
            if len(dm_sids) < 3:
                out["beta_group"][dm_tag] = {"skipped": "too_few_samples"}
                continue
            cat_cols, _ = _eligible_columns(header, rows, types_map, dm_sids)

            # If user specified cols, filter to those
            cols_to_use = list(beta_group_cols or [])
            if cols_to_use:
                cols_to_use = [c for c in cols_to_use if c in cat_cols]
            else:
                cols_to_use = cat_cols

            if len(cols_to_use) == 0:
                out["beta_group"][dm_tag] = {"skipped": "no_valid_categorical_cols"}
                continue

            # Run beta-group-significance per requested column
            percol: Dict[str, str] = {}
            for col in cols_to_use:
                viz = core / f"beta-groups-{dm_tag}-{col}.qzv"
                ok = _safe_call(
                    qiime.diversity_beta_group_significance,
                    input_distance_matrix=dm_path,
                    metadata_file=metadata_file,
                    metadata_column=col,
                    method=beta_method,
                    pairwise=pairwise,
                    permutations=permutations,
                    output_visualization=viz,
                    dry_run=dry_run,
                    show_qiime=show_qiime,
                )
                percol[col] = "ok" if ok else "error"
            out["beta_group"][dm_tag] = percol

            # ADONIS (optional, one viz per dm)
            if adonis_formula:
                viz = core / f"adonis-{dm_tag}.qzv"
                ok = _safe_call(
                    qiime.diversity_adonis,
                    input_distance_matrix=dm_path,
                    metadata_file=metadata_file,
                    formula=adonis_formula,
                    permutations=permutations,
                    output_visualization=viz,
                    dry_run=dry_run,
                    show_qiime=show_qiime,
                )
                out["adonis"][dm_tag] = "ok" if ok else "error"

    # Write summary so users see which steps ran
    summary_path = core / "stats_summary.json"
    try:
        summary_path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    except Exception:
        pass
    return out
