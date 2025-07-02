# modules/iterators.py

import json
import logging
import re
import statistics
import zipfile
from pathlib import Path
from typing import Optional, Tuple, List, Sequence, Dict

from modules import qiime_wrapper

logger = logging.getLogger("qiime_pipeline")


_SCRIPT_RE = re.compile(
    r'<script[^>]*id=["\']table-data["\'][^>]*>\s*(\{.*?\})\s*</script>',
    re.I | re.S,
)

# -----------------------------------------------------------------------------
# helper: sample‑frequency median
# -----------------------------------------------------------------------------

def _extract_median(qzv: Path) -> Optional[int]:
    """Return the median Frequency value contained in *qzv* or **None**."""
    if not qzv.exists():
        return None

    with zipfile.ZipFile(qzv) as z:
        try:
            html = z.read(
                next(p for p in z.namelist() if Path(p).name == "sample-frequency-detail.html")
            ).decode()
        except (KeyError, StopIteration):
            return None

    m = _SCRIPT_RE.search(html)
    if not m:
        return None

    try:
        vals = json.loads(m.group(1))["Frequency"].values()
        return int(statistics.median(vals))
    except Exception:
        return None

# -----------------------------------------------------------------------------
# brute‑force trimming search
# -----------------------------------------------------------------------------

def get_optimal_trimming(
    imported_qza: Path,
    *,
    max_len: int = 250,
    min_len: int = 200,
    step: int = 25,
) -> Tuple[int, int]:
    """Return the trunc‑len pair (F, R) that maximises median sample frequency."""

    out = imported_qza.parent / ".work" / ".optimal_trimming"
    out.mkdir(parents=True, exist_ok=True)

    def run_pair(f: int, r: int) -> Optional[int]:
        tag = f"{f}-{r}"
        table_qza = out / f"{tag}_output_table.qza"
        table_qzv = out / f"{tag}_output_table.qzv"

        if table_qzv.exists():
            logger.info("%s exists, skipping.", table_qzv)
        else:
            qiime_wrapper.dada2_denoise_paired(
                input_seqs=imported_qza,
                trunc_len_f=f,
                trunc_len_r=r,
                output_table=table_qza,
                output_rep_seqs=out / f"{tag}_output_rep_seqs.qza",
                output_denoising_stats=out / f"{tag}_output_denoising_stats.qza",
            )
            qiime_wrapper.feature_table_summarize(table_qza, table_qzv)

        return _extract_median(table_qzv)

    # --- coarse search -------------------------------------------------
    coarse = [
        ((f, r), s)
        for f in range(min_len, max_len + 1, step)
        for r in range(min_len, max_len + 1, step)
        if (s := run_pair(f, r)) is not None
    ]
    if not coarse:
        raise RuntimeError("No valid scores produced.")

    best_score = max(s for _, s in coarse)
    seeds = [p for p, s in coarse if s >= 0.9 * best_score]

    # --- fine search ---------------------------------------------------
    fine = [
        ((f, r), s)
        for f_base, r_base in seeds
        for f in range(f_base - 10, f_base + 11, 5)
        for r in range(r_base - 10, r_base + 11, 5)
        if min_len <= f <= max_len
        if min_len <= r <= max_len
        if (s := run_pair(f, r)) is not None
    ]

    best_pair, best_score = max(fine, key=lambda x: x[1])

    with (out / "optimal_trimming.json").open("w") as fh:
        json.dump(
            {
                "trunc_len_f": best_pair[0],
                "trunc_len_r": best_pair[1],
                "median_score": best_score,
            },
            fh,
            indent=2,
        )
    logger.info("Optimal trimming saved: %s", best_pair)

    return best_pair

# -----------------------------------------------------------------------------
# helper: classification metrics
# -----------------------------------------------------------------------------

def _extract_classification_metrics(qza: Path, min_depth: int = 5) -> Optional[dict]:
    """Return summary metrics for the taxonomy assignments in *qza*."""

    if not qza.exists():
        return None

    with zipfile.ZipFile(qza) as z:
        try:
            taxonomy_tsv = z.read(
                next(p for p in z.namelist() if Path(p).name == "taxonomy.tsv")
            ).decode()
        except (KeyError, StopIteration):
            return None

    lines = taxonomy_tsv.splitlines()
    if not lines:
        return None

    depths: List[int] = []
    confidences: List[float] = []
    full7 = 0

    for line in lines[1:]:
        parts = line.strip().split("\t")
        if len(parts) != 3:
            continue

        taxon_depth = len(parts[1].split(";"))
        try:
            conf = float(parts[2])
        except ValueError:
            continue

        depths.append(taxon_depth)
        confidences.append(conf)
        if taxon_depth >= 7:
            full7 += 1

    if not confidences:
        return None

    return {
        "n_records": len(confidences),
        "mean_conf": statistics.mean(confidences),
        "median_conf": statistics.median(confidences),
        "mean_depth": statistics.mean(depths),
        "n_depth≥7": full7,
        "pct_depth≥7": full7 / len(confidences),
        f"n_depth≥{min_depth}": sum(d >= min_depth for d in depths),
    }

# -----------------------------------------------------------------------------
# helper: simple metric picker
# -----------------------------------------------------------------------------

def _score(m: dict, k: str) -> float:
    """Return metric *k* from *m* or 0 if absent."""
    return m.get(k, 0.0)

# -----------------------------------------------------------------------------
# classifier picker
# -----------------------------------------------------------------------------

def get_optimal_classifier(
    imported_qza: Path,
    classifiers_dir: Path,
    *,
    keys: Sequence[str] = ("pct_depth≥7", "median_conf", "mean_conf"),
) -> Dict[str, Path]:
    """Run each classifier and pick the best for every metric in *keys*.

    Returns a mapping ``metric → classifier_path``.
    Writes a report to ``<work>/optimal_classifiers.json`` so interrupted runs
    can resume without recomputing finished classifications.
    """

    out = imported_qza.parent / ".work" / ".optimal_classifier"
    out.mkdir(parents=True, exist_ok=True)

    # collect metrics --------------------------------------------------
    scores: List[Tuple[Path, dict]] = []

    for cls in classifiers_dir.glob("*.qza"):
        tag = cls.stem
        classification = out / f"{tag}_classification.qza"

        if classification.exists():
            logger.info("%s exists, skipping classify‑sklearn.", classification)
        else:
            qiime_wrapper.classify_sklearn(
                input_reads=imported_qza,
                input_classifier=cls,
                output_classification=classification,
                reads_per_batch=5000,
            )

        metrics = _extract_classification_metrics(classification)
        if metrics:
            scores.append((cls, metrics))
        else:
            logger.warning("No metrics for %s – ignored.", cls)

    if not scores:
        raise RuntimeError("No valid classifier metrics produced.")

    metric_map: Dict[Path, dict] = {p: m for p, m in scores}

    # determine leader per metric -------------------------------------
    leaders: Dict[str, Path] = {}
    for k in keys:
        best_cls, best_met = max(scores, key=lambda t: _score(t[1], k))
        leaders[k] = best_cls
        logger.info("Best %s: %s (value=%s)", k, best_cls.name, best_met.get(k))

    # persist winners --------------------------------------------------
    with (out / "optimal_classifiers.json").open("w") as fh:
        json.dump(
            {
                k: {
                    "classifier": leaders[k].name,
                    "value": metric_map[leaders[k]].get(k),
                }
                for k in leaders
            },
            fh,
            indent=2,
        )

    return leaders

# -----------------------------------------------------------------------------
