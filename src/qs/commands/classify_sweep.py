# src/qs/commands/classify_sweep.py
from __future__ import annotations

import json
import statistics
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple

from qs.utils.logger import get_logger
from qs.qiime import commands as qiime

LOG = get_logger("classify")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "classify-sweep", parents=[parent],
        help="Run a classifier sweep and pick the best model by multiple metrics.",
    )
    p.add_argument("--input-reads", type=Path, required=True, help="FeatureData[Sequence] .qza (rep-seqs).")
    p.add_argument("--classifiers-dir", type=Path, required=True, help="Directory of *.qza sklearn classifiers.")
    p.add_argument("--output-dir", type=Path, required=True, help="Directory to write classifications and report.")
    p.set_defaults(func=run)


def _extract_metrics(tax_qza: Path) -> Dict[str, float]:
    """
    Read taxonomy.tsv inside taxonomy.qza and compute summary metrics.
    """
    try:
        with zipfile.ZipFile(tax_qza) as z:
            tsv = z.read(next(p for p in z.namelist() if p.endswith("taxonomy.tsv"))).decode()
    except Exception:
        return {}
    lines = tsv.splitlines()
    if len(lines) <= 1:
        return {}
    depths, conf = [], []
    full7 = 0
    for ln in lines[1:]:
        parts = ln.strip().split("\t")
        if len(parts) != 3:
            continue
        taxon, cstr = parts[1], parts[2]
        try:
            cval = float(cstr)
        except Exception:
            continue
        depth = len(taxon.split(";")) if taxon else 0
        depths.append(depth)
        conf.append(cval)
        if depth >= 7:
            full7 += 1
    if not conf:
        return {}
    return {
        "n_records": float(len(conf)),
        "mean_conf": float(statistics.mean(conf)),
        "median_conf": float(statistics.median(conf)),
        "mean_depth": float(statistics.mean(depths)) if depths else 0.0,
        "n_depth≥7": float(full7),
        "pct_depth≥7": float(full7 / len(conf)),
    }


def sweep_and_pick(
    *,
    input_reads: Path,
    classifiers_dir: Path,
    out_dir: Path,
    dry_run: bool,
    show_qiime: bool,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_list: List[Tuple[str, Dict[str, float]]] = []

    for cls in sorted(classifiers_dir.glob("*.qza")):
        tag = cls.stem
        out_tax = out_dir / f"{tag}_classification.qza"
        if not out_tax.exists():
            qiime.classify_sklearn(
                input_reads=input_reads,
                input_classifier=cls,
                output_classification=out_tax,
                dry_run=dry_run,
                show_stdout=show_qiime,
            )
        m = _extract_metrics(out_tax)
        if m:
            metrics_list.append((tag, m))
        else:
            LOG.warning("No metrics parsed for %s", cls.name)

    if not metrics_list:
        raise RuntimeError("No valid classifications produced.")

    rank_keys = ("pct_depth≥7", "median_conf", "mean_conf")
    winners: Dict[str, List[Dict[str, object]]] = {}
    for key in rank_keys:
        top = sorted(metrics_list, key=lambda kv: kv[1].get(key, 0.0), reverse=True)[:2]
        winners[key] = [{"classifier": t, "value": kv[1].get(key)} for (t, kv) in top]

    report = {"winners": winners, "all_metrics": {t: m for t, m in metrics_list}}
    (out_dir / "optimal_classifiers.json").write_text(json.dumps(report, indent=2) + "\n")
    LOG.info("Classifier sweep report → %s", out_dir / "optimal_classifiers.json")

    # choose a single final winner by priority
    for key in rank_keys:
        pick = winners.get(key, [])
        if pick:
            return {"priority": key, "classifier_tag": pick[0]["classifier"], "report_path": str(out_dir / "optimal_classifiers.json")}
    # fallback
    return {"priority": "mean_conf", "classifier_tag": sorted(metrics_list, key=lambda kv: kv[1].get("mean_conf", 0.0), reverse=True)[0][0]}
    

def run(args) -> None:
    res = sweep_and_pick(
        input_reads=args.input_reads,
        classifiers_dir=args.classifiers_dir,
        out_dir=args.output_dir,
        dry_run=getattr(args, "dry_run", False),
        show_qiime=getattr(args, "show_qiime", True),
    )
    LOG.info("Winner: %s (by %s)", res["classifier_tag"], res["priority"])
    print(f"[ok] winner={res['classifier_tag']} by {res['priority']}")
