# src/qs/commands/classify_sweep.py
from __future__ import annotations
import json, statistics, zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Union

from qs.utils.logger import get_logger
from qs.qiime import commands as qiime

LOG = get_logger("classify")

def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser("classify-sweep", parents=[parent], help="Sweep sklearn classifiers and pick a winner.")
    p.add_argument("--input-reads", type=Path, required=True)
    p.add_argument("--classifiers-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--reads-per-batch", default="1000")
    p.add_argument("--n-jobs", type=int, default=1)
    p.add_argument("--blas-threads", type=int, default=1)
    p.set_defaults(func=run)

def _extract_metrics(qza: Path) -> Dict[str, float]:
    try:
        with zipfile.ZipFile(qza) as z:
            tsv = z.read(next(p for p in z.namelist() if p.endswith("taxonomy.tsv"))).decode()
    except Exception:
        return {}
    lines = tsv.splitlines()
    if len(lines) <= 1: return {}
    depths, conf, full7 = [], [], 0
    for ln in lines[1:]:
        cols = ln.split("\t")
        if len(cols) < 3: continue
        taxon, cstr = cols[1], cols[2]
        try:
            c = float(cstr)
        except Exception:
            continue
        depth = len(taxon.split(";")) if taxon else 0
        depths.append(depth); conf.append(c)
        if depth >= 7: full7 += 1
    if not conf: return {}
    return {
        "n_records": float(len(conf)),
        "mean_conf": float(statistics.mean(conf)),
        "median_conf": float(statistics.median(conf)),
        "mean_depth": float(statistics.mean(depths)) if depths else 0.0,
        "n_depth≥7": float(full7),
        "pct_depth≥7": float(full7 / len(conf)),
    }

def _parse_rpb(x: Union[str, int]) -> Union[str, int]:
    if isinstance(x, int): return max(1, x)
    s = str(x).strip().lower()
    if s == "auto": return "auto"
    try:
        return max(1, int(s))
    except Exception:
        return 1000

def classify_with_single(
    *, input_reads: Path, classifier_qza: Path, out_dir: Path,
    dry_run: bool, show_qiime: bool, reads_per_batch: Union[str, int] = 1000, n_jobs: int = 1, blas_threads: int = 1,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    env = {
        "OMP_NUM_THREADS": str(blas_threads),
        "OPENBLAS_NUM_THREADS": str(blas_threads),
        "MKL_NUM_THREADS": str(blas_threads),
        "NUMEXPR_NUM_THREADS": str(blas_threads),
    }
    tag = classifier_qza.stem
    out_tax = out_dir / f"{tag}_classification.qza"
    qiime.classify_sklearn(
        input_reads=input_reads, input_classifier=classifier_qza, output_classification=out_tax,
        reads_per_batch=_parse_rpb(reads_per_batch), n_jobs=n_jobs, pre_dispatch="1*n_jobs", confidence=0.7,
        read_orientation="auto", dry_run=dry_run, show_stdout=show_qiime, extra_env=env,
    )
    m = _extract_metrics(out_tax)
    winners = {k: [{"classifier": tag, "value": m.get(k)}] for k in ("pct_depth≥7", "median_conf", "mean_conf")}
    rpt = {"winners": winners, "all_metrics": {tag: m}}
    (out_dir / "optimal_classifiers.json").write_text(json.dumps(rpt, indent=2) + "\n")
    return {"priority": "fixed", "classifier_tag": tag, "report_path": str(out_dir / "optimal_classifiers.json")}

def sweep_and_pick(
    *, input_reads: Path, classifiers_dir: Path, out_dir: Path,
    dry_run: bool, show_qiime: bool, n_jobs: int = 1, reads_per_batch: Union[str, int] = 1000, blas_threads: int = 1,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    env = {
        "OMP_NUM_THREADS": str(blas_threads),
        "OPENBLAS_NUM_THREADS": str(blas_threads),
        "MKL_NUM_THREADS": str(blas_threads),
        "NUMEXPR_NUM_THREADS": str(blas_threads),
    }
    rpb = _parse_rpb(reads_per_batch)
    metrics: List[Tuple[str, Dict[str, float]]] = []
    for cls in sorted(classifiers_dir.glob("*.qza")):
        tag = cls.stem
        out_tax = out_dir / f"{tag}_classification.qza"
        if not out_tax.exists():
            qiime.classify_sklearn(
                input_reads=input_reads, input_classifier=cls, output_classification=out_tax,
                reads_per_batch=rpb, n_jobs=n_jobs, pre_dispatch="1*n_jobs", confidence=0.7,
                read_orientation="auto", dry_run=dry_run, show_stdout=show_qiime, extra_env=env,
            )
        m = _extract_metrics(out_tax)
        if m: metrics.append((tag, m))
    if not metrics:
        raise RuntimeError("No valid classifications produced.")
    keys = ("pct_depth≥7", "median_conf", "mean_conf")
    winners: Dict[str, List[Dict[str, object]]] = {}
    for k in keys:
        top = sorted(metrics, key=lambda it: it[1].get(k, 0.0), reverse=True)[:2]
        winners[k] = [{"classifier": tag, "value": m.get(k)} for (tag, m) in top]
    rpt = {"winners": winners, "all_metrics": {t: m for t, m in metrics}}
    (out_dir / "optimal_classifiers.json").write_text(json.dumps(rpt, indent=2) + "\n")
    for k in keys:
        picks = winners.get(k, [])
        if picks:
            return {"priority": k, "classifier_tag": picks[0]["classifier"], "report_path": str(out_dir / "optimal_classifiers.json")}
    return {"priority": "mean_conf", "classifier_tag": sorted(metrics, key=lambda it: it[1].get("mean_conf", 0.0), reverse=True)[0][0],
            "report_path": str(out_dir / "optimal_classifiers.json")}

def run(args) -> None:
    res = sweep_and_pick(
        input_reads=args.input_reads, classifiers_dir=args.classifiers_dir, out_dir=args.output_dir,
        dry_run=getattr(args, "dry_run", False), show_qiime=getattr(args, "show_qiime", True),
        n_jobs=args.n_jobs, reads_per_batch=args.reads_per_batch, blas_threads=args.blas_threads,
    )
    LOG.info("Winner: %s by %s", res["classifier_tag"], res["priority"])
    print(f"[ok] winner={res['classifier_tag']} by {res['priority']}")
