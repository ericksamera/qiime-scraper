# src/qs/optimize/truncation.py
from __future__ import annotations

import gzip
import json
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from qs.utils.logger import get_logger
from qs.qiime import commands as qiime
from qs.utils.samples import ILLUMINA_RE

LOG = get_logger("opt.trunc")


def _guess_read_len_from_fastqs(run_path: Path, read: str) -> int:
    """
    Best-effort: look for an R1/R2 FASTQ in run_path and return first sequence length.
    Falls back to 150 if nothing found.
    """
    patterns = ["*.fastq.gz", "*.fastq"]
    for pat in patterns:
        for fq in run_path.rglob(pat):
            m = ILLUMINA_RE.match(fq.name)
            if not m:
                continue
            if m.group("read") != read:
                continue
            try:
                if fq.suffix == ".gz":
                    with gzip.open(fq, "rt") as f:
                        f.readline()  # @header
                        seq = f.readline().strip()
                else:
                    with fq.open("rt") as f:
                        f.readline()
                        seq = f.readline().strip()
                if seq:
                    return len(seq)
            except Exception:
                continue
    return 150


def _parse_dada2_stats(stats_qza: Path) -> Tuple[int, int]:
    """
    Return (total_input, total_merged) from a DADA2 denoising-stats artifact.
    Robust to column order; looks for columns named like 'input' and 'merged'.
    """
    if not stats_qza.exists():
        return 0, 0
    try:
        with zipfile.ZipFile(stats_qza) as z:
            # Find first .tsv inside
            tsv_name = next(p for p in z.namelist() if p.endswith(".tsv") or p.endswith(".txt"))
            data = z.read(tsv_name).decode("utf-8")
    except Exception:
        return 0, 0

    lines = [ln for ln in data.splitlines() if ln and not ln.startswith("#")]
    if not lines:
        return 0, 0
    header = [c.strip().lower() for c in lines[0].split("\t")]
    try:
        idx_in = header.index("input")
        idx_mg = header.index("merged")
    except ValueError:
        # fallback heuristic: common positions are 1 and 4
        idx_in, idx_mg = 1, 4
    total_in, total_mg = 0, 0
    for ln in lines[1:]:
        cols = ln.split("\t")
        try:
            total_in += int(cols[idx_in])
            total_mg += int(cols[idx_mg])
        except Exception:
            continue
    return total_in, total_mg


def _score_pair(stats_qza: Path) -> Tuple[float, int, int]:
    ti, tm = _parse_dada2_stats(stats_qza)
    frac = (tm / ti) if ti else 0.0
    return frac, ti, tm


def _trial_run(
    *,
    input_seqs: Path,
    f: int,
    r: int,
    outdir: Path,
    cores: int,
    quick_learn: Optional[int],
    dry_run: bool,
    show_qiime: bool,
    other_params: Dict[str, object],
) -> Tuple[float, int, int, Path]:
    """Run DADA2 once for a given (f, r) and return (frac, total_in, total_mg, stats_path)."""
    table = outdir / f"F{f}_R{r}_table.qza"
    rep   = outdir / f"F{f}_R{r}_rep.qza"
    stats = outdir / f"F{f}_R{r}_stats.qza"
    qiime.dada2_denoise_paired(
        input_seqs=input_seqs,
        trunc_len_f=f,
        trunc_len_r=r,
        output_table=table,
        output_rep_seqs=rep,
        output_stats=stats,
        n_threads=cores,
        n_reads_learn=quick_learn,
        dry_run=dry_run,
        show_stdout=show_qiime,
        **other_params,
    )
    frac, ti, tm = _score_pair(stats)
    return frac, ti, tm, stats


def _neighbors(center: int, lo: int, hi: int, step: int) -> List[int]:
    cand = sorted(set([
        max(lo, min(hi, center)),
        max(lo, min(hi, center - step)),
        max(lo, min(hi, center + step)),
    ]))
    return [c for c in cand if lo <= c <= hi]


def find_optimal_truncation(
    *,
    run_path: Path,
    input_seqs: Path,              # demux.qza or trimmed.qza (already cutadapted if used)
    outdir: Path,
    lower_frac: float = 0.80,      # lower bound fraction of read length
    abs_min_len: int = 50,         # hard lower bound
    coarse_step: int = 10,         # coarse step in bp
    refine_step: int = 5,          # refinement step in bp
    cores: int = 0,
    quick_learn: Optional[int] = 250000,
    dry_run: bool = False,
    show_qiime: bool = False,
    # DADA2 knobs that should match the final run
    trim_left_f: int = 0,
    trim_left_r: int = 0,
    max_ee_f: int = 2,
    max_ee_r: int = 2,
    trunc_q: int = 2,
    min_overlap: int = 12,
    pooling_method: str = "independent",
    chimera_method: str = "consensus",
) -> Dict[str, object]:
    """
    Coarse-to-fine search maximizing merged / input.
    Returns dict with best_f, best_r, best_fraction, trials, bounds.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Bounds from observed read lengths in run FASTQs.
    f_len_max = _guess_read_len_from_fastqs(run_path, "1")
    r_len_max = _guess_read_len_from_fastqs(run_path, "2")
    f_hi, r_hi = f_len_max, r_len_max
    f_lo = max(abs_min_len, int(round(f_hi * lower_frac)))
    r_lo = max(abs_min_len, int(round(r_hi * lower_frac)))

    other = dict(
        trim_left_f=trim_left_f, trim_left_r=trim_left_r,
        max_ee_f=max_ee_f, max_ee_r=max_ee_r,
        trunc_q=trunc_q, min_overlap=min_overlap,
        pooling_method=pooling_method, chimera_method=chimera_method,
    )

    # 2) Coarse grid
    f_vals = list(range(f_hi, f_lo - 1, -coarse_step)) or [f_hi, f_lo]
    r_vals = list(range(r_hi, r_lo - 1, -coarse_step)) or [r_hi, r_lo]

    trials: List[Dict[str, object]] = []
    best = {"f": f_hi, "r": r_hi, "frac": -1.0, "ti": 0, "tm": 0, "stats": None}

    for f in f_vals:
        for r in r_vals:
            frac, ti, tm, sp = _trial_run(
                input_seqs=input_seqs, f=f, r=r, outdir=outdir,
                cores=cores, quick_learn=quick_learn, dry_run=dry_run, show_qiime=show_qiime, other_params=other,
            )
            trials.append({"f": f, "r": r, "fraction": frac, "input": ti, "merged": tm, "stats": str(sp)})
            if (frac > best["frac"]) or (abs(frac - best["frac"]) < 1e-9 and (f + r) > (best["f"] + best["r"])):
                best = {"f": f, "r": r, "frac": frac, "ti": ti, "tm": tm, "stats": str(sp)}

    # 3) Refine around the coarse best
    f_neigh = _neighbors(best["f"], f_lo, f_hi, refine_step)
    r_neigh = _neighbors(best["r"], r_lo, r_hi, refine_step)
    for f in f_neigh:
        for r in r_neigh:
            # Skip already-evaluated pairs
            if any(t["f"] == f and t["r"] == r for t in trials):
                continue
            frac, ti, tm, sp = _trial_run(
                input_seqs=input_seqs, f=f, r=r, outdir=outdir,
                cores=cores, quick_learn=quick_learn, dry_run=dry_run, show_qiime=show_qiime, other_params=other,
            )
            trials.append({"f": f, "r": r, "fraction": frac, "input": ti, "merged": tm, "stats": str(sp)})
            if (frac > best["frac"]) or (abs(frac - best["frac"]) < 1e-9 and (f + r) > (best["f"] + best["r"])):
                best = {"f": f, "r": r, "frac": frac, "ti": ti, "tm": tm, "stats": str(sp)}

    result = {
        "bounds": {"f": {"lo": f_lo, "hi": f_hi}, "r": {"lo": r_lo, "hi": r_hi}},
        "best": {"trunc_len_f": best["f"], "trunc_len_r": best["r"], "merged_fraction": best["frac"], "input": best["ti"], "merged": best["tm"]},
        "trials": trials,
    }
    (outdir / "optimize_trunc.json").write_text(json.dumps(result, indent=2) + "\n")
    LOG.info("Best truncation: F=%d, R=%d (merged=%.1f%%)", best["f"], best["r"], best["frac"] * 100.0)
    return result
