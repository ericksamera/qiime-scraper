# src/qs/primers/detect.py
from __future__ import annotations
import gzip
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from qs.utils.samples import discover_fastqs, collect_pairs, normalize_id

def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def _count_prefixes(fp: Path, *, kmin: int, kmax: int, max_reads: int) -> tuple[Dict[int, Counter], int]:
    counters: Dict[int, Counter] = {k: Counter() for k in range(kmin, kmax + 1)}
    n = 0
    try:
        with _open_text(fp) as fh:
            while n < max_reads:
                if not fh.readline(): break  # header
                seq = (fh.readline() or "").strip().upper()  # sequence
                fh.readline()  # +
                fh.readline()  # qual
                if not seq: break
                n += 1
                # tolerate a couple Ns in the window
                if seq[:kmax].count("N") > 2:
                    continue
                win = seq[:kmax]
                for k in range(kmin, kmax + 1):
                    if len(win) >= k:
                        counters[k][win[:k]] += 1
    except Exception:
        pass
    return counters, n

def _choose_prefix(counters: Dict[int, Counter], total: int, *, min_frac: float) -> str:
    if total <= 0:
        return ""
    best = ("", 0.0, 0)  # (seq, frac, k)
    for k, ctr in counters.items():
        if not ctr: continue
        seq, c = ctr.most_common(1)[0]
        frac = c / float(total)
        # prefer higher fraction; tie-break longer k
        if (frac > best[1]) or (abs(frac - best[1]) < 1e-9 and k > best[2]):
            best = (seq, frac, k)
    return best[0] if best[1] >= min_frac else ""

def _pairs_by_sid(
    fastq_dir: Path,
    *,
    strip_illumina_suffix: bool,
    id_regex: Optional[str],
) -> Dict[str, tuple[Path, Path]]:
    paths = discover_fastqs(fastq_dir)
    pairs = collect_pairs(paths)
    out: Dict[str, tuple[Path, Path]] = {}
    for root, d in pairs.items():
        if "1" not in d or "2" not in d: continue
        sid = normalize_id(root, strip_illumina_suffix=strip_illumina_suffix, id_regex=id_regex)
        if sid not in out:
            out[sid] = (Path(d["1"]), Path(d["2"]))
    return out

def detect_primers_for_samples(
    fastq_dir: Path,
    sample_ids: Iterable[str],
    *,
    strip_illumina_suffix: bool = True,
    id_regex: Optional[str] = None,
    kmin: int = 16,
    kmax: int = 24,
    max_reads: int = 4000,
    min_fraction: float = 0.30,
) -> Dict[str, tuple[str, str]]:
    """Return {sample-id: (forward_from_R1, reverse_from_R2_as_seen)}"""
    sid2pair = _pairs_by_sid(fastq_dir, strip_illumina_suffix=strip_illumina_suffix, id_regex=id_regex)
    want = {s.lstrip("#") for s in sample_ids}
    out: Dict[str, tuple[str, str]] = {}
    for sid, (r1, r2) in sid2pair.items():
        if sid not in want: continue
        c1, n1 = _count_prefixes(r1, kmin=kmin, kmax=kmax, max_reads=max_reads)
        c2, n2 = _count_prefixes(r2, kmin=kmin, kmax=kmax, max_reads=max_reads)
        f = _choose_prefix(c1, n1, min_frac=min_fraction)
        r = _choose_prefix(c2, n2, min_frac=min_fraction)  # as it appears in R2
        if f and r:
            out[sid] = (f, r)
    return out

def detect_groups_for_samples(
    fastq_dir: Path,
    sample_ids: Iterable[str],
    *,
    strip_illumina_suffix: bool = True,
    id_regex: Optional[str] = None,
    kmin: int = 16,
    kmax: int = 24,
    max_reads: int = 4000,
    min_fraction: float = 0.30,
) -> tuple[Dict[str, List[str]], Dict[str, tuple[str, str]]]:
    """Return (groups_map, per_sample_primers) where groups_map is 'F|R' -> [sample-ids]"""
    per = detect_primers_for_samples(
        fastq_dir, sample_ids,
        strip_illumina_suffix=strip_illumina_suffix, id_regex=id_regex,
        kmin=kmin, kmax=kmax, max_reads=max_reads, min_fraction=min_fraction,
    )
    groups: Dict[str, List[str]] = defaultdict(list)
    for sid, (f, r) in per.items():
        groups[f"{f}|{r}"].append(sid)
    for sids in groups.values():
        sids.sort()
    return groups, per
