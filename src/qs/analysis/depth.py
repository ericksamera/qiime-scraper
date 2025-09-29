# src/qs/analysis/depth.py
from __future__ import annotations

import json
import re
import statistics
import zipfile
from pathlib import Path
from typing import Iterable, List, Optional


def _try_parse_table_json(html: str) -> Optional[List[int]]:
    m = re.search(r'id=["\\\']table-data["\\\'][^>]*>\s*(\{.*?\})\s*<', html, flags=re.S)
    if not m:
        return None
    try:
        data = json.loads(m.group(1))
        freqs = list(map(int, (data.get("Frequency", {}) or {}).values()))
        return freqs or None
    except Exception:
        return None


def _extract_sample_frequencies(table_qzv: Path) -> List[int]:
    """
    Parse the feature-table summary .qzv and return per-sample frequencies.

    Strategy:
      1) Look for sample-frequency-detail.html and parse embedded JSON.
      2) Fallback: scan all HTML files in the archive for a 'table-data' blob.
      3) If nothing works, return [] (caller decides fallback).
    """
    if not table_qzv.exists():
        return []
    try:
        with zipfile.ZipFile(table_qzv) as z:
            # 1) canonical path
            names = z.namelist()
            prefer = [p for p in names if p.endswith("sample-frequency-detail.html")]
            if prefer:
                html = z.read(prefer[0]).decode("utf-8", errors="replace")
                freqs = _try_parse_table_json(html)
                if freqs:
                    return freqs
            # 2) scan all htmls
            for p in names:
                if not p.lower().endswith(".html"):
                    continue
                try:
                    html = z.read(p).decode("utf-8", errors="replace")
                except Exception:
                    continue
                freqs = _try_parse_table_json(html)
                if freqs:
                    return freqs
    except Exception:
        return []
    return []


def mean_sample_depth(table_qzv: Path) -> int:
    freqs = _extract_sample_frequencies(table_qzv)
    if not freqs:
        return 0
    return int(round(statistics.mean(freqs)))


def choose_sampling_depth(
    table_qzv: Path,
    retain_fraction: float = 0.90,
    min_depth: int = 1000,
) -> int:
    freqs = sorted(_extract_sample_frequencies(table_qzv))
    if not freqs:
        return min_depth
    k = int((1.0 - retain_fraction) * len(freqs))
    k = max(0, min(k, len(freqs) - 1))
    return max(min_depth, freqs[k])


# ---------- NEW: resilience helpers ------------------------------------------

def count_samples_at_or_above(table_qzv: Path, depth: int) -> int:
    """How many samples have total counts ≥ depth according to table.qzv?"""
    if depth <= 0:
        return 0
    freqs = _extract_sample_frequencies(table_qzv)
    return sum(1 for f in freqs if f >= depth)


def max_depth_with_at_least(table_qzv: Path, n_samples: int, *, floor: int = 1) -> int:
    """
    Largest depth such that ≥ n_samples have counts ≥ that depth.
    If impossible, returns 'floor'.
    """
    freqs = sorted(_extract_sample_frequencies(table_qzv), reverse=True)
    if not freqs:
        return floor
    if n_samples <= 0:
        return floor
    if n_samples > len(freqs):
        n_samples = len(freqs)
    candidate = freqs[n_samples - 1]  # nth largest
    return max(floor, int(candidate))


def choose_depth_and_count(
    table_qzv: Path,
    retain_fraction: float = 0.90,
    min_depth: int = 1000,
) -> tuple[int, int]:
    """
    Return (depth, n_retained) where 'depth' is the quantile-based choice
    and 'n_retained' is how many samples meet that depth.
    """
    depth = choose_sampling_depth(table_qzv, retain_fraction=retain_fraction, min_depth=min_depth)
    return depth, count_samples_at_or_above(table_qzv, depth)
