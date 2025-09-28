# src/qs/analysis/depth.py
from __future__ import annotations

import json
import re
import statistics
import zipfile
from pathlib import Path
from typing import List


def _extract_sample_frequencies(table_qzv: Path) -> List[int]:
    """
    Parse the feature-table summary .qzv and return per-sample frequencies.
    Robust to minor HTML changes; falls back to [] on errors.
    """
    if not table_qzv.exists():
        return []
    try:
        with zipfile.ZipFile(table_qzv) as z:
            html = z.read(
                next(p for p in z.namelist() if p.endswith("sample-frequency-detail.html"))
            ).decode()
    except Exception:
        return []
    m = re.search(r'id=["\\\']table-data["\\\'][^>]*>\s*(\{.*?\})\s*<', html)
    if not m:
        return []
    try:
        data = json.loads(m.group(1))
        freqs = list(map(int, (data.get("Frequency", {}) or {}).values()))
        return freqs
    except Exception:
        return []


def mean_sample_depth(table_qzv: Path) -> int:
    """
    Integer mean of per-sample depths from a .qzv; 0 on failure.
    """
    freqs = _extract_sample_frequencies(table_qzv)
    if not freqs:
        return 0
    return int(round(statistics.mean(freqs)))


def choose_sampling_depth(
    table_qzv: Path,
    retain_fraction: float = 0.90,
    min_depth: int = 1000,
) -> int:
    """
    Choose a rarefaction depth that retains ~retain_fraction of samples,
    but never less than min_depth.
    """
    freqs = sorted(_extract_sample_frequencies(table_qzv))
    if not freqs:
        return min_depth
    k = int((1.0 - retain_fraction) * len(freqs))
    k = max(0, min(k, len(freqs) - 1))
    return max(min_depth, freqs[k])
