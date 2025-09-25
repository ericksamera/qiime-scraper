# src/qs/analysis/depth.py
from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path


def choose_sampling_depth(
    table_qzv: Path,
    retain_fraction: float = 0.90,
    min_depth: int = 1000,
) -> int:
    """
    Parse a feature-table summary .qzv to choose a depth retaining ~retain_fraction of samples.
    Falls back to min_depth on parse errors.
    """
    if not table_qzv.exists():
        return min_depth
    try:
        with zipfile.ZipFile(table_qzv) as z:
            # QIIME's HTML holds a JSON blob in sample-frequency-detail.html
            html = z.read(
                next(p for p in z.namelist() if p.endswith("sample-frequency-detail.html"))
            ).decode()
    except Exception:
        return min_depth

    m = re.search(r'id=["\\\']table-data["\\\'][^>]*>\s*(\{.*?\})\s*<', html)
    if not m:
        return min_depth
    try:
        data = json.loads(m.group(1))
        freqs = sorted(map(int, (data.get("Frequency", {}) or {}).values()))
        if not freqs:
            return min_depth
    except Exception:
        return min_depth
    k = int((1.0 - retain_fraction) * len(freqs))
    k = max(0, min(k, len(freqs) - 1))
    return max(min_depth, freqs[k])
