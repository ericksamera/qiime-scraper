# src/qs/utils/text.py
from __future__ import annotations

import re
import hashlib


def slugify(s: str, max_len: int = 60) -> str:
    """Filesystem-safe slug; keeps [A-Za-z0-9._-]. Falls back to short hash if needed."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "-", s).strip("-_.")
    if not s:
        return "x"
    if len(s) <= max_len:
        return s
    h = hashlib.sha1(s.encode("utf-8")).hexdigest()[:8]
    return s[: max_len - 9] + "-" + h
