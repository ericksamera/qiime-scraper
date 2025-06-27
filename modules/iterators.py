# modules/iterators.py

from pathlib import Path
from .io_utils import run_command
from .logger import log_success
import logging
from modules import qiime_wrapper
from typing import Optional
import json, re, statistics, zipfile
import zipfile

import statistics, zipfile, json, re
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger("qiime_pipeline")


_SCRIPT_RE = re.compile(
    r'<script[^>]*id=["\']table-data["\'][^>]*>\s*(\{.*?\})\s*</script>',
    re.I | re.S,
)

def _extract_median(qzv: Path) -> Optional[int]:
    """Return the median of values in the Frequency object of *qzv* or None."""
    if not qzv.exists():
        return None
    with zipfile.ZipFile(qzv) as z:
        try:
            html = z.read(
                next(p for p in z.namelist()
                     if Path(p).name == "sample-frequency-detail.html")
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

def get_optimal_trimming(
    imported_qza: Path,
    *,
    max_len: int = 250,
    min_len: int = 200,
    step: int = 25,
) -> Tuple[int, int]:
    """
    Brute-force two-stage search for the (trunc_len_f, trunc_len_r) pair
    that maximises *score*, where *score* = median sample frequency extracted
    from the generated `.qzv`. Returns the best pair.

    The function is restart-safe: if a pair's `.qzv` already exists it is
    reused, letting interrupted runs pick up where they left off.
    """
    out = imported_qza.parent / ".work" / ".optimal_trimming"
    out.mkdir(parents=True, exist_ok=True)

    def run_pair(f: int, r: int) -> Optional[int]:
        tag = f"{f}-{r}"
        table_qza = out / f"{tag}_output_table.qza"
        table_qzv = out / f"{tag}_output_table.qzv"

        if not table_qzv.exists():
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

    # persist result for later reuse
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