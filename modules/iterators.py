# modules/iterators.py

from pathlib import Path
from .io_utils import run_command
from .logger import log_success
import logging

from typing import Optional

logger = logging.getLogger("qiime_pipeline")

def get_optimal_trimming(
        imported_qza: Path,
        max_len: int = 250,
        min_len: int = 200,
        step: int = 25
    ) -> None:
    """
    """

    output_path = imported_qza.parent.joinpath('.work').joinpath('.optimal_trimming')
    output_path.mkdir(exist_ok=True, parents=True)

    # coarse-search
    coarse_results = []
    for trunc_len_f in range(min_len, max_len+1, step):
        for trunc_len_r in range(min_len, max_len+1, step):
            score = evaluate_quality(imported_qza, trunc_len_f, trunc_len_r)  # custom function
            coarse_results.append(((trunc_len_f, trunc_len_r), score))
    
    top_2 = sorted(coarse_results, key=lambda x: x[1], reverse=True)[:2]

    fine_results = []

    for (f_base, r_base), _ in top_2:
        for f in range(f_base - 10, f_base + 11, 5):
            for r in range(r_base - 10, r_base + 11, 5):
                if min_len <= f <= max_len and min_len <= r <= max_len:
                    score = evaluate_quality(imported_qza, f, r)
                    fine_results.append(((f, r), score))

    # Final best pair
    best_pair = max(fine_results, key=lambda x: x[1])