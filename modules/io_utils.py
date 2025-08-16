# modules/io_utils.py

from __future__ import annotations

import csv
import re
import subprocess
import logging
import time

from pathlib import Path
from typing import Iterable, Sequence, List

from modules.logger import log_success

logger = logging.getLogger("qiime_pipeline")

# ---------------------------------------------------------------------
# Manifest generation
# ---------------------------------------------------------------------

# Matches Illumina-style filenames:
#   <root>_R1_001.fastq.gz  /  <root>_R2_001.fastq.gz
# Also accepts: <root>_R1.fastq.gz / <root>_R2.fastq.gz
_ILLUMINA_RE = re.compile(r"^(?P<root>.+)_R(?P<read>[12])(?:_[0-9]{3})?\.fastq\.gz$")

def generate_manifest(project_fastq_dir: Path, manifest_dir: Path) -> Path:
    """
    Create a PairedEndFastqManifestPhred33V2 CSV at <manifest_dir>/fastq.manifest with absolute paths.

    Looks for *.fastq.gz directly under `project_fastq_dir`. If your FASTQs are nested
    in subfolders, change the loop to `project_fastq_dir.rglob("*.fastq.gz")`.
    """
    fastq_dir = project_fastq_dir

    # Collect R1/R2 pairs by parsing actual filenames (no string reconstruction).
    pairs: dict[str, dict[str, str]] = {}
    for fq in fastq_dir.glob("*.fastq.gz"):
        m = _ILLUMINA_RE.match(fq.name)
        if not m:
            continue
        root = m.group("root")
        read = m.group("read")  # "1" or "2"
        pairs.setdefault(root, {})[read] = str(fq.resolve())

    rows = []
    for root, d in sorted(pairs.items()):
        if "1" in d and "2" in d:
            rows.append(
                {
                    "sample-id": root,
                    "forward-absolute-filepath": d["1"],
                    "reverse-absolute-filepath": d["2"],
                }
            )

    if not rows:
        raise RuntimeError(f"No R1/R2 pairs found in {fastq_dir}")

    manifest_path = manifest_dir / "fastq.manifest"
    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sample-id",
                "forward-absolute-filepath",
                "reverse-absolute-filepath",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Wrote manifest: {manifest_path}")
    return manifest_path


def run_command(
        command: List[str],
        dry_run=False,
        capture=False):

    """
    """

    logger.info(f"Running: {' '.join(command)}")
    if dry_run:
        logger.debug("Dry-run mode: command not executed.")
        return None
    
    start = time.time()
    result = subprocess.run(command, check=True, capture_output=capture)
    end = time.time()

    duration = end - start
    log_success(f"Command completed in {duration:.2f}s")

    return result