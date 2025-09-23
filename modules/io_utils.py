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

def _collect_pairs(paths: Iterable[Path]) -> list[dict[str, str]]:
    pairs: dict[str, dict[str, str]] = {}
    for fq in paths:
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
    return rows

def generate_manifest(project_fastq_dir: Path, manifest_dir: Path) -> Path:
    """
    Create a PairedEndFastqManifestPhred33V2 CSV at <manifest_dir>/fastq.manifest with absolute paths.

    First scans *.fastq.gz directly under `project_fastq_dir`.
    If no valid R1/R2 pairs are found, it automatically retries recursively.
    """
    fastq_dir = project_fastq_dir

    # Try flat directory
    rows = _collect_pairs(fastq_dir.glob("*.fastq.gz"))

    # Fallback: recursive search
    if not rows:
        logger.info("No R1/R2 pairs found at top level; retrying recursive search…")
        rows = _collect_pairs(fastq_dir.rglob("*.fastq.gz"))

    if not rows:
        raise RuntimeError(f"No R1/R2 pairs found in {fastq_dir} (flat or recursive)")

    manifest_path = manifest_dir / "fastq.manifest"
    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sample-id",
                "forward-absolute-filepath",
                "reverse-absolute-filepath",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Wrote manifest: {manifest_path}")
    return manifest_path

# modules/io_utils.py — run_command

def run_command(
        command: List[str],
        dry_run: bool = False,
        capture: bool = False,
        env: dict | None = None):
    """
    Run a shell command. On failure, print QIIME's stdout/stderr so the user
    sees the *actual* error message (not a Python stack trace).
    If 'capture' is False, QIIME's output streams live to the terminal.
    """
    logger.info(f"Running: {' '.join(command)}")
    if dry_run:
        logger.debug("Dry-run mode: command not executed.")
        return None

    start = time.time()
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=capture,
            text=True,
            env=env
        )
    except subprocess.CalledProcessError as e:
        if capture:
            if e.stdout:
                logger.error("QIIME stdout:\n%s", e.stdout.strip())
            if e.stderr:
                logger.error("QIIME stderr:\n%s", e.stderr.strip())
        logger.error("Command failed (exit %s). See QIIME error above.", e.returncode)
        raise SystemExit(e.returncode)

    duration = time.time() - start
    log_success(f"Command completed in {duration:.2f}s")
    if capture and result.stdout:
        logger.debug("QIIME stdout:\n%s", result.stdout.strip())
    return result
