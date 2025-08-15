# modules/io_utils.py

from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path
from typing import Iterable, Sequence

# ---------------------------------------------------------------------
# Manifest generation
# ---------------------------------------------------------------------

# Matches Illumina-style filenames:
#   <root>_R1_001.fastq.gz  /  <root>_R2_001.fastq.gz
# Also accepts: <root>_R1.fastq.gz / <root>_R2.fastq.gz
_ILLUMINA_RE = re.compile(r"^(?P<root>.+)_R(?P<read>[12])(?:_[0-9]{3})?\.fastq\.gz$")


def generate_manifest(project_fastq_dir: Path) -> Path:
    """
    Create a PairedEndFastqManifestPhred33V2 CSV at <dir>/fastq.manifest with absolute paths.

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

    manifest_path = fastq_dir / "fastq.manifest"
    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sample-id",
                "forward-absolute-filepath",
                "reverse-absolute-filepath",
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    return manifest_path


# ---------------------------------------------------------------------
# Subprocess helper (kept minimal + useful defaults)
# ---------------------------------------------------------------------

def run_command(
    command: Sequence[str] | str,
    *,
    check: bool = True,
    capture: bool = False,
) -> str | subprocess.CompletedProcess:
    """
    Run a shell command. If `capture` is True, return stdout text; otherwise
    return the CompletedProcess. Raises CalledProcessError on failure if `check` is True.
    """
    result = subprocess.run(
        command,
        check=check,
        capture_output=capture,
        text=capture,
    )
    return result.stdout if capture else result
