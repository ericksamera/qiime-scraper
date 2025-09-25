import csv, re, subprocess, logging, time
from pathlib import Path
from typing import List, Optional

from qiime_pipeline.utils.logger import log_success

logger = logging.getLogger("qiime_pipeline")

# Regex to match Illumina FASTQ filenames (e.g., sample_R1_001.fastq.gz)
_ILLUMINA_RE = re.compile(r"^(?P<root>.+)_R(?P<read>[12])(?:_[0-9]{3})?\.fastq\.gz$")

def _collect_pairs(paths) -> List[dict]:
    """Collect forward/reverse FASTQ file paths into sample pairs."""
    pairs = {}
    for fq in paths:
        m = _ILLUMINA_RE.match(fq.name)
        if not m:
            continue
        root = m.group("root")
        read = m.group("read")
        pairs.setdefault(root, {})[read] = str(fq.resolve())
    rows = []
    for root, d in sorted(pairs.items()):
        if "1" in d and "2" in d:
            rows.append({
                "sample-id": root,
                "forward-absolute-filepath": d["1"],
                "reverse-absolute-filepath": d["2"]
            })
    return rows

def generate_manifest(project_fastq_dir: Path, manifest_dir: Path) -> Path:
    """
    Generate a PairedEndFastqManifestPhred33V2 TSV (fastq.manifest) in manifest_dir using FASTQs in project_fastq_dir.
    """
    rows = _collect_pairs(project_fastq_dir.glob("*.fastq.gz"))
    if not rows:
        logger.info("No FASTQ pairs at top level; searching subdirectories...")
        rows = _collect_pairs(project_fastq_dir.rglob("*.fastq.gz"))
    if not rows:
        raise RuntimeError(f"No R1/R2 FASTQ pairs found in {project_fastq_dir}")
    manifest_path = manifest_dir / "fastq.manifest"
    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Manifest written: %s", manifest_path)
    return manifest_path

def run_command(command: List[str], dry_run: bool = False, capture: bool = False, env: Optional[dict] = None):
    """
    Execute a shell command. If dry_run, log the command without executing.
    If capture=True, capture the output (for internal use), else stream it.
    Raises SystemExit on failure, logging QIIME's output for clarity.
    """
    logger.info("Running: %s", " ".join(map(str, command)))
    if dry_run:
        logger.debug("Dry-run mode: command not executed.")
        return None
    start = time.time()
    try:
        result = subprocess.run(command, check=True, capture_output=capture, text=True, env=env)
    except subprocess.CalledProcessError as e:
        if capture:
            if e.stdout:
                logger.error("QIIME stdout:\n%s", e.stdout.strip())
            if e.stderr:
                logger.error("QIIME stderr:\n%s", e.stderr.strip())
        logger.error("Command failed (exit code %s). See QIIME error above.", e.returncode)
        raise SystemExit(e.returncode)
    duration = time.time() - start
    log_success(f"Command completed in {duration:.2f}s")
    if capture and result.stdout:
        logger.debug("QIIME stdout:\n%s", result.stdout.strip())
    return result
