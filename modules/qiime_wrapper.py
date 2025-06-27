# modules/qiime_wrapper.py

from pathlib import Path
from .io_utils import run_command
import logging

logger = logging.getLogger("qiime_pipeline")


def import_data(
        fastq_dir: Path,
        output_path: Path,
        dry_run=False) -> None:
    
    """
    """

    logger.info("Import completed.")

    run_command([
        "qiime", "tools", "import",
        "--type", "SampleData[PairedEndSequencesWithQuality]",
        "--input-path", str(fastq_dir / "fastq.manifest"),
        "--input-format", "PairedEndFastqManifestPhred33V2",
        "--output-path", str(output_path)
    ], dry_run)

    logger.info("Import completed.")


# ---