# modules/io_utils.py

import subprocess
from pathlib import Path
from logging import Logger
from typing import List, Optional
import logging
import time

import csv

from .logger import log_success

logger = logging.getLogger("qiime_pipeline")

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

def generate_manifest(
        fastq_path: Path
        ) -> None:
    """
    Generate a QIIME 2 manifest file for paired-end FASTQ files.
    """
    base_names = {f.stem.split('.R')[0] for f in fastq_path.glob('*.fastq.gz')}
    manifest_path = fastq_path.joinpath("fastq.manifest")
    header = ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"]

    with manifest_path.open(mode='w', newline='') as output_file:
        writer = csv.DictWriter(output_file, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for name in base_names:
            writer.writerow({
                "sample-id": name,
                "forward-absolute-filepath": str(fastq_path.joinpath(f"{name}.R1.fastq.gz")),
                "reverse-absolute-filepath": str(fastq_path.joinpath(f"{name}.R2.fastq.gz"))
            })

# ---