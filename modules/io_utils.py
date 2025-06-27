# modules/io_utils.py

import subprocess
from pathlib import Path
from logging import Logger
from typing import List, Optional

# Set by main
logger = Optional[Logger]

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
    return subprocess.run(command, check=True, capture_output=capture)

# ---