# main.py

import argparse
from pathlib import Path
from modules import logger
from modules import io_utils
from modules import qiime_wrapper

from modules import iterators

PROJECT_BASE = Path("/mnt/c/Users/erick/Documents/metagenomics-test")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing")

    args = parser.parse_args()

    log = logger.setup_logger()

    io_utils.generate_manifest(PROJECT_BASE.joinpath('fastq'))

    # qiime_wrapper.import_data(
    #     input_path=PROJECT_BASE / 'fastq' / 'fastq.manifest',
    #     output_path=PROJECT_BASE / 'output.qza'
    # )

    iterators.get_optimal_trimming(imported_qza=PROJECT_BASE / 'output.qza')

    log.info("asdasdasd")
    

if __name__ == "__main__":
    main()
# ---