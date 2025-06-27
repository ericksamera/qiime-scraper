# main.py

import argparse
from pathlib import Path
from modules import logger
from modules import io_utils

from modules.qiime_wrapper import import_data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing")

    args = parser.parse_args()

    log = logger.setup_logger()

    log.info("asdasdasd")
    

if __name__ == "__main__":
    main()
# ---