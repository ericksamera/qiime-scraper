# modules/logger.py

import logging
from logging import Logger

from colorama import Fore, Style

logger = logging.getLogger("qiime_pipeline")
logger.setLevel(logging.DEBUG)


def log_success(message: str) -> None:
    logger.info(f"{Fore.GREEN}{message}{Style.RESET_ALL}")

def setup_logger() -> Logger:

    ch: logging.StreamHandler = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    fh: logging.FileHandler = logging.FileHandler("qiime_pipeline.log")
    fh.setLevel(logging.DEBUG)

    formatter: logging.Formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)

    return logger

# ---