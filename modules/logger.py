# modules/logger.py

import logging
from logging import Logger

def setup_logger() -> Logger:
    logger: Logger = logging.getLogger("qiime_pipeline")
    logger.setLevel(logging.DEBUG)

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