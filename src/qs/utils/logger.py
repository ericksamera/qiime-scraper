# src/qs/utils/logger.py
from __future__ import annotations

import logging
from typing import Optional

_LOGGER_NAME = "qs"


def get_logger(name: Optional[str] = None) -> logging.Logger:
    if name:
        return logging.getLogger(f"{_LOGGER_NAME}.{name}")
    return logging.getLogger(_LOGGER_NAME)


def setup_logger() -> logging.Logger:
    """
    Configure a root 'qs' logger:
      - INFO to console
      - DEBUG to file (qs.log)
    Idempotent: safe to call multiple times.
    """
    logger = logging.getLogger(_LOGGER_NAME)
    if getattr(setup_logger, "_configured", False):
        return logger

    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    # Remove any pre-existing handlers (only for our logger)
    for h in list(logger.handlers):
        logger.removeHandler(h)

    # Console handler (INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # File handler (DEBUG)
    fh = logging.FileHandler("qs.log", encoding="utf-8")
    fh.setLevel(logging.DEBUG)

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    ch.setFormatter(fmt)
    fh.setFormatter(fmt)

    logger.addHandler(ch)
    logger.addHandler(fh)

    setup_logger._configured = True  # type: ignore[attr-defined]
    return logger
