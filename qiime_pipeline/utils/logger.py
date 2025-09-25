import logging

# Attempt to import colorama for colored output; fall back if not installed
try:
    from colorama import Fore, Style
except ImportError:
    class _DummyFore: GREEN = ""
    class _DummyStyle: RESET_ALL = ""
    Fore = _DummyFore()
    Style = _DummyStyle()

# Module-level logger object used everywhere
logger = logging.getLogger("qiime_pipeline")
logger.setLevel(logging.DEBUG)

def log_success(message: str) -> None:
    """Log a success message in green."""
    logger.info(f"{Fore.GREEN}{message}{Style.RESET_ALL}")

def setup_logger() -> logging.Logger:
    """
    Configure logging once:
      - INFO to console
      - DEBUG to file (qiime_pipeline.log)
    Guard against duplicate handlers.
    """
    if getattr(setup_logger, "_configured", False):
        return logger

    logger.propagate = False  # prevent double logging via root

    # Remove any existing handlers only on our named logger
    for h in list(logger.handlers):
        logger.removeHandler(h)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # File handler
    fh = logging.FileHandler("qiime_pipeline.log")
    fh.setLevel(logging.DEBUG)

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    ch.setFormatter(fmt)
    fh.setFormatter(fmt)

    logger.addHandler(ch)
    logger.addHandler(fh)

    setup_logger._configured = True
    return logger
