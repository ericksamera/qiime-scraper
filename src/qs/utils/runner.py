# src/qs/utils/runner.py
from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Mapping, Optional, Sequence

from qs.utils.logger import get_logger

LOG = get_logger("runner")


def run_command(
    cmd: Sequence[str],
    *,
    dry_run: bool = False,
    capture: bool = False,
    cwd: Optional[Path] = None,
    env: Optional[Mapping[str, str]] = None,
) -> subprocess.CompletedProcess[str]:
    """
    Execute a subprocess with unified logging and error handling.

    - Logs the exact command line.
    - Respects dry_run (no execution).
    - capture=False streams output; capture=True buffers output.
    - Merges provided env with the current process environment (preserves PATH).
    - Raises CalledProcessError on failure (after logging stdout/stderr).
    """
    LOG.info("Running: %s", " ".join(cmd))
    if dry_run:
        LOG.debug("[dry-run] command not executed")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    # Merge env with current environment so PATH and friends are preserved
    env_dict = os.environ.copy()
    if env:
        env_dict.update({str(k): str(v) for k, v in env.items()})

    try:
        result = subprocess.run(
            list(cmd),
            check=True,
            cwd=str(cwd) if cwd else None,
            env=env_dict,
            text=True,
            capture_output=capture,
        )
    except FileNotFoundError as e:
        # Typically means the executable (e.g., 'qiime') is not on PATH
        LOG.error("Executable not found: %s (PATH=%s)", cmd[0], env_dict.get("PATH", ""))
        raise
    except subprocess.CalledProcessError as e:
        if e.stdout:
            LOG.error("STDOUT:\n%s", e.stdout.strip())
        if e.stderr:
            LOG.error("STDERR:\n%s", e.stderr.strip())
        LOG.error("Command failed with exit code %s", e.returncode)
        raise

    LOG.info("Command completed successfully.")
    if capture and result.stdout:
        LOG.debug("Captured STDOUT:\n%s", result.stdout.strip())
    return result
