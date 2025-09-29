# src/qs/commands/doctor.py
from __future__ import annotations

import shutil
import sys
from pathlib import Path

from qs.utils.logger import get_logger
from qs.utils.runner import run_command

LOG = get_logger("doctor")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "doctor", parents=[parent],
        help="Preflight checks for QIIME2, plugins, and environment.",
    )
    p.set_defaults(func=run)


def _ok(x: bool) -> str:
    return "OK" if x else "MISSING"


def run(_args) -> None:
    # 1) qiime presence
    qiime_path = shutil.which("qiime")
    print(f"[check] qiime on PATH: {_ok(bool(qiime_path))} ({qiime_path or 'not found'})")
    if not qiime_path:
        print("error: 'qiime' executable not found on PATH.", file=sys.stderr)
        sys.exit(2)

    # 2) version
    try:
        out = run_command(["qiime", "--version"], capture=True).stdout.strip()
        print(f"[check] qiime --version: {out or 'unknown'}")
    except Exception as e:
        print(f"error: failed to run 'qiime --version': {e}", file=sys.stderr)
        sys.exit(3)

    # 3) core plugins
    plugins = ["cutadapt", "dada2", "feature-classifier", "phylogeny", "diversity", "taxa", "metadata", "tools", "feature-table"]
    bad = []
    for pl in plugins:
        try:
            run_command(["qiime", pl, "--help"], capture=True)
            print(f"[check] plugin '{pl}': OK")
        except Exception:
            print(f"[check] plugin '{pl}': MISSING")
            bad.append(pl)

    if bad:
        print("error: missing plugins: " + ", ".join(bad), file=sys.stderr)
        sys.exit(4)

    print("[ok] environment looks good.")
