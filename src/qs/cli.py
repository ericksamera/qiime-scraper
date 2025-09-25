# src/qs/cli.py
from __future__ import annotations

import argparse
from qs.utils.logger import setup_logger

from qs.commands import init as cmd_init
from qs.commands import import_reads as cmd_import
from qs.commands import manifest as cmd_manifest
from qs.commands import metadata_validate as cmd_mdval
from qs.commands import metadata_tabulate as cmd_mdtab
from qs.commands import trim_primers as cmd_trim
from qs.commands import denoise_runs as cmd_denoise
from qs.commands import classify_sweep as cmd_cls
from qs.commands import auto_run as cmd_auto


def main() -> None:
    logger = setup_logger()
    parser = argparse.ArgumentParser(
        prog="qs",
        description="QIIME2 Pipeline CLI (init, validate, manifest, import, trim, denoise, classify, auto-run).",
    )

    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")
    parent.add_argument("--show-qiime", dest="show_qiime", action="store_true",
                        help="Stream QIIME output live to console (default).")
    parent.add_argument("--no-show-qiime", dest="show_qiime", action="store_false",
                        help="Capture QIIME output (printed on error).")
    parent.set_defaults(show_qiime=True)

    subparsers = parser.add_subparsers(dest="command", required=True)

    cmd_init.setup_parser(subparsers, parent)
    cmd_mdval.setup_parser(subparsers, parent)
    cmd_mdtab.setup_parser(subparsers, parent)
    cmd_manifest.setup_parser(subparsers, parent)
    cmd_import.setup_parser(subparsers, parent)
    cmd_trim.setup_parser(subparsers, parent)
    cmd_denoise.setup_parser(subparsers, parent)
    cmd_cls.setup_parser(subparsers, parent)
    cmd_auto.setup_parser(subparsers, parent)

    args = parser.parse_args()
    logger.debug("Parsed args: %r", args)
    args.func(args)
