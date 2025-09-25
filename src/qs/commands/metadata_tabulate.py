# src/qs/commands/metadata_tabulate.py
from __future__ import annotations

from pathlib import Path

from qs.utils.logger import get_logger
from qs.qiime import commands as qiime

LOG = get_logger("metadata.tabulate")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "metadata-tabulate", parents=[parent],
        help="Create a QIIME2 metadata visualization (.qzv) from a metadata TSV.",
    )
    p.add_argument("--metadata-file", type=Path, required=True, help="Path to metadata.tsv")
    p.add_argument("--output", type=Path, required=True, help="Path to write .qzv")
    p.set_defaults(func=run)


def run(args) -> None:
    qiime.metadata_tabulate(
        input_file=args.metadata_file,
        output_visualization=args.output,
        dry_run=getattr(args, "dry_run", False),
        show_stdout=getattr(args, "show_qiime", True),
    )
    LOG.info("Metadata tabulation complete → %s", args.output)
    print(f"[ok] metadata tabulated → {args.output}")
