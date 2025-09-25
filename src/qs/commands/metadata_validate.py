# src/qs/commands/metadata_validate.py
from __future__ import annotations

import sys
from pathlib import Path

from qs.utils.logger import get_logger
from qs.utils.manifest import read_manifest_sample_ids
from qs.metadata.validate import validate_metadata_file, MetadataError
from qs.metadata import DEFAULT_PROTECTED_COLS

LOG = get_logger("metadata.validate")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "metadata-validate", parents=[parent],
        help="Validate a QIIME2-style metadata TSV (IDs, #q2:types, protected columns).",
    )
    p.add_argument("--metadata-file", type=Path, required=True, help="Path to metadata.tsv")
    p.add_argument("--manifest", type=Path, default=None, help="Optional manifest to cross-check sample IDs.")
    p.add_argument(
        "--protected-cols", type=str, default="__f_primer,__r_primer",
        help="Comma-separated protected columns required."
    )
    p.add_argument("--no-types-row", action="store_true", help="Do not require #q2:types row.")
    p.set_defaults(func=run)


def run(args) -> None:
    metadata_path: Path = args.metadata_file
    against_ids = None
    if args.manifest:
        against_ids = read_manifest_sample_ids(args.manifest)

    required = [c.strip() for c in args.protected_cols.split(",") if c.strip()] or list(DEFAULT_PROTECTED_COLS)

    try:
        warnings = validate_metadata_file(
            metadata_path,
            required_protected=required,
            require_types_row=(not args.no_types_row),
            against_sample_ids=against_ids,
        )
    except MetadataError as e:
        LOG.error("Metadata validation failed: %s", e)
        print(f"[fail] {e}", file=sys.stderr)
        sys.exit(1)

    for w in warnings:
        LOG.warning(w)
        print(f"[warn] {w}")
    print("[ok] metadata validated")
