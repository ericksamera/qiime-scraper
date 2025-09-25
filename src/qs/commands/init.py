# src/qs/commands/init.py
from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Optional

from qs.utils.logger import get_logger
from qs.utils.samples import discover_fastqs, paired_sample_ids
from qs.metadata import DEFAULT_PROTECTED_COLS
from qs.metadata.template import write_metadata_template

LOG = get_logger("init")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "init", parents=[parent],
        help="Scan FASTQs and write a QIIME2-style metadata TSV with protected primer columns.",
        description=(
            "Generate a metadata TSV with '#SampleID' and protected columns "
            "prefixed by '__' (e.g., '__f_primer', '__r_primer'). "
            "By default, SampleIDs drop Illumina tokens like _S1/_L001."
        ),
    )
    p.add_argument("--fastq-dir", type=Path, required=True, help="Directory containing raw FASTQ(.gz) files.")
    p.add_argument("--output-file", type=Path, required=True, help="Path to write metadata TSV (e.g., ./project/metadata.tsv).")
    p.add_argument("--protected-cols", type=str, default="__f_primer,__r_primer",
                   help="Comma-separated protected columns to include (default: __f_primer,__r_primer).")
    p.add_argument("--keep-illumina-suffix", action="store_true",
                   help="Keep _S#/_L### tokens in SampleIDs (default is to strip them).")
    p.add_argument("--id-regex", type=str, default=None,
                   help="Optional regex to extract SampleID from filename root; use group 'id' or group 1.")
    p.add_argument("--force", action="store_true", help="Overwrite output file if it already exists.")
    p.set_defaults(func=run)


def _parse_protected(s: str) -> List[str]:
    cols = [c.strip() for c in s.split(",") if c.strip()]
    return cols or list(DEFAULT_PROTECTED_COLS)


def run(args) -> None:
    fastq_dir: Path = args.fastq_dir
    output_file: Path = args.output_file
    protected_cols: List[str] = _parse_protected(args.protected_cols)
    id_regex: Optional[str] = args.id_regex
    strip_illumina_suffix: bool = not args.keep_illumina_suffix

    if output_file.exists() and not args.force:
        LOG.error("Refusing to overwrite existing file: %s (use --force)", output_file)
        print(f"error: {output_file} exists (use --force)", file=sys.stderr)
        sys.exit(1)

    fastqs = discover_fastqs(fastq_dir)
    if not fastqs:
        LOG.error("No FASTQ files found under %s", fastq_dir)
        print(f"error: no FASTQ files found under {fastq_dir}", file=sys.stderr)
        sys.exit(3)

    sample_ids = paired_sample_ids(
        fastqs,
        strip_illumina_suffix=strip_illumina_suffix,
        id_regex=id_regex,
    )
    if not sample_ids:
        LOG.error("Found FASTQs but no R1/R2 pairs; check naming conventions.")
        print("error: no paired R1/R2 FASTQs detected; check filenames.", file=sys.stderr)
        sys.exit(4)

    LOG.info("Discovered %d paired samples.", len(sample_ids))
    write_metadata_template(sample_ids, output_file, protected_cols)
    LOG.info("Metadata written → %s", output_file)
    print(f"[ok] Wrote metadata with {len(sample_ids)} samples: {output_file}")
