# src/qs/commands/manifest.py
from __future__ import annotations

from pathlib import Path

from qs.utils.logger import get_logger
from qs.utils.manifest import generate_manifest

LOG = get_logger("manifest.cmd")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "manifest", parents=[parent],
        help="Generate a PairedEndFastqManifestPhred33V2 TSV from FASTQs.",
        description=(
            "Writes one row per lane if necessary, reusing the same cleaned SampleID. "
            "QIIME2 will concatenate reads across rows with the same sample-id."
        ),
    )
    p.add_argument("--fastq-dir", type=Path, required=True, help="Directory with FASTQ(.gz) files.")
    p.add_argument("--output", type=Path, required=True, help="Path to write manifest TSV.")
    p.add_argument("--keep-illumina-suffix", action="store_true",
                   help="Keep _S#/_L### tokens in SampleIDs (default is to strip them).")
    p.add_argument("--id-regex", type=str, default=None,
                   help="Optional regex to extract SampleID from filename root; use group 'id' or group 1.")
    p.set_defaults(func=run)


def run(args) -> None:
    path = generate_manifest(
        args.fastq_dir,
        args.output,
        strip_illumina_suffix=not args.keep_illumina_suffix,
        id_regex=args.id_regex,
    )
    LOG.info("Manifest written → %s", path)
    print(f"[ok] manifest → {path}")
