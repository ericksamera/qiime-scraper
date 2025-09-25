# src/qs/commands/trim_primers.py
from __future__ import annotations

import sys
from pathlib import Path

from qs.utils.logger import get_logger
from qs.qiime import commands as qiime
from qs.metadata.read import load_metadata_table, collect_unique_primers

LOG = get_logger("trim")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "trim-primers", parents=[parent],
        help="Trim primers with cutadapt using primer columns from metadata.",
        description=(
            "Reads forward/reverse primer sequences from the metadata TSV "
            "(default columns: __f_primer, __r_primer), aggregates unique sequences, "
            "and runs 'qiime cutadapt trim-paired'."
        ),
    )
    p.add_argument("--project-dir", type=Path, required=True,
                   help="Project directory (defaults input/output paths).")
    p.add_argument("--metadata-file", type=Path, required=True,
                   help="QIIME-style metadata TSV (must include primer columns).")
    p.add_argument("--input-artifact", type=Path, default=None,
                   help="Input .qza (defaults to <project-dir>/output.qza).")
    p.add_argument("--output-artifact", type=Path, default=None,
                   help="Output .qza (defaults to <project-dir>/trimmed.qza).")
    p.add_argument("--front-f-col", type=str, default="__f_primer",
                   help="Metadata column for forward primers (default: __f_primer).")
    p.add_argument("--front-r-col", type=str, default="__r_primer",
                   help="Metadata column for reverse primers (default: __r_primer).")

    # cutadapt knobs
    p.add_argument("--cores", type=int, default=0, help="Number of cores for cutadapt (0=auto).")
    p.add_argument("--discard-untrimmed", action="store_true",
                   help="Discard reads without primer matches (--p-discard-untrimmed).")
    p.add_argument("--no-indels", action="store_true",
                   help="Disallow indels in matches (--p-no-indels).")

    # strictness
    p.add_argument("--require-any-primer", action="store_true",
                   help="Fail if no primers are found in metadata (default).")
    p.add_argument("--allow-empty-primers", dest="require_any_primer", action="store_false",
                   help="Allow empty primer sets (skip trimming).")

    p.set_defaults(func=run, require_any_primer=True)


def run(args) -> None:
    project_dir: Path = args.project_dir
    meta_path: Path = args.metadata_file
    input_qza: Path = args.input_artifact or (project_dir / "output.qza")
    output_qza: Path = args.output_artifact or (project_dir / "trimmed.qza")

    # Load metadata and collect primers
    header, rows = load_metadata_table(meta_path)
    if args.front_f_col not in header or args.front_r_col not in header:
        LOG.error("Metadata missing primer columns: %s / %s", args.front_f_col, args.front_r_col)
        print(f"error: metadata missing primer columns '{args.front_f_col}' / '{args.front_r_col}'", file=sys.stderr)
        sys.exit(2)

    fwd, rev = collect_unique_primers(rows, args.front_f_col, args.front_r_col)
    LOG.info("Unique primers: F=%d, R=%d", len(fwd), len(rev))
    if fwd:
        LOG.info("Forward: %s", ", ".join(fwd))
    if rev:
        LOG.info("Reverse: %s", ", ".join(rev))

    if (not fwd or not rev) and args.require_any_primer:
        LOG.error("No primers found in metadata; aborting (use --allow-empty-primers to skip).")
        print("error: no primers found; use --allow-empty-primers to proceed without trimming.", file=sys.stderr)
        sys.exit(3)

    # If empty primers and allowed, just copy path (no trimming); but here we'll simply no-op
    if not fwd or not rev:
        LOG.warning("Empty primer set; skipping trimming. Input will be used as-is: %s", input_qza)
        print(f"[skip] no primers → skipping trimming; input remains: {input_qza}")
        return

    # Ensure paths
    project_dir.mkdir(parents=True, exist_ok=True)
    if not input_qza.exists():
        LOG.error("Input artifact not found: %s", input_qza)
        print(f"error: input artifact not found: {input_qza}", file=sys.stderr)
        sys.exit(4)

    LOG.info("Trimming primers with cutadapt → %s", output_qza)
    qiime.cutadapt_trim_paired(
        input_seqs=input_qza,
        forward_primers=fwd,
        reverse_primers=rev,
        output_path=output_qza,
        cores=args.cores,
        discard_untrimmed=args.discard_untrimmed,
        no_indels=args.no_indels,
        dry_run=getattr(args, "dry_run", False),
        show_stdout=getattr(args, "show_qiime", True),
    )
    LOG.info("Primer trimming complete → %s", output_qza)
    print(f"[ok] trimmed → {output_qza}")
