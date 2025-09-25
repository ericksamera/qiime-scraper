# src/qs/commands/import_reads.py
from __future__ import annotations

from pathlib import Path

from qs.utils.logger import get_logger
from qs.utils.manifest import generate_manifest
from qs.qiime import commands as qiime

LOG = get_logger("import")


def setup_parser(subparsers, parent) -> None:
    p = subparsers.add_parser(
        "import", parents=[parent],
        help="Import PairedEnd reads using a QIIME2 PairedEndFastqManifestPhred33V2.",
        description=(
            "Create a PairedEnd manifest from FASTQs (with cleaned SampleIDs) "
            "and import to a .qza artifact."
        ),
    )
    p.add_argument("--fastq-dir", type=Path, required=True, help="Directory with FASTQ(.gz) files.")
    p.add_argument("--project-dir", type=Path, required=True, help="Project directory for outputs.")
    p.add_argument("--manifest", type=Path, default=None,
                   help="Optional path to write manifest TSV (defaults to <project-dir>/fastq.manifest).")
    p.add_argument("--output-artifact", type=Path, default=None,
                   help="Optional .qza output path (defaults to <project-dir>/output.qza).")
    p.add_argument("--keep-illumina-suffix", action="store_true",
                   help="Keep _S#/_L### tokens in SampleIDs (default is to strip them).")
    p.add_argument("--id-regex", type=str, default=None,
                   help="Optional regex to extract SampleID from filename root; use group 'id' or group 1.")
    p.set_defaults(func=run)


def run(args) -> None:
    fastq_dir: Path = args.fastq_dir
    project_dir: Path = args.project_dir
    project_dir.mkdir(parents=True, exist_ok=True)

    manifest_path: Path = args.manifest or (project_dir / "fastq.manifest")
    output_qza: Path = args.output_artifact or (project_dir / "output.qza")

    LOG.info("Generating manifest from: %s", fastq_dir)
    generate_manifest(
        fastq_dir,
        manifest_path,
        strip_illumina_suffix=not args.keep_illumina_suffix,
        id_regex=args.id_regex,
    )

    LOG.info("Importing to artifact: %s", output_qza)
    qiime.import_data(
        input_path=manifest_path,
        output_path=output_qza,
        import_type="SampleData[PairedEndSequencesWithQuality]",
        input_format="PairedEndFastqManifestPhred33V2",
        dry_run=getattr(args, "dry_run", False),
        show_stdout=getattr(args, "show_qiime", True),
    )

    LOG.info("Import complete → %s", output_qza)
    print(f"[ok] Imported reads → {output_qza}")
