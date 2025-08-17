import argparse
import json
from pathlib import Path

from modules import io_utils, iterators, logger, qiime_wrapper


def parse_args():
    p = argparse.ArgumentParser(prog="pipeline")
    sp = p.add_subparsers(dest="cmd", required=True)

    # shared opts for steps that still touch FASTQs
    def add_fastq_ops(sub):
        sub.add_argument("--fastq-dir", type=Path, required=True, help="Project FASTQ dir")
        sub.add_argument(
            "--project-dir",
            type=Path,
            required=True,
            help="Named output directory to create/use (pipeline writes here)",
        )
        sub.add_argument("--dry-run", action="store_true", help="Print commands without executing")

    # full run
    run = sp.add_parser("run", help="Run full pipeline")
    add_fastq_ops(run)
    run.add_argument("--classifiers-dir", type=Path, default=Path("/home/erick/qiime"))

    # import-only
    imp = sp.add_parser("import", help="Only import FASTQs to a QIIME artifact")
    add_fastq_ops(imp)

    # trimming-only
    trm = sp.add_parser("trim", help="Only search optimal trimming")
    add_fastq_ops(trm)

    # classify-only (NO fastq-dir required)
    cls = sp.add_parser("classify", help="Run classifier sweep on rep_seqs in <project-dir>")
    cls.add_argument(
        "--project-dir",
        type=Path,
        required=True,
        help="Directory containing work/optimal_trimming and rep_seqs artifact",
    )
    cls.add_argument(
        "--input-reads",
        type=Path,
        default=None,
        help="Optional FeatureData[Sequence] .qza to classify (overrides best rep_seqs)",
    )
    cls.add_argument("--classifiers-dir", type=Path, default=Path("/home/erick/qiime"))
    cls.add_argument("--dry-run", action="store_true")

    return p.parse_args()


def _paths(project_fastq_dir: Path, project_dir: Path):
    work = project_dir / "work"
    imported_qza = project_dir / "output.qza"
    trim_dir = work / "optimal_trimming"
    cls_dir = work / "optimal_classifier"
    return work, imported_qza, trim_dir, cls_dir


def cmd_import(project_fastq_dir: Path, project_dir: Path, dry_run: bool):
    io_utils.generate_manifest(project_fastq_dir, project_dir)
    qiime_wrapper.import_data(
        input_path=project_dir / "fastq.manifest",
        output_path=project_dir / "output.qza",
        dry_run=dry_run,
    )


def cmd_trim(imported_qza: Path):
    # Expected to return (trunc_len_f, trunc_len_r) from your iterator
    return iterators.get_optimal_trimming(imported_qza=imported_qza)


def _load_best_trim(trim_dir: Path) -> tuple[int, int]:
    meta = json.loads((trim_dir / "optimal_trimming.json").read_text())
    return int(meta["trunc_len_f"]), int(meta["trunc_len_r"])


def _best_rep_seqs(project_dir: Path) -> Path:
    trim_dir = project_dir / "work" / "optimal_trimming"
    f, r = _load_best_trim(trim_dir)
    rep = trim_dir / f"{f}-{r}_output_rep_seqs.qza"
    if not rep.exists():
        raise FileNotFoundError(f"rep_seqs not found: {rep}")
    return rep


def cmd_classify(project_dir: Path, classifiers_dir: Path, input_reads: Path | None = None):
    rep_seqs = input_reads if input_reads else _best_rep_seqs(project_dir)
    return iterators.get_optimal_classifier(imported_qza=rep_seqs, classifiers_dir=classifiers_dir)


def main():
    args = parse_args()
    logger.setup_logger()

    # Steps that still need FASTQs
    if args.cmd in {"run", "import", "trim"}:
        project_fastq_dir = args.fastq_dir.resolve()
        project_dir = args.project_dir.resolve()
        project_dir.mkdir(parents=True, exist_ok=True)
        (project_dir / "work").mkdir(parents=True, exist_ok=True)
        work, imported_qza, trim_dir, _ = _paths(project_fastq_dir, project_dir)

        if args.cmd == "import":
            cmd_import(project_fastq_dir, project_dir, args.dry_run)
            return

        if args.cmd == "trim":
            cmd_trim(imported_qza)
            return

        if args.cmd == "run":
            cmd_import(project_fastq_dir, project_dir, args.dry_run)
            cmd_trim(imported_qza)
            cmd_classify(project_dir, args.classifiers_dir)
            return

    # Classify does NOT require --fastq-dir
    if args.cmd == "classify":
        project_dir = args.project_dir.resolve()
        cmd_classify(project_dir, args.classifiers_dir, args.input_reads)
        return


if __name__ == "__main__":
    main()
