import re
from pathlib import Path
from qiime_pipeline.utils import logger, qiime

class Preprocessor:
    """Handles primer removal (Cutadapt trimming) before denoising."""
    def __init__(self, config):
        self.config = config
        self.primers_per_sample = {}
        meta_path = getattr(config, 'metadata_file', None)
        if not (meta_path and Path(meta_path).is_file()):
            if meta_path:
                logger.warning(f"Metadata file not found: {meta_path}. Skipping primer trimming.")
            return
        with open(meta_path, 'r') as fh:
            headers = None
            for line in fh:
                if line.strip() and not line.startswith('#'):
                    headers = line.strip().split('\t')
                    break
            if not headers:
                return
            # Accept a few common column names
            names_f = [c for c in ("primer_f", "forward_primer", "f_primer") if c in headers]
            names_r = [c for c in ("primer_r", "reverse_primer", "r_primer") if c in headers]
            if not (names_f and names_r):
                logger.warning("Metadata provided but missing primer columns. Expected one of primer_f/forward_primer and primer_r/reverse_primer.")
                return
            pf_idx = headers.index(names_f[0]); pr_idx = headers.index(names_r[0])
            sample_idx = 0  # assume first column is sample ID
            fh.seek(0)
            for line in fh:
                if line.strip() and not line.startswith('#'):
                    cols = line.rstrip('\n').split('\t')
                    sid = cols[sample_idx].lstrip('#')
                    self.primers_per_sample[sid] = {
                        'primer_f': cols[pf_idx].strip(),
                        'primer_r': cols[pr_idx].strip()
                    }

    def trim_paired(self, input_artifact: Path, output_artifact: Path,
                    forward_primers: list[str], reverse_primers: list[str], cores: int = 0):
        """Trim given forward/reverse primer sequences from the input artifact using QIIME2 Cutadapt."""
        if not forward_primers or not reverse_primers:
            logger.warning("No primers provided for trimming; skipping.")
            return
        logger.info(f"Trimming primers from {input_artifact} -> {output_artifact}")
        qiime.cutadapt_trim_paired(
            input_seqs=input_artifact,
            forward_primers=forward_primers,
            reverse_primers=reverse_primers,
            output_path=output_artifact,
            cores=cores,
            dry_run=getattr(self.config, 'dry_run', False),
            show_qiime=getattr(self.config, 'show_qiime', False)
        )
        logger.info(f"Primer trimming complete: {output_artifact}")

    def trim_if_needed(self, input_artifact: Path, run_dir: Path) -> Path:
        """Trim primers from the demux artifact if primers are specified for this run. Returns the (possibly new) artifact path."""
        if not self.primers_per_sample:
            return input_artifact  # No primer info; nothing to do
        # Collect all sample IDs in this run's FASTQs
        sample_ids = set()
        for fq in run_dir.glob("*.fastq.gz"):
            m = re.match(r"^(?P<sid>.+)_R[12](?:_[0-9]{3})?\.fastq\.gz$", fq.name)
            if m:
                sample_ids.add(m.group('sid'))
        # Gather unique primers present in this run
        fwd_primers = set(); rev_primers = set()
        for sid in sample_ids:
            if sid in self.primers_per_sample:
                pf = self.primers_per_sample[sid]['primer_f']
                pr = self.primers_per_sample[sid]['primer_r']
                if pf: fwd_primers.add(pf)
                if pr: rev_primers.add(pr)
        if fwd_primers and rev_primers:
            # Determine output artifact name (use run-specific name if multi-run)
            run_name = run_dir.name
            output_artifact = Path(self.config.project_dir) / (
                f"{run_name}_trimmed.qza" if run_dir.resolve() != Path(self.config.fastq_dir).resolve()
                else "trimmed.qza"
            )
            self.trim_paired(input_artifact, output_artifact,
                             forward_primers=list(fwd_primers), reverse_primers=list(rev_primers))
            return output_artifact
        else:
            logger.info(f"No primer sequences to trim for run '{run_dir.name}'. Skipping primer trimming.")
            return input_artifact
