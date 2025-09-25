import logging, json, gzip, zipfile
from pathlib import Path
from qiime_pipeline.utils import qiime

logger = logging.getLogger("qiime_pipeline")

class TrimmingOptimizer:
    """Finds optimal truncation lengths for DADA2 by trying multiple trim parameters."""
    def __init__(self, config):
        self.config = config

    def find_optimal_trimming(self):
        fastq_dir = Path(self.config.fastq_dir)
        read_length = None
        for pattern in ["*_R1*.fastq.gz", "*.fastq.gz"]:
            try:
                fq_path = next(fastq_dir.glob(pattern))
                with gzip.open(fq_path, 'rt') as f:
                    f.readline()
                    seq = f.readline().strip()
                    read_length = len(seq)
                    break
            except StopIteration:
                continue
            except Exception:
                break
        if read_length is None:
            read_length = 150

        if read_length <= 150:
            candidates = [0, max(0, read_length - 20), max(0, read_length - 50)]
        else:
            candidates = [0, max(0, read_length - 50), max(0, read_length - 100)]
        candidates = sorted(set(candidates), reverse=True)
        f_candidates = candidates
        r_candidates = candidates

        trim_dir = Path(self.config.trim_output_dir)
        trim_dir.mkdir(parents=True, exist_ok=True)
        best_fraction = -1.0
        best_pair = None
        best_stats = (0, 0)

        for f_len in f_candidates:
            for r_len in r_candidates:
                if f_len < r_len:
                    continue
                table_qza = trim_dir / f"{f_len}-{r_len}_output_table.qza"
                rep_qza   = trim_dir / f"{f_len}-{r_len}_output_rep_seqs.qza"
                stats_qza = trim_dir / f"{f_len}-{r_len}_output_stats.qza"
                if not (table_qza.exists() and rep_qza.exists() and stats_qza.exists()):
                    logger.info(f"[trim] Testing truncation lengths F={f_len}, R={r_len}")
                    qiime.dada2_denoise_paired(
                        input_seqs=self.config.demux_artifact,
                        trunc_len_f=f_len,
                        trunc_len_r=r_len,
                        output_table=table_qza,
                        output_rep_seqs=rep_qza,
                        output_denoising_stats=stats_qza,
                        dry_run=self.config.dry_run,
                        show_qiime=self.config.show_qiime
                    )
                total_input, total_merged = self._parse_stats(stats_qza)
                if total_input == 0:
                    continue
                fraction = total_merged / total_input
                logger.info(f"[trim] F={f_len}, R={r_len} → {fraction:.1%} reads merged ( {total_merged}/{total_input} )")
                if fraction > best_fraction or (abs(fraction - best_fraction) < 1e-9 and best_pair and f_len + r_len > sum(best_pair)):
                    best_fraction = fraction
                    best_pair = (f_len, r_len)
                    best_stats = (total_input, total_merged)

        if not best_pair:
            raise RuntimeError("Optimal trimming failed: no successful DADA2 runs.")
        best_f, best_r = best_pair
        result = {
            "trunc_len_f": best_f,
            "trunc_len_r": best_r,
            "merged_fraction": best_stats[1] / best_stats[0] if best_stats[0] else 0
        }
        with open(trim_dir / "optimal_trimming.json", 'w') as fh:
            json.dump(result, fh, indent=2)
        self.config.best_trim_lengths = (best_f, best_r)
        self.config.best_table_artifact = trim_dir / f"{best_f}-{best_r}_output_table.qza"
        self.config.best_rep_seqs_artifact = trim_dir / f"{best_f}-{best_r}_output_rep_seqs.qza"
        logger.info(f"Optimal truncation lengths found: trunc_len_f={best_f}, trunc_len_r={best_r}")
        return best_pair

    def _parse_stats(self, stats_qza: Path):
        if not stats_qza.exists():
            return (0, 0)
        try:
            with zipfile.ZipFile(stats_qza) as z:
                stats_file = next(p for p in z.namelist() if p.endswith('.tsv') or p.endswith('.txt'))
                data = z.read(stats_file).decode('utf-8')
        except Exception as e:
            logger.error(f"Could not read {stats_qza}: {e}")
            return (0, 0)
        total_input = total_merged = 0
        input_idx = merged_idx = None
        for line in data.splitlines():
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if 'sample-id' in cols[0].lower():
                try:
                    input_idx = cols.index('input')
                    merged_idx = cols.index('merged')
                except ValueError:
                    input_idx, merged_idx = 1, 4
                continue
            try:
                i_idx = input_idx or 1
                m_idx = merged_idx or 4
                total_input += int(cols[i_idx])
                total_merged += int(cols[m_idx])
            except Exception:
                continue
        return (total_input, total_merged)
