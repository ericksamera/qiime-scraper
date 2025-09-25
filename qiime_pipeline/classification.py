import logging, json, zipfile, statistics, heapq
from pathlib import Path
from typing import Optional, Sequence
from qiime_pipeline.models import ProjectConfig
from qiime_pipeline.utils import qiime

logger = logging.getLogger("qiime_pipeline")

class Classifier:
    """Runs taxonomy classification using provided classifiers and selects the best."""
    def __init__(self, config: ProjectConfig):
        self.config = config
        self.rank_keys: Sequence[str] = ("pct_depth≥7", "median_conf", "mean_conf")

    def _choose_reads(self, input_reads: Optional[Path]) -> Path:
        """Robustly determine FeatureData[Sequence] to classify."""
        if input_reads:
            return input_reads
        if getattr(self.config, "input_reads", None):
            return self.config.input_reads
        # Prefer best trimming output if present
        for p in [
            getattr(self.config, "best_rep_seqs_artifact", None),
            # plain rep-seqs.qza from pipeline run
            (self.config.project_dir / "rep-seqs.qza"),
        ]:
            if p and Path(p).exists():
                return Path(p)
        # Try reading via ProjectConfig helper (may fail if trimming not run)
        try:
            return self.config.get_best_rep_seqs()
        except FileNotFoundError as e:
            raise FileNotFoundError(
                "No representative sequences available for classification. "
                "Expected either work/optimal_trimming/*_rep_seqs.qza or project_dir/rep-seqs.qza."
            ) from e

    def find_optimal_classifier(self, input_reads: Optional[Path] = None):
        reads = self._choose_reads(input_reads)
        classifiers_dir = self.config.classifiers_dir
        if classifiers_dir is None or not classifiers_dir.exists():
            raise FileNotFoundError(f"Classifiers directory not found: {classifiers_dir}")
        output_dir = self.config.classifier_output_dir
        logger.info("Running classification for each classifier in %s", classifiers_dir)

        metrics_list = []
        for classifier_path in sorted(classifiers_dir.glob("*.qza")):
            tag = classifier_path.stem
            output_taxonomy = output_dir / f"{tag}_classification.qza"
            if not output_taxonomy.exists():
                qiime.classify_sklearn(
                    input_reads=reads,
                    input_classifier=classifier_path,
                    output_classification=output_taxonomy,
                    dry_run=self.config.dry_run,
                    show_qiime=self.config.show_qiime
                )
            metrics = self._extract_classification_metrics(output_taxonomy)
            if metrics:
                metrics_list.append((classifier_path, metrics))
            else:
                logger.warning("No metrics extracted for %s; skipping.", classifier_path.name)

        if not metrics_list:
            raise RuntimeError("No valid classification metrics were produced.")

        top_classifiers = {}
        for key in self.rank_keys:
            top_two = heapq.nlargest(2, metrics_list, key=lambda item: item[1].get(key, 0))
            top_classifiers[key] = [path for path, _ in top_two]
            for rank, (cls_path, cls_metrics) in enumerate(top_two, start=1):
                logger.info("Top %d for %s: %s (value=%s)", rank, key, cls_path.name, cls_metrics.get(key))

        report = {
            k: [
                {"classifier": p.name, "value": next(m for c, m in metrics_list if c == p).get(k)}
                for p in top_classifiers[k]
            ] for k in top_classifiers
        }
        with open(output_dir / "optimal_classifiers.json", "w") as fh:
            json.dump(report, fh, indent=2)
        logger.info("Optimal classifiers report saved.")
        return top_classifiers

    def _extract_classification_metrics(self, taxonomy_qza: Path, min_depth: int = 5):
        if not taxonomy_qza.exists():
            return None
        try:
            with zipfile.ZipFile(taxonomy_qza) as z:
                taxonomy_tsv = z.read(next(p for p in z.namelist() if p.endswith("taxonomy.tsv"))).decode()
        except (KeyError, StopIteration):
            return None
        lines = taxonomy_tsv.splitlines()
        if not lines:
            return None
        depths, confidences, full_depth_count = [], [], 0
        for line in lines[1:]:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                continue
            taxon, confidence = parts[1], parts[2]
            try:
                conf_val = float(confidence)
            except ValueError:
                continue
            depth = len(taxon.split(";")) if taxon else 0
            depths.append(depth); confidences.append(conf_val)
            if depth >= 7:
                full_depth_count += 1
        if not confidences:
            return None
        return {
            "n_records": len(confidences),
            "mean_conf": statistics.mean(confidences),
            "median_conf": statistics.median(confidences),
            "mean_depth": statistics.mean(depths) if depths else 0,
            "n_depth≥7": full_depth_count,
            "pct_depth≥7": full_depth_count / len(confidences),
            f"n_depth≥{min_depth}": sum(1 for d in depths if d >= min_depth)
        }
