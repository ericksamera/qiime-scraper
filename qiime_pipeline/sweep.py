import logging, json, time, re, shutil
from pathlib import Path
from collections import Counter
from typing import Iterable, Optional, Dict, List, Tuple
from qiime_pipeline.models import ProjectConfig
from qiime_pipeline.analysis import AnalysisRunner
from qiime_pipeline.utils import qiime

logger = logging.getLogger("qiime_pipeline")

class DiversitySweeper:
    """
    Iterates core diversity analysis over subsets of metadata (supports nested and one-level sweeps).
    """
    SLUG_RE = re.compile(r"[^A-Za-z0-9._-]+")  # regex for safe file names

    def __init__(self, config: ProjectConfig):
        self.config = config

    def run_diversity_sweep(self, by: Optional[Iterable[str]] = None, within: Optional[Iterable[str]] = None,
                             min_samples: int = 5, retain_fraction: float = 0.90,
                             beta_cols: Optional[Iterable[str]] = None, time_column: Optional[str] = None,
                             by_only: bool = False, if_exists: str = "skip"):
        """
        Run diversity analysis on all subsets defined by metadata columns.
        If by_only is True, perform a one-level sweep (only 'by' grouping); otherwise nested by → within loop.
        """
        ana_dir = self.config.analysis_dir
        table = ana_dir / "table.qza"
        if not table.exists():
            raise FileNotFoundError(f"Missing {table}. Run 'stage' to prepare analysis data first.")
        # Ensure phylogenetic tree is present
        AnalysisRunner(self.config).build_phylogeny()
        if not self.config.metadata_file or not self.config.metadata_file.exists():
            raise FileNotFoundError("Metadata file is required for diversity-sweep.")
        headers, types_map, rows = self._parse_metadata_table(self.config.metadata_file)
        cols_in_meta = set(headers)
        categorical_cols = self._categorical_columns(headers, types_map)
        # Determine iteration columns
        outer_cols = list(by) if by else list(categorical_cols)
        inner_cols = [] if by_only else (list(within) if within else list(categorical_cols))
        # Validate specified columns exist
        self._require_columns(cols_in_meta, outer_cols, context="--by")
        if not by_only:
            self._require_columns(cols_in_meta, inner_cols, context="--within")
        if time_column:
            self._require_columns(cols_in_meta, [time_column], context="--time-column")
        if beta_cols:
            self._require_columns(cols_in_meta, beta_cols, context="--beta-cols")
        # Record run parameters and environment
        self._write_run_provenance("diversity-sweep", {
            "project_dir": str(self.config.project_dir),
            "metadata_file": str(self.config.metadata_file),
            "by": outer_cols,
            "within": None if by_only else inner_cols,
            "by_only": by_only,
            "min_samples": min_samples,
            "retain_fraction": retain_fraction,
            "beta_cols": list(beta_cols or []),
            "time_column": time_column,
            "if_exists": if_exists,
            "dry_run": self.config.dry_run,
            "show_qiime": self.config.show_qiime
        })
        # Precompute value counts for each categorical column (for global usage)
        value_counts: Dict[str, Counter] = {col: self._value_counts(rows, col) for col in categorical_cols}
        subsets_root = ana_dir / "subsets"
        subsets_root.mkdir(parents=True, exist_ok=True)
        index_records: List[Dict[str, object]] = []
        # Helper to slugify strings for file names
        def slugify(s: str) -> str:
            return DiversitySweeper.SLUG_RE.sub("-", s.strip()).strip("-._")

        if by_only:
            # One-level sweep over outer_cols
            for outer_col in outer_cols:
                for outer_val, count in sorted(value_counts.get(outer_col, Counter()).items()):
                    if count < min_samples:
                        logger.info("[skip] %s=%s has only %d samples (< %d)", outer_col, outer_val, count, min_samples)
                        continue
                    filters = {outer_col: outer_val}
                    sub_dir = subsets_root / slugify(outer_col) / slugify(outer_val)
                    sub_dir.mkdir(parents=True, exist_ok=True)
                    sub_table = sub_dir / "table.qza"
                    target_table = self._resolve_output(sub_table, if_exists)
                    if target_table is None:
                        continue
                    qiime.feature_table_filter_samples(
                        input_table=table,
                        metadata_file=self.config.metadata_file,
                        where=self._filters_to_where(filters),
                        output_table=target_table,
                        dry_run=self.config.dry_run,
                        show_qiime=self.config.show_qiime
                    )
                    sub_qzv = sub_dir / "table.qzv"
                    qiime.feature_table_summarize(
                        input_table=target_table,
                        output=sub_qzv,
                        sample_metadata_file=self.config.metadata_file,
                        dry_run=self.config.dry_run,
                        show_qiime=self.config.show_qiime
                    )
                    depth = AnalysisRunner.choose_sampling_depth(sub_dir, retain_fraction=retain_fraction, min_depth=1000)
                    sub_core = sub_dir / "core-metrics-phylo"
                    # Ensure a tree exists (reuse global tree)
                    if not (ana_dir / "rooted-tree.qza").exists():
                        AnalysisRunner(self.config).build_phylogeny()
                    target_core_dir = sub_core
                    if sub_core.exists():
                        if if_exists == "skip":
                            logger.info("Reusing existing core metrics at: %s", sub_core)
                        elif if_exists == "overwrite":
                            shutil.rmtree(sub_core)
                        elif if_exists == "new":
                            ts = time.strftime("%Y%m%d-%H%M%S")
                            target_core_dir = sub_dir / f"core-metrics-phylo-{ts}"
                            target_core_dir.mkdir(parents=True, exist_ok=True)
                        elif if_exists == "error":
                            raise FileExistsError(f"Subset output exists: {sub_core}")
                        else:
                            raise ValueError("--if-exists must be one of: skip|overwrite|new|error")
                    if not (sub_core.exists() and if_exists == "skip"):
                        qiime.diversity_core_metrics_phylogenetic(
                            input_phylogeny=ana_dir / "rooted-tree.qza",
                            input_table=target_table,
                            sampling_depth=depth,
                            metadata_file=self.config.metadata_file,
                            output_dir=target_core_dir,
                            dry_run=self.config.dry_run,
                            show_qiime=self.config.show_qiime
                        )
                    done_alpha = self._emit_alpha_visuals(target_core_dir, sub_dir)
                    candidate_cols = [c for c in (beta_cols or []) if c != outer_col]
                    ran_cols, ran_dists = self._emit_beta_tests(target_core_dir, sub_dir, filters, candidate_cols)
                    done_emperor = self._emit_emperor(target_core_dir, sub_dir, time_column)
                    subset_record = {
                        "mode": "by-only",
                        "outer": {"column": outer_col, "value": outer_val, "n_samples": count},
                        "inner": None,
                        "where": self._filters_to_where(filters),
                        "sampling_depth": depth,
                        "paths": {"dir": str(sub_dir), "table": str(target_table), "core": str(target_core_dir)},
                        "retain_fraction": retain_fraction,
                        "min_samples": min_samples,
                        "time_column": time_column,
                        "alpha_tests": done_alpha,
                        "beta_tests": {"columns": ran_cols, "dists": ran_dists},
                        "emperor": done_emperor,
                        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")
                    }
                    (sub_dir / "subset.json").write_text(json.dumps(subset_record, indent=2) + "\n")
                    index_records.append(subset_record)
        else:
            # Nested outer (by) and inner (within) loops
            for outer_col in outer_cols:
                for outer_val, count_outer in sorted(value_counts.get(outer_col, Counter()).items()):
                    if count_outer < min_samples:
                        logger.info("[skip] %s=%s has only %d samples (< %d)", outer_col, outer_val, count_outer, min_samples)
                        continue
                    outer_filters = {outer_col: outer_val}
                    for inner_col in inner_cols:
                        for inner_val, count_inner in sorted(value_counts.get(inner_col, Counter()).items()):
                            if count_inner < min_samples:
                                continue
                            filters = {**outer_filters, **{inner_col: inner_val}}
                            sub_dir = subsets_root / slugify(outer_col) / slugify(outer_val) / slugify(inner_col) / slugify(inner_val)
                            sub_dir.mkdir(parents=True, exist_ok=True)
                            sub_table = sub_dir / "table.qza"
                            target_table = self._resolve_output(sub_table, if_exists)
                            if target_table is None:
                                continue
                            qiime.feature_table_filter_samples(
                                input_table=table,
                                metadata_file=self.config.metadata_file,
                                where=self._filters_to_where(filters),
                                output_table=target_table,
                                dry_run=self.config.dry_run,
                                show_qiime=self.config.show_qiime
                            )
                            sub_qzv = sub_dir / "table.qzv"
                            qiime.feature_table_summarize(
                                input_table=target_table,
                                output=sub_qzv,
                                sample_metadata_file=self.config.metadata_file,
                                dry_run=self.config.dry_run,
                                show_qiime=self.config.show_qiime
                            )
                            depth = AnalysisRunner.choose_sampling_depth(sub_dir, retain_fraction=retain_fraction, min_depth=1000)
                            sub_core = sub_dir / "core-metrics-phylo"
                            if not (ana_dir / "rooted-tree.qza").exists():
                                AnalysisRunner(self.config).build_phylogeny()
                            target_core_dir = sub_core
                            if sub_core.exists():
                                if if_exists == "skip":
                                    logger.info("Reusing existing core metrics at: %s", sub_core)
                                elif if_exists == "overwrite":
                                    shutil.rmtree(sub_core)
                                elif if_exists == "new":
                                    ts = time.strftime("%Y%m%d-%H%M%S")
                                    target_core_dir = sub_dir / f"core-metrics-phylo-{ts}"
                                    target_core_dir.mkdir(parents=True, exist_ok=True)
                                elif if_exists == "error":
                                    raise FileExistsError(f"Subset output exists: {sub_core}")
                                else:
                                    raise ValueError("--if-exists must be one of: skip|overwrite|new|error")
                            if not (sub_core.exists() and if_exists == "skip"):
                                qiime.diversity_core_metrics_phylogenetic(
                                    input_phylogeny=ana_dir / "rooted-tree.qza",
                                    input_table=target_table,
                                    sampling_depth=depth,
                                    metadata_file=self.config.metadata_file,
                                    output_dir=target_core_dir,
                                    dry_run=self.config.dry_run,
                                    show_qiime=self.config.show_qiime
                                )
                            done_alpha = self._emit_alpha_visuals(target_core_dir, sub_dir)
                            candidate_cols = beta_cols or []
                            ran_cols, ran_dists = self._emit_beta_tests(target_core_dir, sub_dir, filters, candidate_cols)
                            done_emperor = self._emit_emperor(target_core_dir, sub_dir, time_column)
                            subset_record = {
                                "mode": "nested",
                                "outer": {"column": outer_col, "value": outer_val, "n_samples": count_outer},
                                "inner": {"column": inner_col, "value": inner_val, "n_samples": count_inner},
                                "where": self._filters_to_where(filters),
                                "sampling_depth": depth,
                                "paths": {"dir": str(sub_dir), "table": str(target_table), "core": str(target_core_dir)},
                                "retain_fraction": retain_fraction,
                                "min_samples": min_samples,
                                "time_column": time_column,
                                "alpha_tests": done_alpha,
                                "beta_tests": {"columns": ran_cols, "dists": ran_dists},
                                "emperor": done_emperor,
                                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")
                            }
                            (sub_dir / "subset.json").write_text(json.dumps(subset_record, indent=2) + "\n")
                            index_records.append(subset_record)
        # Write index of all subset runs
        (subsets_root / "index.json").write_text(json.dumps(index_records, indent=2) + "\n")
        logger.info("Diversity sweep completed across %d subset(s).", len(index_records))

    # ---------------- Helper methods ----------------
    def _parse_metadata_table(self, metadata_file: Path):
        lines = [ln.rstrip("\n") for ln in metadata_file.read_text(encoding="utf-8").splitlines() if ln.strip()]
        if not lines:
            raise RuntimeError(f"Empty metadata file: {metadata_file}")
        headers = [h.strip() for h in lines[0].split("\t")]
        types_map = {h: "" for h in headers}
        start = 1
        if len(lines) > 1 and lines[1].startswith("#q2:types"):
            type_vals = [v.strip().lower() for v in lines[1].split("\t")]
            for i, h in enumerate(headers):
                if i < len(type_vals):
                    types_map[h] = type_vals[i]
            start = 2
        rows = []
        for ln in lines[start:]:
            parts = ln.split("\t")
            row = {}
            for i, h in enumerate(headers):
                row[h] = parts[i] if i < len(parts) else ""
            rows.append(row)
        return headers, types_map, rows

    def _categorical_columns(self, headers: List[str], types_map: Dict[str, str]) -> List[str]:
        cats = []
        for h in headers:
            hl = h.lower()
            # Exclude sample id columns from sweeping
            if hl in {"sample-id", "sampleid", "#sampleid"}:
                continue
            t = (types_map.get(h, "") or "").lower()
            # Prefer explicit categorical, otherwise treat as categorical if not an obvious date/time numeric column
            if t == "categorical" or (t == "" and hl not in {"date", "time"}):
                cats.append(h)
        return cats

    def _value_counts(self, rows: List[Dict[str, str]], column: str) -> Counter:
        c = Counter()
        for r in rows:
            v = (r.get(column) or "").strip()
            if v:
                c[v] += 1
        return c

    def _require_columns(self, available_cols: Iterable[str], needed: Iterable[str], context: str):
        missing = [col for col in (needed or []) if col and col not in available_cols]
        if missing:
            raise RuntimeError(f"[metadata check] {context}: missing columns: {', '.join(missing)}. Available columns: {', '.join(sorted(available_cols))}")

    def _filters_to_where(self, filters: Dict[str, str]) -> str:
        def sql_quote(val: str) -> str:
            return "'" + val.replace("'", "''") + "'"
        return " AND ".join(f"\"{col}\" = {sql_quote(val)}" for col, val in filters.items())

    def _resolve_output(self, path: Path, if_exists: str) -> Optional[Path]:
        """Decide output path based on if_exists policy (return None if skipping)."""
        if path.exists():
            if if_exists == "skip":
                logger.info("Exists, skipping: %s", path)
                return None
            if if_exists == "overwrite":
                try:
                    path.unlink()
                except Exception:
                    pass
            elif if_exists == "new":
                ts = time.strftime("%Y%m%d-%H%M%S")
                path = path.with_name(f"{path.stem}-{ts}{path.suffix}")
            elif if_exists == "error":
                raise FileExistsError(f"Output exists: {path}")
            else:
                raise ValueError("--if-exists must be one of: skip|overwrite|new|error")
        return path

    def _emit_alpha_visuals(self, core_dir: Path, sub_dir: Path) -> List[str]:
        metrics = []
        alpha_vectors = {
            "faith-pd": core_dir / "faith_pd_vector.qza",
            "evenness": core_dir / "evenness_vector.qza",
            "shannon": core_dir / "shannon_vector.qza",
            "observed-features": core_dir / "observed_features_vector.qza"
        }
        for key, vec_path in alpha_vectors.items():
            if not vec_path.exists():
                continue
            out_vis = sub_dir / f"{key}-group-significance.qzv"
            if self._resolve_output(out_vis, if_exists="skip") is None:
                continue
            qiime.diversity_alpha_group_significance(
                input_alpha_vector=vec_path,
                metadata_file=self.config.metadata_file, # type: ignore
                output_visualization=out_vis,
                dry_run=self.config.dry_run,
                show_qiime=self.config.show_qiime
            )
            metrics.append(key)
        return metrics

    def _emit_beta_tests(self, core_dir: Path, sub_dir: Path, filters: Dict[str, str], candidate_cols: Iterable[str]) -> Tuple[List[str], List[str]]:
        ran_cols = []
        ran_dists = set()
        beta_mats = {
            "unweighted-unifrac": core_dir / "unweighted_unifrac_distance_matrix.qza",
            "weighted-unifrac": core_dir / "weighted_unifrac_distance_matrix.qza",
            "bray-curtis": core_dir / "bray_curtis_distance_matrix.qza",
            "jaccard": core_dir / "jaccard_distance_matrix.qza"
        }
        cols_in_meta = []
        if self.config.metadata_file:
            try:
                header = self.config.metadata_file.read_text().splitlines()[0]
                cols_in_meta = [h.strip() for h in header.split("\t")]
            except Exception:
                cols_in_meta = []
        for col in candidate_cols or []:
            if col in filters:  # don't test on the grouping column itself
                continue
            if col not in cols_in_meta:
                continue
            for dist_label, dist_path in beta_mats.items():
                if not dist_path.exists():
                    continue
                out_vis = sub_dir / f"{dist_label}-{DiversitySweeper.SLUG_RE.sub('-', col)}-group-significance.qzv"
                if self._resolve_output(out_vis, if_exists="skip") is None:
                    continue
                qiime.diversity_beta_group_significance(
                    input_distance=dist_path,
                    metadata_file=self.config.metadata_file, # type: ignore
                    metadata_column=col,
                    output_visualization=out_vis,
                    pairwise=True,
                    dry_run=self.config.dry_run,
                    show_qiime=self.config.show_qiime
                )
                ran_dists.add(dist_label)
            ran_cols.append(col)
        return ran_cols, sorted(ran_dists)


    def _emit_emperor(self, core_dir: Path, sub_dir: Path, time_col: Optional[str]) -> List[str]:
        emitted = []
        if not (time_col and self.config.metadata_file):
            return emitted
        pcoa_plots = {
            "unweighted-unifrac": core_dir / "unweighted_unifrac_pcoa_results.qza",
            "bray-curtis": core_dir / "bray_curtis_pcoa_results.qza"
        }
        for label, pcoa_path in pcoa_plots.items():
            if not pcoa_path.exists():
                continue
            out_vis = sub_dir / f"{label}-emperor-{DiversitySweeper.SLUG_RE.sub('-', time_col)}.qzv"
            if self._resolve_output(out_vis, if_exists="skip") is None:
                continue
            qiime.emperor_plot(
                input_pcoa=pcoa_path,
                metadata_file=self.config.metadata_file,
                output_visualization=out_vis,
                custom_axes=time_col,
                dry_run=self.config.dry_run,
                show_qiime=self.config.show_qiime
            )
            emitted.append(label)
        return emitted

    def _write_run_provenance(self, command: str, params: dict):
        """
        Write a RUN.json in analysis/ recording the command, parameters, and environment.
        """
        ana_dir = self.config.analysis_dir
        record = {
            "command": command,
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "args": params,
            "env": {
                "project_dir": str(self.config.project_dir),
                "python": None,
                "platform": None,
                "cwd": str(Path.cwd()),
                "qiime_info": None
            },
            "winners_env": {}
        }
        try:
            import sys, platform, subprocess
            record["env"]["python"] = sys.version
            record["env"]["platform"] = platform.platform()
            result = subprocess.run(["qiime", "info"], check=True, capture_output=True, text=True)
            record["env"]["qiime_info"] = result.stdout.strip()
        except Exception:
            record["env"]["qiime_info"] = None
        winners_file = ana_dir / "WINNERS.env"
        if winners_file.exists():
            try:
                lines = winners_file.read_text().splitlines()
                record["winners_env"] = {k: v for (k, v) in (ln.split("=", 1) for ln in lines if "=" in ln)}
            except Exception:
                record["winners_env"] = {}
        (ana_dir / "RUN.json").write_text(json.dumps(record, indent=2) + "\n")
