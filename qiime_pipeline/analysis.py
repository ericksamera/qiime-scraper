import logging, shutil, json, zipfile, re, time
from pathlib import Path
from typing import Optional, Iterable, Dict
from qiime_pipeline.models import ProjectConfig
from qiime_pipeline.utils import qiime

logger = logging.getLogger("qiime_pipeline")

class AnalysisRunner:
    """
    Handles downstream analysis: staging best results, building phylogeny, core diversity metrics, rarefaction, and taxa barplots.
    """
    def __init__(self, config: ProjectConfig):
        self.config = config

    def stage_winners(self, preferred_metric: str = "pct_depth≥7"):
        """
        Link/copy the winning (or best available) rep-seqs, table, and taxonomy into analysis/.
        Fallbacks:
        - If optimal_trimming results are missing, use project_dir/rep-seqs.qza and table.qza
        - If optimal_classifiers.json missing, pick the first *_classification.qza in classifier_output_dir (if any)
        """
        ana_dir = self.config.analysis_dir

        # Try to use optimal trimming artifacts
        rep_src = tbl_src = tbl_viz_src = None
        f = r = "NA"
        try:
            f, r = self.config.get_best_trim()
            rep_src = self.config.trim_output_dir / f"{f}-{r}_output_rep_seqs.qza"
            tbl_src = self.config.trim_output_dir / f"{f}-{r}_output_table.qza"
            tbl_viz_src = self.config.trim_output_dir / f"{f}-{r}_output_table.qzv"
            if not (rep_src.exists() and tbl_src.exists()):
                raise FileNotFoundError
        except Exception:
            # Fallback: project-level outputs
            rep_src = self.config.project_dir / "rep-seqs.qza"
            tbl_src = self.config.project_dir / "table.qza"
            tbl_viz_src = self.config.project_dir / "table.qzv"
            if not rep_src.exists() or not tbl_src.exists():
                raise FileNotFoundError(
                    "Could not locate representative sequences and table to stage. "
                    "Expected either work/optimal_trimming/* or project_dir/{rep-seqs.qza,table.qza}."
                )

        # Determine best classifier result
        tax_src = None
        opt_json = self.config.classifier_output_dir / "optimal_classifiers.json"
        if opt_json.exists():
            data = json.loads(opt_json.read_text())
            def _normalize_key(key: str) -> str:
                return key.strip().lower().replace(" ", "").replace(">=", "≥")
            pref_key = _normalize_key(preferred_metric)
            winner_name = None
            for key in (pref_key, "pct_depth≥7", "median_conf", "mean_conf"):
                if key in data and data[key]:
                    winner_name = data[key][0]["classifier"]
                    break
            if winner_name:
                winner_stem = Path(winner_name).stem
                cand = self.config.classifier_output_dir / f"{winner_stem}_classification.qza"
                if cand.exists():
                    tax_src = cand
        else:
            # Fallback: take the first available classification
            for cand in sorted(self.config.classifier_output_dir.glob("*_classification.qza")):
                tax_src = cand
                break

        # Stage artifacts
        rep_dest = ana_dir / "rep-seqs.qza"
        tbl_dest = ana_dir / "table.qza"
        tax_dest = ana_dir / "taxonomy.qza"
        self._stage(rep_src, rep_dest)
        self._stage(tbl_src, tbl_dest)
        if tax_src:
            self._stage(tax_src, tax_dest)
        else:
            logger.info("No taxonomy classification found to stage.")

        # Carry summary & (re)build if missing
        if tbl_viz_src and tbl_viz_src.exists():
            self._stage(tbl_viz_src, ana_dir / "table.qzv")
        if self.config.metadata_file:
            qiime.feature_table_summarize(
                input_table=tbl_dest,
                output=ana_dir / "table-with-metadata.qzv",
                sample_metadata_file=self.config.metadata_file,
                dry_run=self.config.dry_run,
                show_qiime=self.config.show_qiime
            )
        qiime.feature_table_tabulate_seqs(
            input_data=rep_dest,
            output=ana_dir / "rep-seqs.qzv",
            dry_run=self.config.dry_run,
            show_qiime=self.config.show_qiime
        )
        (ana_dir / "WINNERS.env").write_text(f"WIN_TRUNC_F={f}\nWIN_TRUNC_R={r}\n")
        logger.info("Winners staged to analysis/: %s, %s%s",
                    rep_dest.name, tbl_dest.name, (", " + tax_dest.name) if tax_src else "")
        return {"rep_seqs": rep_dest, "table": tbl_dest, "taxonomy": (tax_dest if tax_src else None), "analysis_dir": ana_dir}

    def build_phylogeny(self) -> Dict[str, Path]:
        """
        Align sequences and construct a phylogenetic tree (saved as rooted-tree.qza in analysis/).
        """
        ana_dir = self.config.analysis_dir
        rooted_tree = ana_dir / "rooted-tree.qza"
        if rooted_tree.exists():
            logger.info("Phylogenetic tree already exists: %s", rooted_tree)
            return {"rooted_tree": rooted_tree}
        rep_seqs = ana_dir / "rep-seqs.qza"
        if not rep_seqs.exists():
            raise FileNotFoundError(f"Missing {rep_seqs}. Run the 'stage' step first.")
        qiime.phylogeny_align_to_tree_mafft_fasttree(
            input_sequences=rep_seqs,
            output_alignment=ana_dir / "aligned-rep-seqs.qza",
            output_masked_alignment=ana_dir / "masked-aligned-rep-seqs.qza",
            output_tree=ana_dir / "unrooted-tree.qza",
            output_rooted_tree=rooted_tree,
            dry_run=self.config.dry_run,
            show_qiime=self.config.show_qiime
        )
        logger.info("Phylogenetic tree built: %s", rooted_tree)
        return {"rooted_tree": rooted_tree}

    def run_core_diversity(self, sampling_depth: Optional[int] = None, beta_cols: Iterable[str] = ("body-site", "subject"),
                            time_column: Optional[str] = None, include_alpha_tests: bool = True,
                            include_beta_tests: bool = True, include_emperor: bool = True, if_exists: str = "skip") -> Dict[str, Path]:
        """
        Run core phylogenetic diversity analysis and optional group significance tests and Emperor plots.
        """
        ana_dir = self.config.analysis_dir
        table = ana_dir / "table.qza"
        tree = ana_dir / "rooted-tree.qza"
        if not table.exists():
            raise FileNotFoundError(f"Missing {table}. Run 'stage' to prepare analysis data.")
        if not tree.exists():
            logger.info("No rooted-tree.qza found; building phylogeny now...")
            self.build_phylogeny()
        if sampling_depth is None:
            sampling_depth = AnalysisRunner.choose_sampling_depth(ana_dir)
        sampling_depth = int(sampling_depth)
        logger.info("Using sampling depth: %d", sampling_depth)
        core_dir = ana_dir / "core-metrics-phylo"
        # Handle existing output directory according to if_exists policy
        if core_dir.exists():
            if if_exists == "skip":
                # Reuse if outputs seem complete
                needed = ["unweighted_unifrac_distance_matrix.qza", "faith_pd_vector.qza", "evenness_vector.qza"]
                if all((core_dir / n).exists() for n in needed):
                    logger.info("Reusing existing core metrics in %s", core_dir)
                else:
                    raise RuntimeError(f"{core_dir} exists but is incomplete. Rerun with --if-exists overwrite or new.")
            elif if_exists == "overwrite":
                shutil.rmtree(core_dir)
                logger.info("Overwriting existing core metrics at %s", core_dir)
            elif if_exists == "new":
                timestamp = time.strftime("%Y%m%d-%H%M%S")
                core_dir = ana_dir / f"core-metrics-phylo-{timestamp}"
                logger.info("Writing core metrics to new directory: %s", core_dir)
            elif if_exists == "error":
                raise FileExistsError(f"Output directory exists: {core_dir}")
            else:
                raise ValueError("--if-exists must be one of skip|overwrite|new|error")
        # Run core metrics (unless skipping existing completed outputs)
        if not core_dir.exists() or if_exists != "skip":
            qiime.diversity_core_metrics_phylogenetic(
                input_phylogeny=tree,
                input_table=table,
                sampling_depth=sampling_depth,
                metadata_file=self.config.metadata_file,
                output_dir=core_dir,
                dry_run=self.config.dry_run,
                show_qiime=self.config.show_qiime
            )
        # Determine metadata presence
        if self.config.metadata_file is None:
            if include_alpha_tests or include_beta_tests or include_emperor:
                logger.info("No metadata provided; skipping alpha/beta significance tests and Emperor plots.")
            include_alpha_tests = include_beta_tests = include_emperor = False
        # Alpha diversity group significance tests
        if include_alpha_tests:
            alpha_vectors = {
                "faith-pd": core_dir / "faith_pd_vector.qza",
                "evenness": core_dir / "evenness_vector.qza",
                "shannon": core_dir / "shannon_vector.qza",
                "observed-features": core_dir / "observed_features_vector.qza"
            }
            for key, vec_path in alpha_vectors.items():
                if vec_path.exists():
                    qiime.diversity_alpha_group_significance(
                        input_alpha_vector=vec_path,
                        metadata_file=self.config.metadata_file,
                        output_visualization=ana_dir / f"{key}-group-significance.qzv",
                        dry_run=self.config.dry_run,
                        show_qiime=self.config.show_qiime
                    )
                else:
                    logger.info("Skipping alpha group significance for %s (file not found)", vec_path.name)
        # Beta diversity group significance tests
        if include_beta_tests:
            beta_mats = {
                "unweighted-unifrac": core_dir / "unweighted_unifrac_distance_matrix.qza",
                "weighted-unifrac": core_dir / "weighted_unifrac_distance_matrix.qza",
                "bray-curtis": core_dir / "bray_curtis_distance_matrix.qza",
                "jaccard": core_dir / "jaccard_distance_matrix.qza"
            }
            # Determine metadata columns present
            cols_in_meta = []
            if self.config.metadata_file:
                try:
                    with open(self.config.metadata_file, 'r') as fh:
                        header_line = None
                        for line in fh:
                            if not line.strip():
                                continue
                            if line.startswith('#'):
                                first_cell = line.strip().split("\t")[0]
                                if first_cell in ("#SampleID", "#Sample ID", "#OTUID", "#OTU ID"):
                                    header_line = line.strip()
                                    break
                                else:
                                    continue
                            else:
                                header_line = line.strip()
                                break
                        if header_line:
                            cols_in_meta = [h.strip() for h in header_line.split("\t")]
                except Exception:
                    cols_in_meta = []
            for dist_label, dist_path in beta_mats.items():
                if not dist_path.exists():
                    logger.info("Skipping beta group significance: missing %s", dist_path.name)
                    continue
                for col in beta_cols or []:
                    if self.config.metadata_file and col in cols_in_meta:
                        qiime.diversity_beta_group_significance(
                            input_distance=dist_path,
                            metadata_file=self.config.metadata_file,
                            metadata_column=col,
                            output_visualization=ana_dir / f"{dist_label}-{col}-group-significance.qzv",
                            pairwise=True,
                            dry_run=self.config.dry_run,
                            show_qiime=self.config.show_qiime
                        )
                    else:
                        logger.info("Skipping beta group significance for %s: %s", col,
                                    ("column not in metadata." if self.config.metadata_file else "no metadata provided"))
        # Emperor plots for time series if requested
        if include_emperor and time_column and self.config.metadata_file:
            pcoa_files = {
                "unweighted-unifrac": core_dir / "unweighted_unifrac_pcoa_results.qza",
                "bray-curtis": core_dir / "bray_curtis_pcoa_results.qza"
            }
            for label, pcoa_path in pcoa_files.items():
                if pcoa_path.exists():
                    qiime.emperor_plot(
                        input_pcoa=pcoa_path,
                        metadata_file=self.config.metadata_file,
                        output_visualization=ana_dir / f"{label}-emperor-{time_column}.qzv",
                        custom_axes=time_column,
                        dry_run=self.config.dry_run,
                        show_qiime=self.config.show_qiime
                    )
        return {"core_dir": core_dir, "sampling_depth": sampling_depth}

    def run_alpha_rarefaction(self, max_depth: Optional[int] = None, sampling_depth_hint: Optional[int] = None) -> Path:
        """
        Generate an alpha rarefaction visualization for the dataset.
        """
        ana_dir = self.config.analysis_dir
        if max_depth is None:
            if sampling_depth_hint is None:
                sampling_depth_hint = AnalysisRunner.choose_sampling_depth(ana_dir)
            max_depth = int(sampling_depth_hint * 2)
        output_vis = ana_dir / "alpha-rarefaction.qzv"
        qiime.diversity_alpha_rarefaction(
            input_table=ana_dir / "table.qza",
            input_phylogeny=ana_dir / "rooted-tree.qza",
            max_depth=max_depth,
            metadata_file=self.config.metadata_file,
            output_visualization=output_vis,
            dry_run=self.config.dry_run,
            show_qiime=self.config.show_qiime
        )
        logger.info("Alpha rarefaction completed: %s", output_vis)
        return output_vis

    def run_taxa_barplots(self) -> Path:
        """Create taxonomy bar plots visualization."""
        ana_dir = self.config.analysis_dir
        output_vis = ana_dir / "taxa-bar-plots.qzv"
        qiime.taxa_barplot(
            input_table=ana_dir / "table.qza",
            input_taxonomy=ana_dir / "taxonomy.qza",
            metadata_file=self.config.metadata_file,
            output_visualization=output_vis,
            dry_run=self.config.dry_run,
            show_qiime=self.config.show_qiime
        )
        logger.info("Taxonomy bar plots generated: %s", output_vis)
        return output_vis

    def run_full_analysis(self, preferred_metric: str = "pct_depth≥7", sampling_depth: Optional[int] = None,
                          beta_cols: str = "body-site,subject", time_column: Optional[str] = None,
                          include_taxa_barplots: bool = True):
        """
        Execute the complete downstream analysis: stage winners, build phylogeny, core diversity, rarefaction, (optional) barplots.
        """
        self.stage_winners(preferred_metric=preferred_metric)
        self.build_phylogeny()
        beta_list = [c.strip() for c in beta_cols.split(",") if c.strip()]
        diversity_result = self.run_core_diversity(sampling_depth=sampling_depth, beta_cols=beta_list, time_column=time_column)
        depth = diversity_result.get("sampling_depth")
        self.run_alpha_rarefaction(sampling_depth_hint=depth)
        if include_taxa_barplots:
            self.run_taxa_barplots()
        logger.info("Downstream analysis completed in %s", self.config.analysis_dir)
        return {"analysis_dir": self.config.analysis_dir}

    @staticmethod
    def choose_sampling_depth(analysis_dir: Path, retain_fraction: float = 0.90, min_depth: int = 1000) -> int:
        """
        Choose a rarefaction sampling depth that retains ~retain_fraction of samples (default 90%).
        """
        qzv = analysis_dir / "table-with-metadata.qzv"
        if not qzv.exists():
            qzv = analysis_dir / "table.qzv"
        if not qzv.exists():
            logger.warning("No table summary found; using min_depth=%d as default sampling depth.", min_depth)
            return min_depth
        try:
            with zipfile.ZipFile(qzv) as z:
                html = z.read(next(p for p in z.namelist() if p.endswith("sample-frequency-detail.html"))).decode()
        except Exception:
            return min_depth
        m = re.search(r'id=["\\\']table-data["\\\'][^>]*>\s*(\{.*?\})\s*<', html)
        if not m:
            return min_depth
        try:
            data = json.loads(m.group(1))
            freqs = sorted(map(int, data.get("Frequency", {}).values()))
            if not freqs:
                return min_depth
        except Exception:
            return min_depth
        k = int((1.0 - retain_fraction) * len(freqs))
        k = max(0, min(k, len(freqs) - 1))
        return max(min_depth, freqs[k])

    def _stage(self, src: Path, dest: Path):
        """Symlink or copy a file into the analysis directory."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        try:
            dest.symlink_to(src)
            how = "symlinked"
        except Exception:
            shutil.copy2(src, dest)
            how = "copied"
        logger.info("%s -> %s", how.capitalize(), dest)
