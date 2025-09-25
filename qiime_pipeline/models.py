import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, List

import logging
logger = logging.getLogger("qiime_pipeline")

@dataclass
class ProjectConfig:
    fastq_dir: Optional[Path] = None
    project_dir: Path = field(default_factory=Path)
    classifiers_dir: Optional[Path] = None
    input_reads: Optional[Path] = None
    metadata_file: Optional[Path] = None
    dry_run: bool = False
    show_qiime: bool = True
    no_trim_primers: bool = False
    # Internal fields (initialized in __post_init__)
    work_dir: Path = field(init=False)
    demux_artifact: Path = field(init=False)       # Imported .qza path
    trim_output_dir: Path = field(init=False)      # Directory for trimming outputs
    classifier_output_dir: Path = field(init=False) # Directory for classifier outputs
    best_trim_lengths: Optional[Tuple[int, int]] = field(default=None, init=False)
    best_rep_seqs_artifact: Optional[Path] = field(default=None, init=False)
    best_table_artifact: Optional[Path] = field(default=None, init=False)

    def __post_init__(self):
        # Ensure project directories exist
        self.project_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir = self.project_dir / "work"
        self.work_dir.mkdir(parents=True, exist_ok=True)
        # Define standard artifact paths
        self.demux_artifact = self.project_dir / "output.qza"
        self.trim_output_dir = self.work_dir / "optimal_trimming"
        self.trim_output_dir.mkdir(parents=True, exist_ok=True)
        self.classifier_output_dir = self.work_dir / "optimal_classifier"
        self.classifier_output_dir.mkdir(parents=True, exist_ok=True)
        # Resolve optional paths
        if self.classifiers_dir:
            self.classifiers_dir = self.classifiers_dir.resolve()
        if self.fastq_dir:
            self.fastq_dir = self.fastq_dir.resolve()
        if self.metadata_file:
            self.metadata_file = self.metadata_file.resolve()

    @property
    def analysis_dir(self) -> Path:
        """Return (and create if needed) the analysis/ directory for this project."""
        ana = self.project_dir / "analysis"
        ana.mkdir(parents=True, exist_ok=True)
        return ana

    def get_best_trim(self) -> Tuple[int, int]:
        """Load the optimal truncation lengths (F, R) from the trimming results JSON."""
        optimal_json = self.trim_output_dir / "optimal_trimming.json"
        if not optimal_json.exists():
            raise FileNotFoundError(f"Missing {optimal_json}. Run the trimming stage first.")
        data = json.loads(optimal_json.read_text())
        return int(data.get("trunc_len_f")), int(data.get("trunc_len_r"))

    def get_best_rep_seqs(self) -> Path:
        """Get path to the representative sequences artifact of the best trimming run."""
        f, r = self.get_best_trim()
        rep_path = self.trim_output_dir / f"{f}-{r}_output_rep_seqs.qza"
        if not rep_path.exists():
            raise FileNotFoundError(f"Representative sequences not found: {rep_path}. Run the denoising stage first.")
        return rep_path

    def get_best_table(self) -> Path:
        """Get path to the feature table artifact of the best trimming run."""
        f, r = self.get_best_trim()
        table_path = self.trim_output_dir / f"{f}-{r}_output_table.qza"
        if not table_path.exists():
            raise FileNotFoundError(f"Feature table not found: {table_path}. Run the denoising stage first.")
        return table_path

@dataclass
class RunArtifact:
    path: Path
    artifact_type: Optional[str] = None  # e.g., QIIME2 semantic type
    name: Optional[str] = None          # optional human-friendly name

@dataclass
class SampleMetadata:
    path: Path
    def load(self) -> str:
        """Return the raw content of the metadata file."""
        return self.path.read_text(encoding="utf-8")
