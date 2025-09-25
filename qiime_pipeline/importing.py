import logging
from pathlib import Path
from qiime_pipeline.utils import io, qiime

logger = logging.getLogger("qiime_pipeline")

class Importer:
    """Handles importing raw FASTQ reads into a QIIME2 artifact."""
    def __init__(self, config):
        self.config = config

    def run(self):
        """Generate a manifest and import FASTQs into a .qza artifact."""
        fastq_dir = Path(self.config.fastq_dir)
        project_dir = Path(self.config.project_dir)
        manifest_path = io.generate_manifest(fastq_dir, project_dir)
        output_path = project_dir / "output.qza"
        logger.info(f"Importing reads from {manifest_path} into {output_path}")
        qiime.import_data(
            input_path=manifest_path,
            output_path=output_path,
            import_type="SampleData[PairedEndSequencesWithQuality]",
            input_format="PairedEndFastqManifestPhred33V2",
            dry_run=self.config.dry_run,
            show_qiime=self.config.show_qiime
        )
        logger.info(f"Imported data to {output_path}")
        self.config.demux_artifact = output_path
