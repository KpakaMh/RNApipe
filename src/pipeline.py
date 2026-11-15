from typing import Optional
import scanpy as sc
from src.config import PipelineConfig
from src.preprocessing import QualityControlProcessor
from src.clustering import DimensionalityReducer
from src.annotation import CellAnnotator
from src.enrichment import FunctionalEnrichment


class ScRNASeqPipeline:
    """Orchestrateur orienté objet exécutant l'ensemble des étapes d'analyse scRNA-seq."""

    def __init__(
        self,
        config_path: str = "config/config.yaml",
        markers_path: str = "config/markers.yaml",
    ):
        self.config = PipelineConfig.from_yaml(config_path, markers_path)
        self.adata: Optional[sc.AnnData] = None

        self.qc_processor = QualityControlProcessor(self.config)
        self.dim_reducer = DimensionalityReducer(self.config)
        self.annotator = CellAnnotator(self.config)
        self.enricher = FunctionalEnrichment(self.config)

    def run(self, adata: Optional[sc.AnnData] = None) -> sc.AnnData:
        self.adata = self.qc_processor.run(adata)
        self.adata = self.dim_reducer.run(self.adata)
        self.adata = self.annotator.run(self.adata)
        self.adata = self.enricher.run(self.adata)
        return self.adata