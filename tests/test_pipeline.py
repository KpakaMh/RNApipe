from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from src.config import PipelineConfig
from src.preprocessing import QualityControlProcessor
from src.clustering import DimensionalityReducer
from src.annotation import CellAnnotator
from src.pipeline import ScRNASeqPipeline


@pytest.fixture
def synthetic_adata() -> sc.AnnData:
    """
    Génère un mini-jeu de données scRNA-seq synthétique (60 cellules x 50 gènes)
    comprenant 2 populations distinctes et des gènes mitochondriaux pour tester le QC.
    """
    rng = np.random.default_rng(seed=42)
    n_cells = 60
    n_genes = 50

    counts = rng.poisson(lam=3, size=(n_cells, n_genes)).astype(np.float32)

    gene_names = [
        "KRT5", "TP63", "KRT15",
        "FOXJ1", "TPPP3", "PIFO",
        "MT-CO1", "MT-ND1",
        "RPS12", "RPL10",
    ] + [f"GENE_{i}" for i in range(10, n_genes)]

    counts[:30, 0:3] += 25
    counts[30:55, 3:6] += 25

    counts[55:, :] = 1
    counts[55:, 6:8] = 80

    obs = pd.DataFrame(index=[f"CELL_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=gene_names)

    return sc.AnnData(X=counts, obs=obs, var=var)


@pytest.fixture
def test_config(tmp_path: Path) -> PipelineConfig:
    """Crée une instance PipelineConfig isolée pointant vers un répertoire temporaire."""
    return PipelineConfig(
        project_name="RNApipe-Test",
        tissue="synthetic_epithelium",
        paths={
            "raw_data": str(tmp_path / "data" / "raw.h5ad"),
            "qc_data": str(tmp_path / "data" / "01_qc_norm.h5ad"),
            "clustered_data": str(tmp_path / "data" / "02_clustered.h5ad"),
            "annotated_data": str(tmp_path / "data" / "03_annotated.h5ad"),
            "figures_dir": str(tmp_path / "results" / "figures"),
            "tables_dir": str(tmp_path / "results" / "tables"),
        },
        qc={
            "min_genes": 5,
            "max_genes": 500,
            "min_cells": 2,
            "max_pct_mt": 20.0,
        },
        preprocessing={
            "target_sum": 1000,
            "n_top_genes": 30,
        },
        dimred_clustering={
            "n_pcs": 10,
            "n_neighbors": 8,
            "leiden_resolution": 0.4,
        },
        markers={
            "Basal_Cells": ["KRT5", "TP63", "KRT15"],
            "Ciliated_Cells": ["FOXJ1", "TPPP3", "PIFO"],
        },
    )


def test_pipeline_config_properties(test_config: PipelineConfig):
    """Vérifie la création automatique des répertoires de sortie via PipelineConfig."""
    assert test_config.figures_dir.exists()
    assert test_config.tables_dir.exists()
    assert "Basal_Cells" in test_config.markers


def test_quality_control_processor(synthetic_adata: sc.AnnData, test_config: PipelineConfig):
    """Vérifie le filtrage des cellules dégradées (%MT), la normalisation et les HVG."""
    processor = QualityControlProcessor(test_config)
    adata_qc = processor.run(synthetic_adata.copy())

    assert adata_qc.n_obs == 55
    assert "counts" in adata_qc.layers
    assert "highly_variable" in adata_qc.var.columns
    assert adata_qc.raw is not None
    assert Path(test_config.paths["qc_data"]).exists()


def test_dimensionality_reducer(synthetic_adata: sc.AnnData, test_config: PipelineConfig):
    """Vérifie le calcul de la PCA, de l'UMAP et du clustering Leiden."""
    adata_qc = QualityControlProcessor(test_config).run(synthetic_adata.copy())
    reducer = DimensionalityReducer(test_config)
    adata_clust = reducer.run(adata_qc)

    assert "X_pca" in adata_clust.obsm
    assert "X_umap" in adata_clust.obsm
    assert "leiden" in adata_clust.obs.columns
    assert adata_clust.obs["leiden"].nunique() >= 2
    assert Path(test_config.paths["clustered_data"]).exists()


def test_cell_annotator(synthetic_adata: sc.AnnData, test_config: PipelineConfig):
    """Vérifie l'extraction des marqueurs Wilcoxon et l'annotation des deux populations."""
    adata_qc = QualityControlProcessor(test_config).run(synthetic_adata.copy())
    adata_clust = DimensionalityReducer(test_config).run(adata_qc)

    annotator = CellAnnotator(test_config)
    adata_annot = annotator.run(adata_clust)

    assert "cell_type" in adata_annot.obs.columns
    annotated_types = set(adata_annot.obs["cell_type"].unique())
    assert annotated_types == {"Basal_Cells", "Ciliated_Cells"}
    assert (test_config.tables_dir / "03_cluster_markers.csv").exists()
    assert (test_config.figures_dir / "03_markers_dotplot.png").exists()


def test_full_scrnaseq_pipeline_end_to_end(
    synthetic_adata: sc.AnnData, test_config: PipelineConfig, monkeypatch
):
    """Test d'intégration de bout en bout sur ScRNASeqPipeline (sans appel réseau Enrichr)."""
    monkeypatch.setattr(PipelineConfig, "from_yaml", classmethod(lambda cls, *a, **kw: test_config))

    pipeline = ScRNASeqPipeline()
    with patch("src.enrichment.gp.enrichr") as mock_enrichr:
        mock_enrichr.return_value.results = pd.DataFrame(
            {"Term": ["Epithelial development (GO:0060429)"], "Adjusted P-value": [0.001]}
        )
        final_adata = pipeline.run(synthetic_adata.copy())

    assert final_adata is not None
    assert final_adata.n_obs == 55
    assert "cell_type" in final_adata.obs.columns
    assert (test_config.tables_dir / "04_functional_enrichment_go.csv").exists()
