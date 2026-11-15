from pathlib import Path
from typing import Optional
import scanpy as sc
import matplotlib.pyplot as plt
from src.base import BasePipelineStep


class DimensionalityReducer(BasePipelineStep):
    """Exécute la mise à l'échelle, la PCA, le graphe k-NN, l'UMAP et le clustering Leiden."""

    def run(self, adata: Optional[sc.AnnData] = None) -> sc.AnnData:
        if adata is None:
            adata = sc.read_h5ad(Path(self.config.paths["qc_data"]))

        dim_cfg = self.config.dimred_clustering
        fig_dir = self.config.figures_dir

        sc.pp.scale(adata, max_value=10)
        n_pcs = min(int(dim_cfg["n_pcs"]), adata.n_vars - 1, adata.n_obs - 1)
        sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

        sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs, show=False)
        plt.savefig(fig_dir / "02_pca_variance_ratio.png", dpi=300, bbox_inches="tight")
        plt.close()

        sc.pp.neighbors(adata, n_neighbors=int(dim_cfg["n_neighbors"]), n_pcs=n_pcs)
        sc.tl.umap(adata)

        resolution = float(dim_cfg["leiden_resolution"])
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added="leiden",
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )

        sc.pl.umap(
            adata,
            color=["leiden", "n_genes_by_counts", "total_counts"],
            wspace=0.35,
            show=False,
        )
        plt.savefig(fig_dir / "02_umap_clusters_qc.png", dpi=300, bbox_inches="tight")
        plt.close()

        out_path = Path(self.config.paths["clustered_data"])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_path)
        return adata