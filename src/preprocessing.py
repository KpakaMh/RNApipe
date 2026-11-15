from pathlib import Path
from typing import Optional
import scanpy as sc
import matplotlib.pyplot as plt
from src.base import BasePipelineStep


class QualityControlProcessor(BasePipelineStep):
    """Gère le chargement, le contrôle qualité, le filtrage, la normalisation et les HVG."""

    def load_or_fetch_data(self) -> sc.AnnData:
        raw_path = Path(self.config.paths["raw_data"])
        if raw_path.exists():
            return sc.read_h5ad(raw_path)

        raw_path.parent.mkdir(parents=True, exist_ok=True)
        adata = sc.datasets.pbmc68k_reduced()
        if adata.raw is not None:
            adata = adata.raw.to_adata()
        adata.write_h5ad(raw_path)
        return adata

    def _plot_qc_metrics(self, adata: sc.AnnData) -> None:
        fig_dir = self.config.figures_dir
        grid = sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False,
        )
        grid.fig.savefig(fig_dir / "01_qc_violin_before_filtering.png", dpi=300, bbox_inches="tight")
        plt.close(grid.fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
        plt.tight_layout()
        fig.savefig(fig_dir / "01_qc_scatter_metrics.png", dpi=300)
        plt.close(fig)

    def run(self, adata: Optional[sc.AnnData] = None) -> sc.AnnData:
        if adata is None:
            adata = self.load_or_fetch_data()

        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

        obs_qc, var_qc = sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False
        )
        adata.obs[obs_qc.columns] = obs_qc
        adata.var[var_qc.columns] = var_qc

        self._plot_qc_metrics(adata)

        qc_cfg = self.config.qc
        sc.pp.filter_cells(adata, min_genes=int(qc_cfg["min_genes"]))
        adata = adata[adata.obs["n_genes_by_counts"] < qc_cfg["max_genes"], :].copy()
        adata = adata[adata.obs["pct_counts_mt"] < qc_cfg["max_pct_mt"], :].copy()
        sc.pp.filter_genes(adata, min_cells=int(qc_cfg["min_cells"]))

        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=self.config.preprocessing["target_sum"])
        sc.pp.log1p(adata)
        adata.raw = adata

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=int(self.config.preprocessing["n_top_genes"]),
            flavor="seurat",
        )
        sc.pl.highly_variable_genes(adata, show=False)
        plt.savefig(self.config.figures_dir / "01_hvg_dispersion.png", dpi=300, bbox_inches="tight")
        plt.close()

        out_path = Path(self.config.paths["qc_data"])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_path)
        return adata