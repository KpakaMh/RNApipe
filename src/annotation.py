from pathlib import Path
from typing import Optional
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from src.base import BasePipelineStep


class CellAnnotator(BasePipelineStep):
    """Identifie les gènes marqueurs (Wilcoxon) et annote les populations cellulaires."""

    def _extract_top_markers(self, adata: sc.AnnData, n_genes: int = 10) -> pd.DataFrame:
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
        result = adata.uns["rank_genes_groups"]
        groups = result["names"].dtype.names

        records = []
        for group in groups:
            for i in range(min(n_genes, len(result["names"][group]))):
                records.append({
                    "cluster": group,
                    "rank": i + 1,
                    "gene": result["names"][group][i],
                    "logfoldchanges": result["logfoldchanges"][group][i],
                    "pvals_adj": result["pvals_adj"][group][i],
                })
        df_markers = pd.DataFrame(records)
        df_markers.to_csv(self.config.tables_dir / "03_cluster_markers.csv", index=False)
        return df_markers

    def _score_and_annotate(self, adata: sc.AnnData) -> None:
        valid_markers = {}
        available_genes = set(adata.raw.var_names if adata.raw is not None else adata.var_names)

        for cell_type, genes in self.config.markers.items():
            present = [g for g in genes if g in available_genes]
            if present:
                valid_markers[cell_type] = present
                sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{cell_type}")

        if not valid_markers:
            adata.obs["cell_type"] = "Cluster_" + adata.obs["leiden"].astype(str)
            return

        score_cols = [f"score_{ct}" for ct in valid_markers.keys()]
        cluster_means = adata.obs.groupby("leiden", observed=False)[score_cols].mean()
        best_type_per_cluster = cluster_means.idxmax(axis=1).str.replace("score_", "", regex=False)

        adata.obs["cell_type"] = adata.obs["leiden"].map(best_type_per_cluster).astype("category")

        dp = sc.pl.dotplot(adata, var_names=valid_markers, groupby="cell_type", show=False)
        dp["mainplot_ax"].figure.savefig(
            self.config.figures_dir / "03_markers_dotplot.png", dpi=300, bbox_inches="tight"
        )
        plt.close("all")

    def run(self, adata: Optional[sc.AnnData] = None) -> sc.AnnData:
        if adata is None:
            adata = sc.read_h5ad(Path(self.config.paths["clustered_data"]))

        self._extract_top_markers(adata)
        self._score_and_annotate(adata)

        sc.pl.umap(adata, color=["leiden", "cell_type"], wspace=0.4, show=False)
        plt.savefig(self.config.figures_dir / "03_umap_annotated.png", dpi=300, bbox_inches="tight")
        plt.close()

        out_path = Path(self.config.paths["annotated_data"])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_path)
        return adata