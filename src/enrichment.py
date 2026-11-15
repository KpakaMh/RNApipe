from pathlib import Path
from typing import Optional
import pandas as pd
import scanpy as sc
import gseapy as gp
from src.base import BasePipelineStep


class FunctionalEnrichment(BasePipelineStep):
    """Réalise l'analyse d'enrichissement fonctionnel (GO Biological Process) par cluster."""

    def run(self, adata: Optional[sc.AnnData] = None) -> sc.AnnData:
        if adata is None:
            adata = sc.read_h5ad(Path(self.config.paths["annotated_data"]))

        markers_csv = self.config.tables_dir / "03_cluster_markers.csv"
        if not markers_csv.exists():
            return adata

        df_markers = pd.read_csv(markers_csv)
        all_enrichments = []

        for cluster, group_df in df_markers.groupby("cluster"):
            gene_list = group_df["gene"].dropna().tolist()
            if len(gene_list) < 3:
                continue
            try:
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=["GO_Biological_Process_2023"],
                    organism="human",
                    outdir=None,
                )
                res = enr.results.head(5).copy()
                res.insert(0, "cluster", cluster)
                all_enrichments.append(res)
            except Exception as exc:
                print(f"[ENRICH] Avertissement pour le cluster {cluster} : {exc}")

        if all_enrichments:
            final_df = pd.concat(all_enrichments, ignore_index=True)
            out_csv = self.config.tables_dir / "04_functional_enrichment_go.csv"
            final_df.to_csv(out_csv, index=False)

        return adata