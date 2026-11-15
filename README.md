# 🧬 RNApipe — Modular & Object-Oriented scRNA-Seq Architecture

[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Engine-Scanpy%20%7C%20AnnData-00A98F)](https://scanpy.readthedocs.io/)
[![Architecture](https://img.shields.io/badge/Design-OOP%20%26%20YAML--Driven-8A2BE2)](#design-philosophy--positioning)
[![Testing](https://img.shields.io/badge/Tested%20with-Pytest-0A9EDC?logo=pytest&logoColor=white)](#automated-testing)

**RNApipe** is a configuration-driven, object-oriented single-cell RNA sequencing (**scRNA-seq**) pipeline built on top of **Scanpy** and **AnnData**.

### 🎯 Design Philosophy & Positioning
The scRNA-seq ecosystem is already well-served by mature frameworks (**Scanpy**, **Seurat**, and Nextflow pipelines like **nf-core/scrnaseq**). **RNApipe does not aim to compete on algorithmic breadth.** Instead, its primary contribution is **software architecture and maintainability**: decoupling analytical logic into testable, single-responsibility Python classes (`BasePipelineStep`), externalizing all biological parameters into typed YAML configurations, and enforcing automated unit/integration testing via `pytest`.

---

## 🖼️ Visual Benchmark: Reference Labels vs. Pipeline Predictions

> *Full interactive cross-validation, confusion matrices, and differential expression tables are available in [`notebooks/biological_report.ipynb`](notebooks/biological_report.ipynb).*

<table>
  <tr>
    <td colspan="2" align="center">
      <b>1. Cross-Validation Benchmark: Author Labels (CELLxGENE) vs. RNApipe Predictions</b><br/>
      <img src="assets/05_benchmark_reference_vs_predicted.png" alt="Reference vs Predicted UMAP and Confusion Matrix" width="100%"/>
      <sub>Side-by-side UMAP topology comparing published author annotations against automated signature-based predictions, paired with the row-normalized concordance matrix.</sub>
    </td>
  </tr>
  <tr>
    <td width="50%" align="center">
      <b>2. Canonical Epithelial Markers (DotPlot)</b><br/>
      <img src="assets/03_markers_dotplot.png" alt="Canonical Markers DotPlot" width="100%"/>
      <sub>Mean log-normalized expression and fraction of expressing cells across predicted lineages.</sub>
    </td>
    <td width="50%" align="center">
      <b>3. Pre-Filtering Quality Control Distributions</b><br/>
      <img src="assets/01_qc_violin_before_filtering.png" alt="QC Violin Plots" width="100%"/>
      <sub>Empirical distributions of detected genes, UMI counts, and mitochondrial transcript fraction used to calibrate QC cutoffs.</sub>
    </td>
  </tr>
</table>

---

## 📊 Quantitative Validation & Parameter Justification

To avoid arbitrary parameter selection, **RNApipe** was benchmarked on a human epithelial scRNA-seq subset retrieved from **CZI CELLxGENE** (see [`scripts/evaluate_benchmark.py`](scripts/evaluate_benchmark.py) and [`notebooks/biological_report.ipynb`](notebooks/biological_report.ipynb)):

### 1. Annotation Concordance Metrics (`reference_cell_type` vs. `cell_type`)
* **Global Clustering & Label Agreement:** Evaluated using **Adjusted Rand Index (ARI)** and **Normalized Mutual Information (NMI)** between published author annotations and RNApipe cluster-level predictions *(exact metrics exported to `results/tables/05_concordance_matrix_pct.csv`)*.
* **Lineage-Level Precision:** Canonical epithelial compartments (**Basal `KRT5+/TP63+/KRT15+`**, **Multiciliated `FOXJ1+/TPPP3+`**, and **Secretory `SCGB1A1+/MUC5AC+`**) exhibit high diagonal concordance, whereas rare or transitional states show higher label overlap (detailed in the methodological limitations below).

### 2. Leiden Resolution Sensitivity Analysis
Rather than assuming a single arbitrary clustering granularity, we evaluated cluster stability and concordance against reference annotations across four resolution regimes (`results/tables/06_leiden_resolution_sensitivity.csv`):

| Leiden Resolution | Biological Granularity | Behavior on Epithelial Dataset |
| :--- | :--- | :--- |
| **`0.2`** | Coarse | Under-segments distinct secretory subtypes (`club` vs. `goblet` merged) |
| **`0.5` (Default)** | **Balanced** | **Optimal trade-off resolving major epithelial lineages without over-splitting** |
| **`0.8`** | Fine | Splits metabolically active vs. quiescent sub-states within major lineages |
| **`1.0`** | Over-clustered | Fragments homogeneous populations into donor- or stress-driven subclusters |

### 3. Empirical Justification of Quality Control Thresholds
Default thresholds in `config/config.yaml` (`min_genes: 200`, `max_genes: 6000`, `max_pct_mt: 15.0`) were calibrated against the empirical quantiles of the raw count matrix:
* **Lower bound (`min_genes = 200`):** Removes empty droplets and ambient RNA debris falling below the 5th percentile of feature complexity.
* **Upper bound (`max_genes = 6000`):** Excludes putative homotypic/heterotypic doublets exceeding the 99th percentile of detected genes.
* **Mitochondrial ceiling (`max_pct_mt = 15.0%`):** Retains metabolically active airway/epithelial cells (which naturally exhibit higher baseline mitochondrial content than PBMCs) while discarding lysed/apoptotic cells.

---

## ⚠️ Methodological Limitations & Scope

1. **Canonical Signature Scoring vs. Rare/Transitional Cell States:** `CellAnnotator` assigns labels at the Leiden cluster level using curated marker gene sets (`config/markers.yaml`). While transparent, fast, and interpretable, this approach **cannot reliably resolve rare cell types** (e.g., pulmonary ionocytes `<1%` frequency) or continuous differentiation intermediates without sub-clustering or probabilistic label transfer (*CellTypist* / *scArches*).
2. **External API Dependency for Functional Enrichment:** `FunctionalEnrichment` queries the **Enrichr REST API** (`GO_Biological_Process_2023`) via `gseapy`. Because Enrichr is an external web service without pinned database snapshots, **offline execution is not supported** in the current release, and pathway results may evolve as upstream gene-set libraries are updated. Adding a local `.gmt` cache fallback is planned for future versions.

---

## 🏗️ Object-Oriented Architecture

```text
                       ┌──────────────────────────┐
                       │      PipelineConfig      │
                       │ (YAML Loader & Paths)    │
                       └────────────┬─────────────┘
                                    │ injected into
                       ┌────────────▼─────────────┐
                       │     BasePipelineStep     │  ◄── Abstract Base Class (ABC)
                       │   + run(adata) -> AnnData│
                       └────────────┬─────────────┘
                                    │ inherited by
        ┌───────────────────────────┼───────────────────────────┐
        ▼                           ▼                           ▼
┌───────────────────┐     ┌───────────────────┐     ┌───────────────────────┐
│ QualityControl    │     │ Dimensionality    │     │ CellAnnotator &       │
│ Processor         │ ──► │ Reducer           │ ──► │ FunctionalEnrichment  │
│ (QC, Norm, HVG)   │     │ (PCA, kNN, UMAP)  │     │ (Wilcoxon, GO Terms)  │
└───────────────────┘     └───────────────────┘     └───────────────────────┘
        ▲                                                       │
        └─────────────── Orchestrated in memory by ─────────────┘
                           ┌───────────────────┐
                           │ ScRNASeqPipeline  │
                           └───────────────────┘
