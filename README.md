### Repository Structure

* **`assets/`** — Portfolio visual previews rendered in the README (`300 DPI` PNG figures).
* **`config/`**
  * `config.yaml` — Biological & technical thresholds (QC filtering, PCA components, k-NN neighbors, Leiden resolution).
  * `markers.yaml` — Curated canonical gene signatures per human epithelial lineage.
* **`notebooks/`**
  * `biological_report.ipynb` — Interactive Jupyter report with cross-validation confusion matrices, DotPlots, and biological interpretation.
* **`scripts/`**
  * `fetch_epithelial_data.py` — Dynamic CZI CELLxGENE REST API dataset downloader (`.h5ad`) with automated gene symbol indexing.
* **`src/`**
  * `base.py` — Abstract base class (`BasePipelineStep`) enforcing the step contract and headless Matplotlib `Agg` backend.
  * `config.py` — Typed dataclass configuration loader (`PipelineConfig`).
  * `preprocessing.py` — `QualityControlProcessor` (cell/gene filtering, mitochondrial QC, `CP10K` normalization, `log1p`, HVG selection).
  * `clustering.py` — `DimensionalityReducer` (feature scaling, PCA, k-NN graph construction, UMAP embedding, Leiden clustering).
  * `annotation.py` — `CellAnnotator` (Wilcoxon rank-sum DGE, canonical signature scoring, automated cluster labeling, DotPlots).
  * `enrichment.py` — `FunctionalEnrichment` (Enrichr API integration for Gene Ontology Biological Process 2023).
  * `pipeline.py` — `ScRNASeqPipeline` end-to-end orchestrator managing in-memory `AnnData` transitions and checkpoints.
* **`tests/`**
  * `test_pipeline.py` — Pytest unit & integration test suite using synthetic multi-population `AnnData` datasets.
* **`pyproject.toml`** — Pytest root path and project configuration.
* **`requirements.txt`** — Reproducible Python dependencies.
* **`run_pipeline.py`** — Clean command-line entry point.

---

## 🔬 Biological Case Study: Human Epithelial Tissue

To validate biological accuracy, **RNApipe** was benchmarked on human epithelial single-cell transcriptomic data retrieved from **CZI CELLxGENE**:

* **Quality Control & Normalization:** Cells with low transcriptional complexity (`< 200` genes), putative doublets (`> 6,000` genes), or apoptotic signatures (`> 15%` mitochondrial counts) are filtered out. Raw counts are preserved in `adata.layers["counts"]`, normalized to `10,000` counts per cell (`CP10K`), log-transformed (`log1p`), and reduced to the top `2,000` highly variable genes (HVGs).
* **Lineage Segregation & Annotation:** Unsupervised Leiden clustering cleanly resolves major epithelial compartments without batch-driven fragmentation. Automated signature scoring (`CellAnnotator`) demonstrates strong concordance with published author labels (`reference_cell_type`):
  * **Basal / Progenitor Cells (`KRT5+`, `TP63+`, `KRT15+`):** Tissue-resident stem/progenitor population driving epithelial regeneration.
  * **Multiciliated Cells (`FOXJ1+`, `TPPP3+`, `PIFO+`):** Enriched in cilium assembly and microtubule-based movement GO terms.
  * **Secretory / Goblet & Club Cells (`SCGB1A1+`, `MUC5AC+`, `MUC5B+`):** Specialized in mucosal barrier defense and surfactant/mucin secretion.

---

## 🚀 Quickstart & Installation

**1. Clone the Repository & Create a Virtual Environment**

* `git clone https://github.com/KpakaMh/RNApipe.git`
* `cd RNApipe`
* `python -m venv .venv`
* Activate environment:
  * Windows (PowerShell): `.\.venv\Scripts\Activate.ps1`
  * Linux / macOS: `source .venv/bin/activate`
* `pip install -r requirements.txt`

**2. Run the Full Pipeline**

* **Option A — Fetch real human epithelial scRNA-seq data from CELLxGENE and run:**
  `python scripts/fetch_epithelial_data.py`
* **Option B — Run on your own `.h5ad` file (`data/raw.h5ad`) or fallback demo dataset:**
  `python run_pipeline.py`

All publication-ready figures (`300 DPI`) are exported to `results/figures/` and statistical tables (Wilcoxon markers & GO enrichment) to `results/tables/`.

---

## 🧪 Automated Testing

The pipeline includes a full **Pytest** suite (`tests/test_pipeline.py`) that generates a synthetic multi-population `AnnData` object in memory (including simulated apoptotic high-MT cells) to verify each processor class independently as well as the full end-to-end orchestrator:

* Run all unit and integration tests: `python -m pytest -v`

---

## ⚙️ Configuration

All biological parameters are decoupled from the source code inside `config/config.yaml` and `config/markers.yaml`, allowing instant adaptation to new tissues (e.g., PBMC, tumor microenvironment, neural organoids) without modifying Python classes:

* **QC thresholds (`config/config.yaml`):** `min_genes: 200`, `max_genes: 6000`, `min_cells: 3`, `max_pct_mt: 15.0`
* **Clustering parameters (`config/config.yaml`):** `n_pcs: 30`, `n_neighbors: 15`, `leiden_resolution: 0.5`

---

## 👤 Author

**Mouhamadou Adesanya KPAKA**  
M1 Student in Health Biology, Biotechnology & Computational Biology  
GitHub: [@KpakaMh](https://github.com/KpakaMh)
