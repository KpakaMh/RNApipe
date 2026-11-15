import sys
from pathlib import Path
import requests
import scanpy as sc

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.pipeline import ScRNASeqPipeline


DATASETS_API = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) RNApipe/1.0"}


def find_epithelial_h5ad_asset(min_cells: int = 1000, max_dataset_cells: int = 20000) -> tuple[str, str]:
    resp = requests.get(DATASETS_API, headers=HEADERS, timeout=30)
    resp.raise_for_status()
    datasets = resp.json()

    candidates = []
    for ds in datasets:
        organisms = [o.get("label", "") for o in ds.get("organism", [])]
        if "Homo sapiens" not in organisms:
            continue

        title = ds.get("title", "")
        cell_types = " ".join(c.get("label", "") for c in ds.get("cell_type", []))
        tissues = " ".join(t.get("label", "") for t in ds.get("tissue", []))
        search_blob = f"{title} {cell_types} {tissues}".lower()

        if "epithelial" not in search_blob and "epithelium" not in search_blob:
            continue

        cell_count = ds.get("cell_count", 0)
        if not (min_cells <= cell_count <= max_dataset_cells):
            continue

        for asset in ds.get("assets", []):
            if asset.get("filetype") == "H5AD" and asset.get("url"):
                candidates.append((cell_count, title, asset["url"]))
                break

    if not candidates:
        raise RuntimeError("Aucun dataset épithélial humain léger trouvé dans le catalogue.")

    candidates.sort(key=lambda x: x[0])
    chosen_count, chosen_title, chosen_url = candidates[0]
    print(f"[FETCH] Dataset sélectionné : '{chosen_title}' ({chosen_count} cellules)")
    return chosen_title, chosen_url


def download_human_epithelial_h5ad(
    output_path: str = "data/raw.h5ad",
    max_cells: int = 3000,
) -> sc.AnnData:
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    temp_file = out_file.parent / "temp_full_epithelial.h5ad"

    if not temp_file.exists():
        print("[FETCH] Recherche d'un dataset épithélial humain actif sur CELLxGENE...")
        _, h5ad_url = find_epithelial_h5ad_asset()

        print("[FETCH] Téléchargement du fichier .h5ad...")
        with requests.get(h5ad_url, headers=HEADERS, stream=True, timeout=180) as r:
            r.raise_for_status()
            with open(temp_file, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
    else:
        print(f"[FETCH] Réutilisation du fichier déjà téléchargé : {temp_file}")

    print("[FETCH] Lecture et normalisation des métadonnées AnnData...")
    adata = sc.read_h5ad(temp_file)

    if adata.n_obs > max_cells:
        sc.pp.subsample(adata, n_obs=max_cells, random_state=42)
        print(f"[FETCH] Sous-échantillonnage reproductible à {max_cells} cellules.")

    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        raw_adata.obs = adata.obs.copy()
        adata = raw_adata

    # Indexation par symboles de gènes et suppression du nom d'index pour éviter tout conflit h5ad
    for col in ["feature_name", "gene_symbols", "gene_short_name"]:
        if col in adata.var.columns:
            adata.var_names = adata.var[col].astype(str).values
            break
    adata.var.index.name = None
    adata.obs.index.name = None
    adata.var_names_make_unique()

    if "cell_type" in adata.obs.columns:
        adata.obs.rename(columns={"cell_type": "reference_cell_type"}, inplace=True)

    adata.write_h5ad(out_file)
    temp_file.unlink(missing_ok=True)
    print(f"[FETCH] Données prêtes dans {out_file} ({adata.n_obs} cellules x {adata.n_vars} gènes)")
    return adata


if __name__ == "__main__":
    download_human_epithelial_h5ad(output_path="data/raw.h5ad", max_cells=3000)

    pipeline = ScRNASeqPipeline(
        config_path="config/config.yaml",
        markers_path="config/markers.yaml",
    )
    pipeline.run()