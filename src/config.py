from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List
import yaml


@dataclass
class PipelineConfig:
    """Encapsule et valide la configuration du pipeline scRNA-seq."""

    project_name: str
    tissue: str
    paths: Dict[str, str]
    qc: Dict[str, float]
    preprocessing: Dict[str, int]
    dimred_clustering: Dict[str, float]
    markers: Dict[str, List[str]] = field(default_factory=dict)

    @classmethod
    def from_yaml(
        cls,
        config_path: str = "config/config.yaml",
        markers_path: str = "config/markers.yaml",
    ) -> "PipelineConfig":
        with open(config_path, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)

        markers_dict = {}
        m_path = Path(markers_path)
        if m_path.exists():
            with open(m_path, "r", encoding="utf-8") as f:
                m_data = yaml.safe_load(f) or {}
                markers_dict = m_data.get("cell_types", {})

        return cls(
            project_name=cfg["project"]["name"],
            tissue=cfg["project"]["tissue"],
            paths=cfg["paths"],
            qc=cfg["qc"],
            preprocessing=cfg["preprocessing"],
            dimred_clustering=cfg["dimred_clustering"],
            markers=markers_dict,
        )

    @property
    def figures_dir(self) -> Path:
        path = Path(self.paths["figures_dir"])
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def tables_dir(self) -> Path:
        path = Path(self.paths["tables_dir"])
        path.mkdir(parents=True, exist_ok=True)
        return path