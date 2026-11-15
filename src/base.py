from abc import ABC, abstractmethod
from typing import Optional
import matplotlib

# Force un backend sans interface graphique (évite les erreurs Tcl/Tk et fonctionne sur serveur/Docker)
matplotlib.use("Agg")

import scanpy as sc
from src.config import PipelineConfig


class BasePipelineStep(ABC):
    """Classe de base abstraite pour toutes les étapes du pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        sc.settings.verbosity = 1

    @abstractmethod
    def run(self, adata: Optional[sc.AnnData] = None) -> sc.AnnData:
        """Exécute le traitement de l'étape et retourne l'objet AnnData mis à jour."""
        raise NotImplementedError