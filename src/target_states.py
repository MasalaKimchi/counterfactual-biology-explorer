"""Build target-state direction vectors from the dataset's own signature tables.

A "target direction" d in gene space encodes the desired transcriptomic shift, e.g.
"become more Th1-like" or "look transcriptionally younger". We align these vectors to
the same gene ordering as the perturbation dictionary so the solver can operate.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def _align_to_genes(sig: pd.DataFrame, genes: np.ndarray, value_col: str) -> np.ndarray:
    """Project a signature table onto the solver's gene ordering (missing -> 0)."""
    lut = dict(zip(sig["gene_name"].astype(str), sig[value_col].astype(float)))
    return np.array([lut.get(str(g), 0.0) for g in genes], dtype=np.float32)


def polarization_target(genes: np.ndarray, direction: str = "toward_Th1") -> np.ndarray:
    """Th1<->Th2 polarization direction from Th2_Th1_polarization_signature table.

    `direction`:
      "toward_Th1" uses the signature as-is; "toward_Th2" flips the sign.
    Value column is the per-gene z-score of the Th1-vs-Th2 contrast.
    """
    path = DATA_DIR / "Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv"
    sig = pd.read_csv(path)
    d = _align_to_genes(sig, genes, value_col="zscore")
    return d if direction == "toward_Th1" else -d


def aging_target(genes: np.ndarray, direction: str = "toward_young") -> np.ndarray:
    """Reverse-aging direction from the CD4+ T-cell aging signature table.

    The aging signature encodes aged-vs-young change; to move *toward young* we invert.
    """
    path = DATA_DIR / "CD4T_aging_signature_DE_results_full.suppl_table.csv"
    sig = pd.read_csv(path)
    d_aging = _align_to_genes(sig, genes, value_col="zscore")
    return -d_aging if direction == "toward_young" else d_aging


def normalize(d: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(d)
    return d / n if n > 0 else d
