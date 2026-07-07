"""Load and reduce the CD4+ T-cell Perturb-seq DE artifacts for CPU-only work.

Core input: GWCD4i.DE_stats.h5ad
  .obs  -> one row per (perturbed gene, culture condition); includes reproducibility
           columns (guide_correlation_*, donor_correlation_*, ontarget_significant,
           distal_offtarget_flag, keep_test_genes, ...)
  .layers['zscore'] -> perturbation-effect matrix (P x G), our dictionary of causal
           effect vectors.

Design goal: never hold the full multi-layer AnnData in RAM. Load one layer as
float32, subset to high-confidence perturbations and highly-variable genes, cache.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DE_STATS = DATA_DIR / "GWCD4i.DE_stats.h5ad"


@dataclass
class PerturbDictionary:
    """A reduced perturbation-effect dictionary ready for the counterfactual solver."""
    effects: np.ndarray          # (P, G) float32 z-scored logFC
    perturbation_ids: np.ndarray # (P,) perturbed gene names
    conditions: np.ndarray       # (P,) culture condition per row
    genes: np.ndarray            # (G,) gene symbols
    obs_meta: "object"           # pandas DataFrame of reproducibility/QC columns


def load_de_dictionary(
    path: Path = DE_STATS,
    condition: str = "Stim8hr",
    layer: str = "zscore",
    high_confidence_only: bool = True,
    n_hvg: int | None = 2000,
) -> PerturbDictionary:
    """Load a single condition slice of the DE effect matrix as a compact dictionary.

    Parameters
    ----------
    condition : one of {"Rest", "Stim8hr", "Stim48hr"}.
    high_confidence_only : keep perturbations flagged `keep_test_genes` and with a
        significant on-target knockdown, so we build hypotheses from trustworthy rows.
    n_hvg : optionally restrict genes to the top-N by variance across perturbations to
        fit a laptop's RAM. Set None to keep all ~10k genes.

    Notes
    -----
    Uses backed mode so only the requested layer is materialized. This is the key to
    running on CPU/laptop without loading the full ~GB object.
    """
    import anndata as ad  # imported lazily so `--check` works without the file

    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs
    mask = (obs["culture_condition"] == condition).to_numpy()
    if high_confidence_only:
        keep = obs.get("keep_test_genes")
        ont = obs.get("ontarget_significant")
        if keep is not None:
            mask &= keep.fillna(False).to_numpy().astype(bool)
        if ont is not None:
            mask &= ont.fillna(False).to_numpy().astype(bool)

    idx = np.where(mask)[0]
    effects = np.asarray(adata.layers[layer][idx, :], dtype=np.float32)
    genes = adata.var["gene_name"].to_numpy()

    if n_hvg is not None and n_hvg < effects.shape[1]:
        var_rank = np.argsort(np.nanvar(effects, axis=0))[::-1][:n_hvg]
        effects = effects[:, var_rank]
        genes = genes[var_rank]

    effects = np.nan_to_num(effects, nan=0.0)
    sub_obs = obs.iloc[idx].copy()
    return PerturbDictionary(
        effects=effects,
        perturbation_ids=sub_obs["target_contrast_gene_name"].to_numpy(),
        conditions=sub_obs["culture_condition"].to_numpy(),
        genes=genes,
        obs_meta=sub_obs,
    )


def _check() -> None:
    """Lightweight environment/data check that does not require the big file."""
    print("data dir:", DATA_DIR)
    print("DE_stats present:", DE_STATS.exists())
    if not DE_STATS.exists():
        print("  -> see data/README.md to fetch GWCD4i.DE_stats.h5ad")
    for pkg in ("numpy", "pandas", "scipy", "sklearn", "anndata"):
        try:
            __import__(pkg)
            print(f"  [ok] {pkg}")
        except ImportError:
            print(f"  [missing] {pkg} (pip install -r requirements.txt)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check", action="store_true", help="verify env + data presence")
    args = p.parse_args()
    if args.check:
        _check()
