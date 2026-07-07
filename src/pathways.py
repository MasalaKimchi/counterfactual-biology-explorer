"""Pathway / gene-set interpretation of a nominated minimal set.

Turns a bare gene list into biology: which pathways / GO terms are enriched, so a
reviewer understands *why* the set might move the target state. Uses gseapy's Enrichr
interface (over-representation), which is CPU-light.
"""
from __future__ import annotations

import pandas as pd

DEFAULT_LIBRARIES = [
    "GO_Biological_Process_2021",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
]


def enrich(gene_list: list[str], libraries: list[str] | None = None) -> pd.DataFrame:
    """Over-representation analysis on a nominated gene set.

    Returns a tidy DataFrame (term, library, adjusted p, overlap genes). Empty and
    safe if gseapy or network is unavailable — interpretation is a nicety, not a gate.
    """
    libraries = libraries or DEFAULT_LIBRARIES
    if len(gene_list) < 2:
        return pd.DataFrame(columns=["Term", "Gene_set", "Adjusted P-value", "Genes"])
    try:
        import gseapy as gp

        res = gp.enrichr(gene_list=list(gene_list), gene_sets=libraries,
                         outdir=None, no_plot=True)
        df = res.results.sort_values("Adjusted P-value").head(25)
        return df[["Term", "Gene_set", "Adjusted P-value", "Combined Score", "Genes"]]
    except Exception as exc:  # pragma: no cover
        return pd.DataFrame({"error": [str(exc)]})
