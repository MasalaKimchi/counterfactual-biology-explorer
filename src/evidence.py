"""Orthogonal external evidence for a nominated gene.

Evidence is SUPPORT, never proof. We surface prior knowledge so a human can judge a
hypothesis — we never let a literature hit inflate a claim of causality or rescue.

Two lightweight, no-heavy-dependency sources:
  - PubMed E-utilities (esearch) for co-mention counts with the target process.
  - Open Targets GraphQL for gene-disease / target associations.

In the hackathon build these can also be routed through the bio-research MCP tools
(PubMed, Open Targets) if available; the REST fallbacks below keep the repo runnable
standalone.
"""
from __future__ import annotations

from dataclasses import dataclass

import requests

PUBMED = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"


@dataclass
class GeneEvidence:
    gene: str
    pubmed_count: int
    example_pmids: list[str]
    opentargets_hits: list[dict]

    def external_score(self) -> float:
        """Squash co-mention counts into [0,1]. Deliberately conservative:
        absence of literature is NOT evidence of absence, so 0 mentions -> 0.3, not 0.
        """
        import math
        base = 0.3 + 0.7 * (1 - math.exp(-self.pubmed_count / 10.0))
        return float(min(base, 1.0))


def pubmed_comention(gene: str, process: str, retmax: int = 5,
                     timeout: int = 20) -> tuple[int, list[str]]:
    """Count PubMed records co-mentioning a gene and the target process."""
    term = f'("{gene}"[Title/Abstract]) AND ("{process}"[Title/Abstract])'
    params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": retmax}
    r = requests.get(PUBMED, params=params, timeout=timeout)
    r.raise_for_status()
    res = r.json().get("esearchresult", {})
    return int(res.get("count", 0)), list(res.get("idlist", []))


def opentargets_associations(gene_symbol: str, timeout: int = 20) -> list[dict]:
    """Top disease associations for a target from Open Targets (best-effort)."""
    query = """
    query assoc($q: String!) {
      search(queryString: $q, entityNames: ["target"]) {
        hits { id name }
      }
    }
    """
    try:
        r = requests.post(OPENTARGETS, json={"query": query,
                          "variables": {"q": gene_symbol}}, timeout=timeout)
        r.raise_for_status()
        return r.json().get("data", {}).get("search", {}).get("hits", [])
    except Exception:
        return []  # evidence layer must never crash the pipeline


def collect(gene: str, process: str = "T cell polarization") -> GeneEvidence:
    try:
        count, pmids = pubmed_comention(gene, process)
    except Exception:
        count, pmids = 0, []
    return GeneEvidence(
        gene=gene, pubmed_count=count, example_pmids=pmids,
        opentargets_hits=opentargets_associations(gene),
    )
