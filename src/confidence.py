"""Confidence scoring for nominated perturbations.

A hypothesis is only as trustworthy as the evidence behind it. We combine four
orthogonal signals into a transparent, decomposable score (never a black box):

  1. Dataset-native reproducibility (from DE_stats .obs):
       - cross-guide agreement (guide_correlation_signif)
       - cross-donor agreement (donor_correlation_hits_mean)
       - on-target knockdown significance (ontarget_significant)
       - off-target penalty (distal_offtarget_flag, neighboring_gene_KD)
  2. Stability selection: bootstrap the solver; how often is this gene picked?
  3. Held-out-donor generalization: does the effect direction hold on unseen donors?
  4. External evidence (see evidence.py): literature / Open Targets support.

The final score is a weighted average of normalized components, and we ALWAYS return
the component breakdown so a reviewer can audit it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class ConfidenceReport:
    gene: str
    reproducibility: float
    stability: float
    heldout: float
    external: float
    components: dict

    @property
    def score(self) -> float:
        # Equal-weight by default; tune with held-out calibration, don't hard-sell.
        vals = [self.reproducibility, self.stability, self.heldout, self.external]
        vals = [v for v in vals if v is not None and not np.isnan(v)]
        return float(np.mean(vals)) if vals else float("nan")


def reproducibility_score(row: pd.Series) -> float:
    """Map a DE_stats .obs row's QC columns to [0, 1]. Missing -> neutral 0.5."""
    guide = row.get("guide_correlation_signif", np.nan)
    donor = row.get("donor_correlation_hits_mean", np.nan)
    ont = bool(row.get("ontarget_significant", False))
    offtarget = bool(row.get("distal_offtarget_flag", False)) or bool(
        row.get("neighboring_gene_KD", False)
    )
    parts = []
    if not np.isnan(guide):
        parts.append(np.clip((guide + 1) / 2, 0, 1))   # corr in [-1,1] -> [0,1]
    if not np.isnan(donor):
        parts.append(np.clip((donor + 1) / 2, 0, 1))
    parts.append(1.0 if ont else 0.3)
    base = float(np.mean(parts)) if parts else 0.5
    if offtarget:
        base *= 0.6                                     # penalize off-target risk
    return float(np.clip(base, 0, 1))


def stability_selection(solver_fn, E: np.ndarray, d: np.ndarray, gene_ids: np.ndarray,
                        n_boot: int = 200, subsample: float = 0.8,
                        seed: int = 0) -> dict[str, float]:
    """Fraction of bootstraps in which each gene is selected. CPU-cheap and honest."""
    rng = np.random.default_rng(seed)
    P = E.shape[0]
    counts: dict[str, int] = {}
    m = max(2, int(P * subsample))
    for _ in range(n_boot):
        rows = rng.choice(P, size=m, replace=False)
        res = solver_fn(E[rows], d)
        for local_i in res.indices:
            g = str(gene_ids[rows[local_i]])
            counts[g] = counts.get(g, 0) + 1
    return {g: c / n_boot for g, c in counts.items()}


def heldout_donor_consistency(effect_by_donor: dict[str, np.ndarray],
                              d: np.ndarray) -> float:
    """Mean cosine alignment of a gene's per-donor effect with the target direction.

    `effect_by_donor` maps donor_id -> that gene's effect vector (from the by-donor
    DE artifact). High, consistent alignment across donors -> high confidence.
    """
    cosims = []
    for _, e in effect_by_donor.items():
        na, nb = np.linalg.norm(e), np.linalg.norm(d)
        if na > 0 and nb > 0:
            cosims.append(float(e @ d / (na * nb)))
    if not cosims:
        return float("nan")
    # reward both magnitude and consistency (penalize sign flips across donors)
    return float(np.clip(np.mean(cosims) * (1 - np.std(cosims)), -1, 1))
