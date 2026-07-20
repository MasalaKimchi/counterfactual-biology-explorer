"""screenloop: certificate-guided sequential design for combinatorial screens.

A thin, fail-closed layer over :mod:`reachability` and :mod:`combicone`. Where
:mod:`combicone` scores and certifies combinations one at a time, this module puts
the certificate inside a *loop* and asks the two questions a screening *campaign*
faces:

1. **Sequential acquisition** — replay a screen as batches revealed over rounds and
   measure how fast each acquisition policy discovers the emergent combinations.
   :func:`replay_campaign`. This exists to state an honest result: on Norman the
   certificate-adaptive residual does **not** beat the cheap training-free triage
   score at *choosing which combination to run next* (both ~consistent with the
   finding that the cone's value is not ranking accuracy). The harness is the
   instrument that makes that comparison reproducible, not a claim that the
   certificate wins it.

2. **Library augmentation** — the capability the certificate *does* uniquely
   provide: given the combinations that leave the current library's cone, aggregate
   their model-relative separators into a single "unmet demand" direction and rank
   which **new single perturbation** to add to the library next. :func:`nominate_atoms`.
   A forward predictor (GEARS, scGPT, STATE, CPA) emits a prediction for every
   input; it never emits "your library is missing an axis, and here is the axis".
   The separator is exactly that object. :func:`held_out_single_recovery` is the
   falsifiable test of the capability: hide a measured single, and check whether the
   aggregated separator of the combinations that needed it points back at it.

Design contract (inherited from reachability.py / combicone.py and NOT relaxed):
  * ``atoms`` / effect matrices are ``(n_perturbations, n_genes)``.
  * All geometry is *model-relative*: a separator certifies a direction outside the
    non-negative cone of THESE measured effects under THIS metric — never a claim of
    biological necessity or that adding the nominated perturbation will "work" at the
    bench. Every public result echoes that scope in a ``scope`` field.
  * Deterministic given a seed. No hidden global state.

Honest performance envelope (Norman combinatorial CRISPRa screen, file label
A549 / canonically K562; 105 single-gene atoms, 131 measured doubles, 40 emergent
under the two-bar noise-robust label):
  * Sequential acquisition, wells to discover 90% of the 40 emergent combinations
    (batch 8): random ~120, magnitude ~104, training-free triage ~96,
    certificate-adaptive residual ~120. The certificate does not win acquisition.
  * Library augmentation, recover-the-held-out-single (53 singles in >=2 doubles):
    the aggregated separator ranks the true held-out single at median rank 1
    (top-1 0.98, top-5 0.98) versus a naive "average the doubles, take the most
    similar single" baseline (top-1 0.55) and a random null (top-1 ~0.01). This
    recovery is the differentiating result and is validated against a permutation
    null and a magnitude-confound control in ``scripts/`` (see repo docs).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np

import combicone as cc
import reachability as rx

__all__ = [
    "CampaignResult",
    "AtomNomination",
    "replay_campaign",
    "nominate_atoms",
    "held_out_single_recovery",
    "ACQUISITION_POLICIES",
]

_EPS = 1e-12
_SCOPE = (
    "model-relative: separators and nominations are directions outside the "
    "non-negative cone of the supplied effect atoms under the chosen metric; a "
    "nominated perturbation is a ranked hypothesis about an unmet library axis, "
    "not a validated intervention"
)

# Built-in acquisition policies. A policy name maps to how the next batch is
# scored; higher score = acquired earlier. Custom callables are also accepted by
# :func:`replay_campaign` (see its docstring).
ACQUISITION_POLICIES = ("random", "magnitude", "triage", "cone_adaptive")


# --------------------------------------------------------------------------- #
# Small shared helpers
# --------------------------------------------------------------------------- #
def _unit(v: np.ndarray) -> np.ndarray:
    n = rx._stable_norm(v)
    return v / n if n > _EPS else v


def _cos_rows(matrix: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """Cosine of each row of ``matrix`` with ``direction`` (vectorised)."""
    d = np.asarray(direction, dtype=float)
    dn = rx._stable_norm(d)
    if dn <= _EPS:
        return np.zeros(matrix.shape[0], dtype=float)
    mn = np.linalg.norm(matrix, axis=1)
    safe = mn > _EPS
    out = np.zeros(matrix.shape[0], dtype=float)
    out[safe] = (matrix[safe] @ d) / (mn[safe] * dn)
    return out


def _validate_atoms(atoms: np.ndarray, name: str = "atoms") -> np.ndarray:
    a = np.asarray(atoms, dtype=float)
    if a.ndim != 2:
        raise rx.InputError(f"{name} must be 2D (n_atoms, n_genes)")
    if a.shape[0] < 1:
        raise rx.InputError(f"{name} must have at least one row")
    if not np.all(np.isfinite(a)):
        raise rx.InputError(f"{name} must be finite (no NaN/Inf)")
    return a


def _validate_combo_idx(
    combo_atom_idx: np.ndarray, n_atoms: int
) -> np.ndarray:
    idx = np.asarray(combo_atom_idx)
    if idx.ndim != 2 or idx.shape[1] < 2:
        raise rx.InputError("combo_atom_idx must be 2D (n_combos, k>=2)")
    if not np.issubdtype(idx.dtype, np.integer):
        raise rx.InputError("combo_atom_idx must be integer indices into atoms")
    if idx.min() < 0 or idx.max() >= n_atoms:
        raise rx.InputError("combo_atom_idx entries out of range for atoms")
    return idx


# --------------------------------------------------------------------------- #
# 1. Sequential-acquisition campaign harness
# --------------------------------------------------------------------------- #
@dataclass
class CampaignResult:
    """Discovery trajectory of one sequential-acquisition campaign.

    Attributes
    ----------
    policy : str
        Acquisition policy name (or ``"custom"`` for a supplied callable).
    wells : (n_rounds + 1,) int array
        Cumulative number of combinations measured after each round, starting at 0.
    found : (n_rounds + 1,) int array
        Cumulative number of *emergent* combinations discovered, aligned to ``wells``.
    pick_order : (n_combos,) int array
        The order in which combinations were acquired (indices into the input combos).
    n_emergent : int
        Total emergent combinations in the screen (the discoverable maximum).
    grow_cone : bool
        Whether measured combinations were appended to the cone between rounds.
    scope : str
        Model-relative disclaimer.
    """

    policy: str
    wells: np.ndarray
    found: np.ndarray
    pick_order: np.ndarray
    n_emergent: int
    grow_cone: bool
    scope: str = _SCOPE

    def wells_to_fraction(self, frac: float) -> float:
        """Wells needed to discover ``frac`` of all emergent combos (NaN if never)."""
        if not 0.0 < frac <= 1.0:
            raise rx.InputError("frac must be in (0, 1]")
        target = frac * self.n_emergent
        hit = np.asarray(self.found) >= target
        if not hit.any():
            return float("nan")
        return float(self.wells[int(np.argmax(hit))])


def _acquisition_scores(
    policy: str,
    *,
    cone: np.ndarray,
    add_pred: np.ndarray,
    triage_static: np.ndarray,
    mag_static: np.ndarray,
    unmeasured: list[int],
    rng: np.random.Generator,
    gene_weights: np.ndarray | None,
) -> np.ndarray:
    """Score unmeasured combos under a built-in policy (higher = acquire first)."""
    if policy == "random":
        return rng.random(len(unmeasured))
    if policy == "magnitude":
        return mag_static[unmeasured]
    if policy == "triage":
        return triage_static[unmeasured]
    if policy == "cone_adaptive":
        # Residual of each combo's ADDITIVE prediction vs the CURRENT (growing)
        # cone. Prospective: uses only the singles' additive prediction, never the
        # measured combo effect. One NNLS per unmeasured combo per round.
        return np.array(
            [
                rx.project_cone(
                    cone, add_pred[k], gene_weights=gene_weights
                ).residual_fraction
                for k in unmeasured
            ],
            dtype=float,
        )
    raise rx.InputError(
        f"unknown policy {policy!r}; use one of {ACQUISITION_POLICIES} or a callable"
    )


def replay_campaign(
    atoms: np.ndarray,
    combo_effects: np.ndarray,
    combo_atom_idx: np.ndarray,
    labels: np.ndarray,
    *,
    policy: str | Callable = "triage",
    batch_size: int = 8,
    grow_cone: bool = True,
    gene_weights: np.ndarray | None = None,
    pairwise: str = "mean",
    seed: int = 0,
) -> CampaignResult:
    """Replay a combinatorial screen as a sequential-acquisition campaign.

    Starting from a library of single-gene ``atoms``, reveal combinations in
    batches. Each round an acquisition policy ranks the still-unmeasured
    combinations, the top ``batch_size`` are "measured" (their emergence label is
    revealed and, if ``grow_cone``, their measured effect is appended to the cone so
    later rounds see a richer library), and the cumulative discovery of *emergent*
    combinations is recorded. The result's discovery curve is what different
    policies are compared on.

    This harness makes no claim that any policy wins; on Norman the
    certificate-adaptive residual ties or trails the training-free triage score
    (see module docstring). It exists to measure that honestly and reproducibly.

    Parameters
    ----------
    atoms : (n_atoms, n_genes) array
        Single-gene effect vectors — the starting library / cone seed.
    combo_effects : (n_combos, n_genes) array
        Measured effect vectors of the combinations being triaged.
    combo_atom_idx : (n_combos, k) int array
        For each combination, the row indices of its ``k`` constituent atoms. Used
        for the additive prediction (a+b+...) and the triage cosine; the campaign
        never uses ``combo_effects`` to *rank* (only to reveal the label and, when
        growing, to enrich the cone).
    labels : (n_combos,) bool array
        Ground-truth emergence label (e.g. the two-bar noise-robust label from
        :func:`combicone.certify_emergence`). The discovery curve counts these.
    policy : str or callable
        One of :data:`ACQUISITION_POLICIES`, or a callable
        ``policy(cone, add_pred, unmeasured, rng) -> scores`` returning one score
        per index in ``unmeasured`` (higher = acquired first).
    batch_size : int
        Combinations measured per round.
    grow_cone : bool
        If True (default), append each measured combo's effect to the cone between
        rounds, so ``cone_adaptive`` sees a library that grows as the screen runs.
    gene_weights : (n_genes,) array, optional
        Per-gene metric weights passed through to cone projections.
    pairwise : str
        Aggregation for the triage cosine feature ("mean" or "max"); at k=2 both
        equal ``-cos(a, b)``.
    seed : int
        RNG seed (determinism; only the ``random`` policy consumes it).

    Returns
    -------
    CampaignResult
    """
    atoms = _validate_atoms(atoms)
    combo_effects = _validate_atoms(combo_effects, "combo_effects")
    idx = _validate_combo_idx(combo_atom_idx, atoms.shape[0])
    labels = np.asarray(labels).astype(bool)
    n = combo_effects.shape[0]
    if not (idx.shape[0] == labels.shape[0] == n):
        raise rx.InputError(
            "combo_effects, combo_atom_idx and labels must share n_combos"
        )
    if combo_effects.shape[1] != atoms.shape[1]:
        raise rx.InputError("combo_effects and atoms must share n_genes")
    if batch_size < 1:
        raise rx.InputError("batch_size must be >= 1")
    if pairwise not in ("mean", "max"):
        raise rx.InputError("pairwise must be 'mean' or 'max'")

    rng = np.random.default_rng(seed)

    # Additive prediction per combo (singles only) and the two static policy scores.
    add_pred = atoms[idx].sum(axis=1)  # (n_combos, n_genes)
    mag_static = np.linalg.norm(add_pred, axis=1)
    triage_static = np.array(
        [-cc.combo_cosine(list(atoms[row]), agg=pairwise) for row in idx],
        dtype=float,
    )

    unmeasured = list(range(n))
    cone = atoms.copy()
    wells = [0]
    found = [0]
    pick_order: list[int] = []
    n_found = 0

    custom = callable(policy)
    while unmeasured:
        if custom:
            scores = np.asarray(
                policy(cone, add_pred, list(unmeasured), rng), dtype=float
            )
            if scores.shape[0] != len(unmeasured):
                raise rx.InputError(
                    "custom policy must return one score per unmeasured combo"
                )
        else:
            scores = _acquisition_scores(
                policy,
                cone=cone,
                add_pred=add_pred,
                triage_static=triage_static,
                mag_static=mag_static,
                unmeasured=unmeasured,
                rng=rng,
                gene_weights=gene_weights,
            )
        # Deterministic tie-break by combo index so results are reproducible.
        order = sorted(
            range(len(unmeasured)), key=lambda i: (-scores[i], unmeasured[i])
        )
        take = order[:batch_size]
        picked = [unmeasured[i] for i in take]
        for k in picked:
            n_found += int(labels[k])
            pick_order.append(k)
            if grow_cone:
                cone = np.vstack([cone, combo_effects[k][None, :]])
        for k in picked:
            unmeasured.remove(k)
        wells.append(wells[-1] + len(picked))
        found.append(n_found)

    return CampaignResult(
        policy=policy if isinstance(policy, str) else "custom",
        wells=np.array(wells, dtype=int),
        found=np.array(found, dtype=int),
        pick_order=np.array(pick_order, dtype=int),
        n_emergent=int(labels.sum()),
        grow_cone=grow_cone,
    )


# --------------------------------------------------------------------------- #
# 2. Library augmentation — nominate the next single perturbation to add
# --------------------------------------------------------------------------- #
@dataclass
class AtomNomination:
    """Ranking of candidate perturbations to add to a library.

    Attributes
    ----------
    ranking : (n_candidates,) int array
        Candidate indices sorted best-first (most aligned with unmet demand).
    scores : (n_candidates,) float array
        Alignment of each candidate with the aggregated separator direction
        (cosine), in candidate order (NOT sorted).
    aggregate_direction : (n_genes,) float array
        The residual-weighted sum of unit separators — the library's "unmet
        demand" direction. Zero vector if no supplied combination separates.
    n_separating : int
        How many of the supplied combinations were outside the cone and thus
        contributed a separator.
    scope : str
        Model-relative disclaimer.
    """

    ranking: np.ndarray
    scores: np.ndarray
    aggregate_direction: np.ndarray
    n_separating: int
    scope: str = _SCOPE

    def top(self, k: int = 1) -> np.ndarray:
        """Indices of the ``k`` best-ranked candidates."""
        return self.ranking[:k]


def nominate_atoms(
    cone_atoms: np.ndarray,
    measured_combos: np.ndarray,
    candidate_atoms: np.ndarray,
    *,
    weight: str = "residual",
    gene_weights: np.ndarray | None = None,
    separator_tolerance: float | None = None,
) -> AtomNomination:
    """Rank which new perturbation to add to a library, from infeasibility certificates.

    Given a library (``cone_atoms``) and a set of ``measured_combos`` that the
    library cannot additively reach, each combo's cone projection yields a
    model-relative separator: a direction the library is missing. This aggregates
    those separators (residual-weighted by default) into one "unmet demand"
    direction, then ranks ``candidate_atoms`` by how well each aligns with it.

    The nominated candidate is the perturbation whose effect most closely supplies
    the axis the library lacks. This is the one operation a forward predictor cannot
    perform: it is a statement about the *library*, derived from the certificate of
    what the library cannot represent, not a prediction about any single input.

    Parameters
    ----------
    cone_atoms : (n_atoms, n_genes) array
        The current library (the cone whose reach is being tested).
    measured_combos : (m, n_genes) array
        Measured effect vectors tested against the cone. Combinations that fall
        inside the cone contribute no separator and are skipped.
    candidate_atoms : (n_candidates, n_genes) array
        Pool of candidate perturbation effects to rank.
    weight : str
        Separator aggregation weight: ``"residual"`` (weight each unit separator by
        its unreachable fraction; default) or ``"uniform"``.
    gene_weights : (n_genes,) array, optional
        Per-gene metric weights passed through to the cone projection.
    separator_tolerance : float, optional
        Passed through to :func:`reachability.project_cone`.

    Returns
    -------
    AtomNomination
    """
    cone_atoms = _validate_atoms(cone_atoms, "cone_atoms")
    measured_combos = _validate_atoms(measured_combos, "measured_combos")
    candidate_atoms = _validate_atoms(candidate_atoms, "candidate_atoms")
    n_genes = cone_atoms.shape[1]
    if measured_combos.shape[1] != n_genes or candidate_atoms.shape[1] != n_genes:
        raise rx.InputError(
            "cone_atoms, measured_combos and candidate_atoms must share n_genes"
        )
    if weight not in ("residual", "uniform"):
        raise rx.InputError("weight must be 'residual' or 'uniform'")

    aggregate = np.zeros(n_genes, dtype=float)
    n_sep = 0
    for combo in measured_combos:
        pr = rx.project_cone(
            cone_atoms,
            combo,
            gene_weights=gene_weights,
            separator_tolerance=separator_tolerance,
        )
        if pr.dual_separator is None:
            continue  # inside the cone — no unmet demand from this combo
        w = pr.residual_fraction if weight == "residual" else 1.0
        aggregate += w * _unit(np.asarray(pr.dual_separator, dtype=float))
        n_sep += 1

    scores = _cos_rows(candidate_atoms, aggregate)
    # Deterministic tie-break by candidate index.
    ranking = np.array(
        sorted(range(candidate_atoms.shape[0]), key=lambda i: (-scores[i], i)),
        dtype=int,
    )
    return AtomNomination(
        ranking=ranking,
        scores=scores,
        aggregate_direction=aggregate,
        n_separating=n_sep,
    )


def held_out_single_recovery(
    atoms: np.ndarray,
    combo_atom_idx: np.ndarray,
    combo_effects: np.ndarray,
    *,
    min_involved: int = 2,
    weight: str = "residual",
    gene_weights: np.ndarray | None = None,
) -> dict:
    """Falsifiable test of :func:`nominate_atoms`: recover a hidden library axis.

    For each single-gene atom that participates in at least ``min_involved``
    measured combinations, remove it from the library, and use the separators of
    exactly the combinations that involved it to nominate a replacement from the
    full atom pool. If the certificate captures the missing axis, the true held-out
    atom should rank at or near the top.

    Two rankers are reported for every held-out atom:
      * ``sep`` — :func:`nominate_atoms` (the aggregated separator).
      * ``base`` — a naive control: rank candidates by cosine to the *mean* of the
        combos that involved the held-out atom (no cone, no separator). This is the
        obvious "the missing gene's effect dominates its own combinations" shortcut;
        beating it is what shows the separator does non-trivial work.

    Parameters
    ----------
    atoms : (n_atoms, n_genes) array
        Single-gene effect vectors (the full candidate pool).
    combo_atom_idx : (n_combos, k) int array
        Constituent atom indices for each combination.
    combo_effects : (n_combos, n_genes) array
        Measured combination effects.
    min_involved : int
        Minimum number of combinations an atom must appear in to be eligible.
    weight : str
        Passed to :func:`nominate_atoms`.
    gene_weights : (n_genes,) array, optional
        Per-gene metric weights.

    Returns
    -------
    dict
        ``eligible`` (atom indices tested), ``sep_ranks`` / ``base_ranks`` (1-based
        rank of the true atom under each ranker, aligned to ``eligible``), and
        summary top-1 / top-5 / median-rank for both, plus ``scope``.
    """
    atoms = _validate_atoms(atoms)
    idx = _validate_combo_idx(combo_atom_idx, atoms.shape[0])
    combo_effects = _validate_atoms(combo_effects, "combo_effects")
    if idx.shape[0] != combo_effects.shape[0]:
        raise rx.InputError("combo_atom_idx and combo_effects must share n_combos")
    if atoms.shape[1] != combo_effects.shape[1]:
        raise rx.InputError("atoms and combo_effects must share n_genes")

    n_atoms = atoms.shape[0]
    involve: dict[int, list[int]] = {g: [] for g in range(n_atoms)}
    for k, row in enumerate(idx):
        for g in row:
            involve[int(g)].append(k)
    eligible = [g for g in range(n_atoms) if len(involve[g]) >= min_involved]
    if not eligible:
        raise rx.InputError(
            f"no atom appears in >= {min_involved} combinations"
        )

    def rank_of(scores: np.ndarray, g: int) -> int:
        order = np.argsort(-scores, kind="stable")
        return int(np.where(order == g)[0][0]) + 1

    sep_ranks, base_ranks = [], []
    for g in eligible:
        keep = np.array([i for i in range(n_atoms) if i != g])
        cone_g = atoms[keep]
        combos_g = combo_effects[involve[g]]
        nom = nominate_atoms(
            cone_g, combos_g, atoms, weight=weight, gene_weights=gene_weights
        )
        sep_ranks.append(rank_of(nom.scores, g))
        base_dir = combo_effects[involve[g]].sum(axis=0)
        base_scores = _cos_rows(atoms, base_dir)
        base_ranks.append(rank_of(base_scores, g))

    sep = np.array(sep_ranks, dtype=int)
    base = np.array(base_ranks, dtype=int)

    def summ(r: np.ndarray) -> dict:
        return {
            "median_rank": float(np.median(r)),
            "mean_rank": float(r.mean()),
            "top1": float((r == 1).mean()),
            "top5": float((r <= 5).mean()),
        }

    return {
        "eligible": np.array(eligible, dtype=int),
        "sep_ranks": sep,
        "base_ranks": base,
        "n_candidates": n_atoms,
        "separator": summ(sep),
        "baseline": summ(base),
        "scope": _SCOPE,
    }


def _demo() -> None:
    """Tiny synthetic self-check (runs under ``python screenloop.py``)."""
    rng = np.random.default_rng(0)
    n_genes = 40
    atoms = rng.normal(size=(8, n_genes))
    # Build doubles: mostly additive, a few with a planted emergent axis.
    idx = np.array([(i, j) for i in range(8) for j in range(i + 1, 8)])[:15]
    combos = atoms[idx].sum(axis=1)
    labels = np.zeros(len(idx), dtype=bool)
    extra = np.zeros(n_genes)
    extra[0] = 8.0  # an axis no single carries
    for k in (2, 5, 9):
        combos[k] = combos[k] + extra
        labels[k] = True

    camp = replay_campaign(atoms, combos, idx, labels, policy="triage", batch_size=3)
    print("campaign wells:", camp.wells.tolist())
    print("campaign found:", camp.found.tolist())
    print("wells to 100%:", camp.wells_to_fraction(1.0))

    rec = held_out_single_recovery(atoms, idx, combos, min_involved=2)
    print("held-out recovery separator top1:", rec["separator"]["top1"])
    print("scope:", _SCOPE)


if __name__ == "__main__":
    _demo()
