"""combicone: certified triage for combinatorial perturbation screens.

A thin, fail-closed layer over :mod:`reachability`. It answers two questions a
combinatorial screener actually pays for:

1. **Prospective triage** — given only the *single*-perturbation effects, which of
   the O(N^2) unmeasured combinations are most likely to be *emergent* (produce a
   cell state the singles cannot reach)?  :func:`triage_combinations`.

2. **Certified emergence** — given a *measured* combination, is its departure from
   the single-gene cone real, or is it within measurement noise?  Returns a
   certificate: the model-relative separator, a bootstrap CI on the unreachable
   fraction, and a noise-injection p-value.  :func:`certify_emergence`.

Design contract (inherited from reachability.py and NOT relaxed here):
  * ``effects`` / ``singles`` matrices are ``(n_perturbations, n_genes)``.
  * All geometry is *model-relative*: "unreachable" means "outside the non-negative
    cone of THESE measured effects under THIS metric", never "biologically
    impossible". Every public result echoes that scope.
  * Deterministic given a seed. No hidden global state.

Honest performance envelope (Norman combinatorial CRISPRa screen, file label
A549 / canonically K562; 105 single-gene atoms, 131 measured doubles):
  * The default training-free triage score is ``-cos(effA, effB)``. It enriches the
    top-20 picks 2.4x over base rate against the raw unreachable-fraction label,
    but only 1.4x against the noise-robust (magnitude-controlled) emergence label.
    The difference is exactly the signal-to-noise inflation that
    :func:`certify_emergence` guards against. (The optional ``use_gap`` variant
    scored slightly worse: 2.0x raw / 1.2x noise-robust.)
  * A ridge model fitted on a pilot labeled screen (:func:`fit_triage_model`)
    reaches ~2.4x against the noise-robust label (LOO-CV Spearman 0.43, perm
    p=0.002). Use it when you can afford a labeled pilot; otherwise the
    training-free ``-cos`` score is the honest default.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np

import reachability as rx

__all__ = [
    "EmergenceCertificate",
    "TriageResult",
    "triage_combinations",
    "certify_emergence",
    "fit_triage_model",
    "single_effect_cosine",
    "additive_gap",
]

_EPS = 1e-12
_SCOPE = (
    "model-relative: unreachable = outside the non-negative cone of the supplied "
    "single-gene effects under the chosen metric; not a claim of biological "
    "impossibility"
)


# --------------------------------------------------------------------------- #
# Small, individually testable singles-only features
# --------------------------------------------------------------------------- #
def single_effect_cosine(effect_a: np.ndarray, effect_b: np.ndarray) -> float:
    """Cosine between two single-gene effect vectors.

    Low cosine (the two singles push transcription in different directions) is the
    strongest training-free predictor of combinatorial emergence found on Norman
    (Spearman -0.47 vs raw unreachable fraction, -0.37 vs the noise-robust label).
    """
    a = np.asarray(effect_a, dtype=float)
    b = np.asarray(effect_b, dtype=float)
    na = rx._stable_norm(a)
    nb = rx._stable_norm(b)
    if na == 0.0 or nb == 0.0:
        return 0.0
    return float(np.dot(a / na, b / nb))


def additive_gap(
    singles: np.ndarray,
    index_a: int,
    index_b: int,
    *,
    gene_weights: np.ndarray | None = None,
) -> float:
    """Certified leave-the-pair-out residual of the additive prediction A+B.

    Projects the *additive* vector (effect_a + effect_b) onto the cone of all OTHER
    single-gene effects (the pair itself removed) and returns the unreachable
    fraction. High gap = even the naive additive combination already leaves the
    rest-of-library cone, a cheap prospective flag. Spearman +0.35 vs raw label.
    """
    singles = np.asarray(singles, dtype=float)
    n = singles.shape[0]
    if not (0 <= index_a < n and 0 <= index_b < n) or index_a == index_b:
        raise rx.InputError("index_a and index_b must be distinct valid rows")
    keep = [i for i in range(n) if i not in (index_a, index_b)]
    if not keep:
        raise rx.InputError(
            "additive_gap needs at least 3 singles (the pair plus one other atom "
            "to form the leave-pair-out cone)"
        )
    additive = singles[index_a] + singles[index_b]
    result = rx.project_cone(singles[keep], additive, gene_weights=gene_weights)
    return float(result.residual_fraction)


# --------------------------------------------------------------------------- #
# Prospective triage
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TriageResult:
    """Ranked prospective triage of candidate combinations.

    Attributes
    ----------
    pairs : list[tuple[str, str]]
        Candidate combinations, in the input order.
    score : np.ndarray
        Prospective emergence score, higher = more likely emergent. Rank-only;
        not calibrated to a probability.
    rank : np.ndarray
        1-based rank of each pair by descending score (1 = run first).
    cos_ab : np.ndarray
        Single-single cosine per pair (diagnostic).
    gap : np.ndarray
        Leave-pair-out additive gap per pair (diagnostic).
    model : str
        "training-free" or "ridge".
    scope : str
        The model-relative scope disclaimer.
    """

    pairs: list[tuple[str, str]]
    score: np.ndarray
    rank: np.ndarray
    cos_ab: np.ndarray
    gap: np.ndarray
    model: str
    scope: str = _SCOPE

    def top(self, k: int) -> list[tuple[str, str]]:
        """The k pairs to run first (highest score)."""
        order = np.argsort(-self.score, kind="stable")[:k]
        return [self.pairs[i] for i in order]


def _resolve_singles(
    singles: np.ndarray,
    single_names: Sequence[str],
) -> tuple[np.ndarray, dict[str, int]]:
    singles = np.asarray(singles, dtype=float)
    if singles.ndim != 2:
        raise rx.InputError("singles must be 2D (n_singles, n_genes)")
    if len(single_names) != singles.shape[0]:
        raise rx.InputError("single_names length must match singles rows")
    index = {name: i for i, name in enumerate(single_names)}
    if len(index) != len(single_names):
        raise rx.InputError("single_names must be unique")
    return singles, index


def triage_combinations(
    singles: np.ndarray,
    single_names: Sequence[str],
    candidate_pairs: Iterable[tuple[str, str]] | None = None,
    *,
    gene_weights: np.ndarray | None = None,
    model: "TriageModel | None" = None,
    use_gap: bool = False,
) -> TriageResult:
    """Rank unmeasured combinations by predicted emergence, using singles only.

    Parameters
    ----------
    singles : (n_singles, n_genes) array
        Measured single-gene effect vectors (e.g. pseudobulk mean minus control).
    single_names : sequence of str
        Row labels for ``singles``; used to name candidate pairs.
    candidate_pairs : iterable of (name_a, name_b), optional
        Combinations to score. Defaults to ALL distinct unordered pairs.
    gene_weights : (n_genes,) array, optional
        Per-gene metric weights passed through to the cone projection.
    model : TriageModel, optional
        A fitted ridge model from :func:`fit_triage_model`. If given, its learned
        score is used; otherwise a training-free score is used.
    use_gap : bool
        If False (default), the training-free score is ``-cos(effA, effB)`` alone
        — the strongest single validated feature on Norman. If True, the score is
        a rank-average of low cosine and high leave-pair-out additive gap; this is
        more expensive (an O(n) projection per pair) and scored slightly *worse* on
        Norman (2.0x vs 2.4x raw-label enrichment), so it is off by default.
        Always computed (and required) when ``model`` is supplied.

    Returns
    -------
    TriageResult

    Notes
    -----
    The default training-free score is ``-cos(effA, effB)``. On Norman it enriches
    the top-20 picks 2.4x over base rate against a raw unreachable-fraction label
    and 1.4x against the noise-robust label. Enabling ``use_gap`` scored slightly
    worse (2.0x raw / 1.2x noise-robust). For the noise-robust target specifically,
    prefer a fitted :class:`TriageModel` (LOO-CV Spearman 0.43, ~2.4x).
    """
    singles, index = _resolve_singles(singles, single_names)
    if candidate_pairs is None:
        names = list(single_names)
        candidate_pairs = [
            (names[i], names[j])
            for i in range(len(names))
            for j in range(i + 1, len(names))
        ]
    pairs: list[tuple[str, str]] = []
    cos_ab: list[float] = []
    gap: list[float] = []
    for a, b in candidate_pairs:
        if a not in index or b not in index:
            raise rx.InputError(f"unknown single in pair ({a!r}, {b!r})")
        ia, ib = index[a], index[b]
        pairs.append((a, b))
        cos_ab.append(single_effect_cosine(singles[ia], singles[ib]))
        if use_gap or model is not None:
            gap.append(additive_gap(singles, ia, ib, gene_weights=gene_weights))
        else:
            gap.append(np.nan)
    cos_ab_arr = np.asarray(cos_ab, dtype=float)
    gap_arr = np.asarray(gap, dtype=float)

    if model is not None:
        feats = np.column_stack([cos_ab_arr, gap_arr])
        score = model.predict(feats)
        model_name = "ridge"
    elif use_gap:
        # rank-average of (low cosine) and (high gap); rank-based so the two
        # differently-scaled features combine without a fitted weight.
        r_cos = _rankdata(-cos_ab_arr)
        r_gap = _rankdata(gap_arr)
        score = (r_cos + r_gap) / 2.0
        model_name = "training-free"
    else:
        score = -cos_ab_arr
        model_name = "training-free"

    order = np.argsort(-score, kind="stable")
    rank = np.empty(len(score), dtype=int)
    rank[order] = np.arange(1, len(score) + 1)
    return TriageResult(
        pairs=pairs,
        score=np.asarray(score, dtype=float),
        rank=rank,
        cos_ab=cos_ab_arr,
        gap=gap_arr,
        model=model_name,
    )


def _rankdata(x: np.ndarray) -> np.ndarray:
    """Average-rank of x (ties averaged), dependency-free."""
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="stable")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1, dtype=float)
    # average ties
    _, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    return (sums / counts)[inv]


# --------------------------------------------------------------------------- #
# Certified emergence of a MEASURED combination
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class EmergenceCertificate:
    """Certificate that a measured combination is (or is not) emergent.

    Attributes
    ----------
    unreachable_fraction : float
        Observed residual fraction of the measured combination against the
        single-gene cone (rx.ProjectionResult.residual_fraction).
    geometry_status : str
        "outside_model_cone" or "inside_tolerance", from the certified projection.
    separator : np.ndarray | None
        The model-relative dual separator (None when inside the cone). Certifies
        the unreachable direction.
    ci_low, ci_high : float
        Bootstrap confidence interval on ``unreachable_fraction`` under resampled
        measurement noise (percentile CI at ``ci_level``).
    noise_null_mean : float
        Mean residual a *reachable* target of the same magnitude would show under
        the supplied measurement noise (the null the observation is tested against).
    z : float
        (unreachable_fraction - noise_null_mean) / noise_null_sd. Noise-aware
        emergence strength; magnitude-robust.
    p_value : float
        Noise-injection p-value: P(noise-only residual >= observed). Small = the
        unreachable fraction exceeds what measurement noise alone produces.
    floor_ratio : float
        unreachable_fraction / noise_floor, where noise_floor is the residual
        expected from noise. > ~1.9 clears the stricter effect-size bar; a value
        near 1 means the raw residual is essentially noise (the SNR trap).
    verdict : str
        Human-readable summary combining the p-value and the floor ratio.
    scope : str
        Model-relative scope disclaimer.
    """

    unreachable_fraction: float
    geometry_status: str
    separator: np.ndarray | None
    ci_low: float
    ci_high: float
    noise_null_mean: float
    noise_null_sd: float
    z: float
    p_value: float
    floor_ratio: float
    verdict: str
    scope: str = _SCOPE


def certify_emergence(
    singles: np.ndarray,
    measured_combo: np.ndarray,
    *,
    noise_sd: np.ndarray | float | None = None,
    gene_weights: np.ndarray | None = None,
    n_boot: int = 200,
    ci_level: float = 0.9,
    floor_threshold: float = 1.9,
    alpha: float = 0.05,
    seed: int = 0,
) -> EmergenceCertificate:
    """Certify whether a measured combination departs from the single-gene cone.

    The observed departure is tested against the honest null "the target IS
    reachable, and the residual we see is measurement noise". We form that null by
    projecting the measured combo onto the cone (its best reachable approximation
    ``f0``), then repeatedly adding fresh per-gene Gaussian noise of scale
    ``noise_sd`` to ``f0`` and re-projecting. If the observed residual sits well
    above this noise-only distribution, the emergence is certified.

    Parameters
    ----------
    singles : (n_singles, n_genes) array
        Measured single-gene effects (the cone).
    measured_combo : (n_genes,) array
        The measured combination effect vector to test.
    noise_sd : (n_genes,) array or float, optional
        Per-gene standard error of ``measured_combo``. The recommended estimate is
        ``|t1 - t2| / 2`` from a random split-half of the cells. If a scalar, used
        for all genes. If None, the noise test is skipped (p_value / z / floor are
        NaN) and only the point certificate + geometry are returned.
    gene_weights : (n_genes,) array, optional
        Per-gene metric weights passed to the projection.
    n_boot : int
        Number of noise draws for the null and the CI.
    ci_level : float
        Central mass of the percentile bootstrap CI (e.g. 0.9 -> 5th..95th pct).
    floor_threshold : float
        floor_ratio above which the stricter effect-size bar is considered cleared.
    alpha : float
        Significance level for the p-value verdict.
    seed : int
        RNG seed (determinism).

    Returns
    -------
    EmergenceCertificate
    """
    singles = np.asarray(singles, dtype=float)
    target = np.asarray(measured_combo, dtype=float)
    if singles.ndim != 2 or target.ndim != 1:
        raise rx.InputError("singles must be 2D and measured_combo 1D")
    if singles.shape[1] != target.shape[0]:
        raise rx.InputError("gene axis mismatch between singles and measured_combo")

    obs = rx.project_cone(singles, target, gene_weights=gene_weights)
    obs_resid = float(obs.residual_fraction)

    if noise_sd is None:
        return EmergenceCertificate(
            unreachable_fraction=obs_resid,
            geometry_status=obs.geometry_status,
            separator=obs.dual_separator,
            ci_low=float("nan"),
            ci_high=float("nan"),
            noise_null_mean=float("nan"),
            noise_null_sd=float("nan"),
            z=float("nan"),
            p_value=float("nan"),
            floor_ratio=float("nan"),
            verdict="point certificate only (no noise model supplied)",
        )

    se = np.asarray(noise_sd, dtype=float)
    if se.ndim == 0:
        se = np.full(target.shape[0], float(se))
    if se.shape != target.shape:
        raise rx.InputError("noise_sd must be scalar or match the gene axis")
    if np.any(se < 0) or not np.all(np.isfinite(se)):
        raise rx.InputError("noise_sd must be finite and non-negative")

    rng = np.random.default_rng(seed)
    f0 = obs.fitted  # reachable null truth
    null = np.empty(n_boot)
    boot_obs = np.empty(n_boot)
    for i in range(n_boot):
        null[i] = rx.project_cone(
            singles, f0 + rng.normal(0.0, se), gene_weights=gene_weights
        ).residual_fraction
        # CI: resample the OBSERVED target's noise around itself
        boot_obs[i] = rx.project_cone(
            singles, target + rng.normal(0.0, se), gene_weights=gene_weights
        ).residual_fraction

    null_mean = float(null.mean())
    null_sd = float(null.std())
    z = (obs_resid - null_mean) / (null_sd + _EPS)
    # conservative plus-one p-value (reuses the engine's convention)
    p = rx.empirical_p(obs_resid, null)
    lo = float(np.quantile(boot_obs, (1 - ci_level) / 2))
    hi = float(np.quantile(boot_obs, 1 - (1 - ci_level) / 2))
    floor_ratio = obs_resid / (null_mean + _EPS)

    sig = p < alpha
    clears_floor = floor_ratio >= floor_threshold
    if sig and clears_floor:
        verdict = f"certified emergent (p={p:.3g}, {floor_ratio:.1f}x noise floor)"
    elif sig:
        verdict = (
            f"emergent above noise but modest effect size "
            f"(p={p:.3g}, only {floor_ratio:.1f}x noise floor)"
        )
    else:
        verdict = f"within measurement noise (p={p:.3g}); do not call emergent"

    return EmergenceCertificate(
        unreachable_fraction=obs_resid,
        geometry_status=obs.geometry_status,
        separator=obs.dual_separator,
        ci_low=lo,
        ci_high=hi,
        noise_null_mean=null_mean,
        noise_null_sd=null_sd,
        z=float(z),
        p_value=float(p),
        floor_ratio=float(floor_ratio),
        verdict=verdict,
    )


# --------------------------------------------------------------------------- #
# Optional trained triage model (needs a labeled pilot screen)
# --------------------------------------------------------------------------- #
@dataclass
class TriageModel:
    """A tiny ridge model over singles-only features (cos_ab, gap).

    Fitted on a pilot screen where some combinations have been measured and
    labeled with an emergence score (ideally the noise-robust ``z``). Interpretable
    by design: two features, closed-form ridge. Not a deep model.
    """

    coef_: np.ndarray
    intercept_: float
    mean_: np.ndarray
    std_: np.ndarray
    alpha: float

    def predict(self, features: np.ndarray) -> np.ndarray:
        features = np.asarray(features, dtype=float)
        z = (features - self.mean_) / self.std_
        return z @ self.coef_ + self.intercept_


def fit_triage_model(
    singles: np.ndarray,
    single_names: Sequence[str],
    labeled_pairs: Mapping[tuple[str, str], float],
    *,
    gene_weights: np.ndarray | None = None,
    alpha: float = 10.0,
) -> TriageModel:
    """Fit the ridge triage model on a labeled pilot screen.

    Parameters
    ----------
    singles, single_names : as in :func:`triage_combinations`.
    labeled_pairs : mapping (name_a, name_b) -> emergence label
        Measured combinations with a known emergence score (prefer the noise-robust
        ``z`` from :func:`certify_emergence`).
    gene_weights : optional per-gene metric weights.
    alpha : ridge penalty (stable across [1, 100] on Norman).

    Returns
    -------
    TriageModel

    Notes
    -----
    On Norman this reaches LOO-CV Spearman ~0.43 against the noise-robust label
    (permutation p=0.002) and top-20 precision ~2.4x base rate. Use it only when a
    labeled pilot is available; otherwise the training-free score in
    :func:`triage_combinations` is the honest default.
    """
    singles, index = _resolve_singles(singles, single_names)
    feats = []
    y = []
    for (a, b), label in labeled_pairs.items():
        if a not in index or b not in index:
            raise rx.InputError(f"unknown single in labeled pair ({a!r}, {b!r})")
        ia, ib = index[a], index[b]
        feats.append(
            [
                single_effect_cosine(singles[ia], singles[ib]),
                additive_gap(singles, ia, ib, gene_weights=gene_weights),
            ]
        )
        y.append(float(label))
    X = np.asarray(feats, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(y) < 3:
        raise rx.InputError("need at least 3 labeled pairs to fit")
    mean = X.mean(axis=0)
    std = X.std(axis=0)
    std[std == 0] = 1.0
    Xz = (X - mean) / std
    # closed-form ridge on standardized features, intercept = mean(y)
    b = y - y.mean()
    G = Xz.T @ Xz + alpha * np.eye(Xz.shape[1])
    coef = np.linalg.solve(G, Xz.T @ b)
    return TriageModel(
        coef_=coef,
        intercept_=float(y.mean()),
        mean_=mean,
        std_=std,
        alpha=alpha,
    )


def _demo() -> None:
    """Deterministic, data-free smoke run of both public entry points."""
    # A small 4-atom library. Triage ranks all C(4,2) candidate combinations from
    # the singles alone; the certificate then tests two measured "combinations".
    singles = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],  # A
            [0.0, 1.0, 0.0, 0.0],  # B
            [0.0, 0.0, 1.0, 0.0],  # C
            [0.9, 0.1, 0.0, 0.0],  # D (near-collinear with A)
        ]
    )
    names = ["A", "B", "C", "D"]

    triage = triage_combinations(singles, names)
    print("triage model:", triage.model)
    print("run-first pair:", triage.top(1)[0])
    print("candidate pairs:", triage.pairs)
    print("triage score:", np.round(triage.score, 3).tolist())

    emergent = np.array([1.0, 1.0, 0.0, 3.0])  # big unreachable 4th-axis component
    additive = np.array([1.0, 1.0, 0.5, 0.0])  # inside the cone of A,B,C
    for label, combo in [("emergent", emergent), ("additive", additive)]:
        cert = certify_emergence(singles, combo, noise_sd=0.05, n_boot=200, seed=0)
        print(
            f"[{label}] unreachable={cert.unreachable_fraction:.3f} "
            f"z={cert.z:.1f} p={cert.p_value:.3g} "
            f"floor_ratio={cert.floor_ratio:.2f} -> {cert.verdict}"
        )
    print("scope:", _SCOPE)


if __name__ == "__main__":
    _demo()
