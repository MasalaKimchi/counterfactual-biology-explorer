"""neural_baseline: a self-contained learned interaction baseline for CombiCone.

This module answers one question honestly: *does a genuine nonlinear learned model
predict a measured double-perturbation effect better than the trivial additive
baseline, and how does it relate to the reachability cone's reachable fit?*

It is a **baseline for comparison**, not part of the CombiCone engine. It imports
:mod:`reachability` and :mod:`combicone` as READ-ONLY dependencies and never
mutates them.

What it builds
--------------
For each measured double A+B the inputs are the two constituent single-gene effect
atoms ``a`` and ``b``; the target is the measured double effect ``t``. Genes are
reduced to a compact space by an SVD fit on the single-gene atom library (default
50 components, ~95.6% of atom-span variance on Norman). A small 2-layer MLP with a
``tanh`` hidden layer maps interaction features of ``(a, b)`` — the concatenation
``[a, b]`` plus the elementwise product ``a*b`` so pairwise interactions are
representable — to the double effect. The MLP is trained with manual numpy
backprop optimised by ``scipy.optimize.minimize(method="L-BFGS-B")``. It is
deterministic given a seed.

The learner is given its **best fair shot**: by default it is *additive-anchored*
(``anchor=True``), i.e. it predicts ``(a+b) + correction`` and only has to learn a
nonlinear interaction *correction* on top of the additive prediction. This can
only help the learner relative to additivity, so if it still fails to beat the
additive baseline, the conclusion "learned nonlinearity buys nothing here" is
conservative.

Three predictions are compared per double, all scored by cosine to the measured
double in full gene space:

  * **learned MLP**  — blind: never sees the held-out double's target.
    Evaluated leave-pairs-out (train on the other doubles, predict the held-out).
  * **additive** ``a+b`` — blind, training-free.
  * **cone reachable fit** ``rx.project_cone(atoms, t).fitted`` — an in-sample
    RECONSTRUCTION that *requires having measured* ``t`` (it fits non-negative
    coefficients to the observed target). It is NOT a blind predictor and cannot
    rank unmeasured combinations; it is included only to show the reachable
    approximation quality.

Honest scope / caveats (audited)
--------------------------------
  * This module measures **prediction accuracy only**. Accuracy is NOT the
    cone's contribution. On Norman a learned MLP given its best fair shot ties
    the additive baseline for blind prediction, reproducing the field result
    (Ahlmann-Eltze, Huber & Anders, *Nat. Methods* 2025,
    doi:10.1038/s41592-025-02772-6, "Deep-learning-based gene perturbation effect
    prediction does not yet outperform simple linear baselines") that learned models
    ~= simple linear/additive on combinatorial perturbation prediction.
  * The cone's higher cosine is an in-sample fit, not a forecast; comparing it
    to the blind MLP as if both were predictors would be dishonest, so the cone
    number is always labelled a reconstruction here.
  * The defensible, unique capability of the reachability cone is the
    **CERTIFICATE** — a model-relative dual separator plus a noise-injection
    test (:func:`combicone.certify_emergence`) — NOT prediction accuracy. Neither
    the MLP nor the additive baseline emits any certificate of unreachability.
  * All geometry is model-relative (outside the non-negative cone of THESE
    measured single effects under THIS metric), never a claim of biological
    impossibility.
  * Data provenance: the Norman combinatorial CRISPRa screen file is labelled
    cell_type A549 but Norman 2019 is canonically K562; reported as
    "file-label-A549 / canonically-K562".

Dependencies: numpy + scipy only (no torch, sklearn). Runs in the pinned CI env
(numpy==1.26.4, scipy==1.13.1).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import minimize

import reachability as rx  # READ-ONLY dependency

__all__ = [
    "cosine_full",
    "svd_reduce",
    "interaction_features",
    "MLPParams",
    "InteractionMLP",
    "additive_prediction",
    "prediction_error_score",
    "cone_reachable_fit_cosine",
    "leave_pairs_out",
    "LPOResult",
]

_SCOPE = (
    "prediction-accuracy comparison only; the reachability cone's contribution is "
    "the emergence CERTIFICATE (separator + noise test), NOT prediction accuracy. "
    "'cone reachable fit' is an in-sample reconstruction that requires the "
    "measurement, not a blind predictor. Model-relative geometry throughout."
)


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #
def cosine_full(x: np.ndarray, y: np.ndarray) -> float:
    """Plain (unweighted) cosine in full gene space; 0 if either vector is zero."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx = rx._stable_norm(x)
    ny = rx._stable_norm(y)
    if nx == 0.0 or ny == 0.0:
        return 0.0
    return float(np.dot(x / nx, y / ny))


def svd_reduce(atoms: np.ndarray, n_comp: int = 50) -> tuple[np.ndarray, dict]:
    """Fit an orthonormal SVD basis on the single-gene atom library.

    Parameters
    ----------
    atoms : (n_atoms, n_genes) array
        Single-gene effect atoms. The reduction is defined by the SINGLES only, so
        it is identical across leave-pairs-out folds (holding out a double never
        changes the atom library) — no target leakage.
    n_comp : int
        Number of right-singular vectors to keep.

    Returns
    -------
    basis : (n_comp, n_genes) array
        Orthonormal rows; project a gene vector ``v`` via ``v @ basis.T`` and
        reconstruct via ``coords @ basis``.
    info : dict
        ``{"n_comp", "var_retained"}`` where ``var_retained`` is the fraction of
        atom-span squared singular value mass captured.
    """
    atoms = np.asarray(atoms, dtype=float)
    if atoms.ndim != 2:
        raise rx.InputError("atoms must be 2D (n_atoms, n_genes)")
    max_rank = min(atoms.shape)
    if not (1 <= n_comp <= max_rank):
        raise rx.InputError(f"n_comp must be in [1, {max_rank}]")
    _, s, vt = np.linalg.svd(atoms, full_matrices=False)
    var = s ** 2
    var_retained = float(np.sum(var[:n_comp]) / np.sum(var))
    return vt[:n_comp].copy(), {"n_comp": int(n_comp), "var_retained": var_retained}


def interaction_features(a_red: np.ndarray, b_red: np.ndarray) -> np.ndarray:
    """Concatenate ``[a, b, a*b]`` in reduced space so interactions are representable.

    Accepts 1D (single pair) or 2D (batch, n_comp) inputs; returns matching rank.
    """
    a_red = np.asarray(a_red, dtype=float)
    b_red = np.asarray(b_red, dtype=float)
    if a_red.shape != b_red.shape:
        raise rx.InputError("a_red and b_red must have the same shape")
    return np.concatenate([a_red, b_red, a_red * b_red], axis=-1)


# --------------------------------------------------------------------------- #
# The 2-layer interaction MLP (manual numpy backprop, L-BFGS-B optimiser)
# --------------------------------------------------------------------------- #
@dataclass
class MLPParams:
    """Flat learned parameters of a 2-layer tanh MLP plus its normalisation."""

    W1: np.ndarray
    b1: np.ndarray
    W2: np.ndarray
    b2: np.ndarray
    feat_mean: np.ndarray
    feat_std: np.ndarray


def _pack(w1, b1, w2, b2) -> np.ndarray:
    return np.concatenate([w1.ravel(), b1, w2.ravel(), b2])


def _unpack(theta, d_in, hidden, d_out):
    i = 0
    w1 = theta[i:i + d_in * hidden].reshape(d_in, hidden); i += d_in * hidden
    b1 = theta[i:i + hidden]; i += hidden
    w2 = theta[i:i + hidden * d_out].reshape(hidden, d_out); i += hidden * d_out
    b2 = theta[i:i + d_out]
    return w1, b1, w2, b2


def _loss_grad(theta, X, Y, hidden, l2):
    """Mean-squared-error loss + L2, with exact analytic gradient (tanh hidden)."""
    n, d_in = X.shape
    d_out = Y.shape[1]
    w1, b1, w2, b2 = _unpack(theta, d_in, hidden, d_out)
    z1 = X @ w1 + b1
    h1 = np.tanh(z1)
    yh = h1 @ w2 + b2
    r = yh - Y
    loss = 0.5 * np.sum(r * r) / n + 0.5 * l2 * (np.sum(w1 * w1) + np.sum(w2 * w2))
    dyh = r / n
    dw2 = h1.T @ dyh + l2 * w2
    db2 = dyh.sum(0)
    dh1 = dyh @ w2.T
    dz1 = dh1 * (1.0 - h1 * h1)
    dw1 = X.T @ dz1 + l2 * w1
    db1 = dz1.sum(0)
    return loss, _pack(dw1, db1, dw2, db2)


class InteractionMLP:
    """A small, deterministic 2-layer tanh MLP for double-effect prediction.

    Works entirely in the reduced SVD space. If ``anchor=True`` (default) the
    network predicts a *correction* added to the additive baseline ``a+b`` (a skip
    connection), so it can only improve on additivity — the learner's best fair
    shot. If ``anchor=False`` it predicts the double effect directly.

    numpy/scipy only, deterministic given ``seed``.
    """

    def __init__(self, hidden: int = 32, l2: float = 0.5, anchor: bool = True,
                 seed: int = 0, maxiter: int = 500):
        if hidden < 1:
            raise rx.InputError("hidden must be >= 1")
        if l2 < 0:
            raise rx.InputError("l2 must be non-negative")
        self.hidden = int(hidden)
        self.l2 = float(l2)
        self.anchor = bool(anchor)
        self.seed = int(seed)
        self.maxiter = int(maxiter)
        self.params_: MLPParams | None = None
        self.n_iter_: int | None = None

    def fit(self, A_red: np.ndarray, B_red: np.ndarray, T_red: np.ndarray) -> "InteractionMLP":
        """Train on reduced atom pairs -> reduced double effects.

        Parameters use reduced (SVD-space) coordinates: ``A_red``/``B_red`` are
        ``(n, n_comp)`` atom coords, ``T_red`` is ``(n, n_comp)`` double coords.
        """
        A_red = np.asarray(A_red, dtype=float)
        B_red = np.asarray(B_red, dtype=float)
        T_red = np.asarray(T_red, dtype=float)
        X = interaction_features(A_red, B_red)
        base = (A_red + B_red) if self.anchor else np.zeros_like(T_red)
        Y = T_red - base
        mean = X.mean(0)
        std = X.std(0)
        std[std == 0] = 1.0
        Xz = (X - mean) / std
        d_in, d_out = Xz.shape[1], Y.shape[1]
        rng = np.random.default_rng(self.seed)
        theta0 = _pack(
            rng.standard_normal((d_in, self.hidden)) * 0.05, np.zeros(self.hidden),
            rng.standard_normal((self.hidden, d_out)) * 0.05, np.zeros(d_out),
        )
        res = minimize(_loss_grad, theta0, args=(Xz, Y, self.hidden, self.l2),
                       jac=True, method="L-BFGS-B", options={"maxiter": self.maxiter})
        w1, b1, w2, b2 = _unpack(res.x, d_in, self.hidden, d_out)
        self.params_ = MLPParams(w1, b1, w2, b2, mean, std)
        self.n_iter_ = int(res.nit)
        return self

    def predict_reduced(self, A_red: np.ndarray, B_red: np.ndarray) -> np.ndarray:
        """Predict double effect(s) in reduced space. Returns (n, n_comp)."""
        if self.params_ is None:
            raise rx.InputError("model is not fitted")
        A_red = np.atleast_2d(np.asarray(A_red, dtype=float))
        B_red = np.atleast_2d(np.asarray(B_red, dtype=float))
        p = self.params_
        X = interaction_features(A_red, B_red)
        Xz = (X - p.feat_mean) / p.feat_std
        corr = np.tanh(Xz @ p.W1 + p.b1) @ p.W2 + p.b2
        base = (A_red + B_red) if self.anchor else 0.0
        return corr + base


# --------------------------------------------------------------------------- #
# Baseline predictions
# --------------------------------------------------------------------------- #
def additive_prediction(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The trivial additive baseline ``a + b`` (blind, training-free)."""
    return np.asarray(a, dtype=float) + np.asarray(b, dtype=float)


def prediction_error_score(held_out_cosine: np.ndarray) -> np.ndarray:
    """Turn a held-out prediction cosine into a candidate EMERGENCE ranker score.

    A natural bridge between the prediction track and the emergence-ranking track:
    if a combination is emergent, a model trained on the *other* combinations
    should predict it worse, so ``1 - held_out_cosine`` (higher = more surprising)
    is a candidate emergence signal, directly comparable to the columns in the
    project's ``benchmark_rankers`` table.

    HONEST CAVEAT (measured on Norman, 131 doubles): this score is essentially a
    re-expression of the raw unreachable fraction — Spearman +0.87 vs ``cone_raw``,
    +0.68 vs ``nonadditivity`` — and it inherits the SAME magnitude confound
    (Spearman -0.58 vs effect_norm, matching cone_raw's -0.56). It does NOT recover
    the noise-aware emergence label (Spearman only +0.11 vs ``cone_z``). In other
    words a learned model's "surprise" falls into the same signal-to-noise trap as
    every raw ranker and canNOT substitute for the noise-calibrated certificate
    (:func:`combicone.certify_emergence`). Use it as a ranker only with that caveat.
    """
    hc = np.asarray(held_out_cosine, dtype=float)
    return 1.0 - hc


def cone_reachable_fit_cosine(atoms: np.ndarray, target: np.ndarray,
                              *, gene_weights: np.ndarray | None = None) -> float:
    """Cosine of the cone's reachable fit to the target (IN-SAMPLE reconstruction).

    This calls ``rx.project_cone(atoms, target)`` which SEES ``target`` — it is a
    reconstruction quality, NOT a blind forecast, and cannot rank unmeasured
    combinations. Labelled as such everywhere it is reported.
    """
    pr = rx.project_cone(np.asarray(atoms, dtype=float),
                         np.asarray(target, dtype=float), gene_weights=gene_weights)
    return float(pr.cosine)


# --------------------------------------------------------------------------- #
# Leave-pairs-out evaluation
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class LPOResult:
    """Per-double held-out cosines for each method, plus summaries."""

    names: list[str]
    mlp: np.ndarray             # blind, leave-pairs-out
    additive: np.ndarray        # blind, training-free
    cone_fit: np.ndarray        # IN-SAMPLE reconstruction (sees target)
    n_comp: int
    var_retained: float
    hidden: int
    l2: float
    anchor: bool
    seed: int
    scope: str = _SCOPE

    def means(self) -> dict:
        return {
            "mlp_heldout": float(self.mlp.mean()),
            "additive": float(self.additive.mean()),
            "cone_fit_insample": float(self.cone_fit.mean()),
            "mlp_beats_additive_frac": float(np.mean(self.mlp > self.additive)),
        }


def leave_pairs_out(
    atoms: np.ndarray,
    a_index: Sequence[int],
    b_index: Sequence[int],
    targets: np.ndarray,
    names: Sequence[str],
    *,
    n_comp: int = 50,
    hidden: int = 32,
    l2: float = 0.5,
    anchor: bool = True,
    seed: int = 0,
    maxiter: int = 500,
    gene_weights: np.ndarray | None = None,
) -> LPOResult:
    """Leave-pairs-out evaluation of MLP vs additive vs cone reachable fit.

    For each double k: train the MLP on the other doubles (blind to k's target)
    and predict k; the additive baseline needs no training; the cone reachable fit
    is the in-sample reconstruction of k's own target. All cosines are full-gene.

    Parameters
    ----------
    atoms : (n_atoms, n_genes)
        Single-gene atom library (defines the SVD basis and the cone).
    a_index, b_index : sequences of int
        Row indices into ``atoms`` of each double's two constituent singles.
    targets : (n_doubles, n_genes)
        Measured double effects.
    names : sequence of str
        Double labels.
    """
    atoms = np.asarray(atoms, dtype=float)
    targets = np.asarray(targets, dtype=float)
    a_index = np.asarray(a_index, dtype=int)
    b_index = np.asarray(b_index, dtype=int)
    names = list(names)
    n = len(names)
    if not (len(a_index) == len(b_index) == targets.shape[0] == n):
        raise rx.InputError("a_index, b_index, targets, names must align")
    if atoms.shape[1] != targets.shape[1]:
        raise rx.InputError("atoms and targets gene axes must match")

    basis, info = svd_reduce(atoms, n_comp=n_comp)
    A = atoms[a_index]
    B = atoms[b_index]
    A_red = A @ basis.T
    B_red = B @ basis.T
    T_red = targets @ basis.T

    mlp_cos = np.empty(n)
    add_cos = np.empty(n)
    cone_cos = np.empty(n)
    for k in range(n):
        train = np.ones(n, dtype=bool)
        train[k] = False
        model = InteractionMLP(hidden=hidden, l2=l2, anchor=anchor, seed=seed,
                               maxiter=maxiter)
        model.fit(A_red[train], B_red[train], T_red[train])
        pred_red = model.predict_reduced(A_red[k], B_red[k])[0]
        pred_full = pred_red @ basis
        mlp_cos[k] = cosine_full(pred_full, targets[k])
        add_cos[k] = cosine_full(additive_prediction(A[k], B[k]), targets[k])
        cone_cos[k] = cone_reachable_fit_cosine(atoms, targets[k],
                                                gene_weights=gene_weights)
    return LPOResult(
        names=names, mlp=mlp_cos, additive=add_cos, cone_fit=cone_cos,
        n_comp=info["n_comp"], var_retained=info["var_retained"],
        hidden=hidden, l2=l2, anchor=anchor, seed=seed,
    )


# --------------------------------------------------------------------------- #
# Demo (data-free, deterministic) — mirrors reachability.py / combicone.py style
# --------------------------------------------------------------------------- #
def _demo() -> None:
    """Tiny synthetic smoke run: 3 atoms, 2 'doubles'. No data files."""
    rng = np.random.default_rng(0)
    n_genes = 40
    atoms = rng.standard_normal((3, n_genes))
    # two synthetic doubles as noisy additive combos of atom pairs
    a_idx = [0, 1]
    b_idx = [1, 2]
    targets = np.stack([
        atoms[0] + atoms[1] + 0.1 * rng.standard_normal(n_genes),
        atoms[1] + atoms[2] + 0.1 * rng.standard_normal(n_genes),
    ])
    names = ["A+B", "B+C"]
    res = leave_pairs_out(atoms, a_idx, b_idx, targets, names,
                          n_comp=3, hidden=8, l2=0.5, seed=0, maxiter=100)
    m = res.means()
    print(f"reduced to {res.n_comp} comps ({100*res.var_retained:.1f}% atom-span var)")
    print(f"held-out cosine  MLP={m['mlp_heldout']:.3f}  "
          f"additive={m['additive']:.3f}  cone-fit(in-sample)={m['cone_fit_insample']:.3f}")
    print("scope:", _SCOPE)


if __name__ == "__main__":
    _demo()
