"""Convex-cone reachability for cell-state engineering.

The headline method for the cell-state-reachability project. Given a dictionary of
*measured* CRISPRi knockdown effect vectors `E` (P perturbations x G genes) and a target
transcriptomic shift `d` (G,), it answers a feasibility question, not a ranking one:

    Is `d` reachable by some NON-NEGATIVE combination of knockdown effects
    (i.e. is `d` inside the convex cone spanned by the rows of E), and if not,
    which genes make it unreachable (the constructive infeasibility certificate)?

Why a cone (the whole method in two facts):
  1. CRISPRi is loss-of-function: you can apply a knockdown or not, you cannot apply a
     "negative knockdown". So combination weights are non-negative: w >= 0.
  2. You can co-perturb: several knockdowns at once => the effect is (assumed) additive.
Together, the reachable set is C = { E^T w : w >= 0 } — a finitely generated convex
cone. Testing d in C is a non-negative least squares (NNLS) problem.

What makes the verdict falsifiable and honest:
  * reachable  -> residual ~ 0; the weights w are the ranked minimal recipe.
  * outside    -> residual has a component NO non-negative mix can supply. We return the
                  Farkas separating direction rho = d - E^T w*, whose positive, high-|z|
                  coordinates are genes the target wants UP that knockdown cannot deliver
                  => concrete CRISPRa hypotheses. This is the output no ranking method has.
  * "meaningfully outside" is decided against a SHUFFLED-TARGET null, never a hardcoded
    threshold, because with P<<G most random targets are partly outside by construction.

Dependencies: numpy, scipy.optimize.nnls. CPU-cheap (seconds at 34k x 2k). No GPU, no torch.

Author: cell-state-reachability. MIT.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import nnls


# --------------------------------------------------------------------------------------
# Result containers
# --------------------------------------------------------------------------------------
@dataclass
class ReachResult:
    """Geometry of fitting target `d` inside the non-negative cone spanned by rows of E."""
    reachable_cosine: float            # cos(E^T w*, d) at the NNLS optimum, in [-1, 1]
    residual_norm: float               # ||E^T w* - d|| / ||d||  (0 = exactly reachable)
    in_cone_fraction: float            # ||proj onto cone|| / ||d||  (fraction of d explained)
    weights: np.ndarray                # w* >= 0, length P (the ranked knockdown recipe)
    support: np.ndarray                # indices of w* > tol, sorted by weight desc
    fitted: np.ndarray                 # E^T w*  (the closest reachable point), length G
    residual: np.ndarray               # rho = d - E^T w*  (the unmet-demand / certificate direction)
    cert_max_violation: float          # max_p <e_p, rho>  — should be <= tol if w* is optimal
    n_generators: int
    n_genes: int

    def verdict(self, null: "NullResult | None" = None,
                cos_hi: float = 0.9, resid_lo: float = 0.35) -> str:
        """reachable | partially-reachable | outside.

        If a shuffled-target `null` is supplied the call is data-driven: 'outside' means
        the observed cosine is NOT above the null (target is indistinguishable from a
        shuffled one), 'reachable' means it clears the null's upper tail AND the residual
        is small. Without a null we fall back to geometry thresholds and you MUST report
        that the verdict was not null-calibrated.
        """
        if null is not None:
            if self.reachable_cosine <= null.p95:
                return "outside"                       # not distinguishable from noise
            if self.reachable_cosine >= null.p99 and self.residual_norm <= resid_lo:
                return "reachable"
            return "partially-reachable"
        # geometry-only fallback (uncalibrated — say so when you report it)
        if self.residual_norm <= resid_lo and self.reachable_cosine >= cos_hi:
            return "reachable"
        if self.reachable_cosine < 0.5:
            return "outside"
        return "partially-reachable"


@dataclass
class NullResult:
    observed: float
    null_cosines: np.ndarray
    p50: float
    p95: float
    p99: float
    percentile_of_observed: float      # where the observed value sits in the null (0-100)
    z: float                           # (observed - null_mean) / null_sd


@dataclass
class ActivationCertificate:
    """Constructive infeasibility certificate: what a knockdown screen CANNOT reach."""
    gene_index: np.ndarray             # indices into the gene axis, ranked most-unmet first
    residual_value: np.ndarray         # rho_j at those genes (positive = wants up, unmet)
    target_value: np.ndarray           # d_j at those genes
    explanation: str = field(default="")


# --------------------------------------------------------------------------------------
# Core solver
# --------------------------------------------------------------------------------------
def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0.0 or nb == 0.0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def reachability(E: np.ndarray, d: np.ndarray, *,
                 hvg_mask: Optional[np.ndarray] = None,
                 weight_tol: float = 1e-8,
                 max_iter: Optional[int] = None) -> ReachResult:
    """Fit target `d` inside the non-negative cone spanned by the rows of `E`.

    Parameters
    ----------
    E : (P, G) array
        Perturbation dictionary; row p is the measured effect vector of knockdown p.
    d : (G,) array
        Target transcriptomic shift (e.g. z-scored `toward_Th1` signature).
    hvg_mask : (G,) bool array, optional
        Restrict the fit to these genes (e.g. highly-variable genes) on BOTH E and d.
    weight_tol : float
        Weights below this are treated as zero when forming the support.

    Returns
    -------
    ReachResult  (see fields; call .verdict(null) for the classification)
    """
    E = np.asarray(E, dtype=float)
    d = np.asarray(d, dtype=float)
    if E.ndim != 2:
        raise ValueError(f"E must be 2-D (P,G); got {E.shape}")
    if d.shape[0] != E.shape[1]:
        raise ValueError(f"d length {d.shape[0]} != n_genes {E.shape[1]}")

    if hvg_mask is not None:
        hvg_mask = np.asarray(hvg_mask, dtype=bool)
        E = E[:, hvg_mask]
        d = d[hvg_mask]

    A = E.T                                   # (G, P): columns are the cone generators
    kw = {} if max_iter is None else {"maxiter": max_iter}
    w, _ = nnls(A, d, **kw)                    # min_{w>=0} ||A w - d||_2

    fitted = A @ w                             # closest reachable point E^T w*
    residual = d - fitted                      # rho: unmet demand / certificate direction
    dn = np.linalg.norm(d)
    resid_norm = float(np.linalg.norm(residual) / dn) if dn else 0.0
    in_cone = float(np.linalg.norm(fitted) / dn) if dn else 0.0

    # KKT / Farkas check: at the optimum, <e_p, rho> <= 0 for every generator p.
    # (equivalently E @ rho <= 0). The max positive violation certifies optimality.
    cert_violation = float(np.max(E @ residual)) if E.shape[0] else 0.0

    support = np.where(w > weight_tol)[0]
    support = support[np.argsort(-w[support])]

    return ReachResult(
        reachable_cosine=_cosine(fitted, d),
        residual_norm=resid_norm,
        in_cone_fraction=in_cone,
        weights=w,
        support=support,
        fitted=fitted,
        residual=residual,
        cert_max_violation=cert_violation,
        n_generators=E.shape[0],
        n_genes=E.shape[1],
    )


def activation_certificate(res: ReachResult, d: np.ndarray, *,
                           gene_names: Optional[list] = None,
                           top: int = 25,
                           hvg_mask: Optional[np.ndarray] = None) -> ActivationCertificate:
    """The genes that make the target unreachable by knockdown — CRISPRa hypotheses.

    A gene is an 'activation candidate' when the target wants it UP (d_j > 0) yet the
    closest reachable point still under-delivers there (residual rho_j > 0). Those are
    exactly the coordinates of the Farkas separating direction where no non-negative
    knockdown mix can push the transcriptome the way the target demands.
    """
    d = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        d = d[np.asarray(hvg_mask, dtype=bool)]
    rho = res.residual
    # score = unmet upward demand; only where the target itself wants the gene up
    score = np.where((rho > 0) & (d > 0), rho, 0.0)
    order = np.argsort(-score)
    order = order[score[order] > 0][:top]
    expl = (f"{len(order)} genes carry positive unmet demand: the target wants them "
            f"higher than any non-negative knockdown combination can deliver. These are "
            f"the falsifiable CRISPRa (gain-of-function) hypotheses the knockdown assay "
            f"structurally cannot test.")
    return ActivationCertificate(
        gene_index=order,
        residual_value=rho[order],
        target_value=d[order],
        explanation=expl,
    )


def directional_decomposition(res: ReachResult, d: np.ndarray, *,
                              hvg_mask: Optional[np.ndarray] = None) -> dict:
    """Split the target norm into knockdown-reachable vs activation-requiring fractions.

    Returns fractions that sum to ~1 (Pythagorean split of ||d||^2 into the fitted
    component inside the cone and the orthogonal-ish residual outside it).
    """
    d = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        d = d[np.asarray(hvg_mask, dtype=bool)]
    dn2 = float(np.dot(d, d))
    if dn2 == 0.0:
        return {"reachable_fraction": 0.0, "residual_fraction": 0.0}
    reachable = float(np.dot(res.fitted, res.fitted)) / dn2
    residual = float(np.dot(res.residual, res.residual)) / dn2
    # unmet UPWARD demand specifically (the activation story)
    up = res.residual[(res.residual > 0) & (d > 0)]
    up_frac = float(np.dot(up, up)) / dn2 if up.size else 0.0
    return {
        "reachable_fraction": reachable,
        "residual_fraction": residual,
        "activation_required_fraction": up_frac,
    }


# --------------------------------------------------------------------------------------
# Honest nulls
# --------------------------------------------------------------------------------------
def shuffled_target_null(E: np.ndarray, d: np.ndarray, *,
                         n_iter: int = 1000, seed: int = 0,
                         hvg_mask: Optional[np.ndarray] = None) -> NullResult:
    """Permute the gene labels of `d` and refit; the null band is the achievable-by-chance
    reachability cosine. The observed cosine is meaningful only above this band."""
    rng = np.random.default_rng(seed)
    observed = reachability(E, d, hvg_mask=hvg_mask).reachable_cosine
    dd = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, dtype=bool)
        E_use, d_use = E[:, m], dd[m]
    else:
        E_use, d_use = E, dd
    A = E_use.T
    cos = np.empty(n_iter)
    for i in range(n_iter):
        dp = d_use[rng.permutation(d_use.size)]
        w, _ = nnls(A, dp)
        cos[i] = _cosine(A @ w, dp)
    p50, p95, p99 = np.percentile(cos, [50, 95, 99])
    pct = float((cos < observed).mean() * 100.0)
    z = float((observed - cos.mean()) / cos.std()) if cos.std() else np.inf
    return NullResult(observed=observed, null_cosines=cos,
                      p50=float(p50), p95=float(p95), p99=float(p99),
                      percentile_of_observed=pct, z=z)


def reachability_spectrum(E: np.ndarray, d: np.ndarray, *, k_max: int = 12,
                          hvg_mask: Optional[np.ndarray] = None) -> dict:
    """Greedy forward selection under non-negativity: best reachable cosine vs. sparsity k.
    The 'minimal set' is the knee of this curve. Returns arrays k, cosine, residual."""
    E = np.asarray(E, dtype=float)
    d = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, dtype=bool)
        E, d = E[:, m], d[m]
    A = E.T
    P = E.shape[0]
    chosen: list[int] = []
    remaining = set(range(P))
    ks, coss, resids = [], [], []
    dn = np.linalg.norm(d)
    for _ in range(min(k_max, P)):
        best_j, best_cos, best_fit = None, -np.inf, None
        for j in remaining:
            idx = chosen + [j]
            w, _ = nnls(A[:, idx], d)
            fit = A[:, idx] @ w
            c = _cosine(fit, d)
            if c > best_cos:
                best_j, best_cos, best_fit = j, c, fit
        chosen.append(best_j)
        remaining.discard(best_j)
        ks.append(len(chosen))
        coss.append(best_cos)
        resids.append(float(np.linalg.norm(best_fit - d) / dn) if dn else 0.0)
    return {"k": np.array(ks), "cosine": np.array(coss),
            "residual": np.array(resids), "order": np.array(chosen)}


# --------------------------------------------------------------------------------------
# Self-test on synthetic data (runnable now, no Tier-2 download needed)
# --------------------------------------------------------------------------------------
def _selftest(seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    P, G = 300, 800
    E = rng.standard_normal((P, G))                       # mixed-sign effect vectors

    # (1) an IN-cone target: a known non-negative mix of a few generators
    true_idx = rng.choice(P, size=4, replace=False)
    true_w = rng.uniform(0.5, 2.0, size=4)
    d_in = (true_w[:, None] * E[true_idx]).sum(0)
    r_in = reachability(E, d_in)
    null_in = shuffled_target_null(E, d_in, n_iter=200, seed=1)

    # (2) an OUT-of-cone target: add a direction orthogonal to the row space is hard in
    #     P<G, so instead demand the NEGATIVE of an in-cone vector (needs "un-mixing").
    d_out = -d_in + 0.3 * rng.standard_normal(G)
    r_out = reachability(E, d_out)
    null_out = shuffled_target_null(E, d_out, n_iter=200, seed=2)
    cert = activation_certificate(r_out, d_out, top=10)
    dec_in = directional_decomposition(r_in, d_in)
    dec_out = directional_decomposition(r_out, d_out)

    print("=== IN-cone target (should be reachable) ===")
    print(f"  cosine={r_in.reachable_cosine:.4f}  residual={r_in.residual_norm:.4f}  "
          f"|support|={r_in.support.size}  recovered_true={set(true_idx) <= set(r_in.support.tolist())}")
    print(f"  KKT/Farkas max violation (should be <= ~1e-6): {r_in.cert_max_violation:.2e}")
    print(f"  null p95={null_in.p95:.3f} p99={null_in.p99:.3f}  observed pct={null_in.percentile_of_observed:.1f}")
    print(f"  verdict(null)={r_in.verdict(null_in)}")
    print(f"  decomposition reachable={dec_in['reachable_fraction']:.3f} "
          f"residual={dec_in['residual_fraction']:.3f}")

    print("=== OUT-of-cone target (should be less reachable) ===")
    print(f"  cosine={r_out.reachable_cosine:.4f}  residual={r_out.residual_norm:.4f}")
    print(f"  KKT/Farkas max violation (should be <= ~1e-6): {r_out.cert_max_violation:.2e}")
    print(f"  verdict(null)={r_out.verdict(null_out)}")
    print(f"  decomposition reachable={dec_out['reachable_fraction']:.3f} "
          f"activation_required={dec_out['activation_required_fraction']:.3f}")
    print(f"  activation certificate: {cert.gene_index.size} genes; top rho="
          f"{np.round(cert.residual_value[:5], 2).tolist()}")

    # assertions that must hold if the math is right
    assert r_in.reachable_cosine > 0.999, "in-cone target must be essentially fully reachable"
    assert r_in.residual_norm < 1e-3, "in-cone residual must be ~0"
    assert r_in.cert_max_violation < 1e-5, "KKT violation must be ~0 at optimum"
    assert r_out.cert_max_violation < 1e-5, "KKT violation must be ~0 at optimum"
    assert r_out.reachable_cosine < r_in.reachable_cosine, "out-of-cone must fit worse"
    print("\nALL SELF-TESTS PASSED")


if __name__ == "__main__":
    _selftest()
