"""Strict, minimal geometry for screen-relative cell-state direction tests.

The public surface deliberately stops at projection, numerical diagnostics, and
held-out alignment. It does not emit biological verdicts, activation recipes,
minimal panels, or signed-modality claims.

Matrices use ``(perturbations, genes)`` orientation throughout.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.optimize import nnls


class InputError(ValueError):
    """Raised when a numerical problem is invalid or silently ambiguous."""


@dataclass(frozen=True)
class ProjectionResult:
    coefficients: np.ndarray
    fitted: np.ndarray
    residual: np.ndarray
    transformed_residual: np.ndarray
    dual_separator: np.ndarray | None
    cosine: float
    residual_fraction: float
    relative_objective: float
    kkt_violation: float
    geometry_status: str
    polarity_violation: float | None
    orthogonality_error: float | None
    separation_margin: float | None


@dataclass(frozen=True)
class HeldOutResult:
    coefficients: np.ndarray
    fit_cosine: float
    held_out_cosine: float
    fit_indices: np.ndarray
    score_indices: np.ndarray


def _stable_norm(x: np.ndarray) -> float:
    """Overflow-safe Euclidean norm."""
    scale = float(np.max(np.abs(x), initial=0.0))
    if scale == 0.0:
        return 0.0
    return scale * float(np.linalg.norm(x / scale))


def _weighted_cosine(a: np.ndarray, b: np.ndarray, q: np.ndarray) -> float:
    normalized_q = q / np.max(q)
    qa = np.sqrt(normalized_q) * a
    qb = np.sqrt(normalized_q) * b
    norm_a = _stable_norm(qa)
    norm_b = _stable_norm(qb)
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    return float(np.dot(qa / norm_a, qb / norm_b))


def _canonical_problem(
    effects: np.ndarray,
    target: np.ndarray,
    gene_weights: np.ndarray | None,
    gene_mask: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    effects = np.asarray(effects, dtype=float)
    target = np.asarray(target, dtype=float)
    if effects.ndim != 2 or target.ndim != 1:
        raise InputError("effects must be 2D and target must be 1D")
    if effects.shape[1] != target.size:
        raise InputError("effects and target gene axes do not align")
    if effects.shape[0] == 0 or target.size == 0:
        raise InputError("effects and target must be non-empty")
    if not np.all(np.isfinite(effects)) or not np.all(np.isfinite(target)):
        raise InputError("effects and target must contain only finite values")

    if gene_weights is None:
        weights = np.ones(target.size, dtype=float)
    else:
        weights = np.asarray(gene_weights, dtype=float)
        if weights.shape != target.shape:
            raise InputError("gene_weights must match the target axis")
        if not np.all(np.isfinite(weights)) or np.any(weights < 0):
            raise InputError("gene_weights must be finite and non-negative")

    if gene_mask is None:
        mask = np.ones(target.size, dtype=bool)
    else:
        mask = np.asarray(gene_mask)
        if mask.shape != target.shape:
            raise InputError("gene_mask must match the target axis")
        if not np.issubdtype(mask.dtype, np.bool_):
            raise InputError("gene_mask must contain booleans")

    # Zero-weight coordinates have no metric meaning; make their exclusion explicit.
    mask = mask & (weights > 0)
    if not np.any(mask):
        raise InputError("no positive-weight genes remain")

    kept_effects = effects[:, mask]
    kept_target = target[mask]
    kept_weights = weights[mask]
    # A positive global rescaling of W leaves the minimizer and all normalized
    # geometry unchanged. Normalize it to prevent avoidable derived overflow.
    kept_weights = kept_weights / np.max(kept_weights)
    if np.any(kept_weights == 0):
        raise InputError("gene-weight dynamic range is not representable")
    weighted_target = np.sqrt(kept_weights) * kept_target
    if _stable_norm(weighted_target) == 0.0:
        raise InputError("target is zero or below floating-point resolution")
    return kept_effects, kept_target, kept_weights, mask


def project_cone(
    effects: np.ndarray,
    target: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
    gene_mask: np.ndarray | None = None,
    separator_tolerance: float | None = None,
) -> ProjectionResult:
    """Project a target onto the unbounded non-negative cone of measured effects.

    ``dual_separator`` is model-relative. It certifies separation from this cone
    under this metric; it is not evidence of biological impossibility.
    """
    effects_full = np.asarray(effects, dtype=float)
    target_full = np.asarray(target, dtype=float)
    if separator_tolerance is not None and (
        not np.isfinite(separator_tolerance) or separator_tolerance < 0
    ):
        raise InputError("separator_tolerance must be finite and non-negative")
    effects_kept, target_kept, q, mask = _canonical_problem(
        effects_full, target_full, gene_weights, gene_mask
    )

    sqrt_q = np.sqrt(q)
    a_all = effects_kept.T
    aw_all = sqrt_q[:, None] * a_all
    bw = sqrt_q * target_kept

    atom_norms = np.array([_stable_norm(aw_all[:, i]) for i in range(aw_all.shape[1])])
    active_atoms = atom_norms > 0.0
    coefficients = np.zeros(effects_full.shape[0], dtype=float)
    nonzero_norms = atom_norms[active_atoms]
    unit_atoms = aw_all[:, active_atoms] / nonzero_norms
    norm_b = _stable_norm(bw)
    if np.any(active_atoms):
        # Normalize columns so positive atom rescaling cannot change conditioning
        # or the fitted cone point. Normalize the target as well so valid very small
        # or very large directions do not inherit avoidable solver-scale failures.
        try:
            solved, _ = nnls(unit_atoms, bw / norm_b)
        except RuntimeError as exc:
            raise InputError("NNLS solver failed numerical certification") from exc
        coefficients[active_atoms] = solved * norm_b / nonzero_norms

    fitted_kept = coefficients @ effects_kept
    residual_kept = target_kept - fitted_kept
    if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(fitted_kept)):
        raise InputError("projection is not representable at the supplied numerical scale")
    rw = sqrt_q * residual_kept
    fitw = sqrt_q * fitted_kept
    norm_fit = _stable_norm(fitw)
    norm_r = _stable_norm(rw)
    scale = max(norm_b, norm_fit)
    eps = np.finfo(float).eps

    # Atom-scale-invariant KKT diagnostics.
    gamma = unit_atoms.T @ (rw / scale)
    contributions = coefficients[active_atoms] * nonzero_norms / scale
    active = contributions > 100 * np.sqrt(eps)
    primal = float(np.max(np.maximum(-coefficients[active_atoms], 0) * nonzero_norms / scale, initial=0))
    dual = float(np.max(np.maximum(gamma, 0), initial=0))
    stationarity = float(np.max(np.abs(gamma[active]), initial=0))
    complementarity = float(np.max(np.abs(contributions * gamma), initial=0))
    relative_b = norm_b / scale
    relative_fit = norm_fit / scale
    relative_residual = norm_r / scale
    projection_identity = abs(
        relative_b * relative_b
        - relative_fit * relative_fit
        - relative_residual * relative_residual
    ) / (relative_b * relative_b + eps)
    kkt = max(primal, dual, stationarity, complementarity, projection_identity)

    fitted = np.zeros_like(target_full, dtype=float)
    residual = np.zeros_like(target_full, dtype=float)
    transformed_residual = np.zeros_like(target_full, dtype=float)
    fitted[mask] = fitted_kept
    residual[mask] = residual_kept
    transformed_residual[mask] = rw

    relative_separation = norm_r / scale
    geometry_tolerance = separator_tolerance
    if geometry_tolerance is None:
        geometry_tolerance = max(100 * np.sqrt(eps), 1e-10)
    certification_tolerance = 1e-8

    if not np.isfinite(kkt) or kkt > certification_tolerance:
        raise InputError("projection failed KKT numerical certification")

    if relative_separation <= geometry_tolerance:
        geometry_status = "inside_tolerance"
        separator = None
        polarity = orthogonality = separation = None
    else:
        geometry_status = "outside_model_cone"
        separator = np.zeros_like(target_full, dtype=float)
        # q has canonical max(q)=1. A dual separator is only defined up to
        # positive scale, so this avoids avoidable overflow without changing it.
        separator_kept = q * residual_kept
        if not np.all(np.isfinite(separator_kept)):
            raise InputError("separator is not representable at the supplied numerical scale")
        separator[mask] = separator_kept
        polarity = float(
            np.max(
                np.maximum(unit_atoms.T @ (rw / norm_r), 0),
                initial=0,
            )
        )
        orthogonality = (
            0.0
            if norm_fit == 0.0
            else abs(float(np.dot(fitw / norm_fit, rw / norm_r)))
        )
        separation = float(np.dot(bw / norm_b, rw / norm_r))
        if (
            not np.isfinite(polarity)
            or not np.isfinite(orthogonality)
            or not np.isfinite(separation)
            or polarity > certification_tolerance
            or orthogonality > certification_tolerance
            or separation <= 0
        ):
            raise InputError("separator failed numerical certification")

    return ProjectionResult(
        coefficients=coefficients,
        fitted=fitted,
        residual=residual,
        transformed_residual=transformed_residual,
        dual_separator=separator,
        cosine=_weighted_cosine(fitted_kept, target_kept, q),
        residual_fraction=relative_separation,
        relative_objective=0.5 * (norm_r / norm_b) * (norm_r / norm_b),
        kkt_violation=kkt,
        geometry_status=geometry_status,
        polarity_violation=polarity,
        orthogonality_error=orthogonality,
        separation_margin=separation,
    )


def held_out_alignment(
    effects: np.ndarray,
    target: np.ndarray,
    fit_indices: Iterable[int],
    score_indices: Iterable[int],
    *,
    gene_weights: np.ndarray | None = None,
) -> HeldOutResult:
    """Fit coefficients on one gene set and score the frozen fit on another."""
    effects = np.asarray(effects, dtype=float)
    target = np.asarray(target, dtype=float)
    if effects.ndim != 2 or target.ndim != 1:
        raise InputError("effects must be 2D and target must be 1D")
    if effects.shape[1] != target.size or effects.shape[0] == 0 or target.size == 0:
        raise InputError("effects and target must be non-empty with aligned gene axes")
    if not np.all(np.isfinite(effects)) or not np.all(np.isfinite(target)):
        raise InputError("effects and target must contain only finite values")
    try:
        fit_idx = np.asarray(tuple(fit_indices))
        score_idx = np.asarray(tuple(score_indices))
    except (TypeError, ValueError) as exc:
        raise InputError("gene indices must be integer iterables") from exc
    if fit_idx.ndim != 1 or score_idx.ndim != 1:
        raise InputError("gene indices must be one-dimensional")
    if fit_idx.size == 0 or score_idx.size == 0:
        raise InputError("fit_indices and score_indices must be non-empty")
    for indices in (fit_idx, score_idx):
        if np.issubdtype(indices.dtype, np.bool_) or not np.issubdtype(
            indices.dtype, np.integer
        ):
            raise InputError("gene indices must contain integers, not booleans")
    fit_idx = fit_idx.astype(np.intp, copy=False)
    score_idx = score_idx.astype(np.intp, copy=False)
    if np.unique(fit_idx).size != fit_idx.size or np.unique(score_idx).size != score_idx.size:
        raise InputError("gene indices must not contain duplicates")
    if np.intersect1d(fit_idx, score_idx).size:
        raise InputError("fit and score indices must be disjoint")
    if np.any(fit_idx < 0) or np.any(score_idx < 0):
        raise InputError("gene indices must be non-negative")
    if np.any(fit_idx >= target.size) or np.any(score_idx >= target.size):
        raise InputError("gene index is out of range")

    weights = np.ones(target.size) if gene_weights is None else np.asarray(gene_weights, dtype=float)
    if weights.shape != target.shape:
        raise InputError("gene_weights must match the target axis")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise InputError("gene_weights must be finite and non-negative")
    if np.any(weights[fit_idx] == 0) or np.any(weights[score_idx] == 0):
        raise InputError("held-out indices must have positive gene weights")
    fit = project_cone(
        effects[:, fit_idx], target[fit_idx], gene_weights=weights[fit_idx]
    )
    with np.errstate(over="ignore", invalid="ignore"):
        frozen = fit.coefficients @ effects
    if not np.all(np.isfinite(frozen)):
        raise InputError(
            "held-out prediction is not representable at the supplied numerical scale"
        )
    return HeldOutResult(
        coefficients=fit.coefficients,
        fit_cosine=_weighted_cosine(frozen[fit_idx], target[fit_idx], weights[fit_idx]),
        held_out_cosine=_weighted_cosine(
            frozen[score_idx], target[score_idx], weights[score_idx]
        ),
        fit_indices=fit_idx,
        score_indices=score_idx,
    )


def empirical_p(observed: float, null_values: Iterable[float]) -> float:
    """Conservative plus-one empirical p-value with ties counted as exceedances."""
    null = np.asarray(tuple(null_values), dtype=float)
    if null.ndim != 1 or null.size == 0 or not np.all(np.isfinite(null)):
        raise InputError("null_values must be a non-empty finite vector")
    if not np.isfinite(observed):
        raise InputError("observed must be finite")
    return float((1 + np.count_nonzero(null >= observed)) / (null.size + 1))


def _demo() -> None:
    effects = np.eye(4)
    target = np.array([1.0, 0.0, -1.0, 0.0])
    result = project_cone(effects, target)
    print("geometry:", result.geometry_status)
    print("cosine:", f"{result.cosine:.3f}")
    print("residual fraction:", f"{result.residual_fraction:.3f}")
    print("KKT violation:", f"{result.kkt_violation:.3e}")


if __name__ == "__main__":
    _demo()
