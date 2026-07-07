"""Correctness tests for the counterfactual solver on synthetic data.

We build a synthetic dictionary where the target is an exact combination of two known
perturbations plus noise, then check each solver recovers a well-aligned small set.
These tests run in <1s on CPU and guard against regressions.
"""
import numpy as np

from src import counterfactual as cf


def _synthetic(seed=0, P=60, G=200, true_k=2):
    rng = np.random.default_rng(seed)
    E = rng.normal(size=(P, G)).astype(np.float32)
    true_idx = rng.choice(P, size=true_k, replace=False)
    w = np.abs(rng.normal(size=true_k)) + 0.5     # positive -> knockdown-compatible
    d = (E[true_idx].T @ w).astype(np.float32)
    d += 0.01 * rng.normal(size=G).astype(np.float32)
    return E, d, set(true_idx.tolist())


def test_greedy_recovers_alignment():
    E, d, _ = _synthetic()
    res = cf.greedy_minimal_set(E, d, k_max=6)
    assert res.cosine > 0.9
    assert 1 <= len(res.indices) <= 6


def test_omp_recovers_true_support():
    E, d, true_idx = _synthetic()
    res = cf.omp_minimal_set(E, d, k_max=2)
    assert res.cosine > 0.9
    assert set(res.indices) & true_idx        # overlaps the planted perturbations


def test_knockdown_only_forbids_negative_weights():
    E, d, _ = _synthetic()
    res = cf.greedy_minimal_set(E, d, k_max=6, knockdown_only=True)
    assert all(w >= -1e-9 for w in res.weights)


def test_beats_random_null():
    E, d, _ = _synthetic()
    res = cf.greedy_minimal_set(E, d, k_max=4)
    null = cf.random_null(E, d, k=len(res.indices), n_iter=200)
    assert res.cosine > np.quantile(null, 0.95)   # clearly above chance
