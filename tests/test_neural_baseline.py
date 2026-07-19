"""Deterministic tests for neural_baseline.

Run in the pinned CI env (numpy==1.26.4, scipy==1.13.1). No data files, no torch,
no sklearn. The gradient check is the load-bearing correctness test: it certifies
the manual numpy backprop against a finite-difference reference so the L-BFGS-B
optimiser is descending the true objective.
"""

import numpy as np
import pytest

import neural_baseline as nb
import reachability as rx


# --------------------------------------------------------------------------- #
# cosine_full
# --------------------------------------------------------------------------- #
def test_cosine_full_known_values():
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    assert nb.cosine_full(a, a) == pytest.approx(1.0)
    assert nb.cosine_full(a, b) == pytest.approx(0.0)
    assert nb.cosine_full(a, -a) == pytest.approx(-1.0)


def test_cosine_full_zero_vector_is_safe():
    assert nb.cosine_full(np.zeros(3), np.array([1.0, 2.0, 3.0])) == 0.0
    assert nb.cosine_full(np.array([1.0, 2.0, 3.0]), np.zeros(3)) == 0.0


# --------------------------------------------------------------------------- #
# svd_reduce
# --------------------------------------------------------------------------- #
def test_svd_reduce_shape_and_orthonormal():
    rng = np.random.default_rng(0)
    atoms = rng.standard_normal((12, 40))
    basis, info = nb.svd_reduce(atoms, n_comp=8)
    assert basis.shape == (8, 40)
    assert info["n_comp"] == 8
    # rows are orthonormal
    gram = basis @ basis.T
    assert np.allclose(gram, np.eye(8), atol=1e-10)
    assert 0.0 < info["var_retained"] <= 1.0


def test_svd_reduce_variance_monotone_and_full():
    rng = np.random.default_rng(1)
    atoms = rng.standard_normal((10, 30))
    v5 = nb.svd_reduce(atoms, n_comp=5)[1]["var_retained"]
    v9 = nb.svd_reduce(atoms, n_comp=9)[1]["var_retained"]
    assert v9 >= v5
    # full rank keeps ~all variance
    v_full = nb.svd_reduce(atoms, n_comp=10)[1]["var_retained"]
    assert v_full == pytest.approx(1.0, abs=1e-9)


def test_svd_reduce_rejects_bad_n_comp():
    atoms = np.eye(4)
    with pytest.raises(rx.InputError):
        nb.svd_reduce(atoms, n_comp=0)
    with pytest.raises(rx.InputError):
        nb.svd_reduce(atoms, n_comp=5)  # > min(shape)


# --------------------------------------------------------------------------- #
# interaction_features
# --------------------------------------------------------------------------- #
def test_interaction_features_concat_and_product():
    a = np.array([1.0, 2.0])
    b = np.array([3.0, 4.0])
    feat = nb.interaction_features(a, b)
    # [a, b, a*b]
    assert np.allclose(feat, [1.0, 2.0, 3.0, 4.0, 3.0, 8.0])
    assert feat.shape == (6,)


def test_interaction_features_batched():
    A = np.ones((5, 3))
    B = 2 * np.ones((5, 3))
    feat = nb.interaction_features(A, B)
    assert feat.shape == (5, 9)
    # product block is 2.0 everywhere
    assert np.allclose(feat[:, 6:], 2.0)


def test_interaction_features_shape_mismatch():
    with pytest.raises(rx.InputError):
        nb.interaction_features(np.ones(3), np.ones(4))


# --------------------------------------------------------------------------- #
# Gradient correctness (the load-bearing manual-backprop check)
# --------------------------------------------------------------------------- #
def test_loss_grad_matches_finite_difference():
    rng = np.random.default_rng(3)
    n, d_in, hidden, d_out = 7, 6, 5, 4
    X = rng.standard_normal((n, d_in))
    Y = rng.standard_normal((n, d_out))
    theta = rng.standard_normal(d_in * hidden + hidden + hidden * d_out + d_out) * 0.3
    l2 = 0.1
    loss0, grad = nb._loss_grad(theta, X, Y, hidden, l2)
    eps = 1e-6
    num = np.empty_like(grad)
    for k in range(theta.size):
        tp = theta.copy(); tp[k] += eps
        tm = theta.copy(); tm[k] -= eps
        lp = nb._loss_grad(tp, X, Y, hidden, l2)[0]
        lm = nb._loss_grad(tm, X, Y, hidden, l2)[0]
        num[k] = (lp - lm) / (2 * eps)
    # analytic gradient must match finite difference to high precision
    assert np.max(np.abs(grad - num)) < 1e-6


def test_pack_unpack_roundtrip():
    rng = np.random.default_rng(4)
    d_in, hidden, d_out = 6, 5, 4
    w1 = rng.standard_normal((d_in, hidden)); b1 = rng.standard_normal(hidden)
    w2 = rng.standard_normal((hidden, d_out)); b2 = rng.standard_normal(d_out)
    theta = nb._pack(w1, b1, w2, b2)
    w1b, b1b, w2b, b2b = nb._unpack(theta, d_in, hidden, d_out)
    assert np.allclose(w1, w1b) and np.allclose(b1, b1b)
    assert np.allclose(w2, w2b) and np.allclose(b2, b2b)


# --------------------------------------------------------------------------- #
# InteractionMLP
# --------------------------------------------------------------------------- #
def _toy_reduced(seed=0, n=15, k=4):
    """A toy reduced-space dataset: doubles ~ (a+b) + small nonlinear term."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, k))
    B = rng.standard_normal((n, k))
    T = A + B + 0.15 * (A * B)  # additive + genuine interaction
    return A, B, T


def test_mlp_is_deterministic_given_seed():
    A, B, T = _toy_reduced()
    m1 = nb.InteractionMLP(hidden=6, l2=0.3, seed=0, maxiter=100).fit(A, B, T)
    m2 = nb.InteractionMLP(hidden=6, l2=0.3, seed=0, maxiter=100).fit(A, B, T)
    p1 = m1.predict_reduced(A, B)
    p2 = m2.predict_reduced(A, B)
    assert np.allclose(p1, p2)
    assert np.allclose(m1.params_.W1, m2.params_.W1)


def test_mlp_fits_training_data():
    A, B, T = _toy_reduced()
    m = nb.InteractionMLP(hidden=12, l2=1e-3, seed=0, maxiter=500).fit(A, B, T)
    pred = m.predict_reduced(A, B)
    # in-sample cosine per row should be high (it can fit the interaction)
    cos = [nb.cosine_full(pred[i], T[i]) for i in range(len(T))]
    assert np.mean(cos) > 0.95


def test_mlp_anchor_correction_collapses_to_constant_under_huge_l2():
    # L2 regularizes the WEIGHTS but not the biases (standard). With a huge L2 the
    # weights W1,W2 are crushed to ~0, so the learned correction becomes a
    # row-INDEPENDENT constant (the unregularized output bias absorbs the mean
    # training residual). The anchored prediction is therefore (a+b) + const,
    # i.e. it collapses toward additive up to a single constant offset.
    A, B, T = _toy_reduced()
    m = nb.InteractionMLP(hidden=6, l2=1e6, anchor=True, seed=0, maxiter=200).fit(A, B, T)
    assert np.abs(m.params_.W1).max() < 1e-4
    assert np.abs(m.params_.W2).max() < 1e-4
    correction = m.predict_reduced(A, B) - (A + B)
    # row-independent: every row's correction equals the column mean
    assert correction.std(axis=0).max() < 1e-6
    # and that constant is the mean training residual T-(A+B)
    assert np.allclose(correction.mean(0), (T - (A + B)).mean(0), atol=1e-3)


def test_mlp_predict_before_fit_raises():
    m = nb.InteractionMLP()
    with pytest.raises(rx.InputError):
        m.predict_reduced(np.ones((1, 3)), np.ones((1, 3)))


def test_mlp_rejects_bad_hyperparams():
    with pytest.raises(rx.InputError):
        nb.InteractionMLP(hidden=0)
    with pytest.raises(rx.InputError):
        nb.InteractionMLP(l2=-1.0)


# --------------------------------------------------------------------------- #
# additive + cone helpers
# --------------------------------------------------------------------------- #
def test_additive_prediction_is_exact_sum():
    a = np.array([1.0, -2.0, 3.0])
    b = np.array([0.5, 0.5, 0.5])
    assert np.allclose(nb.additive_prediction(a, b), a + b)


def test_prediction_error_score_is_one_minus_cosine():
    hc = np.array([0.9, 0.5, 1.0, 0.0])
    got = nb.prediction_error_score(hc)
    assert np.allclose(got, [0.1, 0.5, 0.0, 1.0])
    # monotone-decreasing in cosine: better prediction -> lower emergence score
    assert np.all(np.diff(nb.prediction_error_score(np.sort(hc))) <= 0)


def test_cone_reachable_fit_matches_project_cone():
    # cone_reachable_fit_cosine must equal rx.project_cone(...).cosine exactly
    # (it is a thin READ-ONLY wrapper, not a reimplementation).
    rng = np.random.default_rng(5)
    atoms = np.abs(rng.standard_normal((6, 20)))
    target = atoms[0] + atoms[1] + 0.3 * np.abs(rng.standard_normal(20))
    got = nb.cone_reachable_fit_cosine(atoms, target)
    ref = rx.project_cone(atoms, target).cosine
    assert got == pytest.approx(ref)


# --------------------------------------------------------------------------- #
# leave_pairs_out
# --------------------------------------------------------------------------- #
def test_leave_pairs_out_contract():
    # 4 atoms, 3 doubles; MLP is blind (LPO), additive is training-free,
    # cone_fit is in-sample and must be >= the blind MLP per-double (it sees target).
    rng = np.random.default_rng(6)
    n_genes = 30
    atoms = np.abs(rng.standard_normal((4, n_genes)))
    a_idx = [0, 1, 2]
    b_idx = [1, 2, 3]
    targets = np.stack([atoms[a] + atoms[b] + 0.1 * np.abs(rng.standard_normal(n_genes))
                        for a, b in zip(a_idx, b_idx)])
    names = ["A+B", "B+C", "C+D"]
    res = nb.leave_pairs_out(atoms, a_idx, b_idx, targets, names,
                             n_comp=4, hidden=6, l2=0.3, seed=0, maxiter=200)
    assert res.names == names
    assert res.mlp.shape == res.additive.shape == res.cone_fit.shape == (3,)
    assert res.n_comp == 4
    m = res.means()
    assert set(m) == {"mlp_heldout", "additive", "cone_fit_insample",
                      "mlp_beats_additive_frac"}
    # cone fit sees the target -> at least as good as blind additive, per double
    assert np.all(res.cone_fit >= res.additive - 1e-9)
    # scope disclaimer present and mentions the certificate
    assert "CERTIFICATE" in res.scope


def test_leave_pairs_out_input_alignment_checked():
    atoms = np.eye(4)
    with pytest.raises(rx.InputError):
        nb.leave_pairs_out(atoms, [0, 1], [1, 2], np.eye(3)[:, :4], ["x", "y"],
                           n_comp=2)


def test_leave_pairs_out_is_deterministic():
    rng = np.random.default_rng(7)
    atoms = np.abs(rng.standard_normal((4, 20)))
    a_idx, b_idx = [0, 1], [2, 3]
    targets = np.stack([atoms[a] + atoms[b] for a, b in zip(a_idx, b_idx)])
    names = ["A+C", "B+D"]
    r1 = nb.leave_pairs_out(atoms, a_idx, b_idx, targets, names, n_comp=4, seed=0, maxiter=100)
    r2 = nb.leave_pairs_out(atoms, a_idx, b_idx, targets, names, n_comp=4, seed=0, maxiter=100)
    assert np.allclose(r1.mlp, r2.mlp)
    assert np.allclose(r1.additive, r2.additive)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
