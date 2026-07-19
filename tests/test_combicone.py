"""Tests for the combicone triage + certify layer.

These use small synthetic cones with known geometry (no data files) so they run in
the pinned CI, plus a couple of invariants that must hold on any input.
"""

import numpy as np
import pytest

import combicone as cc
import reachability as rx


# --------------------------------------------------------------------------- #
# Feature primitives
# --------------------------------------------------------------------------- #
def test_single_effect_cosine_known_values():
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    assert cc.single_effect_cosine(a, a) == pytest.approx(1.0)
    assert cc.single_effect_cosine(a, b) == pytest.approx(0.0)
    assert cc.single_effect_cosine(a, -a) == pytest.approx(-1.0)


def test_single_effect_cosine_zero_vector_is_safe():
    a = np.array([0.0, 0.0])
    b = np.array([1.0, 1.0])
    assert cc.single_effect_cosine(a, b) == 0.0


def test_additive_gap_reachable_is_zero():
    # Additive vector of two axes lies inside the cone of the remaining axes ONLY
    # if those axes span it; here e0+e1 with the pair removed leaves axes 2,3 which
    # cannot represent it -> gap > 0. Use a case that IS reachable instead:
    singles = np.eye(4)
    # additive of rows 0,1 = [1,1,0,0]; remove rows 0,1 -> cone of e2,e3 cannot reach
    gap_unreachable = cc.additive_gap(singles, 0, 1)
    assert gap_unreachable > 0.5
    # a redundant library where the pair is duplicated elsewhere is reachable
    singles2 = np.array(
        [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]  # row2 == row0+row1
    )
    gap_reachable = cc.additive_gap(singles2, 0, 1)
    assert gap_reachable == pytest.approx(0.0, abs=1e-9)


def test_additive_gap_rejects_bad_indices():
    singles = np.eye(3)
    with pytest.raises(rx.InputError):
        cc.additive_gap(singles, 0, 0)
    with pytest.raises(rx.InputError):
        cc.additive_gap(singles, 0, 5)


def test_additive_gap_needs_a_third_atom():
    # Only two atoms: leaving the pair out empties the cone -> clear error.
    singles = np.eye(2)
    with pytest.raises(rx.InputError):
        cc.additive_gap(singles, 0, 1)


# --------------------------------------------------------------------------- #
# triage_combinations
# --------------------------------------------------------------------------- #
def test_triage_default_all_pairs_and_ranks():
    singles = np.eye(4)
    names = ["g0", "g1", "g2", "g3"]
    res = cc.triage_combinations(singles, names)
    assert len(res.pairs) == 6  # C(4,2)
    # ranks are a permutation of 1..6
    assert sorted(res.rank.tolist()) == [1, 2, 3, 4, 5, 6]
    assert res.model == "training-free"
    assert "model-relative" in res.scope
    # top(k) returns k pairs, highest score first
    top2 = res.top(2)
    assert len(top2) == 2
    order = np.argsort(-res.score, kind="stable")
    assert top2 == [res.pairs[order[0]], res.pairs[order[1]]]


def test_triage_orthogonal_singles_score_above_collinear():
    # Two orthogonal singles (low cosine) should be flagged as MORE emergent than
    # two near-collinear singles.
    singles = np.array(
        [
            [1.0, 0.0, 0.0],  # A
            [0.0, 1.0, 0.0],  # B (orthogonal to A)
            [1.0, 0.0, 0.0],  # C
            [0.999, 0.045, 0.0],  # D (near-collinear with C)
        ]
    )
    names = ["A", "B", "C", "D"]
    res = cc.triage_combinations(
        singles, names, [("A", "B"), ("C", "D")], use_gap=False
    )
    # score = -cosine; A,B orthogonal (cos 0) > C,D collinear (cos ~1)
    assert res.score[0] > res.score[1]


def test_triage_rejects_unknown_pair():
    singles = np.eye(3)
    with pytest.raises(rx.InputError):
        cc.triage_combinations(singles, ["a", "b", "c"], [("a", "zzz")])


def test_triage_requires_unique_names():
    with pytest.raises(rx.InputError):
        cc.triage_combinations(np.eye(3), ["a", "a", "b"])


def test_triage_is_deterministic():
    rng = np.random.default_rng(0)
    singles = rng.normal(size=(6, 10))
    names = [f"g{i}" for i in range(6)]
    r1 = cc.triage_combinations(singles, names)
    r2 = cc.triage_combinations(singles, names)
    np.testing.assert_array_equal(r1.score, r2.score)
    np.testing.assert_array_equal(r1.rank, r2.rank)


# --------------------------------------------------------------------------- #
# certify_emergence
# --------------------------------------------------------------------------- #
def test_certify_reachable_target_is_not_emergent():
    # Target exactly inside the cone: residual ~0, must not be called emergent.
    singles = np.eye(5)
    combo = np.array([1.0, 2.0, 0.5, 0.0, 0.0])  # non-negative combo of axes
    cert = cc.certify_emergence(singles, combo, noise_sd=0.01, n_boot=100, seed=0)
    assert cert.geometry_status == "inside_tolerance"
    assert cert.unreachable_fraction < 1e-6
    assert "within measurement noise" in cert.verdict or cert.floor_ratio < 1.9


def test_certify_clearly_outside_is_emergent():
    # Target with a large component orthogonal to the whole cone, tiny noise.
    singles = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])  # spans x,y only
    combo = np.array([1.0, 1.0, 5.0])  # big z component is unreachable
    cert = cc.certify_emergence(singles, combo, noise_sd=0.01, n_boot=200, seed=0)
    assert cert.geometry_status == "outside_model_cone"
    assert cert.separator is not None
    assert cert.p_value < 0.05
    assert cert.floor_ratio > 1.9
    assert "certified emergent" in cert.verdict


def test_certify_high_noise_swamps_signal():
    # Same geometry, but noise as large as the signal -> not certifiable.
    singles = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    combo = np.array([1.0, 1.0, 0.3])
    cert = cc.certify_emergence(singles, combo, noise_sd=5.0, n_boot=200, seed=0)
    # floor ratio should collapse toward 1 (residual is mostly noise)
    assert cert.floor_ratio < 1.9


def test_certify_point_only_without_noise():
    singles = np.array([[1.0, 0.0], [0.0, 1.0]])
    combo = np.array([1.0, -1.0])  # outside the non-negative cone
    cert = cc.certify_emergence(singles, combo, noise_sd=None)
    assert cert.geometry_status == "outside_model_cone"
    assert cert.separator is not None
    assert np.isnan(cert.z)
    assert np.isnan(cert.p_value)
    assert "point certificate only" in cert.verdict


def test_certify_is_deterministic_given_seed():
    singles = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    combo = np.array([1.0, 1.0, 2.0])
    c1 = cc.certify_emergence(singles, combo, noise_sd=0.1, n_boot=100, seed=7)
    c2 = cc.certify_emergence(singles, combo, noise_sd=0.1, n_boot=100, seed=7)
    assert c1.p_value == c2.p_value
    assert c1.z == c2.z
    assert c1.ci_low == c2.ci_low


def test_certify_rejects_bad_noise():
    singles = np.eye(3)
    combo = np.array([1.0, 1.0, -1.0])
    with pytest.raises(rx.InputError):
        cc.certify_emergence(singles, combo, noise_sd=np.array([1.0, 1.0]))  # wrong len
    with pytest.raises(rx.InputError):
        cc.certify_emergence(singles, combo, noise_sd=-1.0)


def test_certify_gene_axis_mismatch():
    singles = np.eye(3)
    with pytest.raises(rx.InputError):
        cc.certify_emergence(singles, np.array([1.0, 2.0]))


# --------------------------------------------------------------------------- #
# fit_triage_model
# --------------------------------------------------------------------------- #
def test_fit_triage_model_learns_and_predicts():
    rng = np.random.default_rng(1)
    singles = rng.normal(size=(8, 12))
    names = [f"g{i}" for i in range(8)]
    # synthetic labels: emergence ~ -cosine so the model should learn a negative
    # weight on cos_ab
    labeled = {}
    for i in range(8):
        for j in range(i + 1, 8):
            c = cc.single_effect_cosine(singles[i], singles[j])
            labeled[(names[i], names[j])] = -c + 0.01 * rng.normal()
    model = cc.fit_triage_model(singles, names, labeled)
    assert model.coef_.shape == (2,)
    # applying the model through triage should reproduce a ranking that correlates
    # with -cosine
    res = cc.triage_combinations(singles, names, list(labeled.keys()), model=model)
    assert res.model == "ridge"
    # higher score should mean lower cosine (Spearman sign check via a coarse test)
    hi = res.score > np.median(res.score)
    assert res.cos_ab[hi].mean() < res.cos_ab[~hi].mean()


def test_fit_triage_model_needs_enough_pairs():
    singles = np.eye(3)
    with pytest.raises(rx.InputError):
        cc.fit_triage_model(singles, ["a", "b", "c"], {("a", "b"): 1.0})
