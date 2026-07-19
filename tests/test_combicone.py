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


# =========================================================================== #
# k-way (order >= 3) extension
#
# The engine's projection is already order-agnostic (cone + one target). These
# tests exercise the triage/certify layer at order 3 and, critically, the
# reachable-from-lower-order flip that is the whole point of the k-way machinery.
# =========================================================================== #


# --------------------------------------------------------------------------- #
# combo_cosine: the k-way generalization of single_effect_cosine
# --------------------------------------------------------------------------- #
def test_combo_cosine_reduces_to_pairwise_at_k2():
    # At k=2, both mean and max reduce to the single pairwise cosine, so a k=2
    # triage built on combo_cosine is identical to -cos(a, b).
    rng = np.random.default_rng(3)
    a = rng.normal(size=20)
    b = rng.normal(size=20)
    pairwise = cc.single_effect_cosine(a, b)
    assert cc.combo_cosine([a, b], agg="mean") == pytest.approx(pairwise)
    assert cc.combo_cosine([a, b], agg="max") == pytest.approx(pairwise)


def test_combo_cosine_mean_and_max_semantics():
    # Two orthogonal axes and one axis collinear with the first.
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    c = np.array([1.0, 0.0, 0.0])  # collinear with a
    # pairwise cosines: (a,b)=0, (a,c)=1, (b,c)=0 -> mean=1/3, max=1
    assert cc.combo_cosine([a, b, c], agg="mean") == pytest.approx(1.0 / 3.0)
    assert cc.combo_cosine([a, b, c], agg="max") == pytest.approx(1.0)


def test_combo_cosine_validates():
    with pytest.raises(rx.InputError):
        cc.combo_cosine([np.ones(3)])  # needs >= 2
    with pytest.raises(rx.InputError):
        cc.combo_cosine([np.ones(3), np.ones(3)], agg="median")  # bad agg


# --------------------------------------------------------------------------- #
# additive_gap: variadic (pairs unchanged, triples supported)
# --------------------------------------------------------------------------- #
def test_additive_gap_pair_call_unchanged():
    # The historical two-index call must behave exactly as before.
    singles = np.eye(4)
    assert cc.additive_gap(singles, 0, 1) > 0.5


def test_additive_gap_triple_reachable_and_unreachable():
    # Unreachable: e0+e1+e2 with those three removed leaves only e3 -> big gap.
    singles = np.eye(5)
    assert cc.additive_gap(singles, 0, 1, 2) > 0.5
    # Reachable: a redundant atom equals the additive triple exactly -> gap 0.
    singles2 = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],  # == row0+row1+row2
        ]
    )
    assert cc.additive_gap(singles2, 0, 1, 2) == pytest.approx(0.0, abs=1e-9)


def test_additive_gap_triple_rejects_bad_input():
    singles = np.eye(4)
    with pytest.raises(rx.InputError):
        cc.additive_gap(singles, 0)  # < 2 indices
    with pytest.raises(rx.InputError):
        cc.additive_gap(singles, 0, 1, 1)  # repeated
    with pytest.raises(rx.InputError):
        cc.additive_gap(singles, 0, 1, 9)  # out of range
    # only three atoms, all in the combo -> no leave-combo-out cone
    with pytest.raises(rx.InputError):
        cc.additive_gap(np.eye(3), 0, 1, 2)


# --------------------------------------------------------------------------- #
# triage_combinations at order 3
# --------------------------------------------------------------------------- #
def test_triage_order3_generates_all_triples():
    singles = np.eye(5)
    names = [f"g{i}" for i in range(5)]
    res = cc.triage_combinations(singles, names, order=3)
    assert len(res.pairs) == 10  # C(5,3)
    assert all(len(c) == 3 for c in res.pairs)
    assert res.order == 3
    assert res.pairwise_agg == "mean"
    assert res.combos == res.pairs  # order-neutral alias
    assert sorted(res.rank.tolist()) == list(range(1, 11))
    assert "model-relative" in res.scope


def test_triage_k2_identical_to_legacy_score():
    # The generalized path at order 2 must reproduce -cos(a, b) exactly.
    rng = np.random.default_rng(5)
    singles = rng.normal(size=(6, 15))
    names = [f"g{i}" for i in range(6)]
    res = cc.triage_combinations(singles, names)  # default order=2
    assert res.order == 2
    # score is exactly the negated pairwise cosine, combo by combo
    for combo, score in zip(res.pairs, res.score):
        ia, ib = names.index(combo[0]), names.index(combo[1])
        expected = -cc.single_effect_cosine(singles[ia], singles[ib])
        assert score == pytest.approx(expected)


def test_triage_order3_orthogonal_scores_above_collinear():
    # A spread-out triple (all mutually orthogonal) should outrank a triple with a
    # collinear pair, under the default mean aggregation.
    singles = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],  # A
            [0.0, 1.0, 0.0, 0.0],  # B
            [0.0, 0.0, 1.0, 0.0],  # C  (A,B,C mutually orthogonal)
            [1.0, 0.0, 0.0, 0.0],  # D  (collinear with A)
        ]
    )
    names = ["A", "B", "C", "D"]
    res = cc.triage_combinations(
        singles, names, [("A", "B", "C"), ("A", "B", "D")]
    )
    assert res.order == 3
    score = {c: s for c, s in zip(res.pairs, res.score)}
    # (A,B,C): mean cos 0 -> score 0; (A,B,D): pair (A,D)=1 -> mean cos 1/3 -> score -1/3
    assert score[("A", "B", "C")] > score[("A", "B", "D")]


def test_triage_max_aggregation_is_more_pessimistic():
    # With one collinear pair inside the triple, max aggregation drives the score
    # to -1 (a single collinear pair dominates) while mean only to -1/3.
    singles = np.array(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
    )
    names = ["A", "B", "C"]  # (A,C) collinear
    mean_res = cc.triage_combinations(singles, names, order=3, pairwise="mean")
    max_res = cc.triage_combinations(singles, names, order=3, pairwise="max")
    assert mean_res.score[0] == pytest.approx(-1.0 / 3.0)
    assert max_res.score[0] == pytest.approx(-1.0)


def test_triage_order3_use_gap_runs():
    singles = np.eye(5)
    names = [f"g{i}" for i in range(5)]
    res = cc.triage_combinations(singles, names, order=3, use_gap=True)
    assert np.all(np.isfinite(res.gap))  # gap computed for every triple
    assert res.order == 3


def test_triage_explicit_combos_infer_order():
    singles = np.eye(5)
    names = [f"g{i}" for i in range(5)]
    res = cc.triage_combinations(singles, names, [("g0", "g1", "g2")])
    assert res.order == 3


def test_triage_candidate_pairs_alias_still_works():
    # The old keyword name must remain accepted.
    singles = np.eye(4)
    names = ["a", "b", "c", "d"]
    res = cc.triage_combinations(singles, names, candidate_pairs=[("a", "b")])
    assert res.pairs == [("a", "b")]
    # passing both aliases is an error
    with pytest.raises(rx.InputError):
        cc.triage_combinations(
            singles, names, [("a", "b")], candidate_pairs=[("a", "c")]
        )


def test_triage_rejects_mixed_order_and_repeats():
    singles = np.eye(4)
    names = ["a", "b", "c", "d"]
    with pytest.raises(rx.InputError):
        cc.triage_combinations(singles, names, [("a", "b"), ("a", "b", "c")])
    with pytest.raises(rx.InputError):
        cc.triage_combinations(singles, names, [("a", "a", "b")])  # repeated member
    with pytest.raises(rx.InputError):
        cc.triage_combinations(singles, names, order=3, pairwise="nope")


# --------------------------------------------------------------------------- #
# certify_emergence at order 3 (cone = any-order atoms)
# --------------------------------------------------------------------------- #
def test_certify_triple_emergent_from_singles():
    # Three singles spanning x,y,z; a measured "triple" with a big 4th-axis
    # component is unreachable from the single-gene cone.
    singles = np.eye(4)[:3]  # e0,e1,e2 (spans first three axes)
    triple = np.array([1.0, 1.0, 1.0, 5.0])  # 4th axis unreachable
    cert = cc.certify_emergence(singles, triple, noise_sd=0.01, n_boot=200, seed=0)
    assert cert.geometry_status == "outside_model_cone"
    assert cert.separator is not None
    assert cert.p_value < 0.05
    assert cert.floor_ratio > 1.9
    assert "certified emergent" in cert.verdict


def test_certify_cone_atoms_alias_matches_singles():
    singles = np.eye(4)[:3]
    triple = np.array([1.0, 1.0, 1.0, 5.0])
    a = cc.certify_emergence(singles, triple, noise_sd=0.01, n_boot=100, seed=1)
    b = cc.certify_emergence(
        cone_atoms=singles, measured_combo=triple, noise_sd=0.01, n_boot=100, seed=1
    )
    assert a.unreachable_fraction == pytest.approx(b.unreachable_fraction)
    assert a.p_value == pytest.approx(b.p_value)
    # passing both, or neither, is an error
    with pytest.raises(rx.InputError):
        cc.certify_emergence(singles, triple, cone_atoms=singles)
    with pytest.raises(rx.InputError):
        cc.certify_emergence(measured_combo=triple)


def test_certify_reachable_from_lower_order_flips_verdict():
    # THE k-way test. A measured triple carries a component that is off the cone of
    # its three singles, but IS present in one measured lower-order combination.
    # Adding that lower-order effect to the cone must flip the verdict from emergent
    # to reachable -- the triple has no epistasis beyond that lower-order part.
    e0 = np.array([1.0, 0.0, 0.0, 0.0])
    e1 = np.array([0.0, 1.0, 0.0, 0.0])
    e2 = np.array([0.0, 0.0, 1.0, 0.0])
    singles = np.vstack([e0, e1, e2])  # spans first three axes only
    # measured triple = additive part + an off-axis (4th) component
    triple = e0 + e1 + e2 + np.array([0.0, 0.0, 0.0, 3.0])
    # a measured lower-order combo that already expresses the 4th-axis direction
    lower = np.array([0.0, 0.0, 0.0, 3.0])

    from_singles = cc.certify_emergence(
        singles, triple, noise_sd=0.01, n_boot=200, seed=0
    )
    enriched = np.vstack([singles, lower])
    from_lower = cc.certify_emergence(
        cone_atoms=enriched, measured_combo=triple, noise_sd=0.01, n_boot=200, seed=0
    )
    # emergent from singles ...
    assert from_singles.geometry_status == "outside_model_cone"
    assert "certified emergent" in from_singles.verdict
    # ... but reachable once the lower-order atom is in the cone
    assert from_lower.geometry_status == "inside_tolerance"
    assert from_lower.unreachable_fraction < 1e-6
    assert "certified emergent" not in from_lower.verdict


# --------------------------------------------------------------------------- #
# fit_triage_model at order 3
# --------------------------------------------------------------------------- #
def test_fit_triage_model_order3_learns_and_predicts():
    rng = np.random.default_rng(11)
    singles = rng.normal(size=(9, 14))
    names = [f"g{i}" for i in range(9)]
    triples = [(names[i], names[j], names[k])
               for i in range(9) for j in range(i + 1, 9) for k in range(j + 1, 9)]
    labeled = {}
    for combo in triples:
        idxs = [names.index(n) for n in combo]
        mc = cc.combo_cosine([singles[i] for i in idxs], agg="mean")
        labeled[combo] = -mc + 0.01 * rng.normal()
    model = cc.fit_triage_model(singles, names, labeled)
    assert model.coef_.shape == (2,)
    assert model.pairwise_agg == "mean"
    res = cc.triage_combinations(singles, names, list(labeled.keys()), model=model)
    assert res.model == "ridge"
    assert res.order == 3
    # higher predicted score should track lower aggregated cosine
    hi = res.score > np.median(res.score)
    assert res.cos_ab[hi].mean() < res.cos_ab[~hi].mean()


def test_fit_triage_model_rejects_mixed_order():
    singles = np.eye(5)
    names = [f"g{i}" for i in range(5)]
    with pytest.raises(rx.InputError):
        cc.fit_triage_model(
            singles, names, {("g0", "g1"): 1.0, ("g0", "g1", "g2"): 0.5, ("g2", "g3"): 0.2}
        )
