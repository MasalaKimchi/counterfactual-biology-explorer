"""Tests for the screenloop sequential-design + library-augmentation layer.

These use small synthetic screens with a *planted* missing axis, so the correct
answer is known without any data file and the tests run in the pinned CI. The
planted-truth construction is the same idea used across this repo's k-way tests: a
combination (or a held-out single) carries an effect direction no other atom
supplies, so a correct certificate must point back at it.

Structure:
  * property invariants that must hold on any input (curve monotonicity, bounds,
    determinism),
  * fail-closed input validation (corrupted inputs raise rx.InputError),
  * regression tests on a planted screen where the separator MUST recover the
    hidden axis and MUST beat the naive baseline / random null.
"""

import numpy as np
import pytest

import reachability as rx
import screenloop as sl


# --------------------------------------------------------------------------- #
# Fixtures: a planted synthetic screen
# --------------------------------------------------------------------------- #
def _planted_screen(seed=0, n_atoms=8, n_genes=40, n_emergent=3):
    """Build atoms + doubles with a few planted-emergent combinations.

    Emergent doubles get an extra axis (gene 0) that no single carries, so their
    measured effect leaves the single-gene cone; additive doubles do not.
    Returns (atoms, combo_effects, combo_atom_idx, labels).
    """
    rng = np.random.default_rng(seed)
    atoms = rng.normal(size=(n_atoms, n_genes))
    atoms[:, 0] = 0.0  # reserve gene 0 as the "unreachable" planted axis
    idx = np.array([(i, j) for i in range(n_atoms) for j in range(i + 1, n_atoms)])
    combos = atoms[idx].sum(axis=1)
    labels = np.zeros(len(idx), dtype=bool)
    emergent_rows = rng.choice(len(idx), size=n_emergent, replace=False)
    for k in emergent_rows:
        combos[k] = combos[k].copy()
        combos[k][0] += 6.0  # inject the planted axis
        labels[k] = True
    return atoms, combos, idx, labels


# --------------------------------------------------------------------------- #
# replay_campaign — property invariants
# --------------------------------------------------------------------------- #
def test_campaign_curve_is_monotone_and_bounded():
    atoms, combos, idx, labels = _planted_screen()
    for policy in sl.ACQUISITION_POLICIES:
        res = sl.replay_campaign(atoms, combos, idx, labels, policy=policy, batch_size=3)
        # wells strictly increase; found is non-decreasing and never exceeds total
        assert np.all(np.diff(res.wells) > 0)
        assert np.all(np.diff(res.found) >= 0)
        assert res.found[0] == 0 and res.wells[0] == 0
        assert res.found[-1] == labels.sum()   # every combo eventually measured
        assert res.wells[-1] == len(labels)
        assert res.n_emergent == int(labels.sum())


def test_campaign_measures_every_combo_exactly_once():
    atoms, combos, idx, labels = _planted_screen()
    res = sl.replay_campaign(atoms, combos, idx, labels, policy="triage", batch_size=4)
    assert sorted(res.pick_order.tolist()) == list(range(len(labels)))


def test_campaign_is_deterministic():
    atoms, combos, idx, labels = _planted_screen()
    a = sl.replay_campaign(atoms, combos, idx, labels, policy="cone_adaptive", batch_size=3, seed=1)
    b = sl.replay_campaign(atoms, combos, idx, labels, policy="cone_adaptive", batch_size=3, seed=1)
    assert np.array_equal(a.pick_order, b.pick_order)
    assert np.array_equal(a.found, b.found)


def test_campaign_random_seed_changes_order_but_not_total():
    atoms, combos, idx, labels = _planted_screen()
    a = sl.replay_campaign(atoms, combos, idx, labels, policy="random", batch_size=3, seed=0)
    b = sl.replay_campaign(atoms, combos, idx, labels, policy="random", batch_size=3, seed=7)
    assert a.found[-1] == b.found[-1] == labels.sum()
    # different seed -> generally different acquisition order
    assert not np.array_equal(a.pick_order, b.pick_order)


def test_campaign_wells_to_fraction():
    atoms, combos, idx, labels = _planted_screen()
    res = sl.replay_campaign(atoms, combos, idx, labels, policy="triage", batch_size=3)
    w = res.wells_to_fraction(1.0)
    assert np.isfinite(w) and 0 < w <= len(labels)
    with pytest.raises(rx.InputError):
        res.wells_to_fraction(0.0)
    with pytest.raises(rx.InputError):
        res.wells_to_fraction(1.5)


def test_campaign_custom_policy_callable():
    atoms, combos, idx, labels = _planted_screen()

    def always_first(cone, add_pred, unmeasured, rng):
        # acquire in index order (score = -index)
        return -np.asarray(unmeasured, dtype=float)

    res = sl.replay_campaign(atoms, combos, idx, labels, policy=always_first, batch_size=2)
    assert res.policy == "custom"
    assert res.pick_order[0] == 0  # lowest index acquired first


def test_campaign_grow_cone_flag_runs_both_ways():
    atoms, combos, idx, labels = _planted_screen()
    for grow in (True, False):
        res = sl.replay_campaign(atoms, combos, idx, labels, policy="cone_adaptive",
                                 batch_size=3, grow_cone=grow)
        assert res.grow_cone is grow
        assert res.found[-1] == labels.sum()


# --------------------------------------------------------------------------- #
# replay_campaign — fail closed
# --------------------------------------------------------------------------- #
def test_campaign_rejects_shape_mismatch():
    atoms, combos, idx, labels = _planted_screen()
    with pytest.raises(rx.InputError):
        sl.replay_campaign(atoms, combos[:-1], idx, labels)          # combos/labels mismatch
    with pytest.raises(rx.InputError):
        sl.replay_campaign(atoms, combos, idx[:, :1], labels)        # k<2
    with pytest.raises(rx.InputError):
        sl.replay_campaign(atoms, combos[:, :-1], idx, labels)       # gene mismatch


def test_campaign_rejects_bad_atom_index():
    atoms, combos, idx, labels = _planted_screen()
    bad = idx.copy(); bad[0, 0] = atoms.shape[0] + 5
    with pytest.raises(rx.InputError):
        sl.replay_campaign(atoms, combos, bad, labels)


def test_campaign_rejects_bad_policy_and_batch():
    atoms, combos, idx, labels = _planted_screen()
    with pytest.raises(rx.InputError):
        sl.replay_campaign(atoms, combos, idx, labels, policy="not_a_policy")
    with pytest.raises(rx.InputError):
        sl.replay_campaign(atoms, combos, idx, labels, batch_size=0)


def test_campaign_rejects_nonfinite():
    atoms, combos, idx, labels = _planted_screen()
    bad = atoms.copy(); bad[0, 0] = np.nan
    with pytest.raises(rx.InputError):
        sl.replay_campaign(bad, combos, idx, labels)


# --------------------------------------------------------------------------- #
# nominate_atoms — the library-augmentation engine
# --------------------------------------------------------------------------- #
def test_nominate_recovers_planted_axis():
    # A library missing the gene-0 axis; combos that need it must nominate the one
    # candidate atom that supplies it.
    n_genes = 30
    rng = np.random.default_rng(0)
    cone = rng.normal(size=(5, n_genes)); cone[:, 0] = 0.0     # library lacks gene 0
    supplier = rng.normal(size=n_genes); supplier[0] = 5.0     # the missing axis
    candidates = np.vstack([cone, supplier[None, :]])          # candidate 5 supplies it
    # measured combos that leave the cone along gene 0
    combos = np.array([cone[0] + cone[1], cone[2] + cone[3]])
    combos = combos + np.array([1.0] + [0.0] * (n_genes - 1)) * 4.0   # push along gene 0
    nom = sl.nominate_atoms(cone, combos, candidates)
    assert nom.top(1)[0] == 5           # the supplier is nominated first
    assert nom.n_separating == 2
    assert nom.aggregate_direction.shape == (n_genes,)


def test_nominate_inside_cone_gives_no_separator():
    # A combo that IS inside the cone contributes no separator.
    atoms = np.eye(4)
    inside = atoms[0] + atoms[1]        # exactly reachable
    nom = sl.nominate_atoms(atoms, inside[None, :], atoms)
    assert nom.n_separating == 0
    assert np.allclose(nom.aggregate_direction, 0.0)


def test_nominate_weight_options():
    atoms, combos, idx, labels = _planted_screen()
    for w in ("residual", "uniform"):
        nom = sl.nominate_atoms(atoms, combos[labels], atoms, weight=w)
        assert nom.ranking.shape[0] == atoms.shape[0]
    with pytest.raises(rx.InputError):
        sl.nominate_atoms(atoms, combos, atoms, weight="bogus")


def test_nominate_rejects_shape_mismatch():
    atoms, combos, idx, labels = _planted_screen()
    with pytest.raises(rx.InputError):
        sl.nominate_atoms(atoms, combos[:, :-1], atoms)     # combo gene mismatch
    with pytest.raises(rx.InputError):
        sl.nominate_atoms(atoms, combos, atoms[:, :-1])     # candidate gene mismatch


def test_nominate_is_deterministic():
    atoms, combos, idx, labels = _planted_screen()
    a = sl.nominate_atoms(atoms, combos[labels], atoms)
    b = sl.nominate_atoms(atoms, combos[labels], atoms)
    assert np.array_equal(a.ranking, b.ranking)
    assert np.allclose(a.scores, b.scores)


# --------------------------------------------------------------------------- #
# held_out_single_recovery — the falsifiable test
# --------------------------------------------------------------------------- #
def test_held_out_recovery_beats_random_null():
    # On a planted screen the separator should recover the held-out single far
    # better than chance (random null top-1 ~ 1/n_atoms).
    atoms, combos, idx, labels = _planted_screen(seed=2, n_atoms=8)
    rec = sl.held_out_single_recovery(atoms, idx, combos, min_involved=2)
    n = rec["n_candidates"]
    random_top1 = 1.0 / n
    assert rec["separator"]["top1"] > 5 * random_top1
    # separator ranks the true atom at least as well as the naive baseline on average
    assert rec["separator"]["mean_rank"] <= rec["baseline"]["mean_rank"] + 1e-9
    # ranks are 1-based and within range
    assert rec["sep_ranks"].min() >= 1 and rec["sep_ranks"].max() <= n


def test_held_out_recovery_requires_eligible_atoms():
    # Two atoms, one double -> no atom appears in >=2 combos.
    atoms = np.eye(3)[:2]
    idx = np.array([[0, 1]])
    combos = atoms.sum(axis=0)[None, :]
    with pytest.raises(rx.InputError):
        sl.held_out_single_recovery(atoms, idx, combos, min_involved=2)


def test_held_out_recovery_summary_shape():
    atoms, combos, idx, labels = _planted_screen()
    rec = sl.held_out_single_recovery(atoms, idx, combos, min_involved=2)
    assert set(rec["separator"]) == {"median_rank", "mean_rank", "top1", "top5"}
    assert len(rec["sep_ranks"]) == len(rec["base_ranks"]) == len(rec["eligible"])
    assert "model-relative" in rec["scope"]


# --------------------------------------------------------------------------- #
# scope discipline (inherited contract)
# --------------------------------------------------------------------------- #
def test_results_carry_scope():
    atoms, combos, idx, labels = _planted_screen()
    camp = sl.replay_campaign(atoms, combos, idx, labels, policy="triage")
    nom = sl.nominate_atoms(atoms, combos[labels], atoms)
    assert "model-relative" in camp.scope
    assert "model-relative" in nom.scope
