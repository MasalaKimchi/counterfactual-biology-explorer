from dataclasses import replace

import numpy as np
import pytest

from reachability import InputError, project_cone
from validation import (
    LabeledEffects,
    LabeledTarget,
    Provenance,
    active_set_oracle,
    align_labeled_problem,
    grouped_gene_splits,
    max_t_empirical_p,
)


@pytest.fixture
def labelled_problem():
    provenance = Provenance(
        dataset_id="effects-v1",
        source_sha256="a" * 64,
        gene_namespace="Ensembl",
        units="z_score",
        modality="CRISPRi",
        context="Rest",
        timepoint="48h",
        orientation="perturbations_by_genes",
    )
    effects = LabeledEffects(
        values=np.eye(3),
        perturbations=("P1", "P2", "P3"),
        genes=("G1", "G2", "G3"),
        provenance=provenance,
    )
    target = LabeledTarget(
        values=np.array([3.0, 2.0, 1.0]),
        genes=("G3", "G2", "G1"),
        provenance=replace(
            provenance,
            dataset_id="target-v1",
            source_sha256="b" * 64,
            modality="bulk-RNA contrast",
            orientation="genes",
        ),
    )
    return effects, target


def test_label_alignment_uses_names_not_positions(labelled_problem):
    effects, target = labelled_problem
    aligned = align_labeled_problem(effects, target)
    assert aligned.genes == effects.genes
    np.testing.assert_array_equal(aligned.target, [1.0, 2.0, 3.0])
    assert aligned.effects_gene_fraction == aligned.target_gene_fraction == 1.0


def test_missing_genes_require_explicit_coverage_gate(labelled_problem):
    effects, target = labelled_problem
    reduced = replace(target, values=target.values[:2], genes=target.genes[:2])
    with pytest.raises(InputError, match="intersection"):
        align_labeled_problem(effects, reduced)
    aligned = align_labeled_problem(
        effects,
        reduced,
        allow_intersection=True,
        min_shared_genes=2,
        min_shared_fraction=2 / 3,
    )
    assert aligned.genes == ("G2", "G3")
    assert aligned.effects_gene_fraction == pytest.approx(2 / 3)
    assert aligned.target_gene_fraction == 1.0


@pytest.mark.parametrize("field,value", [("gene_namespace", "HGNC"), ("units", "log2_fc")])
def test_incompatible_provenance_fails(labelled_problem, field, value):
    effects, target = labelled_problem
    with pytest.raises(InputError):
        align_labeled_problem(
            effects, replace(target, provenance=replace(target.provenance, **{field: value}))
        )


def test_active_set_oracle_matches_production_fitted_point():
    effects = np.array([[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 2.0]])
    target = np.array([1.0, -0.5, 0.3])
    production = project_cone(effects, target)
    _, fitted, relative_objective = active_set_oracle(effects, target)
    np.testing.assert_allclose(fitted, production.fitted, atol=1e-10)
    assert relative_objective == pytest.approx(production.relative_objective, abs=1e-12)


def test_active_set_oracle_handles_scaled_rank_deficient_atoms():
    rng = np.random.default_rng(17)
    effects = rng.normal(size=(5, 12)) * np.geomspace(1e-12, 1e12, 5)[:, None]
    effects[-1] = effects[-2]
    target = rng.normal(size=12)
    production = project_cone(effects, target)
    _, fitted, relative_objective = active_set_oracle(effects, target)
    np.testing.assert_allclose(fitted, production.fitted, rtol=1e-8, atol=1e-10)
    assert relative_objective == pytest.approx(production.relative_objective, abs=1e-10)


def test_grouped_splits_never_split_a_group():
    groups = ("a", "a", "b", "b", "c", "c", "d", "d")
    first = grouped_gene_splits(groups, n_splits=5, seed=9)
    second = grouped_gene_splits(groups, n_splits=5, seed=9)
    for (fit, score), (fit_again, score_again) in zip(first, second):
        np.testing.assert_array_equal(fit, fit_again)
        np.testing.assert_array_equal(score, score_again)
        assert {groups[index] for index in fit}.isdisjoint(
            {groups[index] for index in score}
        )


def test_max_t_is_plus_one_and_no_less_conservative_than_marginal():
    observed = np.array([2.0, 1.0])
    null = np.array([[0.0, 3.0], [2.0, 0.0], [1.0, 1.0]])
    adjusted = max_t_empirical_p(observed, null)
    marginal = np.array([2 / 4, 3 / 4])
    np.testing.assert_allclose(adjusted, [3 / 4, 1.0])
    assert np.all(adjusted >= marginal)
