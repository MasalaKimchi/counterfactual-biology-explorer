#!/usr/bin/env python3
"""Run the deterministic, data-free systemic validation harness."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError, held_out_alignment, project_cone
from validation import (
    LabeledEffects,
    LabeledTarget,
    Provenance,
    active_set_oracle,
    align_labeled_problem,
    grouped_gene_splits,
    max_t_empirical_p,
)


DEFAULT_CONFIG = ROOT / "configs" / "validation_harness.json"
DEFAULT_REPORT = ROOT / "results" / "validation_harness.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _gate(value: float, threshold: float) -> bool:
    return bool(np.isfinite(value) and value <= threshold)


def _oracle_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    fitted_errors = []
    objective_errors = []
    for trial in range(config["oracle_trials"]):
        n_atoms = int(rng.integers(2, 8))
        n_genes = int(rng.integers(n_atoms + 2, 22))
        effects = rng.normal(size=(n_atoms, n_genes))
        if trial % 3 == 0:
            effects *= np.geomspace(1e-12, 1e12, n_atoms)[:, None]
        if trial % 7 == 0 and n_atoms > 2:
            effects[-1] = effects[-2]
        weights = np.exp(rng.uniform(-2, 2, size=n_genes))
        if trial % 2:
            target = rng.uniform(0, 2, size=n_atoms) @ effects
        else:
            target = rng.normal(size=n_genes)
        production = project_cone(effects, target, gene_weights=weights)
        _, fitted, relative_objective = active_set_oracle(
            effects, target, gene_weights=weights, max_atoms=7
        )
        scale = max(float(np.linalg.norm(target)), 1.0)
        fitted_errors.append(float(np.linalg.norm(production.fitted - fitted) / scale))
        objective_errors.append(
            float(abs(production.relative_objective - relative_objective))
        )
    max_fitted = max(fitted_errors)
    max_objective = max(objective_errors)
    thresholds = config["thresholds"]
    passed = _gate(max_fitted, thresholds["max_oracle_fitted_error"]) and _gate(
        max_objective, thresholds["max_oracle_objective_error"]
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "trials": config["oracle_trials"],
        "max_relative_fitted_error": max_fitted,
        "max_relative_objective_error": max_objective,
    }


def _degeneracy_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    effects = rng.normal(size=(6, 30))
    target = rng.normal(size=30)
    base = project_cone(effects, target).fitted
    scales = np.geomspace(1e-8, 1e8, effects.shape[0])
    scaled = project_cone(effects * scales[:, None], target).fitted
    augmented = np.vstack([effects, effects[2], effects[2] * (1 + 1e-12), np.zeros(30)])
    duplicated = project_cone(augmented, target).fitted
    scale = max(float(np.linalg.norm(target)), 1.0)
    errors = [
        float(np.linalg.norm(base - scaled) / scale),
        float(np.linalg.norm(base - duplicated) / scale),
    ]
    maximum = max(errors)
    passed = _gate(maximum, config["thresholds"]["max_degenerate_fitted_error"])
    return {
        "status": "PASS" if passed else "FAIL",
        "max_relative_fitted_error": maximum,
        "challenges": ["atom_rescaling_1e-8_to_1e8", "duplicate_near_duplicate_zero"],
    }


def _label_scenario() -> dict[str, Any]:
    provenance = Provenance(
        dataset_id="synthetic-effects-v1",
        source_sha256="0" * 64,
        gene_namespace="HGNC-symbol",
        units="z_score",
        modality="synthetic-CRISPRi",
        context="synthetic-rest",
        timepoint="synthetic",
        orientation="perturbations_by_genes",
    )
    target_provenance = replace(
        provenance,
        dataset_id="synthetic-target-v1",
        source_sha256="1" * 64,
        modality="RNA-contrast",
        orientation="genes",
    )
    effects = LabeledEffects(
        values=np.eye(4),
        perturbations=("p1", "p2", "p3", "p4"),
        genes=("G1", "G2", "G3", "G4"),
        provenance=provenance,
    )
    target = LabeledTarget(
        values=np.array([1.0, 2.0, 3.0, 4.0]),
        genes=("G1", "G2", "G3", "G4"),
        provenance=target_provenance,
    )
    reordered = replace(target, values=target.values[::-1], genes=target.genes[::-1])
    aligned = align_labeled_problem(effects, reordered)
    reorder_safe = bool(np.array_equal(aligned.target, target.values))

    corruptions: list[Callable[[], object]] = [
        lambda: align_labeled_problem(replace(effects, genes=("G1", "G1", "G3", "G4")), target),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, gene_namespace="Ensembl"))),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, units="log2_fc"))),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, orientation="genes_by_one"))),
        lambda: align_labeled_problem(effects, replace(target, genes=target.genes[:-1], values=target.values[:-1])),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, source_sha256="unknown"))),
        lambda: align_labeled_problem(replace(effects, values=np.where(np.eye(4) == 1, np.nan, 0)), target),
        lambda: align_labeled_problem(replace(effects, values=np.eye(3)), target),
    ]
    caught = 0
    for corruption in corruptions:
        try:
            corruption()
        except InputError:
            caught += 1
    passed = reorder_safe and caught == len(corruptions)
    return {
        "status": "PASS" if passed else "FAIL",
        "safe_reordering": reorder_safe,
        "faults_caught": caught,
        "faults_injected": len(corruptions),
    }


def _grouped_split_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    n_groups = 14
    genes_per_group = 6
    groups = tuple(f"module_{group}" for group in range(n_groups) for _ in range(genes_per_group))
    module_effects = rng.normal(size=(7, n_groups))
    effects = np.repeat(module_effects, genes_per_group, axis=1) + rng.normal(
        scale=0.05, size=(7, n_groups * genes_per_group)
    )
    target = rng.uniform(0.2, 1.5, size=7) @ effects + rng.normal(
        scale=0.1, size=n_groups * genes_per_group
    )
    scores = []
    leakage = 0
    for fit, score in grouped_gene_splits(
        groups, n_splits=config["grouped_splits"], seed=config["seed"] + 1
    ):
        fit_groups = {groups[index] for index in fit}
        score_groups = {groups[index] for index in score}
        leakage += len(fit_groups & score_groups)
        scores.append(held_out_alignment(effects, target, fit, score).held_out_cosine)
    passed = leakage <= config["thresholds"]["max_group_leakage"] and np.all(
        np.isfinite(scores)
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "splits": config["grouped_splits"],
        "group_leakage_count": leakage,
        "held_out_cosine_mean": float(np.mean(scores)),
        "held_out_cosine_sd": float(np.std(scores, ddof=1)),
        "held_out_cosine_min": float(np.min(scores)),
    }


def _null_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    rejections = 0
    adjusted_not_below_marginal = True
    for _ in range(config["null_trials"]):
        observed = rng.normal(size=config["n_hypotheses"])
        null = rng.normal(
            size=(config["null_resamples"], config["n_hypotheses"])
        )
        adjusted = max_t_empirical_p(observed, null)
        marginal = np.asarray(
            [
                (1 + np.count_nonzero(null[:, index] >= value))
                / (config["null_resamples"] + 1)
                for index, value in enumerate(observed)
            ]
        )
        adjusted_not_below_marginal &= bool(np.all(adjusted >= marginal))
        rejections += int(np.any(adjusted <= 0.05))
    rate = rejections / config["null_trials"]
    passed = adjusted_not_below_marginal and _gate(
        rate, config["thresholds"]["max_null_familywise_error_rate"]
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "trials": config["null_trials"],
        "hypotheses_per_trial": config["n_hypotheses"],
        "resamples_per_trial": config["null_resamples"],
        "familywise_error_rate_at_0.05": rate,
        "adjusted_p_never_below_marginal": adjusted_not_below_marginal,
    }


def run(config_path: Path) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "1.0.0":
        raise ValueError("unsupported harness config schema")
    rng = np.random.default_rng(config["seed"])
    scenarios = {
        "active_set_oracle": _oracle_scenario(config, rng),
        "degenerate_cones": _degeneracy_scenario(config, rng),
        "label_and_provenance_faults": _label_scenario(),
        "grouped_gene_holdout": _grouped_split_scenario(config, rng),
        "maxT_null_calibration": _null_scenario(config, rng),
    }
    status = "PASS" if all(item["status"] == "PASS" for item in scenarios.values()) else "FAIL"
    return {
        "schema_version": "1.0.0",
        "generated_on": "2026-07-17",
        "config_sha256": _sha256(config_path),
        "status": status,
        "scope": "deterministic synthetic software/statistical contract; not biological validation",
        "scenarios": scenarios,
    }


def _assert_matches(actual: Any, expected: Any, path: str = "report") -> None:
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise AssertionError(f"{path} keys differ")
        for key in expected:
            _assert_matches(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, float):
        if not np.isclose(actual, expected, rtol=1e-9, atol=1e-12):
            raise AssertionError(f"{path}: {actual} != {expected}")
    elif actual != expected:
        raise AssertionError(f"{path}: {actual} != {expected}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = run(args.config)
    if args.write:
        args.write.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {args.write}")
    elif args.check:
        expected = json.loads(args.check.read_text(encoding="utf-8"))
        _assert_matches(report, expected)
        print(f"validation harness: {report['status']} (frozen report matches)")
    else:
        print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
