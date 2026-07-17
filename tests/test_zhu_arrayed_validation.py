import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from reachability import InputError
from scripts.run_zhu_arrayed_validation import (
    align_profiles,
    build_flow_effects,
    conditional_permutations,
    load_bulk,
    load_flow,
    load_screen,
    profile_metrics,
    summarize_flow,
    verify_input,
)


TARGETS = ("A", "B", "C")


def test_input_identity_enforces_bytes_and_hash(tmp_path, monkeypatch):
    path = tmp_path / "input.csv"
    path.write_bytes(b"frozen\n")
    spec = {
        "path": "input.csv",
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    monkeypatch.setattr("scripts.run_zhu_arrayed_validation.ROOT", tmp_path)
    observed = verify_input(spec)
    assert observed["hash_verified"] is True
    path.write_bytes(b"changed")
    with pytest.raises(InputError, match="SHA-256"):
        verify_input(spec)


def _write_tiny_screen(path: Path):
    h5py = pytest.importorskip("h5py")
    string = h5py.string_dtype("utf-8")
    with h5py.File(path, "w") as handle:
        var = handle.create_group("var")
        var.create_dataset("gene_ids", data=np.asarray(["g1", "g2"], object), dtype=string)
        var.create_dataset("gene_name", data=np.asarray(["A", "X"], object), dtype=string)
        obs = handle.create_group("obs")
        condition = obs.create_group("condition")
        condition.create_dataset("categories", data=np.asarray(["Stim", "Other"], object), dtype=string)
        condition.create_dataset("codes", data=np.asarray([0, 0, 1], dtype=np.int8))
        perturbation = obs.create_group("perturbation")
        perturbation.create_dataset("categories", data=np.asarray(["A", "B", "C"], object), dtype=string)
        perturbation.create_dataset("codes", data=np.asarray([0, 1, 2], dtype=np.int8))
        # An outcome-derived admission field is deliberately absent: the loader must
        # define this panel only from frozen labels and condition.
        layers = handle.create_group("layers")
        layers.create_dataset("log_fc", data=np.asarray([[1.0, 2.0], [3.0, 4.0], [9.0, 9.0]]))


def test_screen_panel_selection_does_not_require_outcome_admission(tmp_path, monkeypatch):
    path = tmp_path / "tiny.h5ad"
    _write_tiny_screen(path)
    config = {
        "inputs": {"screen": {"path": "tiny.h5ad"}},
        "screen": {
            "condition": "Stim",
            "condition_field": "condition",
            "perturbation_field": "perturbation",
            "gene_id_field": "gene_ids",
            "gene_symbol_field": "gene_name",
            "layer": "log_fc",
        },
        "panel": {"perturbations": ["A", "B"]},
        "expected": {"screen_rows": 3, "screen_genes": 2},
    }
    monkeypatch.setattr("scripts.run_zhu_arrayed_validation.ROOT", tmp_path)
    loaded = load_screen(config)
    assert loaded["row_indices"] == [0, 1]
    np.testing.assert_array_equal(loaded["profiles"], [[1.0, 2.0], [3.0, 4.0]])


def test_bulk_loader_rejects_duplicate_target_gene(tmp_path, monkeypatch):
    path = tmp_path / "bulk.csv"
    path.write_text(
        "variable,log_fc,contrast\ng1,1,A\ng1,2,A\n",
        encoding="utf-8",
    )
    config = {
        "inputs": {"bulk_rna": {"path": "bulk.csv"}},
        "panel": {"perturbations": ["A"]},
        "expected": {"bulk_rows": 2, "bulk_genes_per_perturbation": 2},
    }
    monkeypatch.setattr("scripts.run_zhu_arrayed_validation.ROOT", tmp_path)
    with pytest.raises(InputError, match="duplicate"):
        load_bulk(config)


def _flow_config():
    return {
        "inputs": {"flow": {"path": "flow.csv"}},
        "panel": {
            "perturbations": ["A"],
            "cytokines": {"IL10": "IL10_perc", "IL21": "IL21_perc"},
            "control_prefix": "NTC",
        },
        "expected": {"flow_rows": 3, "donors": ["D1"]},
    }


def test_flow_uses_mean_of_within_donor_controls(tmp_path, monkeypatch):
    path = tmp_path / "flow.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Sample", "IL10_perc", "IL21_perc", "Donor", "Perturbation"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {"Sample": "c1", "IL10_perc": 1, "IL21_perc": 2, "Donor": "D1", "Perturbation": "NTC1"},
                {"Sample": "c2", "IL10_perc": 3, "IL21_perc": 6, "Donor": "D1", "Perturbation": "NTC2"},
                {"Sample": "a", "IL10_perc": 4, "IL21_perc": 2, "Donor": "D1", "Perturbation": "A"},
            ]
        )
    monkeypatch.setattr("scripts.run_zhu_arrayed_validation.ROOT", tmp_path)
    rows = load_flow(_flow_config())
    effects = build_flow_effects(rows, _flow_config())
    il10 = next(row for row in effects if row["cytokine"] == "IL10")
    il21 = next(row for row in effects if row["cytokine"] == "IL21")
    assert il10["donor_ntc_mean_percent_positive"] == 2.0
    assert il10["log2_ratio_to_donor_ntc"] == 1.0
    assert il21["donor_ntc_mean_percent_positive"] == 4.0
    assert il21["log2_ratio_to_donor_ntc"] == -1.0


def test_flow_rejects_values_above_one_hundred(tmp_path, monkeypatch):
    path = tmp_path / "flow.csv"
    path.write_text(
        "Sample,IL10_perc,IL21_perc,Donor,Perturbation\n"
        "c1,1,2,D1,NTC1\n"
        "c2,3,6,D1,NTC2\n"
        "a,101,2,D1,A\n",
        encoding="utf-8",
    )
    monkeypatch.setattr("scripts.run_zhu_arrayed_validation.ROOT", tmp_path)
    with pytest.raises(InputError, match=r"\(0, 100\]"):
        load_flow(_flow_config())


def test_profile_metrics_mask_on_target_and_retrieve_identity():
    screen = {
        "targets": TARGETS,
        "genes": ("gA", "gB", "gC", "gX", "gY", "gZ"),
        "symbols": ("A", "B", "C", "X", "Y", "Z"),
        "profiles": np.asarray(
            [
                [10.0, 1.0, 0.0, 1.0, 0.2, 0.1],
                [0.0, 10.0, 1.0, 0.1, 1.0, 0.3],
                [1.0, 0.0, 10.0, 0.4, 0.2, 1.0],
            ]
        ),
    }
    bulk = {
        target: dict(zip(screen["genes"], row))
        for target, row in zip(TARGETS, screen["profiles"])
    }
    config = {
        "analysis": {"mask_all_panel_target_genes": True},
        "expected": {"common_genes": 6, "scored_genes_per_profile": 3},
    }
    aligned = align_profiles(screen, bulk, config)
    rows, _, _ = profile_metrics(aligned, TARGETS, "raw")
    assert all(row["genes_scored"] == 3 for row in rows)
    assert all(row["retrieval_rank"] == 1 for row in rows)
    assert all(row["cosine"] == pytest.approx(1.0) for row in rows)


def test_profile_metrics_reject_retrieval_ties():
    aligned = {
        "source": np.ones((3, 3)),
        "arrayed": np.ones((3, 3)),
        "mask": np.ones((3, 3), dtype=bool),
    }
    with pytest.raises(InputError, match="retrieval similarities contain a tie"):
        profile_metrics(aligned, TARGETS, "raw")


def test_exact_conditional_permutation_is_deterministic():
    similarities = np.eye(3)
    ranks = np.asarray([[1, 2, 3], [2, 1, 3], [2, 3, 1]])
    effects = []
    for cytokine in ("IL10", "IL21"):
        for index, target in enumerate(TARGETS):
            effects.append(
                {
                    "cytokine": cytokine,
                    "target": target,
                    "log2_ratio_to_donor_ntc": float(index),
                }
            )
    flow = summarize_flow(effects, TARGETS)
    rna = {
        "screen_IL10": np.arange(3.0),
        "bulk_IL10": np.arange(3.0),
        "screen_IL21": np.arange(3.0),
        "bulk_IL21": np.arange(3.0),
    }
    first = conditional_permutations({"raw": (similarities, ranks)}, rna, flow, TARGETS)
    second = conditional_permutations({"raw": (similarities, ranks)}, rna, flow, TARGETS)
    assert first == second
    assert first["permutations"] == 6
    assert first["retrieval"]["raw"]["top1"] == 1.0
    assert first["retrieval"]["raw"]["exhaustive_top1_tail_fraction"] == pytest.approx(1 / 6)


def test_profile_alignment_requires_every_panel_target_gene():
    screen = {
        "targets": ("A", "MISSING"),
        "genes": ("gA", "gX"),
        "symbols": ("A", "X"),
        "profiles": np.asarray([[1.0, 0.0], [0.0, 1.0]]),
    }
    bulk = {
        target: {"gA": row[0], "gX": row[1]}
        for target, row in zip(screen["targets"], screen["profiles"])
    }
    config = {
        "analysis": {"mask_all_panel_target_genes": True},
        "expected": {"common_genes": 2, "scored_genes_per_profile": 0},
    }
    with pytest.raises(InputError, match="every panel target gene"):
        align_profiles(screen, bulk, config)


def test_conditional_permutation_rejects_zero_variation():
    similarities = np.eye(3)
    ranks = np.asarray([[1, 2, 3], [2, 1, 3], [2, 3, 1]])
    effects = [
        {"cytokine": cytokine, "target": target, "log2_ratio_to_donor_ntc": 0.0}
        for cytokine in ("IL10", "IL21")
        for target in TARGETS
    ]
    flow = summarize_flow(effects, TARGETS)
    rna = {
        f"{assay}_{cytokine}": np.arange(3.0)
        for assay in ("screen", "bulk")
        for cytokine in ("IL10", "IL21")
    }
    with pytest.raises(InputError, match="nonzero variation"):
        conditional_permutations({"raw": (similarities, ranks)}, rna, flow, TARGETS)
