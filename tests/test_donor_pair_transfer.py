from pathlib import Path

import numpy as np
import pytest

from reachability import InputError
from scripts.run_donor_pair_transfer import (
    fit_training_models,
    load_modality_matrix,
    profile_h5mu,
    run_challenge,
    score_frozen_models,
    validate_partitions,
)


def _write_tiny_h5mu(path: Path, *, include_second=True, nonfinite=False):
    h5py = pytest.importorskip("h5py")
    string = h5py.string_dtype("utf-8")
    with h5py.File(path, "w") as handle:
        handle.attrs["encoding-type"] = "MuData"
        mod = handle.create_group("mod")
        names = ["D1_D2", "D3_D4"] if include_second else ["D1_D2"]
        for index, name in enumerate(names):
            group = mod.create_group(name)
            group.attrs["encoding-type"] = "anndata"
            obs = group.create_group("obs")
            obs.attrs["_index"] = "key"
            key = ["B_Rest", "A_Rest", "C_Stim"] if index == 0 else ["A_Rest", "C_Stim", "B_Rest"]
            obs.create_dataset("key", data=np.asarray(key, object), dtype=string)
            condition = obs.create_group("condition")
            condition.create_dataset(
                "categories", data=np.asarray(["Rest", "Stim"], object), dtype=string
            )
            codes = [0, 0, 1] if index == 0 else [0, 1, 0]
            condition.create_dataset("codes", data=np.asarray(codes, dtype=np.int8))
            var = group.create_group("var")
            var.create_dataset("gene_ids", data=np.asarray(["g1", "g2"], object), dtype=string)
            var.create_dataset("gene_name", data=np.asarray(["G1", "G2"], object), dtype=string)
            layers = group.create_group("layers")
            # Values encode key identity: A=[1,0], B=[0,1], C=[2,2].
            values = (
                np.asarray([[0.0, 1.0], [1.0, 0.0], [2.0, 2.0]])
                if index == 0
                else np.asarray([[1.0, 0.0], [2.0, 2.0], [0.0, 1.0]])
            )
            if nonfinite and index == 0:
                values[0, 0] = np.nan
            for layer in ("log_fc", "lfcSE", "zscore"):
                layers.create_dataset(layer, data=values)


def _tiny_config():
    return {
        "h5mu": {
            "modalities": ["D1_D2", "D3_D4"],
            "condition": "Rest",
            "condition_field": "condition",
            "gene_id_field": "gene_ids",
            "gene_symbol_field": "gene_name",
            "required_layers": ["log_fc", "lfcSE", "zscore"],
        },
        "expected": {
            "genes": 2,
            "rest_atoms": 2,
            "modality_rows": {"D1_D2": 3, "D3_D4": 3},
        },
    }


def test_h5mu_profile_aligns_row_order_by_atom_key(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path)
    profile = profile_h5mu(path, _tiny_config())
    assert profile["rest_atom_keys"] == ("A_Rest", "B_Rest")
    columns = np.asarray([0, 1])
    first = load_modality_matrix(
        path, "D1_D2", profile["modalities"]["D1_D2"]["rest_rows"], columns, "log_fc"
    )
    second = load_modality_matrix(
        path, "D3_D4", profile["modalities"]["D3_D4"]["rest_rows"], columns, "log_fc"
    )
    np.testing.assert_array_equal(first, [[1.0, 0.0], [0.0, 1.0]])
    np.testing.assert_array_equal(first, second)


def test_h5mu_profile_rejects_missing_modality(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, include_second=False)
    with pytest.raises(InputError, match="modalities differ"):
        profile_h5mu(path, _tiny_config())


def test_h5mu_profile_rejects_nonfinite_required_layer(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, nonfinite=True)
    with pytest.raises(InputError, match="non-finite"):
        profile_h5mu(path, _tiny_config())


def test_frozen_models_do_not_depend_on_test_dictionary_or_target():
    train = np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    target = np.asarray([1.0, 0.5])
    keys = ("A", "B", "C")
    first, first_models = run_challenge(
        train, train * 2, target, np.asarray([0.2, -0.1]), keys
    )
    second, second_models = run_challenge(
        train, train[:, ::-1] * -7, target, np.asarray([-9.0, 4.0]), keys
    )
    np.testing.assert_array_equal(
        first_models["cone_coefficients"], second_models["cone_coefficients"]
    )
    assert first["training"] == second["training"]
    assert first["metrics"] != second["metrics"]


def test_training_only_baselines_apply_same_atom_and_scalar_to_test():
    train = np.asarray([[2.0, 0.0], [0.0, 1.0]])
    target = np.asarray([1.0, 0.0])
    models = fit_training_models(train, target, ("A", "B"))
    assert models["best_single_key"] == "A"
    assert models["best_single_alpha"] == pytest.approx(0.5)
    test = np.asarray([[4.0, 2.0], [100.0, -100.0]])
    metrics = score_frozen_models(models, test, np.asarray([2.0, 1.0]))
    assert metrics["training_best_single"]["normalized_rmse"] == pytest.approx(0.0)


def test_partition_labels_are_derived_from_donor_run_metadata():
    config = {
        "donors": {
            "D1": {"id": "A", "rest_run": "R1"},
            "D2": {"id": "B", "rest_run": "R1"},
            "D3": {"id": "C", "rest_run": "R2"},
            "D4": {"id": "D", "rest_run": "R2"},
        },
        "h5mu": {"modalities": ["A_B", "A_C", "A_D", "B_C", "B_D", "C_D"]},
        "partitions": [
            {"id": "ab_cd", "left": "A_B", "right": "C_D", "run_confounded": True},
            {"id": "ac_bd", "left": "A_C", "right": "B_D", "run_confounded": False},
            {"id": "ad_bc", "left": "A_D", "right": "B_C", "run_confounded": False},
        ],
    }
    validate_partitions(config)
    config["partitions"][1]["run_confounded"] = True
    with pytest.raises(InputError, match="run-confounding label"):
        validate_partitions(config)
