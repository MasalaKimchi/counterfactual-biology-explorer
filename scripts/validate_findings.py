#!/usr/bin/env python3
"""Fail-closed validation of canonical findings and artifact lineage."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics

from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
EVIDENCE = RESULTS / "evidence"
MANIFEST = RESULTS / "manifest.json"

MANIFEST_PATHS = [
    ".gitignore",
    "LICENSE",
    "README.md",
    "requirements.txt",
    "requirements-external.txt",
    "reachability.py",
    "validation.py",
    "reproduce.sh",
    "configs/validation_harness.json",
    "configs/source_reconstruction.json",
    "configs/arce_external_validation.json",
    "configs/zhu_arrayed_validation.json",
    "configs/donor_pair_transfer.json",
    "data/README.md",
    "data/fetch_de_stats.sh",
    "tests/test_reachability.py",
    "tests/test_validation.py",
    "tests/test_source_reconstruction.py",
    "tests/test_arce_external_validation.py",
    "tests/test_zhu_arrayed_validation.py",
    "tests/test_donor_pair_transfer.py",
    "scripts/run_validation_harness.py",
    "scripts/run_source_reconstruction.py",
    "scripts/run_arce_external_validation.py",
    "scripts/run_zhu_arrayed_validation.py",
    "scripts/run_donor_pair_transfer.py",
    "scripts/validate_findings.py",
    "docs/FINDINGS.md",
    "docs/METHODS.md",
    "docs/SCIENTIFIC_VALIDATION_PLAN.md",
    "docs/VALIDATION_REPORT.md",
    "docs/figures/make_at_a_glance.py",
    "docs/figures/fig_at_a_glance.png",
    "docs/figures/fig_at_a_glance.pdf",
    "results/findings.json",
    "results/validation_harness.json",
    "results/source_reconstruction.json",
    "results/donor_pair_transfer.json",
    "results/README.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_paths() -> list[str]:
    evidence = [
        str(path.relative_to(ROOT))
        for path in sorted(EVIDENCE.iterdir())
        if path.is_file()
    ]
    return sorted(MANIFEST_PATHS + evidence)


def write_manifest() -> None:
    files = {}
    for relative in artifact_paths():
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(relative)
        files[relative] = {
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
            "executable": bool(path.stat().st_mode & 0o111),
        }
    payload = {
        "schema_version": "1.0.0",
        "generated_on": "2026-07-17",
        "files": files,
    }
    MANIFEST.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {MANIFEST.relative_to(ROOT)} ({len(files)} files)")


def assert_close(actual: float, expected: float, tolerance: float = 1e-12) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{actual} != {expected} within {tolerance}")


def _mean_metric(direction: dict, model: str, metric: str) -> float:
    return statistics.mean(
        split["metrics"][model][metric] for split in direction["splits"]
    )


def validate_activation_csv(ledger: dict, metadata: dict) -> None:
    path = EVIDENCE / "arce_activation_guide_effects.csv"
    expected_fields = [
        "target", "guide", "donor", "context", "n_cells",
        "guide_median_supplied_activation_score",
        "ntc_baseline_median_of_guide_medians",
        "supplied_activation_score_delta", "support_flag_lt_20_cells",
        "generator_available", "generator_admitted",
    ]
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != expected_fields:
            raise AssertionError("Arce activation CSV schema differs")
        rows = list(reader)
    if len(rows) != ledger["perturbation_guide_strata"]:
        raise AssertionError("Arce activation CSV row count differs")
    keys = {
        (row["target"], row["guide"], row["donor"], row["context"])
        for row in rows
    }
    if len(keys) != len(rows):
        raise AssertionError("Arce activation CSV has duplicate stratum keys")
    if len({row["target"] for row in rows}) != 28 or len({row["guide"] for row in rows}) != 56:
        raise AssertionError("Arce activation CSV target/guide axes differ")
    if {row["donor"] for row in rows} != {"A", "B"}:
        raise AssertionError("Arce activation CSV donor axis differs")
    contexts = set(ledger["donor_spearman"])
    if {row["context"] for row in rows} != contexts:
        raise AssertionError("Arce activation CSV context axis differs")
    if Counter(row["context"] for row in rows) != Counter({context: 112 for context in contexts}):
        raise AssertionError("Arce activation CSV context counts differ")
    if Counter(row["donor"] for row in rows) != Counter({"A": 224, "B": 224}):
        raise AssertionError("Arce activation CSV donor counts differ")

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    sparse = []
    for row in rows:
        n_cells = int(row["n_cells"])
        observed_delta = float(row["supplied_activation_score_delta"])
        expected_delta = (
            float(row["guide_median_supplied_activation_score"])
            - float(row["ntc_baseline_median_of_guide_medians"])
        )
        assert_close(observed_delta, expected_delta)
        flag = row["support_flag_lt_20_cells"] == "true"
        if flag != (n_cells < 20):
            raise AssertionError("Arce activation support flag differs from cell count")
        if flag:
            sparse.append((row["target"], row["guide"], row["donor"], row["context"], n_cells))
        grouped[(row["target"], row["donor"], row["context"])].append(row)
    if sparse != [("SOCS3", "SOCS3_2_CRISPRi", "A", "Resting-Treg", 8)]:
        raise AssertionError("Arce activation sparse-stratum identity differs")
    if len(grouped) != 224 or any(len(group) != 2 for group in grouped.values()):
        raise AssertionError("Arce activation CSV does not have two guides per target/donor/context")

    targets = sorted({row["target"] for row in rows})
    for context in sorted(contexts):
        donor_values = {}
        all_four_same = []
        for donor in ("A", "B"):
            donor_values[donor] = [
                statistics.median(
                    float(row["supplied_activation_score_delta"])
                    for row in grouped[(target, donor, context)]
                )
                for target in targets
            ]
        observed_rho = float(spearmanr(donor_values["A"], donor_values["B"]).statistic)
        assert_close(observed_rho, ledger["donor_spearman"][context])
        assert_close(
            observed_rho,
            metadata["contexts"][context]["donor_rank_concordance"]["spearman"],
        )
        for target in targets:
            values = [
                float(row["supplied_activation_score_delta"])
                for donor in ("A", "B")
                for row in grouped[(target, donor, context)]
            ]
            all_four_same.append(all(value > 0 for value in values) or all(value < 0 for value in values))
        fraction = statistics.mean(all_four_same)
        assert_close(fraction, ledger["all_two_guides_two_donors_same_sign_fraction"][context])
        assert_close(
            fraction,
            metadata["contexts"][context][
                "all_two_guides_two_donors_same_nonzero_sign_fraction"
            ],
        )


def validate_zhu_arrayed(ledger: dict) -> None:
    metadata_path = EVIDENCE / "zhu_arrayed_validation_meta.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    config_path = ROOT / "configs" / "zhu_arrayed_validation.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if metadata.get("status") != "PASS" or metadata["benchmark"] != ledger["benchmark"]:
        raise AssertionError("Zhu arrayed benchmark did not pass")
    if metadata["config_sha256"] != sha256(config_path):
        raise AssertionError("Zhu arrayed report/config identity differs")
    for name, identity in metadata["input_verification"].items():
        if (
            not identity["hash_verified"]
            or identity["sha256_actual"] != identity["sha256_expected"]
            or identity["sha256_expected"] != config["inputs"][name]["sha256"]
        ):
            raise AssertionError("Zhu arrayed input identity was not fully verified")
    design = metadata["design"]
    assert len(design["perturbations"]) == ledger["perturbations"]
    assert len(design["followup_donor_labels"]) == ledger["followup_donor_labels"]
    assert design["common_genes"] == ledger["common_genes"]
    if (
        not config["analysis"]["mask_all_panel_target_genes"]
        or not design["all_panel_target_genes_masked"]
        or ledger["scored_genes_per_profile"]
        != ledger["common_genes"] - ledger["perturbations"]
    ):
        raise AssertionError("Zhu arrayed panel-target masking contract differs")

    profile_path = EVIDENCE / "zhu_arrayed_profile_metrics.csv"
    with profile_path.open(newline="", encoding="utf-8") as handle:
        profile_rows = list(csv.DictReader(handle))
    if len(profile_rows) != 2 * ledger["perturbations"]:
        raise AssertionError("Zhu arrayed profile row count differs")
    if {(row["scale"], row["target"]) for row in profile_rows} != {
        (scale, target)
        for scale in ("raw", "panel_centered")
        for target in design["perturbations"]
    }:
        raise AssertionError("Zhu arrayed profile keys differ")
    if {int(row["genes_scored"]) for row in profile_rows} != {
        ledger["scored_genes_per_profile"]
    }:
        raise AssertionError("Zhu arrayed profile gene count differs")
    if any(int(row["retrieval_rank"]) != 1 for row in profile_rows):
        raise AssertionError("Zhu arrayed matching profile is not always rank one")

    profile = ledger["profile_replication"]
    raw = metadata["profile_summary"]["raw"]
    centered = metadata["profile_summary"]["panel_centered"]
    retrieval = metadata["conditional_permutation"]["retrieval"]
    assert_close(raw["cosine"], profile["raw_median_cosine"])
    assert_close(
        raw["median_cosine_gain_over_common_source"],
        profile["raw_median_cosine_gain_over_common_source"],
    )
    assert_close(raw["normalized_rmse"], profile["raw_median_normalized_rmse"])
    assert_close(centered["cosine"], profile["panel_centered_median_cosine"])
    assert_close(
        centered["normalized_rmse"], profile["panel_centered_median_normalized_rmse"]
    )
    assert_close(retrieval["raw"]["top1"], profile["raw_top1_retrieval"])
    assert_close(
        retrieval["panel_centered"]["top1"], profile["panel_centered_top1_retrieval"]
    )
    if metadata["conditional_permutation"]["permutations"] != math.factorial(
        ledger["perturbations"]
    ):
        raise AssertionError("Zhu arrayed permutation enumeration differs")
    for name, expected in ledger["rna_to_donor_median_flow_spearman"].items():
        assert_close(
            metadata["conditional_permutation"]["cytokine_rank_association"][name][
                "spearman"
            ],
            expected,
        )

    flow_path = EVIDENCE / "zhu_arrayed_flow_effects.csv"
    with flow_path.open(newline="", encoding="utf-8") as handle:
        flow_rows = list(csv.DictReader(handle))
    if len(flow_rows) != 78:
        raise AssertionError("Zhu arrayed donor/cytokine row count differs")
    keys = {
        (row["target"], row["donor"], row["cytokine"]) for row in flow_rows
    }
    if len(keys) != len(flow_rows):
        raise AssertionError("Zhu arrayed flow rows have duplicate keys")
    coverage = Counter((row["target"], row["cytokine"]) for row in flow_rows)
    if min(coverage.values()) != 3 or max(coverage.values()) != 6:
        raise AssertionError("Zhu arrayed donor coverage differs")
    for row in flow_rows:
        expected = math.log2(
            float(row["percent_positive"])
            / float(row["donor_ntc_mean_percent_positive"])
        )
        assert_close(float(row["log2_ratio_to_donor_ntc"]), expected)


def validate_donor_pair(ledger: dict) -> None:
    report = json.loads((RESULTS / "donor_pair_transfer.json").read_text(encoding="utf-8"))
    config_path = ROOT / "configs" / "donor_pair_transfer.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if report.get("status") != "PASS" or report["benchmark"] != ledger["benchmark"]:
        raise AssertionError("donor-pair benchmark did not pass")
    if report["config_sha256"] != sha256(config_path):
        raise AssertionError("donor-pair report/config identity differs")
    for name, identity in report["input_verification"].items():
        if (
            not identity["hash_verified"]
            or identity["sha256_actual"] != identity["sha256_expected"]
            or identity["sha256_expected"] != config["inputs"][name]["sha256"]
        ):
            raise AssertionError("donor-pair input identity was not fully verified")

    quality = report["data_quality"]
    if len(quality["modalities"]) != ledger["modalities"]:
        raise AssertionError("donor-pair modality count differs")
    if quality["rest_atoms_complete_all_modalities"] != ledger["rest_atoms_complete_all_modalities"]:
        raise AssertionError("donor-pair complete Rest atom count differs")
    if quality["shared_target_genes"] != ledger["shared_target_genes"]:
        raise AssertionError("donor-pair shared target universe differs")
    if "published two-donor DE-eligible" not in quality["eligibility_warning"]:
        raise AssertionError("donor-pair eligibility-selection warning is missing")

    challenges = report["challenges"]
    balanced = [row for row in challenges if not row["run_confounded"]]
    if len(challenges) != ledger["correlated_challenges"] or len(balanced) != ledger["run_balanced_challenges"]:
        raise AssertionError("donor-pair challenge counts differ")
    for row in challenges:
        for field in ("fit_gene_sha256", "score_gene_sha256"):
            if len(row[field]) != 64:
                raise AssertionError("donor-pair gene split identities are not hash-frozen")
    if "select atom/scalar only on training" not in report["protocol"]["baseline_application"]:
        raise AssertionError("donor-pair baselines are not declared training-only")
    if "no donor-population p-values" not in report["protocol"]["inference"]:
        raise AssertionError("donor-pair inference ceiling is missing")

    observed = report["summary"]["run_balanced_only"]
    expected = ledger["run_balanced"]
    if observed["challenges"] != ledger["run_balanced_challenges"]:
        raise AssertionError("donor-pair run-balanced summary count differs")
    for field in (
        "median_cone_cosine",
        "median_cosine_improvement_over_training_best_single",
        "fraction_cosine_improvement_positive_over_training_best_single",
        "median_cone_normalized_rmse",
        "median_training_best_single_normalized_rmse",
        "fraction_normalized_rmse_improvement_positive_over_training_best_single",
    ):
        assert_close(observed[field], expected[field])


def validate_values(findings: dict) -> None:
    if findings.get("schema_version") != "2.0.0":
        raise AssertionError("unsupported findings schema")

    source = json.loads((RESULTS / "source_reconstruction.json").read_text())
    if source.get("status") != "PASS":
        raise AssertionError("source reconstruction did not pass")
    for status in source["data_quality"]["input_verification"].values():
        if not status["hash_verified"] or status["sha256_actual"] != status["sha256_expected"]:
            raise AssertionError("source input identity was not byte-hash verified")

    split_report = source["source_bound_splits"]
    if split_report["status"] != "PASS":
        raise AssertionError("source-bound split contract did not pass")
    ledger = findings["source_bound_alignment"]
    assert len(split_report["splits"]) == ledger["n_fixed_random_gene_splits"]
    observed_split_values = [row["held_out_cosine"] for row in split_report["splits"]]
    if observed_split_values != ledger["split_values"]:
        raise AssertionError("source-bound split values differ")
    assert_close(split_report["held_out_cosine_mean"], ledger["held_out_cosine_mean"])
    assert_close(split_report["held_out_cosine_sd"], ledger["held_out_cosine_sd"])
    assert_close(split_report["held_out_cosine_min"], ledger["held_out_cosine_min"])
    assert_close(split_report["held_out_cosine_max"], ledger["held_out_cosine_max"])
    if split_report["protocol"] != {
        "gene_order": ledger["gene_order"],
        "rng": ledger["split_rng"],
    }:
        raise AssertionError("source-bound split protocol differs")
    for split in split_report["splits"]:
        for field in ("fit_gene_sha256", "score_gene_sha256"):
            if len(split[field]) != 64:
                raise AssertionError("split gene identities are not hash-frozen")

    lineage = source["data_quality"]["target_lineage"]["zscore"]
    counts = lineage["counts"]
    ledger = findings["target_lineage"]
    assert counts["union"] == ledger["union_genes"]
    assert counts["shared"] == ledger["shared_source_genes"]
    assert counts["sign_concordant"] == ledger["sign_concordant_genes"]
    assert counts["shared_screen"] == ledger["shared_screen_genes"]
    assert counts["final"] == ledger["registered_genes"]
    assert_close(
        source["source_transfer"]["log_fc"]["between_source_screen_target_cosine"],
        ledger["between_source_log_fc_cosine_on_shared_screen"],
    )
    assert_close(
        source["source_transfer"]["zscore"]["between_source_screen_target_cosine"],
        ledger["between_source_zscore_cosine_on_shared_screen"],
    )

    transfer = findings["cross_source_directional_transfer"]
    directions = source["source_transfer"]["log_fc"]["directions"]
    for key in ("ota_to_hollbacker", "hollbacker_to_ota"):
        observed = directions[key]
        expected = transfer[key]
        assert len(observed["splits"]) == transfer["n_correlated_splits"]
        assert_close(_mean_metric(observed, "cone", "cosine"), expected["mean_cone_cosine"])
        assert_close(
            observed["cosine_improvement_mean"],
            expected["mean_cosine_improvement_over_test_selected_better_baseline"],
        )
        assert_close(
            _mean_metric(observed, "cone", "normalized_rmse"),
            expected["mean_cone_normalized_rmse"],
        )
        assert_close(
            _mean_metric(observed, "best_single_atom", "normalized_rmse"),
            expected["mean_best_single_normalized_rmse"],
        )

    arce = json.loads((EVIDENCE / "arce_external_validation_meta.json").read_text())
    arce_config_path = ROOT / "configs" / "arce_external_validation.json"
    arce_config = json.loads(arce_config_path.read_text())
    if arce["config_sha256"] != sha256(arce_config_path):
        raise AssertionError("Arce report/config identity differs")
    arce_input = arce["input_verification"]
    arce_members = arce_input["arce"]["members"]
    if (
        arce_input["arce"]["archive_sha256"]
        != arce_config["dataset"]["archive"]["sha256"]
        or arce_input["zhu_generator"]["predictor_sha256"]
        != arce_config["generator"]["predictor_sha256"]
    ):
        raise AssertionError("Arce input identities differ from the frozen contract")
    for key in (
        "screen_member", "activation_cells_member", "activation_summary_member"
    ):
        if arce_members[key]["sha256"] != arce_config["dataset"][key]["sha256"]:
            raise AssertionError(f"Arce {key} identity differs from the frozen contract")
    source_input = source["data_quality"]["input_verification"]["de_stats"]
    if (
        arce_input["zhu_generator"]["path"] != source_input["path"]
        or arce_input["zhu_generator"]["bytes"] != source_input["bytes"]
    ):
        raise AssertionError("Arce predictor is not cross-bound to the verified Zhu H5AD")
    if "unadjusted exploratory" not in arce["permutation_inference"]:
        raise AssertionError("Arce permutation multiplicity ceiling is missing")
    ledger = findings["arce_external_validation"]
    if arce["benchmark"] != ledger["benchmark"]:
        raise AssertionError("Arce benchmark identity differs")
    assert arce["attrition"]["analysis_eligible"] == ledger["analysis_targets"]
    for context, value in ledger["spearman"].items():
        assert_close(arce["contexts"][context]["point"]["spearman"], value)
    for context, value in ledger["top_25_overlap"].items():
        assert arce["contexts"][context]["point"]["top_k"]["25"]["overlap"] == value
    activation = ledger["activation_score_robustness"]
    data_quality = arce["activation_data_quality"]
    reproduction = arce["activation_summary_reproduction"]
    robustness = arce["activation_robustness"]
    assert data_quality["rows"] == activation["s14_cells"]
    assert data_quality["guide_groups"] == activation["guide_strata"]
    for field in ("min", "median", "max"):
        assert_close(
            data_quality[f"guide_group_cells_{field}"],
            arce_config["expected"][f"activation_guide_group_cells_{field}"],
        )
    assert reproduction["rows"] == activation["s8_aggregate_rows_reproduced"]
    if reproduction["published_p_values_used_for_inference"]:
        raise AssertionError("Arce S8 cell-level p-values entered current inference")
    assert_close(
        max(
            reproduction["maximum_absolute_mean_error"],
            reproduction["maximum_absolute_median_error"],
        ),
        activation["maximum_absolute_s8_error"],
    )
    if not robustness["all_strata_retained"] or robustness["cell_level_inference_emitted"]:
        raise AssertionError("Arce activation robustness violated its claim contract")
    assert (
        robustness["strata_below_20_cells"]
        == arce_config["expected"]["activation_strata_below_20_cells"]
    )
    for context, value in activation["donor_spearman"].items():
        assert_close(
            robustness["contexts"][context]["donor_rank_concordance"]["spearman"],
            value,
        )
    for context, value in activation[
        "all_two_guides_two_donors_same_sign_fraction"
    ].items():
        assert_close(
            robustness["contexts"][context][
                "all_two_guides_two_donors_same_nonzero_sign_fraction"
            ],
            value,
        )
    validate_activation_csv(activation, robustness)

    validate_zhu_arrayed(findings["zhu_arrayed_followup"])
    validate_donor_pair(findings["donor_pair_transfer"])

    harness = json.loads((RESULTS / "validation_harness.json").read_text())
    ledger = findings["systemic_harness"]
    if harness["status"] != ledger["status"]:
        raise AssertionError("systemic harness did not pass")
    assert len(harness["scenarios"]) == ledger["scenarios"]
    max_t = harness["scenarios"]["maxT_null_calibration"]
    assert max_t["familywise_rejections"] == ledger["max_t_false_families"]
    assert max_t["trials"] == ledger["max_t_trials"]
    assert_close(max_t["familywise_error_upper_95"], ledger["max_t_exact_upper_95"])

    retired = split_report["retired_archived_comparison"]
    ledger = findings["retired_provenance"]
    assert_close(retired["mean"], ledger["archived_multisplit_mean"])
    assert_close(retired["sd"], ledger["archived_multisplit_sd"])
    assert_close(retired["fixed_split_cosine"], ledger["archived_fixed_split_cosine"])

    for section in findings.values():
        if not isinstance(section, dict):
            continue
        sources = section.get("sources", [section.get("source")])
        for relative in sources:
            if relative and not (ROOT / relative).is_file():
                raise FileNotFoundError(relative)


def validate_manifest() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.0.0":
        raise AssertionError("unsupported manifest schema")
    if sorted(manifest["files"]) != artifact_paths():
        raise AssertionError("manifest path set is stale")
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        if (
            path.stat().st_size != expected["bytes"]
            or sha256(path) != expected["sha256"]
            or bool(path.stat().st_mode & 0o111) != expected["executable"]
        ):
            raise AssertionError(f"artifact mismatch: {relative}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    findings = json.loads((RESULTS / "findings.json").read_text(encoding="utf-8"))
    validate_values(findings)
    if args.write_manifest:
        write_manifest()
    validate_manifest()
    print("findings and artifact lineage: OK")


if __name__ == "__main__":
    main()
