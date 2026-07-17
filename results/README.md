# Results

This directory is the complete maintained result bundle.

| Path | Role |
|---|---|
| `findings.json` | Canonical machine-readable values, interpretations, and open requirements |
| `manifest.json` | SHA-256 and byte length for every maintained artifact |
| `validation_harness.json` | Frozen data-free systemic stress-test report |
| `evidence/` | Minimal tables that directly support the updated findings |

The public narrative is [`docs/FINDINGS.md`](../docs/FINDINGS.md). Values displayed in
the README and central figure must first appear in `findings.json`.

## Evidence groups

### Headline and baseline

- `headline_heldout_split_stability.csv`
- `historical_fixed_split_null.csv`
- `baseline_comparison.csv`
- `metric_calibration.csv`
- `metric_calibration_provenance.json`

### Context and source construction

- `context_condition_comparison.csv`
- `context_runbalance_caveat.json`
- `reviewer2_ota_hollbacher_meta.json`
- `reviewer2_ota_hollbacher_split.csv`

### Specific robustness challenges

- `confounder_robustness_summary.json`
- `generator_significance_summary.json`
- `ranking_validation_summary.json`
- `reviewer2_deg_survival.csv`

### Combination and transport diagnostics

- `combination_additivity_sensitivity.json`

## Policy

- The synthetic reachable-by-construction calibration is not a biological positive control.
- No recipe, druggability, causal-oracle, clinical, or activation-certificate artifacts are
  maintained.
- Large data and intermediate caches are external to the repository.
- `python scripts/validate_findings.py` fails if values, paths, or hashes drift.
