# Validation Report

**Scope:** consolidated maintained surface, 2026-07-17  
**Canonical facts:** [`results/findings.json`](../results/findings.json)  
**Artifact identity:** [`results/manifest.json`](../results/manifest.json)

## What is technically certified

The numerical tests cover exact inside/boundary/outside projections, weighted separator
polarity and orthogonality, scaling extremes, zero/duplicate/near-duplicate atoms, frozen
held-out coefficients, conservative empirical p-values, and fail-closed malformed inputs.
KKT and separator diagnostics certify at `1e-8`; no biological verdict threshold exists.

The six-scenario systemic harness adds an independent active-set oracle, eight injected
axis/provenance faults, scale/degeneracy challenges, grouped splits, exact uncertainty
around maxT familywise error, and structured specificity. It demonstrates that common
response can yield high raw cosine, random-gene splits can be optimistic relative to
module holdout, and sign selection can inflate scores. The maxT fast check has 24/500
false families and an exact one-sided 95% upper bound of 0.067 under the 0.075 gate. This
is synthetic contract certification, not biological validation.

## What is source reconstructed

`scripts/run_source_reconstruction.py` verifies the full bytes of the 16.8 GB Zhu H5AD
and target CSV, validates grain/schema, rebuilds target lineage, hashes every split's
fit/score genes, and performs bidirectional Ota/Höllbacher transfer on log fold change and
z-score. Any hash or frozen-split mismatch makes the report `FAIL`; the canonical validator
rejects failed reports.

The audit retired the previous 0.446 ± 0.010 table and p=1/61 because their deleted
intermediate did not preserve gene order. The source-bound result is 0.444 ± 0.018. A
separate fixed-split value reproduces within `3e-10`, supporting source alignment without
reconstructing the retired pipeline.

## What is independently evaluated

`scripts/run_donor_pair_transfer.py` verifies the 16.9 GB donor H5MU and target-table
bytes, validates all six two-donor modalities, freezes gene splits, and applies NNLS
weights plus training-selected baselines unchanged to the complementary donor pair and
opposite target source. Tests cover H5 schema corruption, axis alignment, training-only
baseline selection, run-balance classification, and deterministic output. The 24
run-balanced challenges show weak directional gain but worse magnitude error. Released
presence is DE-eligibility-selected, so this certifies a published-pipeline fixed-cohort
sensitivity—not leakage-free donor generality.

`scripts/run_arce_external_validation.py` verifies Arce archive/member hashes and schema,
enforces outcome-independent eligibility, extracts the Zhu `Rest` IL2RA predictor, and
evaluates 480 perturbations in three CRISPR-KO CD25 screen contexts. Tiny workbook/H5
fixtures test failure modes, guide eligibility, selection isolation, orientation,
determinism, and check-mode drift. The same runner streams S14, verifies all 520
donor/guide/context strata, and reproduces 116 S8 aggregates to <5.7×10⁻¹³. Unequal-size
fixtures prove that controls are weighted by guide rather than pooled cell count. This is
modest transfer plus descriptive supplied-score robustness, not donor-population or
whole-state validation.

`scripts/run_zhu_arrayed_validation.py` verifies the aggregate screen and two compact
author-table hashes, masks all nine panel target genes from every profile, and evaluates source-to-arrayed
bulk-RNA identity both raw and after panel centering. It normalizes IL-10/IL-21 flow within
donor against all available NTC measurements and enumerates all 362,880 target-label
permutations with synchronized maxT for four RNA-to-flow rank associations. Tiny fixtures
cover hash drift, outcome-independent screen selection, duplicate bulk keys, donor-control
normalization, all-panel-target masking, retrieval, and exhaustive enumeration determinism. The panel
was source-selected and donor coverage is unbalanced, so this certifies cross-platform
follow-up only.

## Canonical evidence status

| Question | Evidence | Status |
|---|---|---|
| Source-bound random-gene alignment? | 0.444 ± 0.018 across 12 hash-frozen splits | Descriptive; correlated splits |
| Target-source directional transfer? | Positive cosine gain over mean/best-single in 6/6 splits both ways | Directional only; nRMSE not improved |
| Published-eligibility donor-pair transfer? | Run-balanced median cosine gain +0.032; nRMSE 1.153 vs 1.018 | Weak direction; magnitude fails; four fixed donors |
| Independent functional-screen transfer? | Arce Spearman 0.148 / 0.084 / 0.088 | Modest, context dependent |
| Arce selected-panel concordance? | Within-dataset A-vs-B target-rank concordance 0.73–0.93; four-stratum sign agreement 50–64% | Descriptive; preselected 28 regulators, two donors, supplied score |
| Zhu arrayed transcriptome replication? | Matching target retrieves 9/9; median panel-centered cosine 0.580 | Source-selected same-study follow-up; nRMSE 1.052 |
| Zhu cytokine consistency? | Screen/bulk RNA versus donor-median flow Spearman 0.717–0.850 | Six follow-up donor labels with unequal target coverage; conditional panel diagnostic |
| Software/statistical contracts? | All six systemic scenarios pass | Synthetic certification only |
| Leakage-safe primary-model donor/guide generality? | Donor-pair released-object sensitivity only | Not tested |
| Functional state conversion? | No established-state prospective assay | Not tested |

## Artifact consistency

`scripts/validate_findings.py` cross-checks the ledger against all external reports and
the systemic harness, requires report `PASS`, requires full source hashes, checks split-ID
hashes, and validates SHA-256/bytes/executable bits for every canonical file. The central
figure reads only the ledger and is rendered twice to confirm deterministic PNG/PDF output.

## Reproduce

```bash
python -m pip install -r requirements.txt
./reproduce.sh
```

External regeneration additionally requires the gitignored registered inputs and
`requirements-external.txt`. Scientific readiness remains **Needs revision** until the
leakage-safe donor/guide, stronger-baseline, structured-null, measured-outcome, and prospective
gates in the [scientific validation plan](SCIENTIFIC_VALIDATION_PLAN.md) pass.
