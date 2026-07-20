# Validation Report

What is **certified** (math, holds by construction), what is **source
reconstructed** (hash-bound to real data), and what is **independently
evaluated** (descriptive, with a stated ceiling) — for CombiCone. This document
absorbs the trust-dossier argument formerly in
`docs/phase3/certificate_dossier.md`. Numbers are single-sourced from
[`results/findings.json`](../results/findings.json) and the frozen reports under
`results/`; methods are in [`METHODS.md`](METHODS.md), findings in
[`FINDINGS.md`](FINDINGS.md), the forward program in
[`SCIENTIFIC_VALIDATION_PLAN.md`](SCIENTIFIC_VALIDATION_PLAN.md).

---

## 1. What is technically certified (holds by construction)

- **The infeasibility separator.** For every combination outside the cone,
  `reachability.project_cone` returns a Farkas/KKT dual separator at machine
  precision (`1e-8`) or **fails closed**. This is math, not a fit; no forward
  predictor supplies it. On Norman, all **131/131** doubles are certified
  geometrically outside the single-gene cone with a valid separator.
- **The two-bar noise verdict is well-defined and both bars are load-bearing.**
  bar (a) = noise-injection `p<0.05`; bar (b) = `floor_ratio ≥ 1.9×`. The bars
  answer different questions and are reported separately, never merged into one
  significance label (see §2 negative controls for why both are needed).
- **k-way code path.** 38 unit/geometry tests pass (19 pair + 19 k=3); the k=2
  path is byte-identical to the original `−cos`. Fail-closed on mixed-order
  combos, repeated members, unknown genes, degenerate cones, bad aggregation.
- **Learned-baseline gradient.** The MLP's analytic gradient is certified against
  finite differences (max abs error < 1e-6); `combicone.py`/`reachability.py`
  imported read-only, verified byte-identical after the run.
- **Adversarial harness.** 6 data-free scenarios PASS (common-response inflation,
  random-gene optimism, sign-selection inflation); maxT exact one-sided 95% upper
  bound **0.0668** over 24 false families / 500 trials.
  Source: [`results/validation_harness.json`](../results/validation_harness.json).

---

## 2. Certificate trust argument (specificity, sensitivity, construct validity)

A causal-inference-style stress test of the one causal-flavored claim — *a
measured combination lands outside the non-negative cone of its constituent
singles by more than measurement noise can explain* — on the real Norman
substrate (105 atoms, 131 doubles). Reproduce:
`python scripts/certificate_dossier.py`.

### Specificity — negative controls

- **NC-additive (headline specificity).** 150 synthetic `D = a + b` combinations
  (inside the cone by construction) + measurement-scale noise: **0/150
  false-certified** (141 "within noise", 9 "above noise but modest", 0
  emergent). Mean floor ratio 1.01, 95th pct 1.16 — nowhere near the 1.9× bar.
  **False-positive rate 0.000.**
- **NC-sham.** A pure-noise combination scrapes bar (a) but its floor ratio ≈ 1.0
  rejects bar (b) — the concrete demonstration of *why both bars exist*:
  significance alone can be fooled by a vector directionally far from the cone but
  tiny; the noise-floor ratio cannot.
- **NC-reachable (deliberately imperfect).** The 13 least-unreachable real
  doubles re-certified: 3 certified, 10 "above noise but modest", 0 within-noise.
  This does *not* pass cleanly, and that is informative: a small residual (nearly
  additive in size) can still clear both bars when its measurement noise is tiny.
  "Certified emergent" is a noise-relative confidence statement — read the
  residual alongside the verdict.

### Sensitivity radius Γ* — how fragile is the noise assumption

The certificate's one load-bearing assumption is that `SE = |t1 − t2|/2`
faithfully captures measurement noise. Γ* is the smallest inflation factor that
drops a pair below "certified" — the analogue of a Rosenbaum bound / E-value.
**Across 35 certified pairs: median Γ* = 1.25×, range 1.25×–2.0×**; SET+CEBPE
tolerates ~1.9×; only 3% survive a full 2× noise misestimation. Honest reading:
Γ* is a monotone re-expression of the bar-b margin (Γ* ≈ floor_ratio / 1.9), not
an independent axis — the takeaway is calibration, not comfort. "Certified
emergent" means *emergent under the stated noise model, at this margin*; the
strongest pairs carry the widest margins.

### Construct validity — does it measure non-additivity

1. **Recovers classical synergy.** Cone residual vs raw non-additivity Spearman
   **0.736**; partial Spearman **0.62** after controlling for magnitude — the
   synergy signal is not a magnitude artefact.
2. **Removes the magnitude confound.** Raw residual vs magnitude Spearman −0.56;
   noise-robust z vs magnitude +0.14 — removed by design.
3. **Top certified pairs are known interactions.** SET+CEBPE, IRF1+SET,
   MAPK1+PRTG, CEBPE+RUNX1T1 — the myeloid/erythroid master-regulator SET and the
   CEBPE differentiation axis, consistent with Norman et al.'s reported strong
   genetic interactions.

**Key subtlety — size vs. confidence.** The z-score is almost uncorrelated with
the cone residual (ρ = 0.06) and is instead dominated by the per-pair noise floor
(ρ = −0.59). Effect size lives in the residual; z answers the separate question
"how confident are we this is above noise." z must **not** be read as "how large
the synergy is." Reporting both is the honest practice — bar (a) is the
confidence test, bar (b) is an effect-size floor.

---

## 3. What is source reconstructed (hash-bound to real data)

- **Target lineage** from the Zhu Th2→Th1 polarization table: union 25,672
  genes, 11,616 shared source, 7,960 sign-concordant, 6,188 registered;
  between-source log-FC cosine 0.791, z-score cosine 0.698 on shared screen
  genes. Frozen random-gene splits; full-file SHA-256 bound.
  Source: [`results/source_reconstruction.json`](../results/source_reconstruction.json).
- **Substrates** committed at repo root and hash-stable:
  `combicone_substrate.npz` (Norman), `carpool_substrate.npz` (CaRPool).
- **Library-coverage caches** are rebuilt byte-for-byte by a hash-gated
  deterministic builder; the frozen report refuses a cache whose registered
  identity or schema differs. Durable-availability boundary: the Norman/Replogle
  public URLs are mutable and historical S3 versions are not anonymously
  retrievable, so a release archive/DOI is still the manuscript-grade step.

---

## 4. What is independently evaluated (descriptive, with a ceiling)

Each stress benchmark below is a **descriptive** evaluation with an explicit
claim ceiling; none is a calibrated biological or predictive claim. Full ceilings
and provenance: [`SCIENTIFIC_VALIDATION_PLAN.md`](SCIENTIFIC_VALIDATION_PLAN.md)
and [`../data/README.md`](../data/README.md).

| Benchmark | Status | What it shows / ceiling | Source report |
|---|---|---|---|
| Cross-screen (CaRPool) | TRANSFER | Certificate + SNR confound reproduce on an orthogonal screen; triage's magnitude-controlled part does not | `docs/phase2/cross_screen_emergence.csv` |
| Arce Perturb-CITE | STRESS | Modest S1 cross-modality ranking (ρ 0.08–0.15) + S14 supplied-score donor robustness (0.73–0.93) in a preselected panel | `results/evidence/arce_external_validation_meta.json` |
| Schmidt two-donor screens | STRESS | Source-selected top-K concordance; whole-universe agreement limited (0.135–0.332); joint axes not isolated | `results/schmidt_external_validation.json` |
| Zhu arrayed follow-up | SUPPORTING | Same-study cross-platform RNA replication (top-1 1.0) + IL-10/IL-21 flow consistency (0.72–0.85) | `results/evidence/zhu_arrayed_validation_meta.json` |
| Donor-pair transfer | STRESS | Weak directional gain (+0.032, 75% positive), fails magnitude (RMSE worse) | `results/donor_pair_transfer.json` |
| Guide-position transfer | STRESS (negative) | Reciprocal guide-rank transfer negative (0.251 → −0.019); not physical-guide generalization | `results/guide_pair_transfer.json` |
| Goudy triple | STRESS (negative) | Execution PASS, geometric model FAILS, biology INCONCLUSIVE (perfect single-vs-triple confounding) | `results/goudy_combination_validation.json` |
| Library coverage | EXTENSION | Strict membership 0 everywhere; certificate ρ 0.861–0.921 with realized gains but 0/3 top-atom match; retrospective | `results/library_coverage_crossdataset.json` |

---

## 5. Canonical evidence status & artifact consistency

- **`results/findings.json` is the single source of truth** for every headline
  number. Human-readable views ([`FINDINGS.md`](FINDINGS.md), the manuscript, the
  README) must not diverge from it; when a result changes, update the JSON, its
  supporting evidence, the validation rules, and the figure together, then
  regenerate `results/manifest.json`.
- **`results/manifest.json`** records SHA-256, byte length, and executable bit for
  every maintained artifact; **`scripts/validate_findings.py`** fails closed on
  report status, canonical values, paths, hashes, and executable bits.
- **Retired provenance.** The archived 60-shuffle / multisplit pipeline depended
  on an unhashed, deleted `inputs.npz` whose gene order was not preserved; its
  values (`retired_provenance` in findings.json) are provenance, never current
  inference.

---

## 6. Known open item (pre-existing, out of scope for the doc consolidation)

`scripts/validate_findings.py` currently fails at the **guide_pair** check
because `configs/guide_pair_transfer.json` was edited after
`results/guide_pair_transfer.json` was frozen: the config's live SHA-256
(`9f5f1416…`) no longer matches the hash frozen into the report
(`01d69a30…`), and the `claim_ceiling` wording drifted from "positional /
verified physical-guide replication" (frozen report) to "alphanumeric guide-rank
/ named-sgRNA replication" (live config). Honest resolution requires re-running
`scripts/run_guide_pair_transfer.py` against the local 29 GB
`GWCD4i.DE_stats.by_guide.h5mu` (regenerating the frozen report/hash), or a
maintainer decision on the canonical `claim_ceiling` wording. It is **not** a
documentation issue and is tracked in
[`CONSOLIDATION.md`](CONSOLIDATION.md).

---

## 7. Reproduce

```bash
python -m pip install -r requirements.txt
./reproduce.sh          # env gate, test suite, demos, harness, findings + lineage validation
```

External scientific data are gitignored; source routes, hashes, and licenses are
in [`../data/README.md`](../data/README.md).
