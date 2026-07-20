# Methods

CombiCone is a triage layer for combinatorial perturbation screens, built on a
fail-closed convex-cone geometric core. This document is the canonical method
reference; it absorbs the former `docs/phase2/{kway_note,learned_baseline_note,
screenloop_note}.md` and `docs/phase3/{acquisition_note,certificate_dossier}.md`
prose notes. Every number reported here is single-sourced from
[`results/findings.json`](../results/findings.json) (canonical) and the frozen
metric JSON under `docs/phase2/` and `docs/phase3/`. Findings and their claim
ceilings live in [`FINDINGS.md`](FINDINGS.md); what is certified vs. merely
evaluated is in [`VALIDATION_REPORT.md`](VALIDATION_REPORT.md).

**Scope, stated once.** "Unreachable" is *model-relative*: outside the
non-negative cone of **these** measured single-gene effects under **this**
metric — never a claim of biological impossibility. Every result object carries
this disclaimer in a `scope` field.

---

## 1. Geometry and the certificate

Each measured single-gene effect profile is an atom (a row of
perturbations × genes). Their non-negative combinations form a convex cone

```
C = { E^T w : w >= 0 }
```

of everything the singles can additively reach. `reachability.project_cone`
solves the weighted non-negative least-squares projection of a target effect
`d` onto `C`:

- **inside the cone** → the target is representable as a non-negative mix of the
  atoms (numerical membership certified at `1e-8`);
- **outside the cone** → the solver emits a **model-relative separator** (a
  Farkas/KKT dual certificate) proving non-representability, or **fails closed**
  if numerical conditions do not pass.

This is the one object a forward virtual-cell predictor (GEARS, scGPT, STATE,
CPA) structurally cannot supply: a predictor always emits *some* prediction; the
cone emits a certificate of infeasibility when the target is not representable.
`ProjectionResult` carries `coefficients`, `fitted`, `residual`,
`residual_fraction`, `cosine`, `dual_separator`, and the KKT/orthogonality/
polarity diagnostics that back the fail-closed guarantee.

---

## 2. Noise-aware certification (the two-bar verdict)

A raw cone residual is **magnitude-confounded**: low-magnitude combinations show
inflated *normalized* residuals dominated by measurement noise, not biology
(Spearman −0.56 on Norman; §Findings). `combicone.certify_emergence` removes
this confound with a measurement-noise null:

1. **Noise model.** A deterministic split-half of cells gives two independent
   pseudobulk estimates `t1`, `t2` per condition; the per-gene standard error is
   `SE = |t1 − t2| / 2`.
2. **Reachable null.** `f0 = project_cone(singles, target).fitted + N(0, SE)`,
   `B = 200` draws, conservative plus-one empirical p-value.
3. **Floor ratio.** `floor_ratio = residual / noise-floor`.

The verdict has **two distinct bars** (always report which is meant):

- **bar (a)** — noise-injection `p < 0.05` ("emergent at all above noise",
  permissive);
- **bar (b)** — `floor_ratio ≥ 1.9×` ("emergent by a large margin", the
  discriminating gate).

"Certified emergent" requires **both**. The bars are different tests: bar (a)
can be scraped by a vector directionally far from the cone but tiny; bar (b)
cannot (see the NC-sham control in [`VALIDATION_REPORT.md`](VALIDATION_REPORT.md)).

**Deterministic fast path.** `certify_emergence(..., method="analytic")` replaces
the `B` Monte-Carlo re-projections with a closed-form generalized chi-square null
(`reachability.analytic_anisotropy_null`): seed-free, `~B×` fewer solves, and
**conservative by construction** — it can only withhold a certificate, never
inflate one. Across all 289 Norman+CaRPool doubles the analytic certified set is
a strict subset of the Monte-Carlo one (0 analytic-only certifications), and the
analytic-vs-Monte-Carlo gap is a self-reporting certificate-fragility score.
Derivation, guarantees, and cross-screen numbers:
[`analytic_null_note.md`](analytic_null_note.md).

> **Two-bar count reconciliation (documented, not an error).** The Norman
> counts read **129/131 (bar a), 35/131 (bar b)** under the canonical
> `findings.json` pipeline and **128/131 (bar a), 40/131 (bar b)** under the
> shared symmetric cross-screen pipeline. Both are correct for their
> configuration: bar (a) is near-saturated either way (128 vs 129 — nearly every
> double departs the cone *somewhat*), and the bar-(b) tally is the one that
> genuinely wobbles. `certify_emergence` uses a
> Monte-Carlo noise-injection null, and ~38 doubles sit in a dense band of
> `floor_ratio ∈ [1.7, 2.1]` straddling the 1.9 cut, so the *count* clearing
> the bar is intrinsically jittery near the threshold (analytic expected count
> 39.7 ± 0.6 at `n_boot=200` for the symmetric config). The *ordering*
> (CaRPool 48% > Norman ~27–31%) and the z-distribution are stable; the hard
> threshold tally is not. This is exactly the threshold-count fragility the
> noise-aware z was introduced to avoid — the headline compares distributions,
> not bar tallies.

---

## 3. Training-free triage

`combicone.triage_combinations` ranks *unmeasured* combinations from the single
effects alone. The default score is the **negative cosine** between two single
effect vectors, `−cos(effA, effB)` — training-free, no labeled pilot required.
The intuition: two singles pointing in different directions are more likely to
combine into something neither reaches alone. Triage is a cheap ranking
heuristic, **not** the certificate; its defensible value is a modest prospective
enrichment (§Findings), not ranking supremacy.

When a labeled pilot screen is available, `combicone.fit_triage_model` fits a
ridge onto the single-effect features for the magnitude-controlled target,
reported as a ceiling (LOO Spearman 0.435, permutation p=0.002), not the
headline.

---

## 4. k-way generalization (order ≥ 3)

The projection engine is already order-agnostic (it takes any cone and one
target). The k-way work is entirely in the `combicone` triage/certify layer and
is a **backward-compatible superset** — all 19 original pair tests pass
unchanged, and the k=2 path is byte-identical to `−cos(effA, effB)`.

| Function | Generalization | k=2 behavior |
|---|---|---|
| `combo_cosine(effects, agg=)` *(new)* | mean (default) or max of the C(k,2) pairwise cosines | equals the single pairwise cosine |
| `additive_gap(singles, *idx)` | variadic: any number of members | `additive_gap(s, i, j)` unchanged |
| `triage_combinations(..., order=, pairwise=)` | auto-generates C(n, order); score `−agg(pairwise cosines)` | identical to `−cos(a, b)` |
| `certify_emergence(cone_atoms=...)` | order-neutral `cone_atoms` alias | unchanged; `singles=` still works |
| `fit_triage_model(labeled_pairs=...)` | accepts k-tuple keys | unchanged |

For k>2 the triage score defaults to the **mean** of all C(k,2) pairwise cosines
(a `max` option is more pessimistic). All new paths raise
`reachability.InputError` on mixed-order combos, repeated members, unknown
genes, degenerate cones, or bad aggregation names. **38 tests pass** (19 pair +
19 k=3). Validation is two independent ways (see [`FINDINGS.md`](FINDINGS.md)
§k-way): a real reachable-from-lower-order flip on Norman, and a synthetic
planted-3-way-epistasis screen — the latter is a **code-path validation, not
evidence about biological 3-way epistasis** (Norman measures zero triples).

Drivers: `scripts/kway_validate_real.py`, `scripts/synth_triple_screen.py`,
`scripts/kway_validate_synth.py`.

---

## 5. Learned neural baseline (head-to-head)

`neural_baseline.py` is a self-contained learned interaction model built *as a
baseline for comparison*; `combicone.py` and `reachability.py` are imported
read-only (verified byte-identical after the run). For each double `A+B` the
inputs are the two constituent single atoms `a`, `b`; the target is the measured
double `t`. Genes are reduced by an **SVD fit on the 105 single atoms** to 50
components (95.6% of atom-span variance). A 2-layer additive-anchored MLP
(`tanh`, hidden=32) maps interaction features `[a, b, a*b]` to `(a+b) +
correction`, trained by numpy backprop + L-BFGS-B (numpy/scipy only,
deterministic). The analytic gradient is certified against finite differences
(max abs error < 1e-6).

The learner is given its **best fair shot** (additive-anchored, so it only has
to learn a correction), so a failure to beat additive is a conservative
conclusion. Evaluation is **leave-pairs-out** (blind). The SVD reduction does
not handicap the MLP: its ceiling (0.936/0.950) matches the in-sample cone fit
(0.937/0.950), so any accuracy gap reflects modelling, not the basis.

---

## 6. Certificate-guided sequential design (screenloop)

`screenloop.py` puts the certificate inside a *screening campaign* and separates
two questions on two published screens (driver:
`scripts/run_screenloop_campaign.py`; tests: `tests/test_screenloop.py`):

- **Acquisition** — `replay_campaign()` replays a screen as batches (batch 8)
  under a pluggable policy (`random`, `magnitude`, `triage`, `cone_adaptive`)
  and records wells to discover 90% of the two-bar noise-robust emergent
  combinations.
- **Library augmentation** — `nominate_atoms()` aggregates the model-relative
  separators of the combinations a library fails to reach into one
  residual-weighted "unmet demand" direction, then ranks candidate perturbations
  by alignment. The falsifiable test `held_out_single_recovery` removes a single
  that participates in ≥2 combinations and checks whether the separator recovers
  it at the top from the full candidate pool, against a permutation null, a
  magnitude-only ranker, and a dominance-vs-advantage control.

---

## 7. Prospective acquisition recommender

`acquisition.py` turns the per-combination triage score into a **batch
recommendation** for the next round, supporting the closed loop
(recommend → run → observe → refit → recommend). Relevance is the training-free
`−agg_cos` (or the ridge `TriageModel` once a pilot exists); diversity is greedy
max-marginal-relevance (MMR): each pick maximizes
`relevance − diversity_weight · max_similarity_to_batch`, similarity being the
cosine of predicted additive-effect directions (default) or gene-set Jaccard, so
a batch spreads across distinct predicted biology. `diversity_weight` is the
exploit↔explore knob. CLI:
`combicone recommend screen.h5ad --batch-size 10 --strategy diversified`. The
recommender chooses what to **measure**; whether a discovery is genuinely
emergent is decided by `certify_emergence`, not the acquisition score.

---

## 8. Real-data ingestion and CLI

- `effect_dictionary.py` — safe, labeled adapter from cell matrices to portable
  effect dictionaries. Matrix is always `(perturbations, genes)` with unique
  string labels; reads NPZ with `allow_pickle=False`; computes pooled condition
  mean minus pooled control mean. It does **not** normalize raw counts, build
  replicate-aware pseudobulks, correct batch, select features, or estimate
  uncertainty — those upstream choices are the caller's.
- `screen_ingest.py` — ingestion adapter from real screen files to the substrate
  the triage/certify API consumes.
- `combicone_cli.py` — the `combicone` command line (triage, certify, recommend).

---

## 9. Proving-ground / stress methods (retained, demoted)

The single-target reachability work is retained beneath the combinatorial
headline as a proving ground on real genome-scale data. Full provenance
(accessions, SHA-256, licenses, exact object keys, claim ceilings) is in
[`../data/README.md`](../data/README.md); per-benchmark configs are in
`configs/`. Summary of methods:

- **Source-bound reconstruction** (`scripts/run_source_reconstruction.py`) —
  full-file-hash-bound target lineage from the Zhu Th2→Th1 polarization table,
  frozen random-gene splits, and aggregate cross-source transfer (Ota,
  Höllbacker).
- **Independent Arce benchmark** (`scripts/run_arce_external_validation.py`) —
  Perturb-CITE CRISPRi-transcript → CRISPR-KO CD25 transfer and S8/S14 supplied-
  score donor/guide robustness on the authors' 28-regulator panel.
- **Schmidt two-donor functional screens**
  (`scripts/run_schmidt_external_validation.py`) — hash-gated CRISPRa/CRISPRi
  cytokine-sort concordance, two fixed donors, exploratory guide × top-K grid.
- **Zhu arrayed follow-up** (`scripts/run_zhu_arrayed_validation.py`) —
  source-selected bulk-RNA + IL-10/IL-21 flow replication in six follow-up
  donor labels.
- **Donor-pair transfer** (`scripts/run_donor_pair_transfer.py`) — frozen-weight
  complementary donor-pair + target-source transfer sensitivity on released
  H5MU modalities.
- **Guide-position transfer** (`scripts/run_guide_pair_transfer.py`) — negative
  reciprocal transfer across released `guide_1`/`guide_2` alphanumeric guide-rank
  modalities; physical sgRNA IDs are **not** embedded in the H5MU and the exact
  rank-to-ID crosswalk is not hash-verified.
- **Goudy triple cross-experiment stress**
  (`scripts/run_goudy_combination_validation.py`) — a single measured
  FAS+RC3H1+SUV39H1 CRISPRoff triple vs. target-matched singles; a bounded
  negative result (fails the declared geometric model, inconclusive biology).
- **Library coverage** (`library_coverage.py`,
  `scripts/run_library_coverage_crossdataset.py`) — split-first catalog coverage,
  redundancy, gap normals, and retrospective scoring of already-measured
  candidate effects across Zhu/Norman/Replogle.

---

## 10. Reproducibility

- **Pinned environment** `reach-pinned`: numpy 1.26.4 / scipy 1.13.1.
- **`reproduce.sh`** — gates the environment, runs the numerical test suite
  (including `tests/test_combicone.py`), executes the combicone/reachability/
  coverage demos, checks the adversarial harness, and validates the frozen
  findings and artifact lineage via `scripts/validate_findings.py`.
- **External scientific data are gitignored**; source routes, hashes, and
  licenses are in [`../data/README.md`](../data/README.md). Claim-bearing tables
  are regenerated by the external-data runners, never copied from a notebook.
- **Substrates committed at repo root:** `combicone_substrate.npz` (Norman),
  `carpool_substrate.npz` (CaRPool).
