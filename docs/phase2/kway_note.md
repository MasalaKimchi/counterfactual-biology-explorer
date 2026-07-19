# CombiCone k-way certification: extension and validation

**Scope of this note.** CombiCone certifies whether a measured combination effect
lands *outside* the non-negative cone of a supplied set of measured effect atoms,
and backs the geometric separator with a noise-injection null and a two-bar
verdict. Phase 1 shipped this for pairs. This phase generalizes the API to
arbitrary combination order (triples and beyond) and validates the generalization
two independent ways. **Every verdict here is model-relative** (unreachable = outside
the cone of the supplied atoms under the chosen metric; not a claim of biological
impossibility) and, for the triple validation, **synthetic** — the Norman file
(label `cell_type=A549`; canonically K562) contains 131 measured doubles and **zero
measured triples**, so no claim about real biological 3-way epistasis is made or
implied.

---

## 1. What changed in the code (backward-compatible superset)

The projection engine (`reachability.project_cone`) was *already* order-agnostic:
it takes a cone (any set of atoms) and one target vector. The k-way work is
entirely in the `combicone` triage/certify layer, and every change is additive —
**all 19 original pair tests pass unchanged**, and the k=2 code path is
byte-identical in result to the previous `-cos(effA, effB)` score.

| Function | Generalization | k=2 behavior |
|---|---|---|
| `combo_cosine(effects, agg=)` *(new)* | mean (default) or max of the C(k,2) pairwise cosines among the k singles | equals the single pairwise cosine (mean=max) |
| `additive_gap(singles, *idx)` | variadic: any number of member indices | `additive_gap(s, i, j)` unchanged |
| `triage_combinations(..., order=, pairwise=)` | auto-generates C(n, order) combos; `candidate_combos` (alias `candidate_pairs`); score = `-agg(pairwise cosines)` | identical to `-cos(a, b)` |
| `certify_emergence(cone_atoms=...)` | order-neutral `cone_atoms` alias for `singles`; docstring documents any-order cones | unchanged; `singles=` still works |
| `fit_triage_model(labeled_pairs=...)` | accepts k-tuple keys (one fixed order), `pairwise=` agg | unchanged |

**Triage score choice (documented).** For k>2 the single-single-cosine analog is the
**mean** of all C(k,2) pairwise cosines among the k constituent singles (a `max`
option is available and is more pessimistic — a single collinear pair drives the
score to −1). Mean is the default because it degrades gracefully with k and reduces
*exactly* to `-cos(a, b)` at k=2. Triage remains a cheap, training-free *ranking*
heuristic; it is not the certificate.

**Fail-closed and scope preserved.** All new paths raise `reachability.InputError`
on mixed-order combos, repeated members, unknown genes, degenerate cones, or bad
aggregation names; every `TriageResult`/`EmergenceCertificate` still carries the
model-relative `scope` string.

New tests: **19 k=3 tests added** (total 38 pass). They cover `combo_cosine`
semantics and k=2 identity, variadic `additive_gap`, order-3 triage (generation,
aggregation, aliases, rejection paths), order-3 `certify_emergence`, the
`cone_atoms` alias, `fit_triage_model` at order 3, and — the centerpiece — the
**reachable-from-lower-order flip** (below).

---

## 2. Validation A — reachable-from-lower-order on REAL data

Norman has no triples, but it *does* let us exercise the k-way idea between order-1
and order-2 cones. For each of the 131 measured doubles we certify emergence
against two cones:

- **cone (i)** = 105 single-gene atoms;
- **cone (ii)** = 105 singles **+ the other 130 measured doubles** (235 atoms).

A double that is certified emergent from the singles cone but becomes *reachable*
once the other doubles are atoms carries no structure beyond what those lower-order
measurements already express — exactly the order-dependence the k-way machinery is
built to detect.

**Result (n_boot=200, two-bar verdict p<0.05 AND floor_ratio≥1.9):**

| quantity | cone (i) singles | cone (ii) singles+doubles |
|---|---|---|
| certified emergent | **40 / 131** | **16 / 131** |
| geometrically outside cone (point) | 131 / 131 | 131 / 131 |
| bar(a) only (p<0.05, permissive) | 128 / 131 | 125 / 131 |
| median unreachable_fraction | 0.313 | 0.264 |
| median z | 14.8 | 9.0 |

- **24 doubles flip** from certified-emergent to reachable (**60% shrinkage** of the
  certified-emergent set); **0** are newly certified (monotone — a larger cone can
  only shrink the residual; verified: residual strictly non-increasing for all 131).
- The shrinkage is **effect-size driven, not geometric**: all 131 remain point-wise
  *outside* both cones, but among the 40 originally-certified doubles the median
  floor_ratio drops 2.15 → 1.81, pushing 24 below the 1.9 discriminating bar. The
  extra doubles absorb enough of each residual that the leftover no longer clears
  measurement noise by the stricter bar.

This is the honest reading of the two-bar design: bar(a) (does anything stick out
of the cone above noise) barely moves (128→125), because *most* doubles do have
some off-cone component; bar(b) (is that component *large* relative to the floor) is
where order-dependence bites. **Certified emergence is a property of the cone you
certify against, and it is not robust to enriching the cone with other measured
combinations** — which is precisely why higher-order screens need order-aware
certification rather than a single singles-only pass.

*(Flipped doubles include several CEBP-, ETS2-, MAPK1- and UBASH3B-centered pairs;
full per-double table in `kway_real_results.json`.)*

---

## 3. Validation B — synthetic triple screen with planted 3-way epistasis

Because no real triples exist, order-3 recovery is validated on a synthetic screen
(`synth_triple_screen.py`, deterministic, seed 0): K=14 singles, all 91 measured
doubles `D_ij = S_i + S_j + P_ij`, 300 genes, and 120 measured triples in three
planted classes (40 each):

- **additive** `T = S_i + S_j + S_l` — reachable from singles;
- **2-way-reducible** `T =` non-negative mix of the three constituent doubles (+ singles)
  — reachable once doubles are atoms, no genuine 3-way term;
- **3-way-emergent** `T =` reducible base **+ Tt**, where **Tt is orthogonal to the
  span of every atom** (all singles and doubles) — a component present *only* at the
  triple. Planted `‖Tt‖` spans 1.5–7.0 so the problem is non-trivial near the floor.

Only the measured *triples* carry split-half Gaussian noise (τ=0.10); the cone atoms
are the clean true effects. This isolates the certificate's ability to detect a
planted 3-way term above triple-measurement noise from a clean lower-order cone —
the standard planted-truth setup. Each triple is certified at **order 3** against the
cone {singles + its three constituent measured doubles}.

**Result:**

| metric | value |
|---|---|
| AUC, noise-aware z (emergent vs rest) | **1.000** |
| AUC, raw unreachable_fraction | 0.925 |
| two-bar: TP / FN / FP / TN | 36 / 4 / 0 / 80 |
| sensitivity / specificity / precision | 0.900 / 1.000 / 1.000 |
| certified: additive / reducible / emergent | 0/40 / 0/40 / 36/40 |

- The order-3 two-bar verdict certifies **36/40 planted-emergent triples and 0 of the
  80 additive/2-way-reducible ones** (zero false positives).
- The 4 misses are the four smallest planted terms (‖Tt‖ = 1.50–1.92); their
  floor_ratio sits at or just under 1.9 — the certificate correctly *declines* to
  call emergence at the noise floor rather than over-calling.
- **The magnitude confound is visible and is the point.** Raw unreachable_fraction
  for emergent triples (0.096–0.460) *overlaps* the additive range (0.127–0.158), so
  the raw residual only reaches AUC 0.925; the **noise-aware z separates the classes
  perfectly (AUC 1.000)**. This reproduces, at order 3, the same lesson the pair
  analysis taught: report the noise-aware statistic, never the raw residual.

---

## 4. Honest limitations

1. **The triple validation is synthetic.** It validates the order-3 *code path* and
   its discrimination of planted 3-way structure; it is **not** evidence about real
   3-way epistasis, because Norman measures no triples. Real triples would add cone
   jitter (noisy atoms) and biological structure not captured by an orthogonal
   planted term.
2. **No higher-order predictive performance is claimed.** The LOO triage numbers
   (~0.43 Spearman, ~2.4× top-20 precision) are the *pair* results on real data; the
   k-way triage score is a ranking heuristic with no measured-triple benchmark.
3. **Certified emergence is cone-relative and not robust to cone enrichment**
   (Validation A): the same double can be certified against singles and reachable
   against singles+doubles. The certificate answers "outside *this* cone above
   noise", not "irreducibly high-order".
4. **The two-bar verdict is a threshold on a magnitude-confounded quantity made
   noise-aware.** bar(a) (p<0.05) is permissive; bar(b) (floor_ratio≥1.9) is the
   discriminating gate and is where both the real shrinkage and the synthetic misses
   live. Always report both bars, never bar(a) alone.
5. **The defensible unique capability remains the certificate** (a separator + a
   noise test with a stated scope), not prediction accuracy.

---

## 5. Files

- `combicone.py` — k-way superset (backward-compatible; 38 tests pass).
- `test_combicone.py` — 19 original + 19 new k=3 tests.
- `kway_validate_real.py` — Validation A driver (real doubles, two cones).
- `synth_triple_screen.py` — synthetic planted-epistasis generator.
- `kway_validate_synth.py` — Validation B certifier + metrics.
- `fig_kway_validation.png` — panel a (real shrinkage), panel b (synthetic recovery).
- `kway_real_results.json`, `kway_real_summary.json`, `synth_triple_results.json`,
  `synth_triple_screen_config.json` — full numeric results and config.
