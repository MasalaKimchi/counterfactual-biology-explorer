# Learned neural baseline vs the CombiCone certificate — honest verdict

**CombiCone Phase 2, learned-baseline track.** A self-contained learned neural
interaction model built as a *baseline for comparison*, plus an honest head-to-head
against the reachability cone. `combicone.py` and `reachability.py` are imported
**read-only** (verified byte-identical after the run) and never modified.

Data: Norman combinatorial CRISPRa screen (file label **A549 / canonically
K562**), 105 single-gene effect atoms, 131 measured doubles, 5045 genes. All 131
doubles have both constituent singles present as atoms.

---

## What was built

For each double `A+B` the inputs are its two constituent single-gene effect atoms
`a`, `b`; the target is the measured double effect `t = mean(cond) − ctrl`. Genes
are reduced by an **SVD fit on the 105 single atoms** to **50 components (95.6% of
atom-span variance)**. A **2-layer MLP** (`tanh` hidden layer, hidden = 32) maps
interaction features `[a, b, a*b]` (concatenation + elementwise product, so pairwise
interactions are representable) to the double effect. Training is **manual numpy
backprop** optimised by `scipy.optimize.minimize(L-BFGS-B)`; numpy/scipy only,
deterministic given a seed. The analytic gradient is certified against finite
differences in the test suite (max abs error < 1e-6).

The learner is given its **best fair shot**: the MLP is *additive-anchored*, i.e.
it predicts `(a+b) + correction` and only has to learn a nonlinear interaction
*correction* on top of additivity. This can only help it relative to additive, so
a failure to beat additive is a **conservative** conclusion.

Evaluation is **leave-pairs-out**: for each double, train on the other 130 and
predict the held-out one (blind). Seed-stable to the 4th decimal (0.9013–0.9016
over 3 seeds on a 12-fold check).

### Fairness of the SVD reduction
At 50 components the MLP's *ceiling* — cosine of the target's 50-component
projection to the full target — is mean **0.936 / median 0.950**, essentially
identical to the in-sample cone fit (0.937 / 0.950). The reduction does **not**
handicap the MLP relative to the cone; they share the same headroom, so any
accuracy gap reflects modelling, not the basis.

---

## Held-out prediction accuracy (cosine to measured double, full gene space)

| method | what it sees | held-out cosine (mean) | median |
|---|---|---|---|
| **MLP (learned, additive-anchored)** | blind — never sees held-out target | **0.896** | 0.917 |
| **additive `a+b`** | blind, training-free | **0.897** | 0.914 |
| **cone reachable fit** `rx.project_cone` | **in-sample** — *sees* the target | 0.937 | 0.950 |

- The learned MLP given its best fair shot **ties the trivial additive baseline**
  (0.896 vs 0.897) and beats it on a coin-flip **49.6%** of doubles
  (mean Δ = −0.0009).
- Learner variants that were tried and **lost to additive**: an unanchored raw MLP
  (0.838) and a learned linear ridge on the same 150 features (0.79–0.845; overfits
  150 features on 130 samples). The anchored MLP is the learner's genuine best shot.
- The cone reachable fit is higher (0.937) but it is an **in-sample reconstruction
  that requires having measured the target** — it fits 105 non-negative
  coefficients to the observed `t`. It is **not** a blind predictor and cannot rank
  unmeasured combinations. Comparing it to the blind MLP as if both were forecasts
  would be dishonest, so it is always labelled a reconstruction.

**This reproduces the field result** (Ahlmann-Eltze, Huber & Anders, *Nat. Methods*
2025, doi:10.1038/s41592-025-02772-6): benchmarking five foundation models and two
other deep learning models against deliberately simple baselines for single or
double perturbations, none outperformed the baselines — learned / deep models are
≈ simple linear/additive on combinatorial perturbation prediction. (Paper title:
"Deep-learning-based gene perturbation effect prediction does not yet outperform
simple linear baselines".)

---

## The defensible position: the cone's value is the CERTIFICATE, not accuracy

Neither the MLP nor the additive baseline emits any **certificate of
unreachability**. The reachability cone does, and it is verified here on the real
data:

- **On all 131 outside-cone doubles**, the model-relative dual separator `s`
  satisfies the certificate property exactly: `maxᵢ ⟨s, atomᵢ⟩ ≤ 0 < ⟨s, target⟩`
  (131/131). No single-gene atom crosses the separating hyperplane; the measured
  double sits on the far side.
- Worked example **SET+CEBPE** (panel B): every one of the 105 atoms has
  `cos(s, atomᵢ) ≤ 0`; the measured double has `cos(s, t) = +0.50`. Certificate:
  `z = 63.8`, `3.6× noise floor`, `p = 0.005` → **certified emergent (clears both
  bars)** under the split-half noise model.

### Can a learned model's "surprise" substitute for the certificate? No.
A natural bridge: use the learned model's held-out prediction *error*
(`1 − held_out_cosine`) as an emergence ranker — an emergent double should be
predicted worse. Measured against the project's `benchmark_rankers` labels:

| ranker | vs `cone_z` (noise-aware) | vs `cone_raw` | vs `nonadditivity` | vs `effect_norm` (confound) |
|---|---|---|---|---|
| **mlp_pred_error** | **+0.11** | +0.87 | +0.68 | **−0.58** |
| cone_raw (ref) | +0.07 | 1.00 | +0.74 | −0.56 |
| cone_z (ref) | 1.00 | +0.07 | +0.15 | +0.14 |

The learned model's prediction error is **essentially a re-expression of the raw
unreachable fraction** (Spearman +0.87 vs `cone_raw`) and **inherits the same
signal-to-noise confound** (−0.58 vs effect magnitude, matching cone_raw's −0.56).
It does **not** recover the noise-aware emergence signal (+0.11 vs `cone_z`). A
learned "surprise" falls into exactly the SNR trap that
`combicone.certify_emergence` guards against, and cannot replace the
noise-calibrated certificate.

---

## Honest caveats — what this does NOT show

1. **Accuracy is not the cone's claim.** This module measures prediction accuracy
   only. The cone's contribution is the certificate (a separator + a noise test),
   not out-predicting a learner. Accuracy fell as it does in the field: the learner
   ties additive.
2. **The cone fit is in-sample, not a forecast.** Its 0.937 is a reconstruction
   quality that requires the measurement; do not read it as a blind-prediction win.
3. **Modest model, small screen.** 131 doubles, numpy/scipy-only MLP. A larger deep
   model on more data might edge out additive, but the field result and the
   additive-anchored best-shot design make a large blind-prediction gain unlikely,
   and it would still emit no certificate.
4. **Model-relative geometry.** "Unreachable" means outside the non-negative cone
   of *these* measured single effects under *this* metric — never a claim of
   biological impossibility. Every reported result carries that scope.
5. **`cone_z` reproduction, not a new result.** The noise-aware column and the
   confound Spearmans reproduce numbers from the project's prior head-to-head; this
   track's new content is the *learned baseline* and its placement in that picture.

---

## Deliverables

- `neural_baseline.py` — importable, numpy/scipy-only, `__main__` demo, docstrings
  carrying the scope/honesty caveats.
- `test_neural_baseline.py` — 21 deterministic tests (incl. finite-difference
  gradient check), pass in `reach-pinned`.
- `fig_learned_vs_certified.png` — panel A: held-out prediction cosine (MLP vs
  additive vs in-sample cone fit); panel B: the certificate (separator) on
  SET+CEBPE.
- `per_method_heldout_accuracy.csv` — per-double held-out cosines for all three
  methods.
- `benchmark_rankers_with_learned.csv` — the project's `benchmark_rankers` table
  augmented with the learned baseline as an emergence ranker + prediction cosines.
- `learned_ranker_confound_summary.csv` — Spearman of each ranker vs the emergence
  labels and the effect-magnitude confound.
