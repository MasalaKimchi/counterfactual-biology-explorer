# Closed-form anisotropy-corrected noise null

*Methods note for the deterministic emergence certificate. Companion to
[`METHODS.md`](METHODS.md) §2. Every number here is single-sourced from
[`results/analytic_vs_mc_validation.json`](../results/analytic_vs_mc_validation.json)
and the per-double table
[`results/analytic_vs_mc_by_double.csv`](../results/analytic_vs_mc_by_double.csv).*

**Scope, stated once.** "Null" here is the measurement-noise null for a
*model-relative* residual: the residual of a measured combination against the
non-negative cone of **these** measured single-gene effects under **this**
metric. Nothing below is a claim of biological reachability.

---

## 1. The quantity being nulled

`combicone.certify_emergence` asks whether a measured combination effect `d`
departs from the cone of single-gene effects by more than measurement noise. It
projects `d` onto the cone, `f0 = project_cone(singles, d).fitted`, and reports
the **residual fraction** `residual/‖d‖`. The honest null is *"the target is in
fact reachable and the residual we see is noise"*: draw fresh per-gene noise
around `f0`, re-project, and see how large the residual fraction is under noise
alone. The production path (`method="montecarlo"`) does this by `B=200`
re-projections. This note derives the same null in closed form
(`method="analytic"`).

## 2. Derivation (generalized chi-square)

Work in the metric-whitened coordinates the projection uses, `x → W^{1/2} x`,
with `W` the (mask-restricted, canonically scaled) gene weights. Let `S` be the
atoms **active** at `f0` (positive NNLS coefficients) and `P` the orthogonal
projector onto their whitened column span. For a noise draw
`ε ~ N(0, D)`, `D = diag((W^{1/2}·SE)²)`, re-projecting `f0 + ε` keeps the
active facet to first order, so the null residual vector is `(I − P) ε` and

```
Q = ‖(I − P) ε‖²  ,   residual_fraction = √Q / ‖W^{1/2} f0‖ .
```

`Q` is a **generalized chi-square** — a quadratic form in a Gaussian. Its first
two cumulants are exact and need only `|S|×|S|` traces:

```
E[Q]   = tr((I − P) D) = tr(D) − tr(P D)
Var[Q] = 2 · tr( ((I − P) D)² ) .
```

We match a gamma to `(E[Q], Var[Q])` (Satterthwaite/Welch), map the moments of
`√Q` to the residual fraction, and report
`p = P(Q ≥ (observed · ‖W^{1/2} f0‖)²)`. The whole computation is a handful of
small matrix traces — deterministic, seed-free, and roughly `B×` fewer NNLS
solves than the Monte-Carlo path.

The confidence interval on the *observed* residual under its own noise uses the
same active-facet linearization: `‖r(ε)‖² = ‖r0‖² + 2 r0·ε + Q`, whose mean and
variance are exact for zero-mean Gaussian `ε` (the cross term is uncorrelated
with `Q` by Wick/Isserlis). We propagate to the residual fraction by the delta
method and report a symmetric Wald interval — a normal approximation to the
sampling distribution, not an empirical percentile interval.

## 3. Why anisotropy is in the name

Per-gene noise **anisotropy** (heterogeneous `SE`) enters through
`tr(P D) = tr((AᵀA)⁻¹ Aᵀ D A)` — the noise-weighted leverage of the active
facet. A scalar-`SE` (isotropic) surrogate cannot see this term and mis-estimates
the null whenever atom leverage aligns with the noisy genes. The closed form
handles it exactly; that exactness is the point of the "anisotropy-corrected"
name.

**Honest boundary.** On a synthetic case where atom support is deliberately
aligned to the high-noise genes, an isotropic/RMS-`SE` approximation is ~11%
wrong while the anisotropic form stays exact. **On the real screens, though, the
isotropic approximation is small**: across all 289 doubles the null-mean
difference between the anisotropic form and an isotropic RMS-`SE` surrogate has
median 1.7% / max 4.2% on Norman and median 0.03% / max 0.6% on CaRPool. So
anisotropy is a *correctness guarantee* the closed form provides — not a large
empirical mover on Norman/CaRPool.

## 4. Guarantee: conservative direction

`P` projects onto the active-atom **subspace**. That subspace is locally
contained in the cone, so the true cone distance never exceeds the subspace
distance and the analytic null residual is `≥` the Monte-Carlo one. Where the
facet is a genuine vertex/edge of the cone (few active atoms — the biologically
relevant regime) this holds strictly and substantially, so the analytic path can
only **withhold** a certificate, never inflate one. In general position (random
atoms, a full-dimensional facet) the two nulls coincide and the inequality holds
only up to Monte-Carlo estimation error (deviations of a fraction of a percent).
It is therefore conservative *in direction and in the certified regime*, not a
strict pointwise theorem in every geometry.

**Verified.** On all 289 real doubles the analytic null exceeds the Monte-Carlo
null: **0/131** Norman and **0/158** CaRPool violations, ratio 1.01–1.44
(Norman) / 1.01–1.62 (CaRPool).

## 5. Cross-screen agreement (289 doubles)

`method="analytic"` vs `method="montecarlo"` (`n_boot=200`, `seed=0`):

| | Norman CRISPRa | CaRPool-seq |
|---|---:|---:|
| doubles | 131 | 158 |
| conservative-direction violations | 0 | 0 |
| null-mean agreement within 5% | 84% | 53% |
| verdict-tier agreement | 92% | 88% |
| certified (Monte-Carlo) | 40 | 76 |
| certified (analytic) | 34 | 68 |
| analytic-only certifications | 0 | 0 |
| significance (p<0.05) agreement | 97% | 93% |

The decision-level guarantee is the last three rows: **the analytic certified
set is a strict subset of the Monte-Carlo one on both screens**, with zero
analytic-only certifications. Every verdict disagreement is the analytic path
being *stricter*.

## 6. The analytic-vs-Monte-Carlo gap is a fragility diagnostic

The gap between the two nulls is not noise — it is a **certificate-fragility
score** computable from pure geometry. It is large exactly for the low-SNR
certificates the effect-size bar already demotes, and negligible for robust ones:

- **Driven by SNR, not active-set churn.** Gap vs Monte-Carlo floor ratio:
  Spearman −0.53 (Norman) / −0.78 (CaRPool). Gap vs active-set Jaccard churn:
  Spearman ~0 (not significant) — churn does *not* explain the gap. Binned across
  both screens, fragile certificates (<1.5× floor) have a median gap of 7.2% (max
  62%) while strong ones (>3×) have a median of 2.0% (max 2.4%).
- **Self-reporting.** The gap is predictable from *analytic-only* quantities
  (gap vs analytic floor ratio Spearman −0.57 / −0.82), so the deterministic path
  flags its own borderline cases without ever running Monte-Carlo.
- **Decision-safe.** All 14 certificates the analytic path withheld relative to
  Monte-Carlo (6 Norman + 8 CaRPool) sat at a Monte-Carlo floor ratio in
  **[1.90, 1.97]** — right at the 1.9× bar. The fast path diverges only in the
  borderline zone, never on a robustly-certified double.

This subsumes the ad-hoc "bar-count jitter near 1.9×" observation in
[`METHODS.md`](METHODS.md) §2 with a principled, certificate-level fragility
measure.

![Analytic vs Monte-Carlo emergence null across both screens: (a) the analytic
null reproduces the Monte-Carlo null and never falls below it; (b) the gap is a
fragility score, large only for low-SNR certificates; (c) verdicts flip only at
the 1.9× bar.](figures/fig_analytic_null_fragility.png)

## 7. Usage

```python
import combicone as cc
# deterministic closed-form path (seed-free):
cert = cc.certify_emergence(cone_atoms=atoms, measured_combo=effect,
                            noise_sd=se, method="analytic")
# stochastic reference (default, unchanged):
cert = cc.certify_emergence(cone_atoms=atoms, measured_combo=effect,
                            noise_sd=se, method="montecarlo", n_boot=200, seed=0)
```

Both return an `EmergenceCertificate` with identical two-bar verdict semantics.
The low-level null is `reachability.analytic_anisotropy_null`, which returns the
null mean/sd, an analytic p-value, the gamma parameters, an anisotropy-aware
effective degrees of freedom, and the analytic CI.

## 8. Open items

- **Second-order term for the low-SNR regime.** The over-estimate on fragile
  certificates comes from the fixed-facet linearization; a tangent-cone
  correction (re-projecting the noise onto the one-sided cone of inactive atoms)
  was tested and does **not** close the gap, so the divergence is intrinsic
  active-set instability rather than a linearization artifact patchable at second
  order. Documented as a conservative fragility signal, not a defect.
- **Effective-dof reporting.** `analytic_anisotropy_null` returns an
  anisotropy-aware effective dof (`tr((I−P)D)² / tr(((I−P)D)²)`); it is not yet
  surfaced in the `EmergenceCertificate`.
