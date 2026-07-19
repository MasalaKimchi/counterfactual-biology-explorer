# CombiCone Phase 2 — Cross-Screen Emergence Panel

**Question.** Does the CombiCone thesis — *certified* emergence of combinatorial
perturbations, with an explicit guard against the signal-to-noise (SNR) confound —
rest on a peculiarity of the single Norman screen, or does the geometry transfer to
an independent screen measured on a different modality, cell line, and perturbation
direction?

**Verdict: the core claims transfer.** The SNR/magnitude confound and the two-bar
certificate reproduce on a fully orthogonal second screen. The training-free triage
*rank* signal also transfers against the raw label, but its magnitude-controlled
component is screen-specific — reported honestly below, not smoothed over.

---

## The two screens

| | **Norman 2019** | **CaRPool-seq (Wessels 2023)** |
|---|---|---|
| GEO accession | GSE133344 | GSE213957 |
| PMID | 31395745 (*Science*) | 36550277 (*Nat. Biotechnol.*) |
| Modality | CRISPRa (transcriptional **activation**, gain-of-function) | Cas13d (RNA-targeting **knockdown**, loss-of-function) |
| Target molecule | DNA (promoter) | mRNA |
| Cell line | K562 (file-label caveat below) | THP-1 (acute monocytic leukemia) |
| Single-gene atoms | 105 | 28 |
| Measured doubles | 131 (all with both singles) | 158 (all 158 with both singles) |
| Genes (HVG space) | 5045 | 2000 |

These screens are orthogonal on **every** axis that could matter: perturbation
*direction* (up vs down), perturbation *target* (DNA vs RNA), CRISPR *effector*
(dCas9-activator vs Cas13d nuclease), and *cell context* (chronic myeloid vs acute
monocytic leukemia). A property shared by both is unlikely to be a screen artifact.

### Provenance caveat (carried verbatim)
- **Norman:** the CZI-redistributed file this project used is labeled
  `cell_type == A549`, whereas Norman 2019 is canonically **K562** (*Science*
  365:786). We report the finding as **"file-label-A549 / canonically-K562"** and do
  not silently assert one. This does not affect the cross-screen conclusion (the
  analysis is within-screen relative geometry).
- **CaRPool:** GEO GSE213957 sample metadata labels the THP-1 arm
  `cell type = "Monocyte isolated from peripheral blood from an acute monocytic
  leukemia patient"` (THP-1), `tissue = blood`, `Homo sapiens`. Reported as the file
  states it.

---

## Method (identical pipeline, `reach-pinned`: numpy 1.26.4 / scipy 1.13.1)

1. **Substrate.** For CaRPool, cells were assigned to perturbation groups from the
   authors' own `CRISPR.Array` / `Guide.Class` annotation (single = gene paired with
   a non-targeting guide; dual = two genes; NT_NT = control). Raw 10x UMI counts
   (4 lanes, 30,707 QC-annotated cells) were CP10k-normalized and `log1p`-transformed;
   pseudobulk mean per group; **effect = group mean − NT_NT control**. Top 2000
   HVGs by cross-condition variance. A deterministic split-half of cells gives two
   independent pseudobulk estimates per condition for the noise model.
2. **Single-effect cone.** The 28 (Norman: 105) single-gene effect vectors are the
   non-negative cone atoms.
3. **Triage** (`combicone.triage_combinations`, training-free `−cos(effA, effB)`),
   scored over the measured pairs.
4. **Certificate** (`combicone.certify_emergence`): weighted-NNLS projection onto the
   cone → residual (`unreachable_fraction`); noise-injection null
   `f0 = fitted + N(0, noise_sd)` with `noise_sd = |split1 − split2| / 2`, B=200 draws
   → plus-one p-value; and `floor_ratio = residual / noise-floor`. Two-bar verdict:
   **bar(a)** p<0.05 (permissive), **bar(b)** floor_ratio ≥ 1.9 (discriminating);
   "certified emergent" requires **both**.

All geometry is **model-relative**: "unreachable" = outside the non-negative cone of
*these* single-gene effects under *this* metric — never a claim of biological
impossibility.

---

## Findings

### 1. The SNR/magnitude confound is universal — and the fix works in both screens
The raw unreachable fraction is strongly anti-correlated with the combination's
effect magnitude (small effects look most "unreachable" because their residual is
dominated by noise):

- Norman: Spearman ρ(raw, ‖eff‖) = **−0.56** (p = 3.6e-12)
- CaRPool: ρ = **−0.82** (p = 3.7e-39) — *even more severe*

The noise-aware **z** removes it in **both**:
- Norman: ρ(z, ‖eff‖) = **+0.13** (p = 0.13, n.s.)
- CaRPool: ρ(z, ‖eff‖) = **+0.12** (p = 0.12, n.s.)

This is the single most important cross-screen result: **reporting the raw residual
as the emergence headline is a trap in every screen, and the noise-aware z is the
honest statistic in every screen.** The project's hard-honesty rule is not
Norman-specific.

### 2. The two-bar certificate transfers
| bar | Norman | CaRPool |
|---|---|---|
| (a) p<0.05 (permissive) | 128/131 (98%) | 158/158 (100%) |
| (b) floor ≥ 1.9× (discriminating) | 40/131 (31%) | 76/158 (48%) |
| **both (certified emergent)** | **40/131 (31%)** | **76/158 (48%)** |

The permissive bar passes almost everything in both screens (as intended — nearly
every measured double departs the cone *somewhat*). The discriminating bar is what
separates real emergence from noise, and it does so in both. CaRPool shows **more**
certified emergence (48% vs 31%): Cas13d knockdown of chromatin/Mediator/SAGA
regulators in AML produces genuinely more non-additive combinations than Norman's
CRISPRa pairs — a biological difference the certificate surfaces, not smooths.

> **Reconciliation — Norman bar counts differ slightly from an earlier project run.**
> An earlier audited pass on the same Norman substrate reported bar(a) = 129/131 and
> bar(b) = **35**/131 (27%); this run reports 128/131 and **40**/131 (31%). The gap is
> **not** simply RNG seed noise. `certify_emergence` uses a Monte-Carlo noise-injection
> null (`n_boot=200`), so `floor_ratio = residual / null_mean` inherits the sampling
> error of `null_mean`. On Norman, **38 of 131 doubles sit in a dense band of
> floor_ratio ∈ [1.7, 2.1]** straddling the 1.9 cut, so the *count* of doubles clearing
> that bar is intrinsically jittery near the threshold. At `n_boot=200` the analytic
> expected count is **39.7 ± 0.6** (SD from threshold jitter; 38 confidently pass,
> 90 confidently fail, 3 within ±2 SE of the cut) — i.e. 40 is stable for *this*
> configuration, and the earlier 35 reflects a different noise-injection configuration
> (e.g. `n_boot`, noise-SD construction, or HVG/atom set), not run-to-run seed drift.
> **Both estimates put Norman's certified rate at ~27–31%, robustly below CaRPool's 48%**
> — the cross-screen conclusion is unchanged. The absolute count near a hard threshold
> is not a stable statistic; the *ordering* (CaRPool > Norman) and the *z*-distribution
> are. This is exactly the kind of threshold-count fragility the noise-aware z was
> introduced to avoid, and it is why the headline compares distributions, not bar tallies.

### 3. Triage: rank signal transfers; the magnitude-controlled part is screen-specific
The training-free `−cos(effA, effB)` triage score predicts the **raw** emergence
label in both screens, and gives ~2.2× top-quartile enrichment in both:
- Norman: ρ(−cos, raw) = +0.47; CaRPool: +0.64. Enrichment 2.2× / 2.2×.

But against the **noise-robust z** label the two screens diverge:
- Norman: ρ(−cos, z) = **+0.37** (p=1.1e-05, survives) — 1.6× enrichment
- CaRPool: ρ(−cos, z) = **+0.08** (p=0.34, n.s.) — 1.6× enrichment

**Honest reading.** The part of the triage signal that is *just* the SNR confound
transfers trivially (both screens). The part that survives magnitude control —
the genuinely useful prospective signal — is present in Norman but not detectable
in CaRPool. Triage is a *rank heuristic whose robust value is screen-dependent*; the
**certificate**, not the prospective score, is the defensible cross-screen capability.

---

## What this does NOT show
- It does **not** show CombiCone predicts unmeasured combinations well; triage is a
  weak rank prior and its magnitude-controlled signal did not replicate in CaRPool.
- It does **not** establish biological interaction mechanism; "emergent" is strictly
  geometric and model-relative.
- The two HVG spaces (2000 vs 5045 genes) and atom counts (28 vs 105) differ, so the
  absolute z / floor_ratio scales are **not** directly comparable across screens —
  only the *within-screen* relationships (confound sign, bar pass-rates, triage ρ)
  are compared. That is exactly the comparison made here.
- CaRPool cell-per-group counts (single median 67, double median 149) are lower than
  Norman's; the noise model absorbs this, but very small groups have wider CIs.

## Reproducibility & CSV notes
- **`n_boot` dependence.** All certificates use `n_boot=200`, `seed=0`. `p_value` and
  `floor_ratio` are Monte-Carlo estimates; near the 1.9 bar the *count* of passing
  doubles wobbles by ~±1 (see the reconciliation above). Rank-based statistics
  (triage ρ, the z-distribution, CaRPool > Norman ordering) are stable; hard threshold
  tallies are not — treat the ~27–31% / 48% rates as ranges, not exact integers.
- **`triage_rank` column.** Both screens now carry a fully populated 1-based
  `triage_rank` (1 = run first, by descending `−cos(A,B)`): Norman 1–131, CaRPool
  1–158. (An earlier CSV draft populated this column for CaRPool only; fixed — sorting
  or filtering the Norman rows by `triage_rank` is now meaningful.)
- **`ci_low`/`ci_high`.** Bootstrap percentile CI on the unreachable fraction,
  populated for every row in both screens.

## Deliverables
- `fig_cross_screen.png` — 2×2 cross-screen panel.
- `cross_screen_emergence.csv` — per-double certificate + triage for all 289 doubles
  (both screens, symmetric columns).
- `carpool_substrate.npz` — reusable CaRPool substrate (atoms, doubles, split-halves).
