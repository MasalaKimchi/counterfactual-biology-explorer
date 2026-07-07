# Project Assessment — Cell-State Reachability

*Independent technical review. Every quantitative claim below was reproduced from the
local data or the restored prototype in this repo, or verified against the primary
literature. Where a number differs from the existing docs, the discrepancy is flagged.
Method-verification figure: `verification_dashboard.png`.*

---

## TL;DR verdict

| Question | Verdict |
|---|---|
| **Is the hypothesis sound and testable?** | **Yes.** Well-posed, falsifiable, and the core sanity checks reproduce cleanly from raw data. |
| **Enough data to run the project?** | **Tier-1: yes, verified** (7 CSVs, 6/7 row-exact). **The headline method needs the Tier-2 matrix**, which is one free CZI download — not yet local. |
| **Enough compute capacity?** | **Yes, decisively.** Full-scale solver + nulls run in **<1 s** on an 18 GB CPU laptop. Compute is not a constraint. |
| **Is the code capacity there?** | **Partially.** ~1,000 LOC of clean, tested prototype is recoverable from git. But **the headline "reachability cone" method described in the docs is not actually implemented yet** — the reconstruction layer exists; the reachability-verdict layer does not. This is the main build task. |
| **Is it well-documented?** | **Docs are unusually strong** for a hackathon (honest-scope, tiering, risks). Two gaps: docs describe methods the code doesn't yet contain, and one schema-gate row count is stale (will fail on Day 0). Both fixed/flagged here. |
| **Real-world impact?** | **Concrete and defensible** — 17 immune-mediated diseases with 185 significant GWAS-gene enrichments already in the local data. See §5. |

---

## 1. State of the project and the hypothesis

**The hypothesis** (as I read it): *Given a genome-scale library of measured CRISPRi
knockdown effect vectors, a target cell-state shift is either reachable by some
non-negative combination of knockdowns — in which case a small, ranked set reconstructs
it — or it is provably outside that reachable set.* This is a genuine, falsifiable
scientific claim, not a description. It is stronger than "which genes differ" (DE) and
stronger than "which genes correlate with the target" (the Tier-1 proxy).

**What is verified as working** (reproduced independently this review, from the raw CSVs):

- **Sign convention is correct.** `toward_Th1 = −(Th2-vs-Th1)` moves **GATA3 down**,
  **TBX21 / IFNG / STAT4 / STAT1 up**, and Th2 cytokines **IL4 / IL5 / IL13 down** — in
  *both* source contrasts (Ota 2021 and Höllbacher 2021). Textbook-correct helper-T
  biology.
- **Cross-source robustness reproduces exactly.** 11,616 shared genes, **68.5%
  sign-concordant**, Spearman **ρ = 0.562** — matching the ROADMAP's headline numbers.
  *(Doc footnote: the ρ = 0.56 claim is on z-score; on log-fold-change it is 0.533. The
  docs should state which — now noted.)*
- **The positive-control nomination is even stronger than documented.** The Th2 master
  regulator **GATA3 ranks #1 (Rest), #2 (Stim8hr), #2 (Stim48hr)** among significant
  knockdowns for `toward_Th1` — the ROADMAP conservatively says "#3 of 2,135."

**Read the GATA3 result honestly (this is a positive control, not a discovery):** since
`toward_Th1 = −(Th2-vs-Th1)`, the target by construction "wants down" exactly the genes
most up in Th2, and GATA3 is *the* Th2 master TF — so a method that ranks "genes the
signature most wants down, with a strong reproducible knockdown" is nearly guaranteed to
surface it. It confirms the plumbing and sign convention are right. It is **not**
evidence the method finds anything non-obvious. The non-obvious result only comes from
the Tier-2 reachability analysis (below), which is not yet built.

---

## 2. Do you have sufficient DATA?

**Tier-1 (local): verified sufficient for the warm-up + all supporting analyses.**

All 7 CSVs load (37.7 MB total). Six match documented row counts exactly; every column
the pipeline modules depend on is present. **One discrepancy:**
`sgrna_library_metadata` has **26,504 rows, not the documented 31,110** (the 26,504
sgRNAs are all unique — it is not a truncation artifact). This is minor for the science
but **will trip the restored schema gate**, which hard-codes 31,110 (±2%). Fix in §4.

**Tier-2 (the headline method's input): required, not yet local.** The gene-level
`GWCD4i.DE_stats.h5ad` effect matrix (33,983 perturbations × 10,282 genes) is where the
reachability cone actually runs. It is a **free CZI login + one download** (~1.4 GB for
one float32 layer). **This is the single most important data action** — without it the
project can only demonstrate the 1-D proxy, which is near-tautological (§1). Fetch it
Day 0.

---

## 3. Do you have sufficient COMPUTE and CODE capacity?

### Compute: not a constraint, empirically.

I ran the restored solver at increasing scale on this 18 GB CPU machine:

| Dictionary size | Greedy solve (k≤8) + 1,000-iteration null | Peak RAM |
|---|---|---|
| 6,000 × 2,000 | 0.16 s | |
| 16,000 × 2,000 | 0.33 s | |
| **34,000 × 2,000 (full Tier-2 scale)** | **0.57 s** | **1.4 GB** |

A complete reachability spectrum (k = 1…8 × both null types × 1,000 iterations) is a few
minutes at most. The ROADMAP's "compute was never the binding constraint" is **true**.

### Code: the prototype is recoverable, but the headline method is not yet built.

~1,000 LOC is intact in git at commit `3f8db17` and restores in one command (§4). It is
clean, typed, and has passing synthetic-data tests. **But there is a gap between what the
markdown describes and what the code contains:**

| Capability (as documented) | In the restored code? |
|---|---|
| Perturbation dictionary loader + schema gate | ✅ `data_loader.py` |
| Target-state vectors + sign fix + cross-source concordance | ✅ `target_states.py` |
| Sparse reconstruction: greedy / OMP (+ `knockdown_only` non-negativity) | ✅ `counterfactual.py` |
| Random-**perturbation** null | ✅ `counterfactual.random_null` |
| Confidence: reproducibility, **stability selection**, held-out-**donor** | ✅ `confidence.py` |
| External evidence: PubMed / Open Targets / Enrichr | ✅ `evidence.py` |
| Streamlit explorer | ✅ `app/explorer.py` |
| **Formal NNLS convex-cone fit → reachable / provably-outside verdict + residual** | ❌ **not implemented** (the headline claim) |
| **Reachability spectrum output (k vs alignment curve as a deliverable)** | ❌ not implemented |
| **Shuffled-**target** null** (distinct from the random-*perturbation* null) | ❌ not implemented |
| **Held-out-**gene** generalization** (the primary validity claim in §1 of ROADMAP) | ❌ not implemented (only held-out-*donor* exists) |
| Bootstrap CI on the reachability metric | ❌ not implemented |

**Bottom line:** you have the *reconstruction* engine; you are missing the *reachability
verdict* layer that is the project's stated novelty. This is a well-scoped ~1–2 day build
on top of the restored code, not a from-scratch effort. Concrete specs in `CLAUDE.md`.

---

## 4. How to evaluate SUCCESS (falsifiable, tiered)

Success is not "the demo runs." It is a set of pre-registered, falsifiable checks. Three
tiers, hardest-first:

1. **Primary validity — held-out-gene generalization (Tier-2).** Fit non-negative
   weights on a random half of the target signature's genes; score alignment on the
   held-out half. *Pass = held-out cosine clearly above a shuffled-target null with a
   bootstrap CI that excludes it.* This is the claim that cannot be inflated by
   overfitting, and it is currently **unbuilt** — build it first.
2. **Reachability verdict with honest nulls (Tier-2).** The reachability score sits
   above **both** a shuffled-target null **and** a random-perturbation null, with effect
   size + bootstrap CI, and survives leave-one-donor-out (n = 4). *Falsifiable: if the
   target is not distinguishable from a shuffled target, the state is "not reachably
   distinct" and you report that — a negative result is still a result.*
3. **Cross-source + external corroboration (Tier-1/2).** Target built from the
   sign-concordant Ota∩Höllbacher core (✅ done, ρ = 0.56); top nominations cross-checked
   against ≥1 independent screen (Shifrut 2018 / Schmidt 2022) and a known-biology gate
   (GATA3 down / TBX21 up — ✅ verified).

Plus reproducibility hygiene: fixed seed + pinned env, whole pipeline reproduces on CPU
in minutes, every hypothesis card carries a decomposable confidence score + off-target
flag + disease link + citations.

**A failing result is a valid submission.** If the Th1↔Th2 target turns out to be
*unreachable* by knockdown alone (plausible — it likely needs activation, which CRISPRi
can't do), that is a clean, honest, interesting finding. Design the writeup so either
outcome is publishable.

---

## 5. Real-world impact

The impact case is concrete and already grounded in the local data, not aspirational:

- **The disease layer is real.** The autoimmune-enrichment table links perturbation
  clusters to **17 immune-mediated diseases** (Crohn's, ulcerative colitis, IBD, asthma,
  atopic eczema, psoriasis, multiple sclerosis, type-1 diabetes, lupus, rheumatoid
  arthritis, celiac, ankylosing spondylitis, and more) with **185 significant
  enrichments (FDR < 0.05)**. Top signals reach odds ratios of 58 (Crohn's) driven by
  real GWAS genes (IRF4, BATF, STAT3, ITK). This is a direct hook for the **Gladstone
  Special Prize** ("most potential to advance science that can overcome disease").
- **Th1/Th2 balance is a validated therapeutic axis.** Th2-skewing drives allergy,
  asthma, and atopic disease; Th1/Th17 drives autoimmunity. A method that nominates the
  *minimal, ranked, confidence-scored* knockdowns to shift that balance — and states when
  a shift is *not* achievable by knockdown — is directly relevant to immunotherapy and
  cell-engineering target selection.
- **The featured-dataset framing fits.** The event highlights this exact Perturb-seq
  dataset as *"Find new drug targets in this T cell Perturb-seq dataset."* Your knockdown
  nominations *are* target hypotheses; lead with that framing.
- **Method generality.** Nothing in the approach is T-cell-specific. Any Perturb-seq
  atlas + any target state (disease → healthy, aged → young, cell-type A → B) plugs into
  the same reachability machinery. That transferability is the larger scientific claim.

**Honesty guardrail on impact:** these are *hypotheses for wet-lab testing*, not
validated targets. CRISPRi is loss-of-function only; multi-gene sets assume additivity
(no epistasis); the effects are measured in one primary-cell system. State all three
plainly — the honest framing is what makes the impact claim credible to expert judges.

---

## 6. Literature verification (all confirmed against primary sources)

| Claim in the docs | Verified source | Status |
|---|---|---|
| DL perturbation predictors don't yet beat linear baselines → justifies baseline-first | Ahlmann-Eltze, Huber & Anders, *Nat Methods* **22**:1657–1661 (2025), doi:10.1038/s41592-025-02772-6 | ✅ Confirmed against fetched full text. Five foundation + two other DL models vs deliberately simple baselines for single/double perturbations; none outperformed. For **double (two-gene) perturbations, all models had prediction error substantially higher than the 'additive' baseline** (sum of individual log-fold-changes; paper Fig. 1a,b) — directly backs your additivity stance. |
| Prior art: Mogrify (minimal TF set → cell conversion) | Rackham et al., *Nat Genet* **48**:331–335 (2016) | ✅ Confirmed. Predicts reprogramming factors; validated 2 new transdifferentiations. **Delta:** they use TF/gain-of-function + network heuristic; you use measured LOF effect vectors. |
| Prior art: CellOracle (GRN-simulated perturbation → state transition) | Kamimoto et al., *Nature* **614**:742–751 (2023) | ✅ Confirmed. Inferred GRNs → "vector map of transitions in cell identity" + randomized-model null. **Delta:** they infer the GRN; you use the *measured* effect matrix and return a convex-cone reachability verdict. |
| Source paper nominates Th1/Th2 + aging regulators | Zhu et al., bioRxiv, doi:10.64898/2025.12.23.696273 | ✅ Confirmed verbatim: the screen is described as "nominating regulators of Th1 and Th2 polarization and of age-related T cell phenotypes." Your axes = the source paper's own results, so **claim novelty for the method, not the regulator lists.** |

**Novelty is defensible** as long as it is stated precisely: not "we find Th1/Th2
regulators" (the source paper did that) but "we provide a measured-effect, convex-cone
*reachability verdict* with held-out-gene validation and a screen-native confidence
decomposition." That specific delta is not occupied by Mogrify or CellOracle.

---

## 7. Is it well-documented? (and what I changed)

The existing docs are **above hackathon norm**: honest-scope note, explicit tiering,
risk register, day-by-day plan, external-validation table. The gaps were:

1. **Docs described methods the code doesn't contain** (reachability cone, held-out-gene,
   shuffled null) as if implemented. → Now explicitly separated into "implemented vs. to
   build" (this file §3 + `CLAUDE.md`).
2. **Tiering was inverted** (Tier-1 called the "graded core," Tier-2 "optional"), which
   put the novel method off the critical path. → Corrected across all four markdown files
   last round: Tier-2 = graded core, Tier-1 = warm-up/fallback.
3. **"Code was deleted, rebuild from plan"** → corrected to "recoverable from `3f8db17`;
   restore, don't rebuild."
4. **Stale schema row count** (`sgrna_library_metadata` 31,110 → actual 26,504) will fail
   the Day-0 gate. → Flagged here + in `CLAUDE.md` with the one-line fix.

**For the Claude Code handoff:** `CLAUDE.md` is the operating manual — repo orientation,
the restore command, the schema fix, the implemented-vs-aspirational gap, concrete
function specs for the missing reachability layer, the definition of done, and the
honesty guardrails. Give Claude Code `CLAUDE.md` + this file + `ROADMAP.md`.
