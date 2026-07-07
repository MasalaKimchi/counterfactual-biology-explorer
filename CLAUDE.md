# CLAUDE.md — operating manual for Claude Code

*This file orients an autonomous coding agent working in this repo. Read it fully before
touching code. Companion docs: `ASSESSMENT.md` (why the project is sound + the
implemented-vs-missing analysis), `ROADMAP.md` (7-day plan + analysis catalog),
`README.md` (framing), `data/README.md` (how to fetch data).*

---

## 0. What this project is (one paragraph)

Genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells gives a *measured* causal
effect vector for each gene knockdown. This project treats those vectors as a dictionary
and asks a geometric question: **is a target transcriptomic state (e.g. Th2→Th1
polarization) reachable by some non-negative combination of knockdowns, and if so, what
is the smallest ranked set that reaches it?** Because knockdowns are loss-of-function and
weights are non-negative, the reachable shifts form a **convex cone**; the headline
output is a falsifiable **reachable / provably-outside-the-cone** verdict with a residual,
each nomination carrying a decomposable confidence score. This is *not* differential
expression and *not* a similarity ranking — it is a reachability decision.

---

## 1. FIRST ACTIONS (do these before anything else)

```bash
# (a) Restore the validated prototype (~1,000 LOC) from git history — DO NOT REWRITE IT.
#     It was removed from HEAD but is intact at 3f8db17.
git checkout 3f8db17 -- src app tests environment.yml requirements.txt data/fetch_de_stats.sh

# (b) Create the environment (CPU-only; no GPU/torch needed).
conda env create -f environment.yml   # or: pip install -r requirements.txt

# (c) APPLY THE SCHEMA-GATE FIX below BEFORE running the gate, or it will fail.

# (d) Run the gate + tests.
python -m src.data_loader --check      # schema gate over the 7 local CSVs
pytest -q                              # synthetic-data solver tests (<1s)
```

### 1.1 Known Day-0 bug — fix before running the gate

`src/data_loader.py` hard-codes an expected row count for the sgRNA table that no longer
matches the local file. **Verified actual counts** (this repo's `data/*.suppl_table.csv`):

| file | actual rows | what the gate expects |
|---|---|---|
| `sgrna_library_metadata.suppl_table.csv` | **26,504** | 31,110 ← **stale, will fail** |

Fix: in the `LOCAL_CSV_CONTRACT` dict in `src/data_loader.py`, change the
`sgrna_library_metadata` expected count from `31110` to `26504`. (All other 6 tables
match exactly: DE_stats 33,983 · Th2_Th1 37,288 · aging 10,000 · guide_kd 73,765 ·
autoimmune 5,236 · sample_metadata 12.)

---

## 2. Repo map

```
cell-state-reachability/
├── CLAUDE.md            ← you are here (operating manual)
├── ASSESSMENT.md        ← technical review + implemented-vs-missing gap analysis
├── ROADMAP.md           ← 7-day plan, tiering, risks, analysis catalog
├── README.md            ← framing, method, related work, hackathon fit
├── data/
│   ├── README.md        ← Tier-1 (local) + Tier-2 (fetch) instructions
│   └── *.suppl_table.csv (7 files, ~37.7 MB, local, gitignored)
├── src/     (restore from 3f8db17)   app/  (restore)   tests/ (restore)
└── notebooks/README.md
```

### 2.1 Restored modules and their public API (verified from `3f8db17`)

- **`src/data_loader.py`** — `PerturbDictionary` (dataclass), `load_de_dictionary(...)`
  (loads one condition slice of the Tier-2 matrix, subsets to significant on-target
  perturbations + HVGs), `validate_local_data(...)` (the schema gate), `_check()`.
- **`src/target_states.py`** — `polarization_target(direction="toward_Th1")`,
  `aging_target(...)`, `concordance_report(value_col="zscore")`, `normalize(d)`. Reads
  the correct gene column per file (`variable` for polarization, `gene_name` for aging),
  builds the target from the sign-concordant Ota∩Höllbacher core, fixes the sign.
- **`src/counterfactual.py`** — `greedy_minimal_set(E, d, k_max, knockdown_only=False)`,
  `omp_minimal_set(E, d, k_max)`, `random_null(E, d, k, n_iter)`. Returns `MinimalSet`
  (indices, weights, cosine, residual_norm). `knockdown_only=True` forbids negative
  weights (no invented activation).
- **`src/confidence.py`** — `reproducibility_score(row)`, `stability_selection(...)`,
  `heldout_donor_consistency(...)`, `ConfidenceReport`.
- **`src/evidence.py`** — `pubmed_comention(gene, process)`, `opentargets_diseases(gene)`,
  `collect(gene, process)`, `enrich(gene_list)`. All no-auth public APIs (NCBI E-utils,
  Open Targets GraphQL, Enrichr).
- **`app/explorer.py`** — Streamlit hypothesis-card UI.
- **`tests/test_counterfactual.py`** — synthetic-data recovery + knockdown-only +
  beats-random-null tests. Run first; they must stay green.

---

## 3. THE BUILD TASK — what exists vs. what to implement

The restored code is the **reconstruction** engine. The project's stated novelty — the
**reachability verdict layer** — is **not yet implemented**. This is the primary work.
Do NOT re-implement what already exists; build the missing pieces on top.

### 3.1 Already implemented (reuse, don't rewrite)
✅ dictionary loader + schema gate ✅ target vectors + sign fix + concordance
✅ greedy/OMP reconstruction with non-negative (`knockdown_only`) weights
✅ random-**perturbation** null ✅ reproducibility + stability-selection + held-out-**donor**
✅ PubMed/Open Targets/Enrichr evidence ✅ Streamlit explorer

### 3.2 To implement (the headline method — priority order)

1. **NNLS convex-cone reachability verdict** — `src/reachability.py`, new.
   ```python
   def reachability(E, d, *, hvg_mask=None):
       """Non-negative least squares fit of target d inside the cone spanned by
       knockdown effect rows E (P×G). Returns:
         reachable_cosine : cosine(E_S^T w*, d) at the NNLS optimum
         residual_norm    : ||E^T w* - d|| / ||d||   (0 = fully reachable)
         verdict          : 'reachable' | 'partially' | 'outside'   (thresholds
                            calibrated against the shuffled-target null, NOT hardcoded)
         w                : non-negative weights (the ranked nominations)
       Use scipy.optimize.nnls or a non-negative Lasso path. CPU-cheap (<1s at 34k×2k)."""
   ```
   The verdict is the deliverable. "outside" means d has a component outside the
   non-negative span of the (mixed-sign) knockdown vectors — a *falsifiable* negative.

2. **Shuffled-target null** — add to `src/reachability.py`. Distinct from the existing
   random-*perturbation* null: permute the gene labels of `d` (n≥1,000), recompute the
   reachability cosine each time, report the observed value's percentile + bootstrap CI.
   This is what makes the verdict honest.

3. **Held-out-GENE evaluation** — `src/validation.py`, new. Fit non-negative weights on a
   random half of the target's genes; score alignment on the held-out half; repeat over
   folds/seeds. **This is the primary validity claim (ROADMAP §1) and is currently
   missing.** It cannot be inflated by overfitting — build it first after the cone fit.

4. **Reachability spectrum** — sweep k = 1…K (or the NNLS L1 path), record
   cosine-vs-k / residual-vs-k as a saved table + figure. The "minimal set" is the knee.

5. **Bootstrap CI** on the reachability metric (resample genes or donors).

6. **Wire the above into the confidence card + Streamlit explorer.**

### 3.3 Verified quantitative facts (trust these; they were reproduced from the data)

| fact | value | note |
|---|---|---|
| Local CSVs | 7 files, 37.7 MB | 6/7 row-exact; sgRNA is 26,504 (see §1.1) |
| Sign convention | `toward_Th1` → GATA3 ↓, TBX21/IFNG/STAT4/STAT1 ↑, IL4/IL5/IL13 ↓ | both source contrasts |
| Cross-source concordance | 11,616 shared genes · 68.5% sign-concordant · ρ=0.562 (z-score) | ρ=0.533 on log-fc |
| GATA3 rank (`toward_Th1`) | #1 Rest, #2 Stim8hr, #2 Stim48hr | positive control, NOT a discovery |
| Solver speed @ 34k×2k | 0.57 s greedy + 1,000-null; 1.4 GB RAM | compute is not a constraint |
| Autoimmune enrichment | 17 diseases, 185 hits FDR<0.05, top OR 58 (Crohn's) | disease-impact layer |
| Conditions | Rest, Stim8hr, Stim48hr | in DE_stats `culture_condition` |
| Donors | 4 (age 22–34) | for leave-one-donor-out |

---

## 4. DATA — Tier-2 is required for the headline method

- **Tier-1 (local, no auth):** the 7 CSVs. Enough for target vectors, concordance,
  QC, autoimmune enrichment, the 1-D directional proxy, and the whole evidence layer.
- **Tier-2 (fetch early — this is the graded core):** `GWCD4i.DE_stats.h5ad` (33,983
  perturbations × 10,282 genes). Free CZI login; ~1.4 GB for one float32 layer. Use
  `data/fetch_de_stats.sh` (restored from `3f8db17`). Column grain: perturbation ×
  condition. **The reachability cone runs on this matrix; Tier-1 cannot demonstrate the
  method's thesis.** If the download truly fails, Tier-1 is the fallback and you say so.

Key columns in `DE_stats`: `target_contrast_gene_name`, `culture_condition`,
`ontarget_effect_size`, `ontarget_significant` (bool), `offtarget_flag`,
`crossdonor_correlation_mean/min`, `crossguide_correlation`.

---

## 5. DEFINITION OF DONE (falsifiable — see ASSESSMENT.md §4)

A finished submission satisfies:
1. **Held-out-gene**: held-out cosine clearly above shuffled-target null, bootstrap CI excludes it.
2. **Reachability verdict**: above BOTH nulls (shuffled-target + random-perturbation), with effect size + CI, survives leave-one-donor-out (n=4).
3. **Cross-source + external**: target from concordant core (✅), top hits checked vs Shifrut 2018 / Schmidt 2022, GATA3-down/TBX21-up gate passes (✅).
4. Every hypothesis card: decomposable confidence + off-target flag + disease link + citations.
5. Whole pipeline reproduces from **fixed seed + pinned env** on CPU in minutes.

**A negative result is valid.** If Th1↔Th2 is *unreachable by knockdown alone* (plausible
— it may need activation, which CRISPRi can't do), report that cleanly. Design the writeup
so either outcome stands.

---

## 6. GUARDRAILS (non-negotiable — this project's credibility depends on them)

- **Knockdown-only.** CRISPRi is loss-of-function. Never invent activation effects; keep
  weights non-negative (`knockdown_only=True`). If a target needs a gene *up*, the honest
  answer is "not reachable by knockdown" — say it.
- **Additivity is an assumption.** Multi-gene sets assume no epistasis. Flag every
  multi-gene nomination as an extrapolation to be tested, not a prediction.
- **Nominations are hypotheses, not validated targets.** Always framed for wet-lab test.
- **Claim novelty for the METHOD, not the regulator lists.** The source paper (Zhu et al.)
  already reports Th1/Th2 + aging regulators; recovering them = validation, not discovery.
  Novelty = the measured-effect convex-cone reachability verdict + held-out-gene validation.
- **Nulls before claims.** No reachability/alignment number is reported without its null
  and CI. The GATA3 sanity check is a positive control, not evidence of non-obvious signal.
- **No invented judging rubric.** The event page publishes no scored criteria; do not
  optimize for or cite an imagined rubric.
- **Reproducibility hygiene.** Fixed seed everywhere; `fig, ax = plt.subplots()` then
  `fig.savefig()` (never bare `plt.savefig()`); fetches in their own step; pin the env.
- **Do not commit the raw data** (`data/*.csv` is gitignored; force-add only the small
  Tier-1 CSVs if you want a self-contained repo, never the Tier-2 h5ad).

---

## 7. Literature anchors (verified — cite these exactly)

- **Baseline-first justification:** Ahlmann-Eltze, Huber & Anders, *Nat Methods*
  22:1657–1661 (2025), doi:10.1038/s41592-025-02772-6. Five foundation + two other DL
  models vs deliberately simple baselines; none beat them. For double perturbations, all
  models had prediction error *substantially higher* than the additive baseline (Fig. 1a,b)
  — the empirical backing for the additivity assumption. (Verified against full text.)
- **Prior art — Mogrify:** Rackham et al., *Nat Genet* 48:331–335 (2016). Minimal TF set →
  cell conversion (gain-of-function + network). *Your delta:* measured LOF vectors.
- **Prior art — CellOracle:** Kamimoto et al., *Nature* 614:742–751 (2023). Inferred-GRN
  perturbation → cell-identity transition vector + randomized null. *Your delta:* measured
  effect matrix + convex-cone reachability verdict.
- **Source dataset:** Zhu et al., bioRxiv, doi:10.64898/2025.12.23.696273.
- **External validation:** Shifrut et al., *Cell* (2018, PMID 30449619); Schmidt et al.,
  *Science* (2022, PMID 35113687) — CRISPRa arm tests up-regulation hypotheses CRISPRi can't.
