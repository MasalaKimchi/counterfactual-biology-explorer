# Counterfactual Biology Explorer

**Built with Claude — Life Sciences Hackathon (Research Track)**

*What is the minimal set of gene perturbations that shifts a primary human CD4+ T cell from one transcriptional state toward a target "healthier" state — and how confident should we be?*

---

## TL;DR

Rather than reporting *what is different* between two cell states (differential
expression), this project generates *interpretable, falsifiable counterfactual
hypotheses* about **what changes might move a cell toward a target state**. It does
this by treating each genome-scale CRISPRi perturbation as a measured causal
"effect vector" and asking a geometric question: **is a target state reachable by some
combination of knockdowns**, and if so, what is the **smallest set** that gets closest?
The headline output is a falsifiable **reachable / provably-outside-the-cone** verdict,
not just a similarity score. Every hypothesis ships with an explicit **confidence
score** built from the dataset's own reproducibility metrics, held-out validation, and
orthogonal literature/database evidence.

## Why this framing (and an honest scope note)

The hackathon abstract originally targeted Acute Myeloid Leukemia (AML) using the
Marson/Pritchard Perturb-seq dataset. **That dataset is a genome-scale CRISPRi
Perturb-seq screen in primary human CD4+ T cells — not AML/HSPC data.** We
therefore reframe the identical *method* onto a state transition that this dataset
can actually support:

- **Primary axis:** Th1 ↔ Th2 helper-T polarization balance
- **Secondary axis:** CD4+ T-cell aging signature (aged → young-like)

Both target-state signatures are provided directly by the dataset authors. The
polarization signature ships from **two independent source contrasts** (Ota 2021,
Höllbacher 2021), which we use for a cross-source robustness check in lieu of an
arrayed wet-lab validation table (not present in our local data). See
[`ROADMAP.md`](./ROADMAP.md) for the full rationale and the analysis catalog.

## How this maps to the hackathon (from the event page)

*Verified against the [event details page](https://cerebralvalley.ai/e/built-with-claude-life-sciences)
on 2025-07-07. The page publishes **no scored judging rubric** — only that finalists
are chosen by a panel (Anthropic + our partner, Gladstone Institutes) and that those
selections are final. The priorities below are read off the **stated deliverable and
prize structure**, not an official criteria list.*

- **Track — Research ("Build From the Bench").** The mandate is to *start from a
  biological question, use Claude Science to answer it, and submit something discrete —
  "a finding, a trained model, an analysis others can reproduce."* A reproducible
  analysis with figures is exactly the expected shape; this repo is built to that.
- **Featured dataset, verbatim prompt.** This CD4+ T-cell Perturb-seq screen (Marson +
  Pritchard/Stanford) is a highlighted Research-track dataset, framed as **"Find new
  drug targets in this T cell Perturb-seq dataset."** Our knockdown nominations *are*
  target hypotheses — align the framing accordingly.
- **Gladstone Special Prize.** A dedicated award for the entry with **"most potential
  to advance science that can overcome disease"** ($10k). Our autoimmune-enrichment +
  Open Targets disease layer is the direct hook for it — make disease relevance
  prominent, not an appendix.
- **Logistics.** Fully virtual, 7-day build, team size ≤ 2; each participant gets one
  month of Max 20x + $200 API credits. Prize pool $100k in credits across both tracks.
- **What this implies for us.** With no published rubric, we optimize for what the page
  *does* reward: a **discrete, reproducible deliverable** (fixed seed + pinned env), a
  **new-drug-target** reading of the output, and a **disease-impact** narrative — all
  wrapped in the honest-limitations stance below so the claims survive expert scrutiny.

## Data

Marson & Pritchard labs, *Genome-scale perturb-seq in primary human CD4+ T cells*
(Zhu et al., bioRxiv 2025). Hosted on the CZI Virtual Cells Platform under an MIT
license.

- Dataset card: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- Analysis repo & supplementary tables: https://github.com/emdann/GWT_perturbseq_analysis_2025

**Feasibility on a laptop (CPU-only):** the raw dataset is ~22M cells and is *not*
tractable on a laptop. We build on the authors' **precomputed derived artifacts**,
and the work is split into two tiers by what's needed.

**Tier 1 — already local (checked into `data/`, ~36 MB of CSVs):** enough for the
**warm-up + zero-dependency fallback** (see the Tier note under Tier 2) — target
signatures, per-perturbation QC/effect summaries, guide QC, off-target design,
autoimmune-disease enrichment, and donor metadata. *Tier 1 is a 1-D diagonal proxy, not
the headline method; the graded core is the Tier-2 reachability cone.*

| Local file | Rows | What it is |
|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | Per **perturbation × condition** summary (DE-gene counts, on/off-target flags, cross-donor & cross-guide reproducibility). *Summary, not the gene-level matrix.* |
| `Th2_Th1_polarization_signature…csv` | 37,288 | Th2→Th1 target signature, two source contrasts |
| `CD4T_aging_signature…csv` | 10,000 | Aged→young target signature |
| `guide_kd_efficiency.suppl_table.csv` | 73,765 | Per-guide knockdown QC |
| `sgrna_library_metadata.suppl_table.csv` | 31,110 | Guide design / off-target annotation |
| `cluster_autoimmune_enrichment…csv` | 5,236 | Perturbation-cluster × 17 autoimmune diseases |
| `sample_metadata.suppl_table.csv` | 12 | 4 donors × 3 conditions, demographics |

**Tier 2 — one download away (not local yet):** the full gene-level effect matrix,
needed for the **reachability-cone + reconstruction solver** — i.e. the headline method.
*The graded core should run here, not on the Tier-1 proxy (see the Tier note below).*

| Artifact | Shape | What it is |
|---|---|---|
| `GWCD4i.DE_stats.h5ad` | 33,983 pert×cond × 10,282 genes | Per-perturbation effect matrix (logFC, z-score, p) — **input to the reachability + reconstruction solver** |
| `GWCD4i.pseudobulk_merged.h5ad` | guide×donor×cond × 18,129 genes | Pseudobulk expression profiles |

Working from the ~34k × 10k matrix (subset to significant on-target perturbations +
HVGs) is comfortably CPU-tractable once fetched. See [`data/README.md`](./data/README.md).

## Method (baseline-first, honestly benchmarked)

1. **Perturbation dictionary** `E ∈ R^{P×G}`: each row is one perturbation's measured
   causal effect on the transcriptome (z-scored logFC), per stimulation condition.
2. **Target direction** `d ∈ R^G`: the desired transcriptomic shift (e.g. the Th2→Th1
   signature, or the reverse-aging vector).
3. **Reachability first, then minimal set.** Because CRISPRi is loss-of-function and
   weights are **non-negative**, the reachable transcriptome shifts form a **convex
   cone**. We first solve a non-negative least-squares fit of `d` inside that cone to get
   the *reachable-vs-outside* verdict + closest point + residual; then apply sparse
   selection (LASSO / OMP / greedy forward) *within* the cone to find the **smallest**
   perturbation set that reconstructs the reachable component. Sparsity operationalizes
   "minimal"; non-negativity operationalizes "achievable by knockdown."
4. **Confidence** = combination of (a) dataset-native reproducibility (cross-guide and
   cross-donor correlation, on-target knockdown significance, off-target flags),
   (b) bootstrap **stability selection** frequency, (c) held-out-donor generalization,
   and (d) orthogonal **literature/Open Targets** evidence for the nominated gene.

A deliberately simple **linear/additive baseline is the primary model, not a
fallback.** A 2025 *Nature Methods* benchmark (Ahlmann-Eltze, Huber & Anders) shows
current deep-learning perturbation predictors do **not** yet beat simple linear
baselines — so any DL component here is an explicitly-optional comparison, never an
unbenchmarked claim.

## Related work — and our specific novelty delta

"Find a minimal set of perturbations that moves a cell toward a target state" is a
**named, established problem**, so our contribution has to be stated as a precise
*delta*, not a category. The two closest prior methods:

- **Mogrify** (Rackham et al., *Nat Genet* 48:331–335, 2016) "combines gene expression
  data with regulatory network information to predict the reprogramming factors
  necessary to induce cell conversion," applied across 173 human cell types, and
  validated two new transdifferentiations. That is the minimal-set-to-target idea — but
  with **transcription factors / gain-of-function** and a **regulatory-network** basis.
- **CellOracle** (Kamimoto et al., *Nature* 614:742–751, 2023) uses **inferred GRNs** to
  simulate KO/overexpression effects, "converted into a vector map of transitions in
  cell identity," benchmarked against a **randomized-model null**. That is the
  perturbation-to-state-transition idea, with null-model hygiene close to ours.

**Our delta (what is genuinely new here):**

1. **Measured, not inferred, effect vectors.** Our perturbation dictionary is the
   *empirically measured* genome-scale CRISPRi effect matrix — not a GRN inferred from
   wild-type data (CellOracle) nor a network-influence heuristic (Mogrify).
2. **A convex-cone reachability verdict.** Loss-of-function effects + non-negative
   weights span a convex cone; we return a formal **reachable vs. (provably) outside the
   knockdown cone** answer with the closest reachable point and a residual — not just a
   similarity score. *Precise statement: "unreachable" means the target has a component
   outside the non-negative span of the measured knockdown effect vectors (which are
   themselves mixed-sign); this often — not always — corresponds to needing activation.*
3. **Held-out-gene validity + a screen-native confidence decomposition** (see Method §4).

Note also that the biology our pipeline surfaces (Th1/Th2 and aging regulators) is what
the **source paper itself reports** — the Zhu et al. abstract states perturbation
signatures let them nominate "regulators of Th1 and Th2 polarization and of age-related
T cell phenotypes." We therefore treat regulator recovery as **validation that the
method works**, and claim novelty only for the method + decision layer above — never for
the regulator lists themselves.

## Honest limitations (read before trusting any output)

- **CRISPRi is loss-of-function only.** We can directly nominate *knockdowns*.
  Gain-of-function hypotheses are extrapolations that this assay cannot test.
- **Additivity is an assumption.** The screen perturbs single genes; any multi-gene
  "minimal set" is an *untested combinatorial extrapolation* that ignores epistasis.
- **Transcriptome ≠ phenotype.** Matching a transcriptional signature does not prove
  functional rescue.
- **The target state is a proxy.** "Healthier" is operationalized as a transcriptomic
  signature, not a clinical outcome.
- Outputs are **ranked, falsifiable hypotheses for future experimental validation** —
  not conclusions.

## Status: docs live on `HEAD`; the validated prototype is one `git checkout` away

The working tree currently holds the design docs plus the local Tier-1 data. The
runnable pipeline (`src/`, `app/`, `tests/`, build config) is **not deleted — it is
preserved in git history** and does **not** need to be rebuilt from scratch. It was
removed from `HEAD` by the two "Remove legacy pipeline code" commits, so it lives intact
at the commit immediately before them:

```bash
# restore the entire validated prototype (~1,000 LOC) onto the working tree:
git checkout 3f8db17 -- src app tests environment.yml requirements.txt data/fetch_de_stats.sh
python -m src.data_loader --check          # schema gate should pass on all 7 CSVs
pytest -q                                  # smoke tests should be green
```

Recoverable at `3f8db17` ("Add Tier-2 fetch script and update pipeline modules"):

| Path | ~LOC | What it is |
|---|---|---|
| `src/data_loader.py` | 170 | schema gate + CSV loaders |
| `src/target_states.py` | 175 | target vectors, sign fix, cross-source concordance |
| `src/counterfactual.py` | 208 | the ranking / reconstruction engine |
| `src/confidence.py` | 123 | confidence decomposition |
| `src/evidence.py` | 129 | PubMed / Open Targets / Enrichr |
| `app/explorer.py` | 80 | Streamlit explorer |
| `tests/test_counterfactual.py` | 46 | smoke tests |
| `data/fetch_de_stats.sh` | 39 | Tier-2 (`GWCD4i.DE_stats.h5ad`) fetch script |
| `environment.yml`, `requirements.txt` | — | pinned environment |

**So the week's first move is `restore + green tests`, not `rewrite`.** Current tree:

```
cell-state-reachability/
├── README.md             # this file — framing, method, related work, hackathon fit
├── ROADMAP.md            # 7-day plan, reachability reframe, risks, live checklist
├── LICENSE               # MIT
├── data/
│   ├── README.md         # how the data is sourced (no-auth first)
│   └── *.suppl_table.csv # Tier-1 supplementary tables (local, gitignored)
└── notebooks/README.md   # notes for exploratory analysis (no notebooks yet)
```

The method, tiers, and analysis catalog live in [`ROADMAP.md`](./ROADMAP.md); the data
provenance lives in [`data/README.md`](./data/README.md).

## License

MIT. Data: MIT (CZI Virtual Cells Platform). Please cite Zhu et al. 2025.
