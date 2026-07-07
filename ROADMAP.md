# Roadmap — Counterfactual Biology Explorer (7 days, solo, CPU-only)

This is the working roadmap. The polished versions live in the accompanying Word
document and slide deck. This copy is the source of truth for the repo.

## 0. Framing decision (why we pivoted)

The featured Marson/Pritchard dataset is a **CD4+ T-cell** genome-scale CRISPRi
Perturb-seq screen, not AML/HSPC data. Forcing an AML narrative onto it would
require sourcing separate AML scRNA-seq and treating the Perturb-seq only as a weak
prior — more moving parts, weaker validation, higher risk in a one-week solo build.

Instead we keep the *counterfactual method* exactly as proposed and apply it to a
transition the dataset can support **and independently validate**:

- **Primary axis — Th1 ↔ Th2 polarization.** Chosen because the dataset ships an
  arrayed CRISPRi validation table (`Th1Th2_validation_summary.suppl_table.csv`),
  i.e. wet-lab ground truth to check nominated regulators against.
- **Secondary axis — CD4+ T-cell aging (aged → young-like).** The dataset ships an
  aging signature and pre-fit regulator coefficients we can benchmark against.

The counterfactual question becomes: *what is the minimal set of gene knockdowns
whose measured effects best move the CD4+ T-cell transcriptome from state A toward
state B?*

## 1. Success criteria (falsifiable, not vibes)

A submission is "done" if, on a held-out donor:

1. The sparse solver reconstructs the target direction with cosine similarity
   meaningfully above a random-perturbation null (report effect size + CI).
2. Nominated Th1/Th2 regulators overlap the arrayed-validation hits above chance
   (Fisher's exact test, report odds ratio + p).
3. Every hypothesis card renders a confidence score with its component breakdown
   and at least one literature/Open Targets citation.
4. The whole pipeline reproduces from a fixed seed + pinned environment on CPU in
   under ~15 min.

## 2. Day-by-day plan

### Day 1 — Data access, load, sanity
- Register on CZI Virtual Cells Platform; install `vcp-cli`; pull `GWCD4i.DE_stats.h5ad`.
  Fallback with **no auth**: the open supplementary CSVs on the analysis GitHub repo.
- Load DE_stats; subset to `keep_test_genes` perturbations and ~2,000 HVGs to fit RAM.
- **Reproduce a known result** (e.g. a canonical Th regulator's known targets) as a
  correctness gate before building anything new.

### Day 2 — Target-state vectors + QC wiring
- Build `d_polarization` (Th2→Th1) and `d_aging` (aged→young-like) from the
  supplementary signature tables.
- Wire in per-perturbation reproducibility columns (cross-guide, cross-donor,
  on/off-target flags) — these feed confidence later.

### Day 3 — Counterfactual engine
- Implement the sparse minimal-set solver: LASSO, Orthogonal Matching Pursuit, and
  greedy forward selection over the perturbation dictionary.
- Enforce **CRISPRi = knockdown-only** as a hard constraint on admissible directions;
  surface upregulation needs as clearly-labeled non-testable hypotheses.
- Output: ranked minimal sets of size k = 1…N with reconstruction quality.

### Day 4 — Confidence + honest benchmarking
- Confidence module: dataset reproducibility + **bootstrap stability selection** +
  **held-out-donor** generalization.
- Benchmarks: additive linear baseline (primary), random-perturbation null, and
  (optional/stretch) a small scGen-style latent-arithmetic model — reported as a
  comparison only, per the Nature Methods 2025 finding.

### Day 5 — Evidence + interpretation
- `evidence.py`: for each nominated gene, query **PubMed** and **Open Targets** for
  known roles in the target process; attach citations. Evidence is *support*, never proof.
- `pathways.py`: gene-set / GO enrichment on the nominated set for interpretability.
- Cross-check polarization nominations against the arrayed-validation table.

### Day 6 — Explorer UI + hypothesis cards
- Lightweight Streamlit app (`app/explorer.py`): pick an axis and condition → view
  ranked minimal sets → expandable hypothesis cards (confidence breakdown + pathways
  + citations + limitations banner).
- CPU-only friendly; no GPU, no heavy model serving.

### Day 7 — Reproducibility, write-up, demo
- Pin environment, fix seeds, add `tests/`, record a short demo.
- Finalize README, slides, and a one-page limitations statement.

## 3. Risks & mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Dataset registration / S3 access friction | Medium | Open supplementary CSVs on GitHub as no-auth fallback |
| RAM blow-up loading DE_stats | Medium | Load one layer (z-score) as float32; subset to high-confidence perturbations + HVGs |
| Overclaiming causality/rescue | High impact | Hard-coded limitations banner on every output; hypotheses framed as falsifiable |
| Combinatorial predictions untestable | Medium | Flag multi-gene sets as extrapolations; prioritize k=1–2 nominations that map to arrayed validation |
| Scope creep across both axes | Medium | Polarization is the graded deliverable; aging is a stretch |
| MCP auth for some literature tools | Low | PubMed + Open Targets are the primary evidence sources |

## 4. Explicitly out of scope

- Wet-lab validation or any clinical claim.
- Retraining a foundation model or GPU-scale DL.
- Reprocessing the 22M-cell raw matrix from FASTQ/counts.

## 5. Key references

- Zhu et al. *Genome-scale perturb-seq in primary human CD4+ T cells* bioRxiv 2025.
- Ahlmann-Eltze, Huber & Anders. *Deep-learning-based gene perturbation effect
  prediction does not yet outperform simple linear baselines.* Nat Methods 22,
  1657–1661 (2025).
- Lotfollahi et al. scGen (Nat Methods 2019); CPA (Mol Syst Biol 2023);
  Roohani et al. GEARS (Nat Biotechnol 2023).
