# CombiCone manuscript

**Which Combinations Are Worth Running? Certified Triage of Combinatorial
Perturbation Screens with Convex-Cone Emergence Certificates**

An ICML-2026-format manuscript reframed around the CombiCone thesis: use a
convex-cone emergence certificate to decide which combinatorial perturbations
reach cell states that single perturbations cannot, and to triage which
untested combinations are worth measuring.

## Compile

Requires a working TeX installation (TeX Live, MacTeX, or Overleaf). One command:

```bash
bash compile.sh          # pdflatex -> bibtex -> pdflatex x2  ->  main.pdf
```

On Overleaf: upload this folder, set `main.tex` as the root document, compiler
`pdfLaTeX`.

> **Build provenance.** The committed `main.pdf` (20 pages) was compiled with
> TinyTeX v2026.07 (TeX Live 2026) via `bash compile.sh` — final pass clean:
> 0 undefined citations, 0 undefined references, 0 overfull boxes, 65
> bibliography entries resolved. (The conda `texlive-core` package ships
> binaries without the LaTeX macro tree, and the macOS Tectonic build crashes
> in a native system-configuration call, so TinyTeX is the working toolchain
> on this platform.)

## Structure

| File | Content |
|------|---------|
| `main.tex` | Preamble, title, author block, `\input` order |
| `sections/00_abstract.tex` | Abstract — combinatorial-triage thesis |
| `sections/10_introduction.tex` | Two-camps framing; learned-baseline tie as hook |
| `sections/20_methods.tex` | Convex-cone engine (Eq., Farkas/KKT) + noise-aware certification, training-free triage, self-contained learned baseline, k-way generalization |
| `sections/30_results.tex` | Six figures: hook, certificate, triage, head-to-head, k-way, cross-screen |
| `sections/40_related_work.tex` | 91-method survey reframed to combinatorial emergence |
| `sections/50_discussion.tex` | Causal-inference agenda, emergence/triage interpretation |
| `sections/60_ai_declaration.tex` | Honest AI-usage disclosure |
| `sections/90_supplement.tex` | Three supplementary figures + reproducibility note |
| `references.bib` | 106 entries |
| `figures/` | `fig1_hook`, `fig2`–`fig6`, `figS1`–`figS3` (png + pdf) |

## Figures

- **Fig 1 (hook)** — learned MLP ties additive on blind prediction
  (0.896 vs 0.897); only the cone emits a separator. Accuracy is the wrong axis.
- **Fig 2 (certificate)** — the emergence certificate and the measurement-noise
  bar it must clear; recovers classical synergy at partial Spearman 0.62.
- **Fig 3 (triage)** — training-free triage flags emergent combinations 2.4x
  above base rate.
- **Fig 4 (head-to-head)** — learned baseline vs cone; only the cone certifies.
- **Fig 5 (k-way)** — order-k generalization (real reachable-from-lower-order +
  synthetic planted-epistasis triple screen).
- **Fig 6 (cross-screen)** — certificate transfers to the orthogonal CaRPool-seq
  Cas13d screen.

## Provenance note

The Norman screen H5AD used here is labelled cell type **A549**, whereas the
canonical Norman 2019 line is **K562**; this discrepancy is carried verbatim
through the manuscript. All headline numbers are single-sourced from the
Phase-2 deliverable artifacts (see `../docs/metrics/`).
