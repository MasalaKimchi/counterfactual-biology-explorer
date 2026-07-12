# Reproducibility Check — cell-state-reachability

**Date:** 2026-07-11  **Env:** `cellreach` (Python 3.13.14, numpy 2.4.6, scipy 1.18.0)
**Scope:** every `.py` file, all 8 batch scripts, and all 10 Jupyter notebooks, verified before final submission.

## Verdict

**All code compiles, the core method reproduces bit-for-bit, and all 10 notebooks run end-to-end.**
Four latent defects were found and fixed (details below); none had been caught because each sits on a
code path that the frozen result files did not exercise. All fixes are staged on top of the existing
pre-submission working tree without disturbing any of the other staged consolidation changes.

## Core method — GREEN

| Check | Result |
|---|---|
| `pytest tests/test_reachability.py` | **11 passed** (38-assert packaged self-test + 10 independent property tests) |
| `./reproduce.sh` | **EXIT 0** — "reproduction succeeded"; verdicts / KKT-Farkas certification / signed decomposition all reproduce |
| `.py` compile sweep (14 files) | **all pass** |
| notebook JSON (`nbformat.validate`) | **all 10 valid** |

## Defects found and fixed

All four were present in both `HEAD` (`f99930d`) and the staged working tree — i.e. they would have shipped.

1. **`scripts/a2_conditional_reachability_scaffold.py` — `IndentationError` (would not compile).**
   The consolidation's doc-reference edit (`CAUSAL.md` → "Technical Dossier (Part 4 — Trust & Causal
   Inference)") was accompanied by a reformatting pass that flattened all nested indentation to a single
   space, leaving `for` loops with empty bodies. **Fix:** took the known-good `HEAD` version (correct
   4-space indentation) and re-applied only the intended one-line doc substitution.

2. **`scripts/run_iv_compliance.py` — `IndentationError` (would not compile).** Same root cause
   (99/99 lines reformatted). This is a live batch driver that generates `iv_compliance_verdicts.csv`
   and `late_rescaling_invariance.csv`, read by notebook 09. **Fix:** same method. Verified functionally:
   `run_iv_compliance.py --invariance-only` regenerates `late_rescaling_invariance.csv` **bit-identical**
   to the committed file (max abs diff = 0.0 over all numeric cells).

3. **`notebooks/04_experimental_design_toolkit.ipynb` — stale data path.** Loaded
   `../atlas_work/inputs.npz`; the repo reorg (commit `9b7a2b1`) moved this to
   `../analysis_cache/atlas_work/inputs.npz` but missed this notebook. It failed at its first data cell.
   **Fix:** corrected the path; notebook now runs 12/12 cells and reproduces `cache/design_cards.json`
   byte-identically.

4. **`notebooks/05_target_id_showcase.ipynb` — same stale path** (`os.path.join(REPO, "atlas_work",
   "inputs.npz")`). Failed at cell 2. **Fix:** inserted `"analysis_cache"` into the join; runs 8/8.

**Also harmonized:** `notebooks/06_reinforcement_analyses.ipynb` was the lone notebook whose repo-root
fallback was `os.getcwd()` (works only when launched from the repo root), while the other nine assume
CWD = `notebooks/`. Changed to `os.path.abspath("..")` (matching notebook 09) so all notebooks share one
working-directory convention. Runs 8/8 from `notebooks/` and reproduces all `analysis_cache/nb_out/L*.{csv,json}`
byte-identically.

## Notebooks — all 10 GREEN

Notebooks cannot be executed via `jupyter`/`nbconvert` in this sandbox (the kernel's socket bind is denied).
Each was validated by extracting its code cells into a script run sequentially in one namespace, from
`cwd=notebooks/`, with matplotlib forced non-interactive and the IPython-kernel display builtins
(`display`, `Image`, `HTML`, `Markdown`, `clear_output`) injected exactly as a live kernel provides them.

| Notebook | Code cells | Result |
|---|---|---|
| 01_exploratory_data_analysis | 15 | **OK** |
| 02_reachability_on_tier2 (heavy: full-E NNLS + MC nulls) | 10 | **OK** |
| 03_generalizability_and_impact | 13 | **OK** |
| 04_experimental_design_toolkit | 12 | **OK** (after path fix) |
| 05_target_id_showcase | 8 | **OK** (after path fix) |
| 06_reinforcement_analyses | 8 | **OK** (after CWD fix) |
| 07_cross_celltype_transfer | 7 | **OK** |
| 08_deg_weighted_evaluation | 7 | **OK** |
| 09_causal_validation_dossier | 10 | **OK** |
| bring_your_own_target | 6 | **OK** |

No notebook loads the full 16 GB Tier-2 h5ad: nb01/nb02 take only `h5py` row-subset slices (≤280 MB) and
nb02's effect-matrix build short-circuits to the committed `cache/E_*.npz` bundles. All 30 result/cache
files the notebooks reference are present.

## Batch scripts — status

| Script | Status | Note |
|---|---|---|
| `run_iv_compliance.py` | **runs; reproduces** | fixed; `--invariance-only` reproduces its CSV bit-identical |
| `run_atlas.py`, `run_nulls.py` | **runnable; skip-existing verified** | read the present `inputs.npz`; heavy (~44 s per full-E solve, tens of minutes for a cold run). Both were executed: with the 24 committed `point_*.json`/`cell_*.json` present they correctly skip every cell and exit 0, leaving the tree untouched — a full cold recompute was not run. |
| `run_bootstrap.py` | **runnable (inner loop exercised)** | reads the present `inputs.npz`; heavy (~44 s point estimate + ~14 s per subsample refit). Its `signed_reachability` + NNLS subsample-refit inner loop was reproduced directly; the full 20-draw script was not run to completion. |
| `a2_conditional_reachability_scaffold.py` | **runs; SystemExit by design** | fixed; requires raw single-cell counts (not in repo) and stops with a documented message |
| `run_a1_sensitivity.py` | **provenance-only** | requires env vars `A1_INPUTS` + `A1_SE`; `A1_SE` → `atlas_lfcSE.npz`, a ~283 MB measurement-error matrix that is gitignored/absent. Its committed `a1_sensitivity_radius.csv` is frozen. |
| `build_effect_matrices.py` | **provenance-only** | requires `k562_essential.h5ad` / `rpe1_essential.h5ad` (gitignored CZI downloads, absent). Its output `cross_celltype_effects.npz` is already committed. |

The two provenance-only scripts and the a2 scaffold are the documented "heavy inputs are gitignored" cases —
their frozen outputs are consumed downstream, and the downstream notebooks that reproduce those numbers all pass.

## Tree integrity

Every script and notebook was run in place, so this repo's committed outputs served as the reproduction oracle.
After all runs, no regenerated result or figure file was left modified: everything the scripts and notebooks
wrote either matched the committed versions byte-for-byte or went only to gitignored locations
(`notebooks/figures/`, `notebooks/cache/E_*.npz`). The only tracked-tree changes are the five staged fixes above.
The one new untracked file is this report itself (`REPRODUCIBILITY_CHECK.md`, written to the repo root); it is not
a stray regenerated result, and the user can track, keep, or delete it.
