# Fetching the data

**Two tiers, and the tiering was corrected — read this.** The local no-login
supplementary CSVs (this folder) are the **Tier-1 warm-up + zero-dependency fallback**,
and they need no CZI / Synapse / Wiley account. But the **graded core is Tier-2**: the
gene-level `GWCD4i.DE_stats.h5ad` effect matrix, where the novel reachability-cone method
actually runs. It is **required for the headline result** (not "optional") — fetch it
early via a free CZI login (`vcp login` + `vcp data download`, or `fetch_de_stats.sh`).
The raw ~22M-cell dataset is **not** needed (hundreds of GB — only the derived matrix).
Tier 1 stays fully no-auth so there is still a defensible submission if the download
fails. See the ROADMAP "Two tiers" section for the full rationale.

> Note (`.gitignore`): `data/*.csv` is currently ignored, so the 7 CSVs are present
> locally but **not tracked** by git. To ship a self-contained repo, force-add them:
> `git add -f data/*.suppl_table.csv`.

## What is already local (checked in)

The following **supplementary CSV tables (~36 MB total)** are committed to this
repo and are enough to run every analysis in the "CSV-only" tier of `ROADMAP.md`:

| File | Rows | What it is |
|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | Per **perturbation × culture-condition** summary: DE-gene counts, on-target knockdown effect/significance, off-target flag, cross-donor & cross-guide reproducibility. **Not** the full gene-level effect matrix. |
| `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` | 37,288 | Th2-vs-Th1 target signature (per-gene logFC/zscore) from **two** source contrasts (Ota 2021, Höllbacher 2021). |
| `CD4T_aging_signature_DE_results_full.suppl_table.csv` | 10,000 | CD4+ T-cell aging signature (per-gene logFC/zscore), Yaza 2022 discovery contrast. |
| `guide_kd_efficiency.suppl_table.csv` | 73,765 | Per-guide × condition knockdown QC (guide vs NTC expression, t-stat, `signif_knockdown`). |
| `sgrna_library_metadata.suppl_table.csv` | 31,110 | Per-guide design + off-target annotation (TSS distance, nearby/non-target genes, alternate alignments). |
| `cluster_autoimmune_enrichment_results.suppl_table.csv` | 5,236 | Perturbation-cluster × **autoimmune-disease** GWAS-gene enrichment (odds ratio, FDR, intersecting genes) across 17 diseases and 4 gene sets. |
| `sample_metadata.suppl_table.csv` | 12 | Sample sheet: 4 donors × 3 conditions, donor demographics (age, sex, ethnicity). |

**Not local (fetch this — it is the graded core):** the large `GWCD4i.DE_stats.h5ad`
gene×perturbation *effect matrix* (and `GWCD4i.pseudobulk_merged.h5ad`). This is the
input to the **reachability-cone + reconstruction solver** — the headline method — so
it is required for the primary result, not a "nice to have" (see `ROADMAP.md`
§"Two tiers"). Fetch via Option A below early in the week.

## Option A — CZI Virtual Cells Platform CLI (recommended)

1. Register (free) at https://virtualcellmodels.cziscience.com/
2. Install the CLI: see https://chanzuckerberg.github.io/vcp-cli/usage/data.html
3. Search & download:

```bash
vcp data search "Primary Human CD4+ T Cell Perturb-seq" --exact
# then download the artifact(s) you need, e.g. GWCD4i.DE_stats.h5ad
```

Core artifact for this project: **`GWCD4i.DE_stats.h5ad`**
(33,983 perturbation×condition rows × 10,282 genes; layers: `log_fc`, `zscore`,
`p_value`, `adj_p_value`, `baseMean`, `lfcSE`).

## Option B — No-auth fallback (supplementary tables on GitHub)

Open, MIT-licensed supplementary tables (enough to prototype the pipeline):
https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/master/metadata/suppl_tables

The files listed in the table above are the ones we already pulled from here. If
you re-clone the analysis repo, also look for an **arrayed CRISPRi validation
table** (e.g. a `Th1Th2_validation*` file) if the authors publish one — we do
**not** currently have it locally, so our polarization ground-truth check instead
relies on internal reproducibility (`guide_kd_efficiency` + `DE_stats`
cross-donor/cross-guide columns) and orthogonal literature/Open Targets evidence.

## Memory tips (CPU-only laptop)

- Load a single layer (e.g. `zscore`) as `float32`, not the full AnnData with all layers.
- Subset to perturbations with `keep_test_genes == True` and the top ~2,000 HVGs.
- Cache the reduced matrix to `.npz` so you never reload the h5ad during iteration.
