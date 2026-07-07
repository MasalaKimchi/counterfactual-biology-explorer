# Fetching the data

The raw dataset is ~22M cells and is **not** needed. We work with precomputed
derived artifacts.

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

Useful files:
- `DE_stats.suppl_table.csv` — per-perturbation DE metadata + reproducibility columns
- `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` — polarization target
- `CD4T_aging_signature_DE_results_full.suppl_table.csv` — aging target
- `Th1Th2_validation_summary.suppl_table.csv` — arrayed CRISPRi ground truth
- `guide_kd_efficiency.suppl_table.csv` — knockdown QC
- `sgrna_library_metadata.suppl_table.csv` — guide/off-target metadata

## Memory tips (CPU-only laptop)

- Load a single layer (e.g. `zscore`) as `float32`, not the full AnnData with all layers.
- Subset to perturbations with `keep_test_genes == True` and the top ~2,000 HVGs.
- Cache the reduced matrix to `.npz` so you never reload the h5ad during iteration.
