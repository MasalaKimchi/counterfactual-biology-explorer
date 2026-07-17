# External data

Scientific source data stay local and are ignored by Git. The repository ships a compact
frozen evidence bundle and a data-free validation harness, not the large source matrices.

## Primary source: Zhu et al.

Genome-scale CRISPRi Perturb-seq in primary human CD4 cells from four donors and three
assay conditions:

- preprint: https://doi.org/10.64898/2025.12.23.696273
- GEO/SRA: GSE314342 / SRP643211
- source code: https://github.com/emdann/GWT_perturbseq_analysis_2025
- public data prefix: `s3://genome-scale-tcell-perturb-seq/marson2025_data/`
- VCP release: [v1.0 dataset page](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq?access_dataset=true)
- licensing: VCP labels the dataset MIT and applies its Acceptable Use Policy; the
  analysis code is separately MIT; mirror terms must be recorded independently

List the public objects without credentials:

```bash
aws s3 ls --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/
```

Highest-priority objects:

| Object | Bytes | Use |
|---|---:|---|
| `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` | 6,155,771 | Exact target-construction input; fetched from the author repository |
| `IL10IL21bulkRNAseq_DESeq2_results.csv` | 13,952,871 | Arrayed bulk-RNA follow-up for nine source-selected perturbations |
| `IL10_IL21_arrayed_validation.csv` | 2,200 | IL-10/IL-21 flow percentages for Donor5–Donor10 and matched NTC controls |
| `GWCD4i.DE_stats.by_donors.h5mu` | 16,866,278,447 | Disjoint donor-pair transfer |
| `GWCD4i.DE_stats.by_guide.h5mu` | 29,424,424,894 | Guide-held-out replication |
| `GWCD4i.pseudobulk_merged.h5ad` | 44,566,657,140 | Controlled re-aggregation and alternate DE models |
| `GWCD4i.DE_stats.h5ad` | 16,786,240,107 | 33,983 perturbation-condition profiles × 10,282 genes; log FC, z-score, p, adjusted p, base mean, and LFC SE |

Frozen source identities used by `configs/source_reconstruction.json` and
`configs/zhu_arrayed_validation.json`:

| Object | SHA-256 |
|---|---|
| `GWCD4i.DE_stats.h5ad` | `c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62` |
| `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` | `c47d2df21414ca85e7aa255f4148904eec700fbcd9debc2f734ec97049698444` |
| `IL10IL21bulkRNAseq_DESeq2_results.csv` | `c20418a9285b10104dbae362b825971f86f97425800a92269e4433ce780e666d` |
| `IL10_IL21_arrayed_validation.csv` | `f60cdda392d6f29d10a539727ff7324b04d17e35c0512c889b733e00380b83dc` |
| `GWCD4i.DE_stats.by_donors.h5mu` | `2ee3cf90925600eb044619021da2bdd47d661f306a204586652256facf17af64` |

The target CSV is pinned to author commit
[`848d62f`](https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/848d62fc2b7027f7218d6fc5f5b0c37255dc94af),
not the mutable default branch.

The VCP schema also documents `Th1Th2_validation_summary.suppl_table.csv`, but no working
public S3 or author-repository object was verified on 2026-07-17. Treat retrieval as
unresolved; do not present this table as an executable benchmark until a source route is
confirmed.

Download one allow-listed object with:

```bash
./data/fetch_de_stats.sh Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv
./data/fetch_de_stats.sh IL10IL21bulkRNAseq_DESeq2_results.csv
./data/fetch_de_stats.sh IL10_IL21_arrayed_validation.csv
./data/fetch_de_stats.sh GWCD4i.DE_stats.by_donors.h5mu
```

The default remains `GWCD4i.DE_stats.h5ad`. HTTPS downloads resume from a `.part` path,
verify registered byte length, and then rename atomically. Record SHA-256 and retrieval
date in the dataset card before analysis.

With both frozen files present, run the source-bound reconstruction separately from the
small CI suite:

```bash
python -m pip install -r requirements-external.txt
python scripts/run_source_reconstruction.py --profile
python scripts/run_source_reconstruction.py --check results/source_reconstruction.json
```

This reconstructs target lineage and aggregate-screen source-transfer diagnostics. It
does not substitute for donor- or guide-level objects.

With the two compact follow-up tables present, the source-selected arrayed benchmark adds
8,967 transcript coordinates per perturbation after masking all nine panel target genes,
plus donor-normalized IL-10 and
IL-21 flow readouts:

```bash
python scripts/run_zhu_arrayed_validation.py --check
```

The nine targets were chosen upstream using the same source study. Donor coverage ranges
from three to six, and target-label permutations are conditional on that panel. This is
same-study cross-platform follow-up in six additional donor labels—not held-out discovery,
donor-population inference, or state conversion.

The donor-pair H5MU stage is implemented:

```bash
python scripts/run_donor_pair_transfer.py --check
```

It fits on one two-donor modality and target source, then freezes coefficients and
training-selected baselines for the complementary donor pair, opposite target source,
and held-out genes. Across 24 run-balanced correlated challenges, median cosine gain over
the frozen best single is +0.032, while normalized RMSE is worse (1.153 versus 1.018).
This is a published-pipeline sensitivity, not predictive utility. The next stages are
reciprocal-guide H5MU and structural-QC re-estimation from pseudobulk counts. Released
H5MU modalities omit DE-ineligible targets, so presence itself can select on effectiveness.
Only pseudobulk can rebuild a leakage-safe universe without `keep_effective_guides`,
held-out DE, sign, or correlation fields. Four donors support fixed-cohort robustness,
not donor-population significance.

## Target sources

- Ota et al., *Cell* 2021, DOI 10.1016/j.cell.2021.03.056, NBDC E-GEAD-397.
- Höllbacher et al., *ImmunoHorizons* 2020, DOI
  10.4049/immunohorizons.2000037, GEO GSE149090.

The current target uses the Zhu supplementary polarization table. Its registered merged
analysis has 6,188 genes after requiring both sources, concordant signs, and screen
measurement. Source-transfer validation must return to the independent Ota and
Höllbacher inputs and must not select coordinates using the held-out source.

## Compact external benchmarks

These archives are small enough to exercise the external-data harness before downloading
the 17–45 GB donor/guide/pseudobulk objects:

| Dataset | Exact archive | Bytes | Terms and allowed use |
|---|---|---:|---|
| Arce Perturb-CITE, Zenodo 13924126 / GSE278572 | [`data_tables.zip`](https://zenodo.org/api/records/13924126/files/data_tables.zip/content), MD5 `886ed0fea0b9dc0625355c2e4928077c`, SHA-256 `dc9e2efb04d24f1a6d4b8db6a8b1d5cd01c935777c3740088be339de5b5062b4` | 57,967,623 | **Implemented**: S1 CRISPRi-transcript → CRISPR-KO CD25 transfer and S14/S8 donor/guide supplied-score robustness; not state or donor-population validation |
| Schmidt screens, Zenodo 5784651 / GSE174292 | [`Genome-wide-screens.zip`](https://zenodo.org/api/records/5784651/files/Genome-wide-screens.zip/content), upstream MD5 `e0392eb7b2512720bb8cbf705ce9854f` | 26,152,593 | CC-BY-4.0; CRISPRa/CRISPRi orientation and IFNG/IL2 transfer with only two donors |

The Arce runner verifies the archive plus S1/S8/S14 workbook bytes, requires four guides in every
screen context, freezes Zhu `Rest` admission before reading Arce outcomes, and checks the
current deterministic evidence:

```bash
python scripts/run_arce_external_validation.py --check
```

The retained S1 benchmark has 1,347 four-guide genes, 1,259 present in the Zhu `Rest`
dictionary, and 480 source-admitted analysis targets. S1 is aggregate and cannot supply
donor uncertainty. Arce S9 is significant-only (absence is censored, never zero), and S4
contains unlabelled technical replicates that must be aggregated within biological key.
Schmidt's cell-level marker tables are likewise not donor-level inference.

S14 contributes 100,087 singlet cells, 520 complete guide×donor×context strata, and four
contexts. S8's 116 pooled summaries are exactly reproduced but are not independent data;
its pooled-cell tests are excluded. The supplied `activation.score` lacks a frozen local
formula/gene set and is used only for descriptive within-object robustness.

## Next dataset priorities

- Schmidt primary-T-cell CRISPRa/CRISPRi: GSE174292 and subseries GSE174255,
  GSE190604, GSE190846.
- Goudy primary-T-cell CRISPRoff/Cas9: GSE306915; descriptive only until sample mapping
  and reuse terms are resolved.
- Zhu donor/guide objects above, after compact benchmark calibration is stable.

See [`docs/SCIENTIFIC_VALIDATION_PLAN.md`](../docs/SCIENTIFIC_VALIDATION_PLAN.md) for
the exact evaluation and claim ceiling for each resource.

## Data policy

- Never force-add H5AD, H5MU, CSV, NPZ, or raw-count payloads under `data/`.
- Every analysis input needs an accession, exact object key, version/retrieval date,
  byte length, SHA-256, license/terms status, gene namespace, orientation, units,
  donor/guide/context/batch fields, and missingness profile.
- A missing donor or guide estimate remains missing; it is never silently replaced by a
  collapsed average.
- Claim-bearing tables must be regenerated by the external-data runner, not copied from
  an interactive notebook.
