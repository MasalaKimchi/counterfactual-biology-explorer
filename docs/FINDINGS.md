# Updated Findings

This narrative is controlled by [`results/findings.json`](../results/findings.json); the
artifact identities are controlled by [`results/manifest.json`](../results/manifest.json).

## Claim boundary

The project tests directional support under a non-negative linear-combination model. The
primary dictionary contains donor-collapsed primary-CD4 CRISPRi differential-expression
profiles from the source study's post-expansion `Rest` condition. The target is a
source-study-reused, selectively constructed cross-sectional population contrast oriented
toward the reported Th1 centroid. It is not independent, an observed trajectory, or a
measurement from established polarized Th2 cells.

## 1. Source-bound directional alignment is variable across gene splits

With the gene universe ordered lexicographically, split IDs hashed, and NumPy
`default_rng(seed)` frozen for seeds 0–11, held-out cosine is **0.444 ± 0.018** (mean ±
SD; range **0.417–0.473**) over 12 half-gene splits. This is split variability, not donor
uncertainty. Correlated genes make these splits neither independent replicates nor a
substitute for module holdout.

The earlier **0.446 ± 0.010** table is retired: its script consumed an unhashed, deleted
`analysis_cache/atlas_work/inputs.npz`, and that object's gene ordering was not preserved.
The separately archived fixed split (**0.448154**) agrees with the new source-bound seed-0
split within `3e-10`, confirming source alignment but not reconstructing the old multisplit
protocol. Its 60 target shuffles were not regenerated and **p = 1/61 is no longer a current
finding**.

Source: [`source_reconstruction.json`](../results/source_reconstruction.json).

## 2. Target construction materially narrows the estimand

The target table has 25,672 genes in the source union, 11,616 shared by Ota and
Höllbacher, 7,960 with concordant effect signs, and 6,188 remaining after screen
intersection. Across the 8,950 shared, screen-measured genes used for leakage-safe source
transfer, between-source cosine is **0.791** for log fold change and **0.698** for z-score.
These sources differ in cohort, assay, and preprocessing; agreement is not independence
from the source study.

## 3. Cross-source transfer is directional, not magnitude accurate

The cone is fit to one target source and evaluated against the other without selecting
genes by held-out-source sign. On log fold change, cone cosine exceeds the better of the
mean-response ray and best single atom in every one of six correlated random-gene splits:

| Direction | Mean cone cosine | Mean cosine gain | Mean cone nRMSE | Best-single nRMSE |
|---|---:|---:|---:|---:|
| Ota → Höllbacher | 0.255 | +0.087 | 1.185 | 1.003 |
| Höllbacher → Ota | 0.290 | +0.090 | 1.020 | 0.981 |

The directional gain is consistent across these splits, but magnitude error is worse than
the best single and sometimes worse than zero prediction. Baselines remain incomplete
(no nested PCA, ridge, unconstrained least squares, or matched random cone), and selecting
the better baseline on the test coordinates is disclosed. This is a construction-sensitivity
diagnostic, not predictive superiority, donor generalization, or biological validation.
Z-score results are sensitivity-only.

## 4. Donor-pair transfer retains weak direction but fails magnitude

The released donor H5MU contains six two-donor modalities. All six share 1,584 complete
`Rest` perturbation atoms, and 8,949 genes are common to both target sources and the
donor object. For each of three complementary 2-vs-2 donor partitions, the runner fits
NNLS weights on one donor-pair dictionary, one target source, and one random half of the
genes. It applies those weights unchanged to the complementary donor pair, opposite
target source, and held-out genes. The common-ray and best-single identity/scalar are
also selected only from the training side and frozen before testing.

The two mixed-run partitions contribute 24 run-balanced correlated challenges across
both donor directions, both target-source directions, and three fixed gene splits. Median
cone cosine is **0.031**, a **+0.032** gain over the training-selected best single; 18/24
challenges have positive cosine gain. Magnitude fails: median normalized RMSE is **1.153**
for the cone versus **1.018** for the best single, with improvement in only 1/24
challenges and none over the zero/common-response baselines. The 12 fully run-confounded
donors-1–2 versus donors-3–4 challenges are reported separately and not headlined.

This is a published-pipeline sensitivity, not leakage-free donor generalization. The
modalities are two-donor summaries, released presence is conditioned on DE eligibility,
gene splits are correlated, and four fixed donors cannot support donor-population
inference at 0.05. The result therefore narrows the utility claim: some directional
structure transfers, but predictive magnitude does not.

Source: [`donor_pair_transfer.json`](../results/donor_pair_transfer.json).

## 5. Independent Arce transfer is modest and context dependent

The compact Arce benchmark compares the frozen negative Zhu IL2RA log-fold-change score
with an independent CRISPR-KO CD25/IL2RA screen. Selection uses only four-guide coverage
in every Arce context and source-side Rest admission, leaving 480 targets. No Arce outcome
is used for selection or choosing the sign; orientation is fixed from assay semantics.

| Context | Signed Spearman | Signed Kendall | Direction agreement | Magnitude top-25 overlap |
|---|---:|---:|---:|---:|
| Resting Teff | 0.148 | 0.102 | 0.538 | 8 |
| Stimulated Teff | 0.084 | 0.057 | 0.521 | 9 |
| Resting Treg | 0.088 | 0.061 | 0.535 | 5 |

Resting-Teff rank association is strongest; the other contexts are weaker. Top-k overlap
ranks absolute magnitudes, while global Spearman/Kendall retain the prespecified sign.
Permutations exchange target labels, not donors. All 18 context-by-metric permutation
p-values are unadjusted exploratory diagnostics across correlated tests and support no
FWER- or FDR-controlled inference. The result supports only modest cross-study,
cross-modality ranking alignment, not state reachability, donor generalization, or causal
treatment validation.

Sources: [`arce_external_validation_meta.json`](../results/evidence/arce_external_validation_meta.json)
and [`arce_il2ra_context_predictions.csv`](../results/evidence/arce_il2ra_context_predictions.csv).

## 6. Arce donor/guide strata expose both robustness and heterogeneity

The same hash-bound archive contains 100,087 S14 singlet cells across two donors, two
guides for each member of the authors' preselected 28-regulator panel, nine Non-Targeting guides, and four Teff/Treg ×
resting/stimulated contexts. All 520 donor×guide×context strata are present. The runner
first takes each guide's median supplied `activation.score`, then subtracts the median of
the nine Non-Targeting guide medians within donor and context; cells are never counted as
biological replicates.

All 116 published S8 target×context means and medians reproduce from S14 with maximum
absolute error **5.7×10⁻¹³**. Across 28 targets, donor A-versus-B rank concordance is
**0.73–0.93** across contexts, but only **50–64%** of targets have the same nonzero effect
sign across both guides and both donors. One retained stratum has 8 cells; every other
stratum has at least 20. These are descriptive robustness diagnostics. S8 is derived from
S14, its pooled-cell p-values are unused, and two donors cannot estimate population
uncertainty. The score's gene set, formula, normalization, and independence are not frozen,
so it is not interpreted as CD25 protein, functional activation, Th1/Th2 identity, or
state conversion. Panel selection used prior screen/state-specific evidence, so these
metrics are conditional within-panel concordance, not genome-wide generality or independent
validation. The endpoint-mismatched Zhu IL2RA-to-`activation.score` correlation is
intentionally not reported.

Source: [`arce_activation_guide_effects.csv`](../results/evidence/arce_activation_guide_effects.csv).

## 7. Source-selected arrayed follow-up is transcriptomically specific and cytokine-consistent

The source authors arrayed nine screen-selected perturbations and released bulk RNA-seq
plus IL-10/IL-21 flow percentages for six additional donor labels (Donor5–Donor10). The
runner binds both tables to author commit `848d62f`, verifies every input byte, uses the
`Stim8hr` screen profiles, intersects 8,976 Ensembl genes, and masks all nine panel target
genes from every profile before transcriptome scoring (8,967 coordinates per profile).

Every arrayed bulk profile retrieves its matching screen perturbation first among all nine
(top-1 **9/9**, MRR **1.0**). Median raw cosine is **0.549**, versus a median **0.283**
gain over the common source-panel response. After subtracting each assay's across-panel
gene mean—a Systema-style specificity sensitivity—median cosine is **0.580** and retrieval
remains 9/9. Median normalized RMSE is **1.012** raw and **1.052** centered, so magnitude
calibration is not supported despite strong identity and directional structure.

Flow effects are `log2(percent positive / donor mean NTC percent positive)` within donor;
targets then give every observed donor equal weight. Donor coverage is uneven (three to
six). Across all nine targets, Spearman correlation between donor-median flow effect and
screen/bulk RNA cytokine log fold change is **0.817/0.717** for IL-10 and **0.833/0.850**
for IL-21. Synchronized exhaustive target-label diagnostics over these four correlations
are conditional on the selected panel; maximum-statistic tail fractions are 0.020, 0.062,
0.015, and 0.011, respectively. Unequal donor count/composition means target
exchangeability is not established. These fractions are not inferential p-values or
multiplicity-adjusted inference, do not undo upstream selection, and are not
donor-population inference.

This supports same-study cross-platform replication and direct cytokine consistency for
the measured panel. It does not establish held-out target discovery, guide robustness,
population generalization, durable state conversion, function, fitness, chromatin
remodeling, or intervention efficacy.

Sources: [`zhu_arrayed_validation_meta.json`](../results/evidence/zhu_arrayed_validation_meta.json),
[`zhu_arrayed_profile_metrics.csv`](../results/evidence/zhu_arrayed_profile_metrics.csv), and
[`zhu_arrayed_flow_effects.csv`](../results/evidence/zhu_arrayed_flow_effects.csv).

## 8. The systemic harness exposes nonspecific success modes

The data-free harness independently checks the NNLS solver and fails closed on axis and
provenance corruption. Its structured-specificity scenario produces high raw cosine from
a nuisance common response, shows random-gene optimism relative to module holdout, and
shows large sign-selection inflation while retaining power for a sparse true alternative.
The fast maxT check has 24/500 false family rejections; its exact one-sided 95% upper bound
is 0.067 under the declared 0.075 gate. This certifies software and synthetic statistical
contracts only.

## What remains unknown

Claim-bearing biological evidence still requires leakage-safe donor- and guide-held-out
effects rebuilt from pseudobulk, module/pathway holdouts, stronger nested baselines,
calibrated structured nulls, measured
in-domain combinations, paired CRISPRi/CRISPRa, chromatin and protein/function endpoints,
durability and fitness, and prospective testing from an established polarized starting
state. Until then the appropriate output is a model-relative direction test and the next
measurement—not a reachability verdict or intervention recommendation.
