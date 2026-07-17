# Method

## Question and geometry

Given measured perturbation profiles and a target transcriptional direction in a shared
gene space, how closely can a non-negative linear combination of those profiles align with
the target? The method tests one dictionary, target, and metric; it is narrower than
perturbation-response prediction.

Let `E` have shape perturbations × genes, `A = E.T`, and `d` be the target. Coefficients
solve

```text
minimize_w>=0  1/2 ||sqrt(W) (d - A w)||²
```

for diagonal gene metric `W`. Zero-weight coordinates are excluded and the largest
retained weight is normalized to one. The fitted direction is `A w`; the residual is
`rho = d - A w`; and a metric-relative separator is proportional to `W rho`. The separator
describes this measured cone only. It does not prove biological impossibility or identify
a gene that must be activated.

## Numerical contract

The core reports coefficients, fitted direction, weighted cosine, residual fraction,
relative objective, atom-scale-invariant KKT diagnostics, numerical geometry status, and
a separator only when residual energy clears a declared tolerance. KKT, polarity, and
orthogonality certify at `1e-8` or the projection fails closed. Exact-zero, non-finite,
misaligned, all-zero-weight, and unrepresentable-scale inputs fail before a result is
returned. No biological verdict threshold exists.

## Source-bound held-out challenge

`held_out_alignment` fits on one gene subset and scores the frozen coefficients on a
disjoint subset. Source-bound case-study genes are ordered lexicographically by symbol;
NumPy `default_rng(seed)` generates half-gene splits for seeds 0–11; and each split's
fit/score gene identifiers are SHA-256-bound in the report. Random-gene splits treat
correlated genes as exchangeable, so module/pathway holdout remains required.

The former 60-shuffle target diagnostic was retired because its deleted intermediate did
not preserve gene order. No claim-bearing p-value is attached to the current case-study
cosine. The systemic harness separately tests structured common-response, module,
sign-selection, and multiplicity failure modes.

## Case-study inputs and target construction

- Dictionary: donor-collapsed primary-human-CD4 CRISPRi differential-expression profiles
  in the source study's post-engineering/post-expansion `Rest` condition.
- Target: source-study-reused, selectively constructed cross-sectional population
  contrast oriented toward the reported Th1 centroid; not independent or a trajectory.
- Unit: perturbation-condition profile, not donor effect.
- Primary source-transfer scale: log fold change; z-score is sensitivity-only.

The registered merged target retains genes present in both Ota and Höllbacher target
contrasts, requires concordant effect signs, averages their z-scores, reverses the reported
Th2-versus-Th1 orientation, and intersects with screen genes. Counts are 25,672 source
union, 11,616 shared, 7,960 sign-concordant, and 6,188 registered coordinates. Every
availability, sign, and screen filter is a selection step.

The target sources are Ota et al., *Cell* 2021, DOI
`10.1016/j.cell.2021.03.056` (E-GEAD-397), and Höllbacher et al.,
*ImmunoHorizons* 2020, DOI `10.4049/immunohorizons.2000037` (GSE149090). The
dictionary is Zhu et al. 2025, DOI `10.64898/2025.12.23.696273` (SRP643211 / GSE314342).
The exact target CSV and 16.8 GB H5AD bytes are verified before source-bound analysis.

## Cross-source transfer

Ota and Höllbacher are evaluated on the 8,950 genes present in both sources and the
screen. This universe uses identifiers/presence only; held-out-source signs and magnitudes
do not select coordinates. For each direction and split, coefficients are fit to one
source on fit genes and scored against the other on score genes. Mean-response and
best-single baselines are trained on the fit source. Improvement is reported against the
larger test cosine of those comparators, a disclosed test-selected descriptive baseline.
Cosine, normalized RMSE, norm ratio, sign agreement, support size, and split hashes are
retained.

## Independent Arce benchmark

Arce S1 is read directly from its hash-verified Zenodo archive. Genes require exactly
four guides in all three contexts. Intersecting with Zhu `Rest` availability and
source-side `ontarget_significant` admission—without reading Arce outcomes—leaves 480
targets. S1 positive log fold change denotes enrichment of knockouts lowering CD25/IL2RA,
so the regulator score is frozen as negative Zhu IL2RA transcript log fold change.

Per context, the runner reports signed Spearman/Kendall after the prespecified orientation,
direction agreement, absolute-magnitude top-25/50/100 overlap, and target-label permutation
diagnostics. Directional agreement uses a null center fixed analytically from the two sign
margins. Targets, not cells or donors, are the exchangeable units. All 18 permutation
p-values are unadjusted exploratory diagnostics across correlated contexts/metrics and
provide no multiplicity-controlled inference or donor uncertainty. The claim ceiling is
cross-study/cross-modality ranking alignment, not whole-state reachability or causal
validation.

S14 is streamed from the same verified archive with its complete 32-column schema checked,
while retaining only registered QC and score fields. Every cell must have one sgRNA,
`has_sgrna=true`, a unique cell ID, finite supplied `activation.score`, and global HTO
Singlet status. For target *t*, guide *g*, donor *d*, and context *c*, the guide contrast is
the guide-cell median minus the median of the nine Non-Targeting guide-cell medians within
the same donor/context. The two guide contrasts are equally weighted for a target/donor
summary. All strata are retained and their cell counts disclosed.

The runner exactly reproduces S8 pooled target/context means and medians as an archival
provenance gate. It does not use S8 pooled-cell p-values, bootstrap cells as donors, or
emit donor-population inference from two donors. Guide-pair and donor rank/sign concordance
are descriptive. Because the local tables do not freeze the supplied score's genes,
formula, normalization, or independence, it is not equated with functional activation or
with the one-gene Zhu IL2RA predictor.
The 28 regulators were preselected by the source authors using prior screen/state-specific
evidence; target-rank and sign concordance are therefore conditional panel diagnostics,
not genome-wide estimates or independent validation.

## Source-selected arrayed RNA/protein follow-up

The Zhu author repository contributes arrayed bulk-RNA differential-expression summaries
for nine source-selected perturbations and IL-10/IL-21 flow percentages in follow-up donor
labels Donor5–Donor10. The tables and aggregate H5AD are full-file SHA-256-bound. The
screen side uses `Stim8hr` log fold change without consulting aggregate admission,
guide-correlation, donor-correlation, significance, or follow-up outcomes.

Bulk and screen profiles intersect on 8,976 Ensembl IDs. All nine panel target genes are
masked from every profile before computing cosine, Pearson, Spearman, normalized RMSE, and
retrieval rank among all nine source profiles. Raw profiles are compared with the common
source-panel response. A specificity sensitivity subtracts the across-nine-perturbation
mean separately for every gene and assay before rescoring; this panel centering is not a
replacement biological estimand.

Within each donor and cytokine, the flow baseline is the arithmetic mean of every
available NTC-labeled measurement. A perturbation effect is
`log2(percent_positive / donor_NTC_mean)`. Donor measurements are never pooled as cells;
target summaries are the equal-donor median, with mean/range and donor count retained.
Coverage is unbalanced and no missing target/donor is imputed.

All `9! = 362,880` target-label assignments are enumerated. Retrieval top-1/MRR and four
primary Spearman associations (screen/bulk RNA × IL-10/IL-21 donor-median flow) are
recomputed synchronously; the correlation family records a one-sided maximum-statistic
tail fraction. Because donor count, composition, and NTC replication differ by target,
target exchangeability is not established. These exhaustive fractions are conditional
diagnostics, not inferential p-values or multiplicity-adjusted inference. They do not
replay upstream panel selection and provide neither held-out-discovery nor donor-population
inference. The claim ceiling is measured same-study cross-platform replication and
cytokine consistency.

## Staged donor/guide validation contract

The released donor-pair H5MU stage is implemented on `Rest` log fold change. Its six
modalities represent every two-donor combination among four donors. The analysis requires
the 1,584 perturbation atoms present in every modality and constructs an 8,949-gene target
universe from both target sources plus donor-H5MU gene symbols without using held-out
source signs. For each of three complementary 2-vs-2 partitions, both donor directions,
both target-source directions, and seeds 0–2, the runner fits NNLS on one donor pair,
source, and half-gene split. It applies identical coefficients to the complementary donor
pair, other source, and held-out genes.

The common-response ray and best-single atom/scalar are selected solely on the training
dictionary, target source, and fit genes; their identity and scale remain frozen on test.
The two mixed-run donor partitions yield 24 primary descriptive challenges. The donors
1–2 versus 3–4 partition yields 12 fully run-confounded sensitivities reported separately.
All challenges are correlated and no donor-population p-value or confidence interval is
emitted. The released modalities omit DE-ineligible targets, so complete-case presence is
potentially effectiveness-selected.

The next released-object stage is reciprocal `guide_1`/`guide_2` H5MU transfer under the
same frozen-coefficient and training-only-baseline rules.

A leakage-safe donor/guide claim therefore requires structural-QC re-estimation from the
pseudobulk count object. Universe construction must exclude `keep_effective_guides`,
aggregate `ontarget_significant`, held-out DE/sign/correlation fields, and any filter whose
outcome independence is not demonstrated. With four donors, all partitions remain
fixed-cohort descriptive challenges; exact donor-population significance at 0.05 is
mathematically unavailable.

## Removed legacy behavior

The maintained surface intentionally excludes categorical reachability labels, staged
`E`/`-E` modality proxies, activation certificates, greedy recipes, uncalibrated analytic
nulls, automatic experiment recommendations, and evidence tables derived from deleted,
unhashed intermediates. Git history preserves provenance without presenting those outputs
as current evidence.

## Reproducibility

`./reproduce.sh` runs exact numerical tests and the demo, checks the frozen systemic
harness, validates every canonical scalar against the source-bound and Arce reports,
rejects any failed report, and verifies artifact hashes. External runners regenerate their
reports from gitignored registered inputs; acquisition is documented in
[`data/README.md`](../data/README.md).
