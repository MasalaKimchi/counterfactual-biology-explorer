# Method

## Question

Given measured perturbation profiles and a target transcriptional direction in a shared
gene space, how closely can a non-negative linear combination of those profiles align with
the target?

The method is intentionally narrower than perturbation-response prediction. It tests one
dictionary and one target under one metric.

## Geometry

Let `E` have shape perturbations × genes and let `d` be the target vector. With
`A = E.T`, non-negative coefficients solve

```text
minimize_w>=0  1/2 ||sqrt(W) (d - A w)||²
```

where `W` is a diagonal gene metric. Zero-weight coordinates are excluded before solving;
all retained weights are positive. A global positive scale of `W` is immaterial, so the
implementation canonically sets its largest retained value to one and rejects an
unrepresentable dynamic range. The fitted direction is `d_hat = A w`, the raw residual
is `rho = d - d_hat`, and the metric-relative separator is proportional to `W rho` at
that canonical scale.

The separator is a statement about the fitted cone only. It does not prove biological
impossibility or identify a gene that must be activated.

## Numerical contract

The maintained core reports:

- non-negative coefficients and fitted direction;
- weighted cosine, residual fraction, and dimensionless relative objective;
- atom-scale-invariant KKT diagnostics;
- `inside_tolerance` or `outside_model_cone` geometry status;
- a separator only when residual energy clears a numerical threshold;
- dimensionless polarity, orthogonality, and separation diagnostics.

KKT, separator polarity, and separator orthogonality must each certify at `1e-8` or the
projection fails closed. The user-facing separator tolerance controls only whether a
numerically certified residual is treated as negligible; it cannot weaken certification.
The reported relative objective is `0.5 * ||sqrt(W) residual||² / ||sqrt(W) target||²`;
the raw dimensional objective is not returned because its square can underflow or overflow
even when the normalized geometry is representable.

Exact-zero, non-finite, misaligned, or all-zero-weight targets fail before solver dispatch.
The core contains no hard-coded biological verdict threshold.

## Held-out challenge

`held_out_alignment` fits coefficients on one gene subset and scores the frozen fit on a
disjoint subset. The case study uses 12 registered random half-gene splits. This catches
some in-sample flexibility but treats correlated genes as if partitions were exchangeable;
module and pathway holdouts remain required.

The historical target-shuffle diagnostic uses a plus-one empirical p-value with ties
counted as exceedances. Sixty shuffles are descriptive and cannot estimate an extreme
tail. Claim-bearing inference requires a structured null, multiplicity handling, and more
resamples.

## Case-study inputs

- Dictionary: donor-collapsed primary-human-CD4 CRISPRi differential-expression z-score
  profiles in the source study's post-engineering/post-expansion `Rest` condition.
- Target: a sign-concordant direction constructed from two external Th1-vs-Th2 contrasts.
- Unit of analysis: perturbation-condition profile, not an independent donor effect.
- Effect scale: z-score geometry; coefficients are not intervention dose.

The target construction is exact and fixed: pivot the Zhu supplementary polarization
table by gene and source contrast, retain genes with non-missing Ota and Höllbacher Wald
z-scores of the same sign, average those two z-scores, flip the Th2-vs-Th1 sign to obtain
the Th1-like direction, then intersect with genes measured in the CRISPRi screen. The
union of source tables contains 25,672 unique genes and 9,831 overlap the screen; those
are coverage counts, not the final estimand. The two sources share 11,616 genes, 7,960
are sign-concordant, and 6,188 remain in the registered merged analysis after screen
intersection. Thirty-eight of the 50 strongest target-table DE genes overlap the screen.
Every source-sharing, sign, and screen-intersection filter is a selection step, not
neutral preprocessing.

The two target inputs are Ota et al., *Cell* 2021,
doi:10.1016/j.cell.2021.03.056, NBDC E-GEAD-397, and Höllbacher et al.,
*ImmunoHorizons* 2020, doi:10.4049/immunohorizons.2000037, GEO GSE149090. The
Höllbacher contrast was re-estimated from raw counts using a donor- and PC1-adjusted
DESeq2 model. The primary dictionary is Zhu et al. 2025,
doi:10.64898/2025.12.23.696273, SRA SRP643211 / GEO GSE314342. The retained Norman
diagnostic uses Norman et al., *Science* 2019, doi:10.1126/science.aax4438, GEO
GSE133344. Source retrieval and transformations are frozen through the Zhu supplement
and its `4_polarization_signatures` analysis.

## Descriptive calibration anchor

The 0.463 dynamic-range fraction uses a synthetic reachable-by-construction reference,
not an experimental positive control. Forty targets are generated with seed 5;
each combines five effect atoms with independent coefficients sampled uniformly from
0.5–2.0, then adds independent Gaussian noise with per-gene standard deviation
`0.25 * ||target|| / sqrt(n_genes)`. The reported reference is the median held-out cosine.
These fractions describe this metric under that simulation and have no probability or
efficacy interpretation.

## Removed legacy behavior

The consolidated code intentionally removes:

- categorical “reachable” and “weakly reachable” labels;
- staged `E` then `-E` modality decomposition;
- gene-level activation-certificate claims;
- greedy panels described as minimal recipes;
- analytic nulls that could bypass empirical calibration;
- automatic experimental recommendations.

Historical implementations remain available through Git history, not the maintained API.

## Reproducibility

`./reproduce.sh` runs exact synthetic tests, executes a small demo, validates the findings
ledger against selected evidence tables, and verifies every canonical artifact hash.
Large source data stay outside Git; acquisition instructions are in [`data/README.md`](../data/README.md).
