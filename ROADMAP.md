# Roadmap

The detailed gates, datasets, expert handoffs, and prospective assay are in
[`docs/SCIENTIFIC_VALIDATION_PLAN.md`](docs/SCIENTIFIC_VALIDATION_PLAN.md). This page is
the concise execution order.

## Completed in the current revision

- target-lineage correction: 25,672 union → 11,616 shared → 7,960 sign-concordant →
  6,188 analyzed;
- target-scale normalization and fail-closed KKT/separator certification;
- independent active-set oracle and degenerate-cone stress tests;
- labelled axis/provenance contract and injected input faults;
- grouped split and maxT utilities;
- deterministic data-free harness report in pull-request CI.

## Next frozen release: source-data reconstruction

1. Ingest Zhu donor-pair, guide-level, pseudobulk, and target-construction files; resolve
   the currently unavailable arrayed Th1/Th2 validation-table route.
2. Record immutable dataset cards and hashes.
3. Regenerate the real-data evidence from source summaries.
4. Make log fold change primary and z-score geometry sensitivity-only.
5. Reproduce Zhu's own `pert2state_model` as the same-question baseline.

Exit: a clean external-data run regenerates all claim-bearing tables.

## Next inferential release: transfer and calibration

1. Ota→Höllbacher and Höllbacher→Ota source transfer without held-out-source sign
   selection.
2. Complementary donor-pair and guide-held-out transfer.
3. Module/pathway/context/run splits and covariance-preserving nulls.
4. Nested baselines, multiple metrics, maxT control, interval coverage, power, and
   coefficient-equivalence analysis.

Exit: target-specific performance clears the best frozen baseline with calibrated
source/donor/guide-held-out uncertainty.

## Biological validation release

1. Resolve and then score the VCP-documented arrayed polarized-culture RNA and flow
   validation without refitting.
2. Quantify additivity error for one measured four-donor primary-T-cell triple.
3. Test CRISPRa orientation and RNA/protein context transfer.

Exit: frozen predictions agree with measured molecular outcomes across donors and
modalities, within each dataset's stated claim ceiling.

## Prospective biology

Preregister established-memory-Th2 and naïve-derived-Th2 experiments with independent
donors, matched CRISPRi/CRISPRa, singles and selected combinations, RNA, protein/cytokine,
chromatin, fitness, lineage, durability, and functional readouts.

Exit: only coordinated, durable, within-lineage molecular and functional movement permits
cautious partial-reprogramming language.
