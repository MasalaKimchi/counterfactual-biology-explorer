# Notebooks

Exploratory, day-by-day. Keep heavy exploration here; promote stable logic into `src/`.

Suggested sequence (mirrors `ROADMAP.md`):

1. `01_data_access_qc.ipynb` — fetch DE_stats, subset high-confidence perturbations,
   sanity-check a known regulator (Day 1).
2. `02_target_states.ipynb` — build & visualize Th1/Th2 and aging target vectors (Day 2).
3. `03_counterfactual.ipynb` — run greedy/OMP/lasso; plot k-vs-alignment curves (Day 3).
4. `04_confidence_benchmark.ipynb` — reproducibility, stability selection, held-out
   donor, random null, linear baseline (Day 4).
5. `05_evidence_pathways.ipynb` — PubMed/Open Targets + enrichment; compare to the
   arrayed CRISPRi validation table (Day 5).

Kept out of version control by default (see `.gitignore` if you add outputs).
