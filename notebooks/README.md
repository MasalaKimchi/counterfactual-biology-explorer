# Notebooks

Exploratory, day-by-day. Keep heavy exploration here; promote stable logic into `src/`.

> **Day 0 first:** restore the validated pipeline before writing new notebooks —
> `git checkout 3f8db17 -- src app tests environment.yml requirements.txt data/fetch_de_stats.sh`,
> then `python -m src.data_loader --check` + `pytest -q`. Also fetch the Tier-2
> `GWCD4i.DE_stats.h5ad` matrix (the graded core; see `data/README.md`).

Suggested sequence (mirrors `ROADMAP.md`):

1. `01_data_load_qc.ipynb` — load the local Tier-1 CSVs, assert shapes, index guide
   QC, sanity-check a known regulator (Day 1). **Also load the Tier-2 `GWCD4i.DE_stats.h5ad`
   matrix — required for the reachability core, not optional.**
2. `02_target_states.ipynb` — build & visualize Th1/Th2 and aging target vectors,
   plus the Ota-vs-Höllbacher concordance core (Day 2).
3. `03_counterfactual.ipynb` — Tier-1 directional ranking as a warm-up, then the
   **Tier-2 reachability cone + spectrum** (NNLS + greedy/OMP/lasso, k-vs-alignment
   curves with shuffled-target and random-perturbation nulls) as the headline (Day 3).
4. `04_confidence_benchmark.ipynb` — reproducibility, off-target audit, stability
   selection, leave-one-donor-out (n=4), random null, linear baseline (Day 4).
5. `05_disease_evidence_pathways.ipynb` — autoimmune-enrichment linkage +
   PubMed/Open Targets/Consensus evidence + enrichment (Day 5).

Kept out of version control by default (see `.gitignore` if you add outputs).
