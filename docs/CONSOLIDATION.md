# Consolidation changelog (2026-07-20)

Documentation + directions consolidation for CombiCone. This records what merged
into what, what moved, the one canonical number that was corrected, and a
clearly-labelled set of **pre-existing / out-of-scope** items that this pass did
**not** fix.

---

## 1. Canonical number corrected: partial Spearman 0.60/0.599 → 0.62

The synergy-recovery partial Spearman (residual vs classic non-additivity,
controlling for effect magnitude) was **recomputed directly** from
`combicone_substrate.npz` using the repo's own `reachability.project_cone` over
all 131 Norman doubles:

- raw Spearman(residual, non-additivity) = **0.7361** (matches the reported 0.736);
- **partial Spearman = 0.6244 → 0.62** by both the standard three-correlation
  estimator and a double-rank estimator (they agree);
- Spearman(residual, magnitude) = −0.5599 (matches the reported −0.56).

So **0.62 is the reproducible value**; the previously-carried **0.60/0.599 was a
retired non-standard double-rank variant**. Corrected in:

| File | Was | Now |
|---|---|---|
| `results/findings.json` | `0.599` | `0.624` |
| `README.md` (certificate table) | `0.60` | `0.62` |
| `CLAUDE.md` (Frozen headline + reconcile-flag) | `0.60`, inverted flag | `0.62`, flag corrected |
| `manuscript/README.md` | `0.60` | `0.62` |

The manuscript body (`.tex`) and `docs/phase2/manuscript_numbers.json` already
carried 0.62 and were left unchanged. *(Note: `results/findings.json` = `0.624`
was subsequently absorbed into the concurrent commit `f7d94579`; see §5.)*

---

## 2. Documentation consolidation: 6 phase notes → 4 canonical docs

The user asked for a consolidated doc set, mirroring the earlier 19→7 pass. The
six prose notes under `docs/phase2/` and `docs/phase3/` were merged into four
canonical documents and then Trashed (their JSON/CSV data artifacts were **kept**).

**Created (canonical human-readable set):**

- [`docs/METHODS.md`](METHODS.md) — geometry, noise-aware two-bar certificate,
  triage, k-way, learned baseline, screenloop, acquisition, ingestion, proving-ground
  methods, reproducibility.
- [`docs/FINDINGS.md`](FINDINGS.md) — 10 findings + open requirements, tiered
  MAIN / TRANSFER / EXTENSION / PROVING-GROUND.
- [`docs/VALIDATION_REPORT.md`](VALIDATION_REPORT.md) — certified vs. source-reconstructed
  vs. independently-evaluated, plus the certificate trust argument (specificity /
  sensitivity Γ* / construct validity).
- [`docs/SCIENTIFIC_VALIDATION_PLAN.md`](SCIENTIFIC_VALIDATION_PLAN.md) — forward
  program (6 workstreams), next-dataset priorities, per-capability claim ceilings.

**Merge mapping (Trashed → absorbed into):**

| Trashed prose note | Absorbed into |
|---|---|
| `docs/phase2/cross_screen_analysis.md` | FINDINGS §4, VALIDATION_REPORT §4 |
| `docs/phase2/kway_note.md` | METHODS §4, FINDINGS §8 |
| `docs/phase2/learned_baseline_note.md` | METHODS §5, FINDINGS §3 |
| `docs/phase2/screenloop_note.md` | METHODS §6, FINDINGS §6 |
| `docs/phase3/acquisition_note.md` | METHODS §7, FINDINGS §7 |
| `docs/phase3/certificate_dossier.md` | VALIDATION_REPORT §2 |

**Retained (not touched):** all data artifacts under `docs/phase2/` (13 JSON/CSV
files incl. `manuscript_numbers.json`, `emergence_certificate.csv`), and
`docs/analytic_null_note.md` (a separate, committed methods note — cross-referenced
from METHODS §2, not one of the merged six).

The 35-vs-40 bar-(b) count and the 129-vs-128 bar-(a) count are carried as a
**documented pipeline-configuration reconciliation** (METHODS §2, FINDINGS §1),
not silently reconciled.

---

## 3. Cross-reference repointing

- `README.md` — the merged `docs/phase2/screenloop_note.md` link → `docs/FINDINGS.md`
  §6 (+ `docs/METHODS.md` §6).
- `data/README.md` — its `SCIENTIFIC_VALIDATION_PLAN.md` link now resolves.
- Repo-wide sweep: no dangling links to any merged note remain (the one residual
  match is intentional "formerly in …" provenance prose in VALIDATION_REPORT §intro).

---

## 4. Codebase tidy

- **Trashed** (byte-identical to canonical copies, SHA-256 verified, not in the
  manifest): root-level `fig_screenloop_{A_discovery,B_recovery}.{png,pdf}`
  (canonical copies live in `docs/figures/`).
- **`.gitignore`** — added a root-anchored scratch-hygiene block covering
  `/_view_*.png`, `/*_groundtruth.npz`, `/screenloop_poc.*`,
  `/screenloop_validation_*.json`, `/fig_screenloop_*.{png,pdf}`, and the
  **divergent** root regenerations `/main.pdf`, `/certificate_dossier.json` (these
  two are *not* byte-identical to `manuscript/main.pdf` / `results/certificate_dossier.json`
  — ignored so they don't shadow the canonical copies; inspect before deleting).
  Rules are root-anchored so the committed substrates (`combicone_substrate.npz`,
  `carpool_substrate.npz`) and the `docs/figures/` copies are never matched
  (verified with `git check-ignore`).
- **CLAUDE.md** (gitignored, local-only) — `## Directions and next moves` rewritten
  to mark Dirs 1/2/4 shipped and Dir 3 in-progress with pointers; `## Tiers` and
  `## Maintained surface` extended with the Phase-3 modules, the analytic null, and
  the canonical `docs/` set.

---

## 5. Concurrent work by another session (coordination note)

A **sibling analysis session ran in this same working tree during this pass** and
committed `f7d94579` ("Analytic closed-form anisotropy-corrected emergence null")
to branch **`analytic-null`** (the working tree was on `phase2-combicone-rigor` at
the start of this pass and was switched to `analytic-null` by that session). The
two efforts are **complementary**, not conflicting:

- The number fix to `results/findings.json` (0.599→0.624) was cleanly absorbed into
  `f7d94579`; HEAD now carries `0.624` alongside the sibling's analytic-null block.
- The sibling's updated `scripts/validate_findings.py` **already requires the four
  canonical docs** created here (plus `docs/analytic_null_note.md` and the analytic
  artifacts) — so the doc set produced here is exactly what its validator expects.

Because that session was still actively editing `scripts/validate_findings.py` and
regenerating `results/manifest.json`, the **validator/manifest sync (plan step 8)
was deliberately deferred to it** to avoid a concurrent-write collision on those
files. See §6 for the precise remaining edit.

Nothing was committed by this pass (branch/commit discipline is the user's; the
July-14 careful-staging constraint applies).

---

## 6. Pre-existing / out-of-scope (NOT fixed here)

These were present before this pass and are outside "consolidate the docs":

1. **`validate_findings.py` fails at `validate_guide_pair`** (the active
   `reproduce.sh` blocker). Root cause: `configs/guide_pair_transfer.json` was
   edited *after* `results/guide_pair_transfer.json` was frozen, so the config's
   live SHA-256 (`9f5f1416…`) no longer matches the hash frozen into the report
   (`01d69a30…`, which the validator also expects). The `claim_ceiling` wording
   drifted **from** "positional modalities / verified physical-guide replication"
   (the frozen report, matching the original hash) **to** "alphanumeric guide-rank
   modalities / named-sgRNA replication" (the live config). Honest resolution needs
   re-running `scripts/run_guide_pair_transfer.py` against the local 29 GB
   `GWCD4i.DE_stats.by_guide.h5mu` (regenerating the frozen report + hash), **or** a
   maintainer decision on which `claim_ceiling` wording is canonical. Not a
   documentation issue.

2. **3 stale `MANIFEST_PATHS` entries** in `scripts/validate_findings.py` (lines
   83–85): `docs/figures/make_at_a_glance.py`, `fig_at_a_glance.png`,
   `fig_at_a_glance.pdf` — all deleted in the pivot commit `8f16ccab` and superseded
   by `fig_emergence_keystone.png`. These are **latent behind the guide_pair
   blocker** (the path checks run after `validate_values`, so they are never
   reached until guide_pair is resolved). The precise fix is to delete those three
   lines and regenerate `results/manifest.json` — **deferred to the concurrent
   session that owns that file** (§5).

3. **Other uncommitted WIP left untouched** (belongs to the concurrent
   drug-combination / manuscript work, not this consolidation): the modified
   `tutorial/tutorial_combicone.ipynb`; `results/drug_combination_generalization.json`
   and `docs/figures/fig_drug_generalization.*`, `fig_predictor_anchor.*`,
   `fig_efficiency_*.png`; `manuscript/figures/*`, `manuscript/editorial_verdict.json`,
   `manuscript/limitations_and_reinforcement_plan.*`, `manuscript/manuscript_facts.json`.

---

## 7. Verification

- All cross-references in the four canonical docs + `CLAUDE.md` resolve (figure
  links, `../results`, `../data`, inter-doc links); the only forward reference is to
  this file.
- `results/findings.json` parses and carries `0.624`.
- **Change surface is markdown + `.gitignore` + file moves only — zero Python/shell
  touched** — so the Python test suite and the validator's logic are unaffected by
  this pass, and the validator reaches the *same* pre-existing `guide_pair` blocker
  it hit at session start (no regression introduced by the consolidation).
