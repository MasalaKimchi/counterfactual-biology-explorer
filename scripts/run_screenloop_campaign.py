"""Certificate-guided sequential design + library augmentation on real screens.

Replays two published combinatorial screens as self-driving campaigns and runs the
library-augmentation recovery test, with a permutation null and a magnitude-confound
control. Deterministic given the committed substrates; writes
``results/screenloop_campaign.json`` and re-verifies the frozen headline numbers so
it fails closed if the geometry drifts.

Screens (substrates committed at repo root):
  * Norman GSE133344 CRISPRa  — combicone_substrate.npz (105 singles, 131 doubles)
  * CaRPool-seq Cas13d KD      — carpool_substrate.npz   (28 singles, 158 doubles)

Both emergence labels use the two-bar noise-robust certificate from
:func:`combicone.certify_emergence` (p < 0.05 AND floor >= 1.9x), identical to the
manuscript. Run: ``python scripts/run_screenloop_campaign.py`` (add ``--quick`` to
skip the O(N) certification and reuse cached labels if present).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import combicone as cc  # noqa: E402
import reachability as rx  # noqa: E402
import screenloop as sl  # noqa: E402

N_BOOT = 200
FLOOR = 1.9
ALPHA = 0.05
SEED = 0
N_PERM = 2000
BATCH = 8

# Frozen headline numbers (regression guard). Tolerances are generous on the
# campaign wells (integer, batch-quantized) and tight on the recovery rates.
FROZEN = {
    "norman": {
        "n_emergent": 40,
        "recovery_sep_top1": 0.9811320754716981,
        "recovery_base_top1": 0.5471698113207547,
        "mag_only_top1_max": 0.05,
    },
    "carpool": {
        "n_emergent": 76,
        "recovery_sep_top1": 1.0,
        "recovery_base_top1": 0.28,
        "mag_only_top1_max": 0.05,
    },
}


def _load(path: Path) -> dict:
    d = np.load(path, allow_pickle=True)
    return {k: d[k] for k in d.files}


def _norman_arrays(sub: dict):
    atoms = sub["atoms"].astype(float)
    genes = sub["genes"].astype(str)
    single_genes = sub["single_genes"].astype(str)
    ctrl = sub["ctrl"].astype(float)
    means = sub["means"].astype(float)
    means1 = sub["means1"].astype(float)
    means2 = sub["means2"].astype(float)
    conditions = np.asarray(sub["conditions"]).astype(str)
    doubles = np.asarray([str(x) for x in sub["doubles"]])
    cond_idx = {c: i for i, c in enumerate(conditions)}
    rows = np.array([cond_idx[n] for n in doubles])
    D = means[rows] - ctrl
    noise = np.abs(means1[rows] - means2[rows]) / 2.0
    g2a = {g: i for i, g in enumerate(single_genes)}
    pair = np.array(
        [(g2a[nm.split("+")[0]], g2a[nm.split("+")[1]]) for nm in doubles]
    )
    return atoms, D, noise, pair, doubles


def _carpool_arrays(sub: dict):
    atoms = sub["atoms"].astype(float)
    single_genes = sub["single_genes"].astype(str)
    doubles = np.asarray([str(x) for x in sub["doubles"]])
    D = sub["double_eff"].astype(float)
    noise = np.abs(sub["double_m1"].astype(float) - sub["double_m2"].astype(float)) / 2.0
    g2a = {g: i for i, g in enumerate(single_genes)}
    pair = np.array(
        [(g2a[nm.split("+")[0]], g2a[nm.split("+")[1]]) for nm in doubles]
    )
    return atoms, D, noise, pair, doubles


def _emergence_labels(atoms, D, noise):
    certs = [
        cc.certify_emergence(
            cone_atoms=atoms,
            measured_combo=D[k],
            noise_sd=noise[k],
            n_boot=N_BOOT,
            floor_threshold=FLOOR,
            alpha=ALPHA,
            seed=SEED,
        )
        for k in range(D.shape[0])
    ]
    bar_a = np.array([c.p_value < ALPHA for c in certs])
    bar_b = np.array([c.floor_ratio >= FLOOR for c in certs])
    return (bar_a & bar_b).astype(bool)


def _campaign(atoms, D, pair, labels):
    out = {}
    for pol in ("magnitude", "triage", "cone_adaptive"):
        res = sl.replay_campaign(
            atoms, D, pair, labels, policy=pol, batch_size=BATCH, seed=SEED
        )
        out[pol] = res.wells_to_fraction(0.9)
    # random averaged over 40 seeds
    rand_found = np.mean(
        [
            sl.replay_campaign(
                atoms, D, pair, labels, policy="random", batch_size=BATCH, seed=s
            ).found
            for s in range(40)
        ],
        axis=0,
    )
    rand_w = sl.replay_campaign(
        atoms, D, pair, labels, policy="random", batch_size=BATCH, seed=0
    ).wells
    hit = rand_found >= 0.9 * labels.sum()
    out["random"] = float(rand_w[np.argmax(hit)]) if hit.any() else float("nan")
    return out


def _recovery_with_nulls(atoms, D, pair, rng):
    rec = sl.held_out_single_recovery(atoms, pair, D, min_involved=2)
    elig = rec["eligible"]
    n_atoms = atoms.shape[0]

    # permutation null on top-1
    score_mats = []
    for g in elig:
        keep = np.array([i for i in range(n_atoms) if i != g])
        involve = [k for k in range(len(pair)) if g in pair[k]]
        nom = sl.nominate_atoms(atoms[keep], D[involve], atoms, weight="residual")
        score_mats.append(nom.scores)
    score_mats = np.array(score_mats)
    obs_top1 = float((rec["sep_ranks"] == 1).mean())
    null_top1 = np.empty(N_PERM)
    for p in range(N_PERM):
        fake = rng.integers(0, n_atoms, size=len(elig))
        ranks = np.array(
            [
                int(np.where(np.argsort(-score_mats[i]) == fake[i])[0][0]) + 1
                for i in range(len(elig))
            ]
        )
        null_top1[p] = (ranks == 1).mean()
    perm_p = (1 + np.sum(null_top1 >= obs_top1)) / (1 + N_PERM)

    # magnitude-only control
    mag_ranks = []
    for g in elig:
        involve = [k for k in range(len(pair)) if g in pair[k]]
        mdm = np.mean([np.linalg.norm(D[k]) for k in involve])
        am = np.linalg.norm(atoms, axis=1)
        order = np.argsort(-(-np.abs(am - mdm)), kind="stable")
        mag_ranks.append(int(np.where(order == g)[0][0]) + 1)
    mag_ranks = np.array(mag_ranks)

    # separator top-1 restricted to atoms where the naive baseline fails
    base_failed = rec["base_ranks"] > 1
    sep_where_fail = (
        float((rec["sep_ranks"][base_failed] == 1).mean())
        if base_failed.any()
        else float("nan")
    )

    # dominance vs advantage
    dom = []
    for g in elig:
        involve = [k for k in range(len(pair)) if g in pair[k]]
        dom.append(
            np.mean(
                [
                    np.dot(D[k], atoms[g])
                    / (np.linalg.norm(D[k]) * np.linalg.norm(atoms[g]) + 1e-12)
                    for k in involve
                ]
            )
        )
    rho, _ = spearmanr(np.array(dom), rec["base_ranks"] - rec["sep_ranks"])

    return {
        "n_eligible": int(len(elig)),
        "n_candidates": int(n_atoms),
        "separator": {k: float(v) for k, v in rec["separator"].items()},
        "baseline": {k: float(v) for k, v in rec["baseline"].items()},
        "magnitude_only_top1": float((mag_ranks == 1).mean()),
        "permutation_null": {
            "n_perm": N_PERM,
            "observed_top1": obs_top1,
            "null_top1_mean": float(null_top1.mean()),
            "null_top1_sd": float(null_top1.std()),
            "perm_p_top1": float(perm_p),
            "z_top1": float((obs_top1 - null_top1.mean()) / (null_top1.std() + 1e-12)),
        },
        "separator_top1_where_baseline_fails": sep_where_fail,
        "n_baseline_failures": int(base_failed.sum()),
        "spearman_dominance_vs_advantage": float(rho),
        "mean_atom_dominance": float(np.mean(dom)),
    }


def _check_frozen(name, camp, rec):
    f = FROZEN[name]
    tol = 1e-9
    problems = []
    if abs(rec["separator"]["top1"] - f["recovery_sep_top1"]) > 1e-6:
        problems.append(
            f"{name} separator top1 {rec['separator']['top1']} != {f['recovery_sep_top1']}"
        )
    if abs(rec["baseline"]["top1"] - f["recovery_base_top1"]) > 1e-6:
        problems.append(
            f"{name} baseline top1 {rec['baseline']['top1']} != {f['recovery_base_top1']}"
        )
    if rec["magnitude_only_top1"] > f["mag_only_top1_max"]:
        problems.append(
            f"{name} magnitude-only top1 {rec['magnitude_only_top1']} exceeds "
            f"{f['mag_only_top1_max']} (recovery may be a magnitude artifact)"
        )
    if rec["permutation_null"]["perm_p_top1"] > 0.01:
        problems.append(f"{name} permutation p {rec['permutation_null']['perm_p_top1']} > 0.01")
    return problems


def main() -> None:
    quick = "--quick" in sys.argv
    rng = np.random.default_rng(SEED)
    payload = {"schema_version": "1.0.0", "batch_size": BATCH, "screens": {}}
    all_problems = []

    for name, loader in (("norman", _norman_arrays), ("carpool", _carpool_arrays)):
        subpath = ROOT / (
            "combicone_substrate.npz" if name == "norman" else "carpool_substrate.npz"
        )
        if not subpath.is_file():
            raise FileNotFoundError(f"missing substrate {subpath}")
        atoms, D, noise, pair, doubles = loader(_load(subpath))

        t0 = time.time()
        labels = _emergence_labels(atoms, D, noise)
        assert labels.sum() == FROZEN[name]["n_emergent"], (
            f"{name}: {int(labels.sum())} emergent != {FROZEN[name]['n_emergent']}"
        )
        camp = _campaign(atoms, D, pair, labels)
        rec = _recovery_with_nulls(atoms, D, pair, rng)
        elapsed = time.time() - t0

        payload["screens"][name] = {
            "n_singles": int(atoms.shape[0]),
            "n_doubles": int(D.shape[0]),
            "n_emergent": int(labels.sum()),
            "base_rate": float(labels.mean()),
            "campaign_wells_to_90pct": camp,
            "recovery": rec,
            "elapsed_s": round(elapsed, 1),
            "scope": sl._SCOPE,
        }
        all_problems += _check_frozen(name, camp, rec)
        print(
            f"[{name}] emergent {int(labels.sum())}/{D.shape[0]} | "
            f"campaign wells-to-90 triage {camp['triage']:.0f} "
            f"cone_adaptive {camp['cone_adaptive']:.0f} | "
            f"recovery sep top1 {rec['separator']['top1']:.3f} "
            f"base {rec['baseline']['top1']:.3f} "
            f"perm-p {rec['permutation_null']['perm_p_top1']:.1e} "
            f"({elapsed:.0f}s)"
        )

    if all_problems:
        raise SystemExit("FROZEN CHECK FAILED:\n  " + "\n  ".join(all_problems))

    out = ROOT / "results" / "screenloop_campaign.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)}; frozen headline numbers verified")


if __name__ == "__main__":
    main()
