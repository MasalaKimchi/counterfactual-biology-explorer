"""Validation A: reachable-from-lower-order on the REAL Norman doubles.

For each of the 131 measured doubles, certify its emergence against two cones:
  (i)  singles only                       (105 single-gene atoms)
  (ii) singles + all OTHER measured doubles (105 + 130 atoms)

A double that is *certified emergent* from the single-gene cone but becomes
*reachable* once the other measured doubles are atoms in the cone carries no
structure beyond what those other lower-order measurements already express.
This exercises the k-way machinery on real data with NO triples: it is an
order-dependence test between order-1 and order-2 cones.

Deterministic (fixed seed per double). Writes kway_real_results.json.
"""

from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, ".")
import combicone as cc
import reachability as rx

N_BOOT = 200
SEED = 0
FLOOR = 1.9
ALPHA = 0.05
N_WORKERS = 6

# Shared read-only arrays (populated in main). NNLS (scipy/Fortran) and the BLAS
# matmuls in project_cone both release the GIL, so a thread pool gives real
# speedup here; ProcessPoolExecutor is unavailable in the sandbox (semaphores).
_ATOMS = None
_DBL_EFF = None
_NOISE = None


def _init(atoms, dbl_eff, noise):
    global _ATOMS, _DBL_EFF, _NOISE
    _ATOMS, _DBL_EFF, _NOISE = atoms, dbl_eff, noise


def _certify_one(k: int) -> dict:
    """Certify double k from cone (i) singles and cone (ii) singles+other doubles."""
    target = _DBL_EFF[k]
    sd = _NOISE[k]
    out = {"index": k}
    # (i) singles-only cone
    try:
        ci = cc.certify_emergence(
            _ATOMS, target, noise_sd=sd, n_boot=N_BOOT,
            floor_threshold=FLOOR, alpha=ALPHA, seed=SEED,
        )
        out["i_status"] = ci.geometry_status
        out["i_unreach"] = ci.unreachable_fraction
        out["i_z"] = ci.z
        out["i_p"] = ci.p_value
        out["i_floor"] = ci.floor_ratio
        out["i_verdict"] = ci.verdict
        out["i_certified"] = bool(ci.p_value < ALPHA and ci.floor_ratio >= FLOOR)
        out["i_barA"] = bool(ci.p_value < ALPHA)
    except rx.InputError as exc:
        out["i_error"] = str(exc)
    # (ii) singles + all OTHER measured doubles
    other = np.delete(np.arange(_DBL_EFF.shape[0]), k)
    cone_ii = np.vstack([_ATOMS, _DBL_EFF[other]])
    try:
        cii = cc.certify_emergence(
            cone_ii, target, noise_sd=sd, n_boot=N_BOOT,
            floor_threshold=FLOOR, alpha=ALPHA, seed=SEED,
        )
        out["ii_status"] = cii.geometry_status
        out["ii_unreach"] = cii.unreachable_fraction
        out["ii_z"] = cii.z
        out["ii_p"] = cii.p_value
        out["ii_floor"] = cii.floor_ratio
        out["ii_verdict"] = cii.verdict
        out["ii_certified"] = bool(cii.p_value < ALPHA and cii.floor_ratio >= FLOOR)
        out["ii_barA"] = bool(cii.p_value < ALPHA)
    except rx.InputError as exc:
        out["ii_error"] = str(exc)
    return out


def main() -> None:
    d = np.load("combicone_substrate.npz", allow_pickle=True)
    atoms = d["atoms"]
    ctrl = d["ctrl"]
    means = d["means"]
    means1 = d["means1"]
    means2 = d["means2"]
    conditions = np.asarray(d["conditions"]).astype(str)
    doubles = [str(x) for x in d["doubles"]]

    cond_idx = {c: i for i, c in enumerate(conditions)}
    rows = np.array([cond_idx[n] for n in doubles])
    dbl_eff = means[rows] - ctrl
    noise = np.abs(means1[rows] - means2[rows]) / 2.0

    n = len(doubles)
    t0 = time.time()
    _init(atoms, dbl_eff, noise)  # threads share the same address space
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        results = list(ex.map(_certify_one, range(n)))
    results.sort(key=lambda r: r["index"])
    for r, name in zip(results, doubles):
        r["name"] = name

    payload = {
        "n_doubles": n,
        "n_boot": N_BOOT,
        "floor_threshold": FLOOR,
        "alpha": ALPHA,
        "cone_i_atoms": int(atoms.shape[0]),
        "cone_ii_atoms": int(atoms.shape[0] + n - 1),
        "results": results,
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open("kway_real_results.json", "w") as fh:
        json.dump(payload, fh, indent=2)
    # terse liveness summary
    ci_cert = sum(r.get("i_certified", False) for r in results)
    cii_cert = sum(r.get("ii_certified", False) for r in results)
    lost = sum(
        r.get("i_certified", False) and not r.get("ii_certified", False)
        for r in results
    )
    print(f"done {n} doubles in {payload['elapsed_s']}s")
    print(f"certified emergent: cone_i={ci_cert}  cone_ii={cii_cert}  flipped={lost}")


if __name__ == "__main__":
    main()