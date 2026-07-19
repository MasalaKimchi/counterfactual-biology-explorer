"""Validation B: run order-3 certify_emergence on the synthetic triple screen.

Loads synth_triple_screen.npz (from synth_triple_screen.py), certifies every
measured triple at order 3 against the cone {singles + its three constituent
measured doubles}, and scores recovery of the planted 3-way-emergent class against
the additive and 2-way-reducible classes. Writes synth_triple_results.json.

SYNTHETIC: validates the order-3 code path and its discrimination; NOT evidence
about real biological 3-way epistasis (Norman measures no triples).
"""

from __future__ import annotations

import json
import sys

import numpy as np

sys.path.insert(0, ".")
import combicone as cc

N_BOOT = 200
SEED = 0
FLOOR = 1.9
ALPHA = 0.05


def _auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney AUC, ties at 0.5 (no sklearn in the pinned env)."""
    pos, neg = scores[labels], scores[~labels]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def main() -> None:
    d = np.load("synth_triple_screen.npz", allow_pickle=True)
    S = d["singles"]
    dpairs = [tuple(int(x) for x in p) for p in d["double_pairs"]]
    datoms = d["double_atoms"]
    dmap = {tuple(sorted(dpairs[i])): datoms[i] for i in range(len(dpairs))}
    Tm = d["triples_measured"]
    members = [tuple(int(x) for x in m) for m in d["triple_members"]]
    cls = [str(c) for c in d["triple_class"]]
    nsd = d["noise_sd"]
    tt = d["tt_norm"]

    zs, unr, floors, ps, certified, barA = [], [], [], [], [], []
    for idx, (mem, sd) in enumerate(zip(members, nsd)):
        i, j, l = mem
        cone = np.vstack(
            [S,
             dmap[tuple(sorted((i, j)))],
             dmap[tuple(sorted((i, l)))],
             dmap[tuple(sorted((j, l)))]]
        )
        c = cc.certify_emergence(
            cone_atoms=cone, measured_combo=Tm[idx], noise_sd=sd,
            n_boot=N_BOOT, floor_threshold=FLOOR, alpha=ALPHA, seed=SEED,
        )
        zs.append(c.z); unr.append(c.unreachable_fraction); floors.append(c.floor_ratio)
        ps.append(c.p_value)
        certified.append(bool(c.p_value < ALPHA and c.floor_ratio >= FLOOR))
        barA.append(bool(c.p_value < ALPHA))

    zs = np.array(zs); unr = np.array(unr); floors = np.array(floors)
    ps = np.array(ps); certified = np.array(certified); barA = np.array(barA)
    lab = np.array([c == "emergent" for c in cls])

    TP = int((certified & lab).sum()); FN = int((~certified & lab).sum())
    FP = int((certified & ~lab).sum()); TN = int((~certified & ~lab).sum())
    sens = TP / (TP + FN) if (TP + FN) else float("nan")
    spec = TN / (TN + FP) if (TN + FP) else float("nan")
    prec = TP / (TP + FP) if (TP + FP) else float("nan")

    def per_class(name):
        m = np.array([c == name for c in cls])
        return [int(certified[m].sum()), int(m.sum())]

    out = {
        "n_triples": len(cls),
        "n_boot": N_BOOT, "floor_threshold": FLOOR, "alpha": ALPHA,
        "classes": {n: int(np.sum([c == n for c in cls]))
                    for n in ("additive", "reducible", "emergent")},
        "auc_z": _auc(zs, lab),
        "auc_unreachable_fraction": _auc(unr, lab),
        "auc_floor_ratio": _auc(floors, lab),
        "two_bar": {"TP": TP, "FN": FN, "FP": FP, "TN": TN,
                    "sensitivity": sens, "specificity": spec, "precision": prec},
        "bar_a_only": {  # permissive p<alpha bar alone (no effect-size gate)
            "additive": int(barA[[c == "additive" for c in cls]].sum()),
            "reducible": int(barA[[c == "reducible" for c in cls]].sum()),
            "emergent": int(barA[lab].sum()),
        },
        "certified_per_class": {"additive": per_class("additive"),
                                "reducible": per_class("reducible"),
                                "emergent": per_class("emergent")},
        "unreachable_range_per_class": {
            n: [float(unr[[c == n for c in cls]].min()),
                float(unr[[c == n for c in cls]].max())]
            for n in ("additive", "reducible", "emergent")},
        "per_triple": [
            {"members": list(members[i]), "class": cls[i], "tt_norm": float(tt[i]),
             "z": float(zs[i]), "unreachable_fraction": float(unr[i]),
             "floor_ratio": float(floors[i]), "p_value": float(ps[i]),
             "certified": bool(certified[i])}
            for i in range(len(cls))],
        "note": ("SYNTHETIC planted-epistasis screen. Validates the order-3 code "
                 "path and its discrimination of 3-way-emergent triples from "
                 "additive / 2-way-reducible ones. NOT evidence about real "
                 "biological 3-way epistasis: Norman measures 0 triples."),
    }
    with open("synth_triple_results.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"AUC z={out['auc_z']:.3f}  AUC raw-unreach={out['auc_unreachable_fraction']:.3f}")
    print(f"two-bar: TP={TP} FN={FN} FP={FP} TN={TN} "
          f"sens={sens:.3f} spec={spec:.3f} prec={prec:.3f}")
    print(f"certified per class: add={per_class('additive')} "
          f"red={per_class('reducible')} em={per_class('emergent')}")


if __name__ == "__main__":
    main()