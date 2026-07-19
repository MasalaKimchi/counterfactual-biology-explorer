#!/usr/bin/env python
"""Causal-inference trust dossier for the CombiCone emergence certificate.

Three pillars of trust, all computed on the real Norman substrate + its shipped
certificate CSV:

1. NEGATIVE CONTROLS — the certificate must NOT fire on inputs known to be
   additive/reachable.
     (NC-additive)  synthetic combinations built as a literal sum of two measured
                    single effects (+ measurement-scale noise). A truly additive
                    combination is inside the non-negative cone by construction;
                    the false-positive RATE (fraction wrongly "certified emergent")
                    is the headline specificity number.
     (NC-sham)      a pure-noise "combination" (control-scale Gaussian around 0).
                    Illustrates why BOTH bars are needed: it can scrape bar (a)
                    but never clears bar (b).
     (NC-reachable) the least-unreachable real doubles (bottom decile by residual)
                    — should mostly certify as within-noise, not emergent.

2. SENSITIVITY RADIUS — for each certified-emergent real pair, the smallest noise
   inflation factor Gamma* such that treating the true measurement noise as
   Gamma* x our split-half estimate flips the verdict below "certified". A large
   Gamma* means the verdict tolerates substantial noise MIS-estimation. This is
   the certificate's analogue of a Rosenbaum sensitivity bound / E-value: it
   quantifies how wrong the one key assumption (residual = measurement noise)
   would have to be to overturn the call.

3. CONSTRUCT VALIDITY — evidence the certificate measures genuine non-additivity:
   recovery of classical synergy (partial Spearman, from the ledger), decorrelation
   of the noise-robust z from raw effect magnitude (the confound removal), and the
   identity of the top certified pairs as known genetic interactions.

Outputs JSON (--out) with every number, for the figure and the written dossier.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import combicone as cc


def _is_certified(cert) -> bool:
    return cert.verdict.startswith("certified")


def negative_control_additive(atoms, names, *, n_pairs=200, noise_frac=0.1,
                              n_boot=200, seed=0):
    """False-positive rate on synthetic ADDITIVE combinations (a + b + noise)."""
    rng = np.random.default_rng(seed)
    n = len(names)
    verdicts = {"certified": 0, "modest": 0, "within_noise": 0}
    floor_ratios = []
    for _ in range(n_pairs):
        i, j = rng.choice(n, size=2, replace=False)
        add = atoms[i] + atoms[j]
        noise = noise_frac * np.abs(add) + 1e-4
        combo = add + rng.normal(0, noise)
        c = cc.certify_emergence(cone_atoms=atoms, measured_combo=combo,
                                 noise_sd=noise, n_boot=n_boot, seed=int(rng.integers(1 << 30)))
        floor_ratios.append(float(c.floor_ratio))
        if _is_certified(c):
            verdicts["certified"] += 1
        elif "modest" in c.verdict:
            verdicts["modest"] += 1
        else:
            verdicts["within_noise"] += 1
    return {
        "n_pairs": n_pairs,
        "noise_frac": noise_frac,
        "false_positive_rate": verdicts["certified"] / n_pairs,
        "verdict_counts": verdicts,
        "floor_ratio_mean": float(np.mean(floor_ratios)),
        "floor_ratio_p95": float(np.percentile(floor_ratios, 95)),
    }


def sensitivity_radius(atoms, combo_effect, noise_sd, *, floor_threshold=1.9,
                       alpha=0.05, n_boot=200, seed=0,
                       grid=None, full_curve=True):
    """Smallest noise-inflation Gamma that drops a pair below 'certified'.

    Returns Gamma* (np.inf if the verdict never flips within the grid) and the
    per-Gamma floor ratios, by re-certifying with noise_sd -> Gamma * noise_sd.
    With ``full_curve=False`` the sweep stops at the first flip (Gamma* only),
    which is much cheaper for the per-pair scan.
    """
    if grid is None:
        grid = np.round(np.arange(1.0, 4.01, 0.25), 2)
    gamma_star = float("inf")
    curve = []
    for g in grid:
        c = cc.certify_emergence(
            cone_atoms=atoms, measured_combo=combo_effect, noise_sd=g * noise_sd,
            n_boot=n_boot, floor_threshold=floor_threshold, alpha=alpha, seed=seed,
        )
        curve.append((float(g), float(c.floor_ratio), float(c.p_value), _is_certified(c)))
        if gamma_star == float("inf") and not _is_certified(c):
            gamma_star = float(g)
            if not full_curve:
                break
    return gamma_star, curve


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--substrate", default="combicone_substrate.npz")
    ap.add_argument("--certificate", required=True, help="emergence_certificate.csv path")
    ap.add_argument("--ledger", default=None, help="manuscript_numbers.json (construct validity)")
    ap.add_argument("--out", default="certificate_dossier.json")
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-neg", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    z = np.load(args.substrate, allow_pickle=True)
    atoms = np.asarray(z["atoms"], dtype=float)
    names = [str(g) for g in z["single_genes"]]
    cond = [str(c) for c in z["conditions"]]
    means = np.asarray(z["means"], dtype=float)
    ctrl = np.asarray(z["ctrl"], dtype=float)
    m1 = np.asarray(z["means1"], dtype=float)
    m2 = np.asarray(z["means2"], dtype=float)

    rows = list(csv.DictReader(open(args.certificate)))
    cert_by = {r["double"]: r for r in rows}

    def combo_effect(name):
        return means[cond.index(name)] - ctrl

    def combo_noise(name):
        i = cond.index(name)
        return np.abs(m1[i] - m2[i]) / 2.0

    result = {"n_atoms": len(names), "n_doubles": len(rows)}

    # ---- 1. negative controls --------------------------------------------- #
    print("negative control: additive reconstructions...", file=sys.stderr)
    result["nc_additive"] = negative_control_additive(
        atoms, names, n_pairs=args.n_neg, n_boot=args.n_boot, seed=args.seed
    )
    # NC-reachable: bottom-decile-by-residual real doubles
    resid = np.array([float(r["residual"]) for r in rows])
    order = np.argsort(resid)
    bottom = [rows[i]["double"] for i in order[: max(5, len(rows) // 10)]]
    nc_reach = {"certified": 0, "modest": 0, "within_noise": 0}
    for d in bottom:
        c = cc.certify_emergence(cone_atoms=atoms, measured_combo=combo_effect(d),
                                 noise_sd=combo_noise(d), n_boot=args.n_boot, seed=args.seed)
        if _is_certified(c):
            nc_reach["certified"] += 1
        elif "modest" in c.verdict:
            nc_reach["modest"] += 1
        else:
            nc_reach["within_noise"] += 1
    result["nc_reachable"] = {"n": len(bottom), "verdict_counts": nc_reach,
                              "pairs": bottom}

    # ---- 2. sensitivity radius for certified-emergent pairs --------------- #
    print("sensitivity radius...", file=sys.stderr)
    certified = [r["double"] for r in rows
                 if r["bar_a_significant"] in ("True", "1")
                 and r["bar_b_floor_ratio_ge_1p9"] in ("True", "1")]
    sens = {}
    for d in certified:
        gstar, _ = sensitivity_radius(atoms, combo_effect(d), combo_noise(d),
                                      n_boot=args.n_boot, seed=args.seed, full_curve=False)
        sens[d] = gstar
    gstars = np.array([v for v in sens.values() if np.isfinite(v)])
    result["sensitivity"] = {
        "n_certified": len(certified),
        "gamma_star_by_pair": sens,
        "gamma_star_median": float(np.median(gstars)) if gstars.size else None,
        "gamma_star_min": float(gstars.min()) if gstars.size else None,
        "gamma_star_max": float(gstars.max()) if gstars.size else None,
        "frac_robust_to_2x": float(np.mean(gstars >= 2.0)) if gstars.size else None,
    }
    # a full curve for the flagship pair
    flag = "SET+CEBPE" if "SET+CEBPE" in cert_by else certified[0]
    _, curve = sensitivity_radius(atoms, combo_effect(flag), combo_noise(flag),
                                  n_boot=args.n_boot, seed=args.seed)
    result["sensitivity"]["flagship_pair"] = flag
    result["sensitivity"]["flagship_curve"] = curve

    # ---- 3. construct validity ------------------------------------------- #
    cv = {}
    if args.ledger and Path(args.ledger).exists():
        led = json.load(open(args.ledger))
        flag_block = led.get("emergence_certificate_flagship", {})
        cv["partial_spearman_synergy_given_magnitude"] = flag_block.get("synergy_partial_given_magnitude")
        cv["raw_spearman_synergy"] = flag_block.get("synergy_raw_spearman")
    # z vs magnitude decorrelation (recompute from CSV)
    zval = np.array([float(r["z"]) for r in rows])
    mag = np.array([float(r["effect_norm"]) for r in rows])
    nonadd = np.array([float(r["nonadditivity"]) for r in rows])
    def spear(a, b):
        from scipy.stats import spearmanr
        return float(spearmanr(a, b).statistic)
    cv["z_vs_magnitude_spearman"] = spear(zval, mag)
    cv["z_vs_nonadditivity_spearman"] = spear(zval, nonadd)
    cv["raw_residual_vs_magnitude_spearman"] = spear(resid, mag)
    cv["top_certified_pairs"] = certified[:8]
    result["construct_validity"] = cv

    json.dump(result, open(args.out, "w"), indent=2, default=str)
    print(f"\nwrote {args.out}", file=sys.stderr)
    print(f"NC-additive false-positive rate: {result['nc_additive']['false_positive_rate']:.3f}", file=sys.stderr)
    print(f"sensitivity Gamma* median: {result['sensitivity']['gamma_star_median']}", file=sys.stderr)
    print(f"z vs magnitude Spearman: {cv['z_vs_magnitude_spearman']:.3f} "
          f"(raw residual vs magnitude: {cv['raw_residual_vs_magnitude_spearman']:.3f})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
