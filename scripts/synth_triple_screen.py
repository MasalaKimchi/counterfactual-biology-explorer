"""Synthetic combinatorial screen with PLANTED 3-way epistasis.

Norman contains no measured triples, so the only way to validate order-3
certification against a known ground truth is a synthetic screen. This generator
plants three classes of triples and hands them to :func:`combicone.certify_emergence`
at order 3 (cone = singles + the three constituent measured doubles):

  * ADDITIVE       T = S_i + S_j + S_l                       (no interaction)
  * 2-WAY-REDUCIBLE T = non-neg combo of the 3 constituent   (reachable: its
                        doubles (+ singles)                    structure is fully
                                                               expressed by lower
                                                               order measurements)
  * 3-WAY-EMERGENT  T = (reducible base) + Tt, with Tt        (a component that
                        orthogonal to the span of EVERY        exists ONLY at the
                        atom (all singles + all doubles)       triple)

Everything is built to be consistent with the engine's non-negative-cone model:
"reachable" means expressible as a non-negative mixture of the measured
lower-order effects. The additive and reducible classes are placed inside that
cone by construction; the emergent class is pushed outside it by exactly the
planted 3-way term Tt and nothing else. Tt magnitudes span a range so the recovery
problem is non-trivial (some planted triples sit near the measurement-noise floor).

Design choices, disclosed:
  * The cone atoms (singles, doubles) are the NOISE-FREE true effects; only the
    measured TRIPLES carry Gaussian measurement noise (and a split-half estimate
    of it). This isolates the certificate's ability to detect a planted 3-way term
    above triple-measurement noise from a clean lower-order cone -- the standard
    setup for a planted-truth recovery test. Noisy atoms would add cone jitter that
    is not the quantity under test.
  * SYNTHETIC. This validates the order-3 code path and its discrimination; it is
    NOT evidence about real biological 3-way epistasis (none is measured in Norman).

Deterministic given ``seed``. ``python synth_triple_screen.py`` writes
synth_triple_screen.npz and prints a summary.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass, asdict

import numpy as np


@dataclass(frozen=True)
class ScreenConfig:
    n_genes: int = 300
    n_singles: int = 14
    single_norm: float = 5.0          # typical single-effect magnitude
    single_density: float = 0.4       # fraction of genes a single touches
    pair_norm: float = 2.5            # true 2-way interaction magnitude
    pair_density: float = 0.35
    n_additive: int = 40
    n_reducible: int = 40
    n_emergent: int = 40
    tt_norm_lo: float = 1.5           # planted 3-way term magnitude range
    tt_norm_hi: float = 7.0
    noise_tau: float = 0.10           # per-split-half per-gene measurement noise sd
    seed: int = 0


def _sparse_vec(rng, n_genes, density, norm):
    v = np.zeros(n_genes)
    k = max(1, int(round(density * n_genes)))
    idx = rng.choice(n_genes, size=k, replace=False)
    v[idx] = rng.normal(size=k)
    nv = np.linalg.norm(v)
    return v * (norm / nv) if nv > 0 else v


def generate(config: ScreenConfig | None = None) -> dict:
    cfg = config or ScreenConfig()
    rng = np.random.default_rng(cfg.seed)
    G, K = cfg.n_genes, cfg.n_singles

    # --- singles ---------------------------------------------------------- #
    S = np.vstack([_sparse_vec(rng, G, cfg.single_density, cfg.single_norm)
                   for _ in range(K)])
    single_names = [f"S{i}" for i in range(K)]

    # --- true pairwise interactions + measured doubles (ALL C(K,2) pairs) - #
    all_pairs = list(itertools.combinations(range(K), 2))
    P = {p: _sparse_vec(rng, G, cfg.pair_density, cfg.pair_norm) for p in all_pairs}
    D = {p: S[p[0]] + S[p[1]] + P[p] for p in all_pairs}
    double_pairs = all_pairs
    double_atoms = np.vstack([D[p] for p in double_pairs])   # (n_pairs, G)

    # Full atom span for orthogonalizing the planted 3-way terms.
    atom_span = np.vstack([S, double_atoms])                 # (K + n_pairs, G)
    # Orthonormal basis of the atom row space via SVD (row space = span of atoms).
    U, sv, Vt = np.linalg.svd(atom_span, full_matrices=False)
    rank = int(np.sum(sv > 1e-9 * sv[0]))
    span_basis = Vt[:rank]                                   # (rank, G) orthonormal

    def orthogonalize(v):
        return v - span_basis.T @ (span_basis @ v)

    # --- choose triples, one class each, disjoint sets -------------------- #
    all_triples = list(itertools.combinations(range(K), 3))
    rng.shuffle(all_triples)
    n_total = cfg.n_additive + cfg.n_reducible + cfg.n_emergent
    if n_total > len(all_triples):
        raise ValueError("requested more triples than C(K,3) available")
    chosen = all_triples[:n_total]
    classes = (["additive"] * cfg.n_additive
               + ["reducible"] * cfg.n_reducible
               + ["emergent"] * cfg.n_emergent)
    # interleave so classes are not blocked (guards against any index-order bias)
    order = rng.permutation(n_total)
    chosen = [chosen[i] for i in order]
    classes = [classes[i] for i in order]

    tt_targets = np.linspace(cfg.tt_norm_lo, cfg.tt_norm_hi, cfg.n_emergent)
    rng.shuffle(tt_targets)

    triples_true = []
    triple_members = []
    triple_class = []
    tt_norm = []
    emergent_seen = 0
    for (i, j, l), cls in zip(chosen, classes):
        members = (i, j, l)
        cpairs = [(i, j), (i, l), (j, l)]
        if cls == "additive":
            T = S[i] + S[j] + S[l]
            planted = 0.0
        elif cls == "reducible":
            cd = rng.uniform(0.5, 1.5, size=3)     # non-neg coeffs on the 3 doubles
            cs = rng.uniform(0.0, 0.5, size=3)     # non-neg coeffs on the 3 singles
            T = sum(cd[k] * D[cpairs[k]] for k in range(3)) \
                + cs[0] * S[i] + cs[1] * S[j] + cs[2] * S[l]
            planted = 0.0
        else:  # emergent
            cd = rng.uniform(0.5, 1.5, size=3)
            base = sum(cd[k] * D[cpairs[k]] for k in range(3))
            raw = rng.normal(size=G)
            tt = orthogonalize(raw)
            ntt = np.linalg.norm(tt)
            target_norm = float(tt_targets[emergent_seen]); emergent_seen += 1
            tt = tt * (target_norm / ntt)
            T = base + tt
            planted = target_norm
        triples_true.append(T)
        triple_members.append(members)
        triple_class.append(cls)
        tt_norm.append(planted)

    triples_true = np.vstack(triples_true)                   # (n_total, G)

    # --- measure the triples with split-half noise ------------------------ #
    tau = cfg.noise_tau
    h1 = triples_true + rng.normal(0.0, tau, size=triples_true.shape)
    h2 = triples_true + rng.normal(0.0, tau, size=triples_true.shape)
    triples_measured = 0.5 * (h1 + h2)
    noise_sd = np.abs(h1 - h2) / 2.0                         # per-gene SE estimate

    labels = np.array([1 if c == "emergent" else 0 for c in triple_class], dtype=int)

    return {
        "config": asdict(cfg),
        "singles": S,
        "single_names": np.array(single_names),
        "double_atoms": double_atoms,
        "double_pairs": np.array(double_pairs),
        "triples_measured": triples_measured,
        "triples_true": triples_true,
        "triple_members": np.array(triple_members),
        "triple_class": np.array(triple_class),
        "triple_label_emergent": labels,
        "tt_norm": np.array(tt_norm),
        "noise_sd": noise_sd,
    }


def main() -> None:
    data = generate()
    np.savez("synth_triple_screen.npz", **{k: v for k, v in data.items()
                                            if k != "config"})
    with open("synth_triple_screen_config.json", "w") as fh:
        json.dump(data["config"], fh, indent=2)
    cls, cnt = np.unique(data["triple_class"], return_counts=True)
    print("classes:", dict(zip(cls.tolist(), cnt.tolist())))
    print("n_singles:", data["singles"].shape[0],
          "n_double_atoms:", data["double_atoms"].shape[0],
          "n_triples:", data["triples_measured"].shape[0],
          "n_genes:", data["singles"].shape[1])
    print("planted Tt norm range (emergent):",
          round(float(data["tt_norm"][data["tt_norm"] > 0].min()), 2), "..",
          round(float(data["tt_norm"].max()), 2))


if __name__ == "__main__":
    main()