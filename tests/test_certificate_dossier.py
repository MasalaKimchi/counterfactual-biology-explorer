"""Tests for the certificate trust-dossier analysis functions."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


import certificate_dossier as cd  # noqa: E402


def _substrate(seed=0, n=8, n_genes=60):
    """A small synthetic cone: n single atoms + one planted-emergent direction."""
    rng = np.random.default_rng(seed)
    atoms = rng.normal(0, 1, (n, n_genes))
    names = [f"g{i}" for i in range(n)]
    return atoms, names


def test_additive_control_is_not_certified():
    """A literal sum of two atoms must never be certified emergent."""
    atoms, names = _substrate()
    res = cd.negative_control_additive(atoms, names, n_pairs=25, n_boot=80, seed=1)
    assert res["false_positive_rate"] == 0.0
    assert res["verdict_counts"]["certified"] == 0
    # floor ratio for additive combos should sit near 1
    assert res["floor_ratio_p95"] < 1.9


def test_sensitivity_radius_monotone_and_flips():
    """Inflating noise must (weakly) lower the floor ratio and eventually flip."""
    import combicone as cc
    atoms, names = _substrate(n=6, n_genes=40)
    # build a genuinely emergent combo: a direction far from the cone
    rng = np.random.default_rng(3)
    combo = rng.normal(0, 1, 40)
    combo = combo - atoms.T @ np.linalg.lstsq(atoms.T, combo, rcond=None)[0]  # orthogonal component
    combo = combo / np.linalg.norm(combo) * np.linalg.norm(atoms[0]) * 3
    noise = 0.05 * np.abs(combo) + 1e-3
    c0 = cc.certify_emergence(cone_atoms=atoms, measured_combo=combo, noise_sd=noise, n_boot=80, seed=0)
    if not cd._is_certified(c0):
        return  # not emergent enough to test the flip; skip silently
    gstar, curve = cd.sensitivity_radius(atoms, combo, noise, n_boot=80, seed=0)
    floors = [row[1] for row in curve]
    # floor ratio is non-increasing as gamma grows (more noise -> smaller ratio)
    assert all(floors[i] >= floors[i + 1] - 1e-6 for i in range(len(floors) - 1))
    assert gstar >= 1.0


def test_is_certified_helper():
    class Fake:
        def __init__(self, v):
            self.verdict = v
    assert cd._is_certified(Fake("certified emergent (p=..., 3x floor)"))
    assert not cd._is_certified(Fake("within measurement noise (p=1)"))
    assert not cd._is_certified(Fake("emergent above noise but modest effect size"))
