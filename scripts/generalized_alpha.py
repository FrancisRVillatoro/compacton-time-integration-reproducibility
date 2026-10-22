"""Canonical first-order generalized-alpha amplification matrix.

This is the corrected 2 x 2 matrix used in the manuscript's spectral plots.
It is independent of the nonlinear trajectory solver used in the production
experiments.
"""
from __future__ import annotations

import numpy as np


def ga_matrix(z: complex, rho_inf: float) -> np.ndarray:
    r = float(rho_inf)
    hden = 3.0 + 2.0 * r - r * r - 2.0 * z
    return np.array(
        [
            [
                (3.0 + 2.0 * r - r * r + 2.0 * r * z) / hden,
                (1.0 - r * r) / hden,
            ],
            [
                2.0 * z * (1.0 + r) ** 2 / hden,
                (1.0 - 2.0 * r - 3.0 * r * r + 2.0 * r * z) / hden,
            ],
        ],
        dtype=complex,
    )


def ga_rho(z: complex, rho_inf: float) -> float:
    return float(np.max(np.abs(np.linalg.eigvals(ga_matrix(z, rho_inf)))))
