"""Explicit mappings from phenotype distributions to hybrid fitness."""
from __future__ import annotations

import math
import numpy as np
from numpy.typing import ArrayLike

from .linalg import Array, as_symmetric, assert_spd, spd_sqrt


def expected_negative_log_fitness(
    mean: ArrayLike,
    covariance: ArrayLike,
    optimum: ArrayLike,
    curvature: ArrayLike,
) -> float:
    """Half the expected quadratic displacement from a local optimum."""
    mu = np.asarray(mean, dtype=float)
    theta = np.asarray(optimum, dtype=float)
    p = as_symmetric(covariance)
    h = assert_spd(curvature)
    if mu.shape != theta.shape or p.shape != h.shape or p.shape[0] != mu.size:
        raise ValueError("incompatible phenotype dimensions")
    d = mu - theta
    return float(0.5 * (d @ h @ d + np.trace(h @ p)))


def gaussian_mean_fitness(
    mean: ArrayLike,
    covariance: ArrayLike,
    optimum: ArrayLike,
    curvature: ArrayLike,
    *,
    peak_fitness: float = 1.0,
) -> float:
    """Exact arithmetic mean fitness for Gaussian phenotypes and quadratic log fitness."""
    if peak_fitness <= 0:
        raise ValueError("peak_fitness must be positive")
    mu = np.asarray(mean, dtype=float)
    theta = np.asarray(optimum, dtype=float)
    p = as_symmetric(covariance)
    h = assert_spd(curvature)
    root = spd_sqrt(h)
    middle = np.eye(mu.size) + root @ p @ root
    d = root @ (mu - theta)
    sign, logdet = np.linalg.slogdet(middle)
    if sign <= 0:
        raise ValueError("I + H^(1/2) P H^(1/2) must be positive definite")
    exponent = -0.5 * d @ np.linalg.solve(middle, d)
    return float(peak_fitness * math.exp(-0.5 * logdet + exponent))


def covariance_cross_load(
    g1: ArrayLike,
    g2: ArrayLike,
    h1: ArrayLike,
    h2: ArrayLike,
    *,
    segregation_excess: ArrayLike | None = None,
) -> float:
    """Average covariance contribution to hybrid log-load across parental environments.

    P_H = (G1+G2)/2 + R.  The compact identity is
    -tr[(H2-H1)(G2-G1)]/8 + tr[(H1+H2)R]/4.
    """
    gg1, gg2 = as_symmetric(g1), as_symmetric(g2)
    hh1, hh2 = assert_spd(h1), assert_spd(h2)
    if not (gg1.shape == gg2.shape == hh1.shape == hh2.shape):
        raise ValueError("all matrices must have the same shape")
    r = np.zeros_like(gg1) if segregation_excess is None else as_symmetric(segregation_excess)
    return float(
        -0.125 * np.trace((hh2 - hh1) @ (gg2 - gg1))
        + 0.25 * np.trace((hh1 + hh2) @ r)
    )


def direct_covariance_benchmark_load(
    g1: ArrayLike,
    g2: ArrayLike,
    h1: ArrayLike,
    h2: ArrayLike,
    hybrid_covariance: ArrayLike,
) -> float:
    """Direct parental benchmark expression used to verify the compact trace identity."""
    gg1, gg2 = as_symmetric(g1), as_symmetric(g2)
    hh1, hh2 = assert_spd(h1), assert_spd(h2)
    p = as_symmetric(hybrid_covariance)
    return float(0.25 * (np.trace(hh1 @ (p - gg1)) + np.trace(hh2 @ (p - gg2))))


def ridge_covariance_load(delta_g: float, delta_h: float, angles: ArrayLike) -> float:
    """Exact channel-ridge bridge: (delta_g delta_h / 4) sum sin^2(psi_j)."""
    if delta_g < 0 or delta_h < 0:
        raise ValueError("ridge contrasts must be nonnegative")
    psi = np.asarray(angles, dtype=float)
    return float(0.25 * delta_g * delta_h * np.sum(np.sin(psi) ** 2))


def bounded_isolation_from_log_load(load: float) -> float:
    """Map a nonnegative additive log barrier to [0,1)."""
    if load < 0:
        raise ValueError("a bounded isolation index requires a nonnegative load")
    return float(1.0 - math.exp(-load))
