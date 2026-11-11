"""Standing-variation dynamics and an explicitly phenomenological DMI handoff."""
from __future__ import annotations

import math
import numpy as np
from numpy.typing import ArrayLike
import pandas as pd

from .channels import channel_basis, compare_channels
from .linalg import Array, as_symmetric


def trait_mean(frequencies: ArrayLike, effects: ArrayLike, *, reference: ArrayLike | None = None) -> Array:
    """Diploid additive trait mean 2 sum_l (p_l-p_l0) a_l."""
    p = np.asarray(frequencies, dtype=float)
    a = np.asarray(effects, dtype=float)
    p0 = np.zeros_like(p) if reference is None else np.asarray(reference, dtype=float)
    if a.ndim != 2 or a.shape[0] != p.size or p0.shape != p.shape:
        raise ValueError("frequency and effect dimensions do not agree")
    return 2.0 * np.sum((p - p0)[:, None] * a, axis=0)


def additive_covariance(frequencies: ArrayLike, effects: ArrayLike) -> Array:
    """Linkage-equilibrium additive covariance 2 sum p_l(1-p_l) a_l a_l^T."""
    p = np.asarray(frequencies, dtype=float)
    a = np.asarray(effects, dtype=float)
    if a.ndim != 2 or a.shape[0] != p.size:
        raise ValueError("frequency and effect dimensions do not agree")
    weights = 2.0 * p * (1.0 - p)
    return as_symmetric(np.einsum("l,li,lj->ij", weights, a, a))


def allele_frequency_velocity(
    focal: ArrayLike,
    other: ArrayLike,
    effects: ArrayLike,
    beta: ArrayLike,
    *,
    migration: float = 0.0,
) -> Array:
    """Weak-selection allele-frequency dynamics with symmetric migration."""
    p = np.asarray(focal, dtype=float)
    q = np.asarray(other, dtype=float)
    a = np.asarray(effects, dtype=float)
    b = np.asarray(beta, dtype=float)
    if p.shape != q.shape or a.shape != (p.size, b.size):
        raise ValueError("incompatible dimensions")
    if migration < 0:
        raise ValueError("migration must be nonnegative")
    return p * (1.0 - p) * (a @ b) + migration * (q - p)


def mean_response_identity(
    focal: ArrayLike,
    other: ArrayLike,
    effects: ArrayLike,
    beta: ArrayLike,
    *,
    migration: float,
) -> tuple[Array, Array]:
    """Return direct and matrix forms of dz/dt for identity checks."""
    dp = allele_frequency_velocity(focal, other, effects, beta, migration=migration)
    direct = 2.0 * np.sum(dp[:, None] * np.asarray(effects, float), axis=0)
    g = additive_covariance(focal, effects)
    zi = trait_mean(focal, effects)
    zj = trait_mean(other, effects)
    formula = g @ np.asarray(beta, float) + migration * (zj - zi)
    return direct, formula


def simulate_standing_variation(
    *,
    generations: float = 80.0,
    dt: float = 0.1,
    loci: int = 160,
    seed: int = 20260807,
    migration: float = 0.002,
) -> pd.DataFrame:
    """Illustrative deterministic channel rotation using standing variation only."""
    if generations <= 0 or dt <= 0 or loci < 4:
        raise ValueError("invalid simulation settings")
    rng = np.random.default_rng(seed)
    effects = rng.normal(scale=0.075, size=(loci, 3))
    base = rng.uniform(0.12, 0.88, size=loci)
    p1 = base.copy()
    p2 = base.copy()
    theta1 = np.array([1.0, -0.55, 0.25])
    theta2 = np.array([-0.4, 0.95, -0.45])
    h1 = np.array([[1.25, 0.22, 0.0], [0.22, 0.70, 0.08], [0.0, 0.08, 0.45]])
    # Both lineages begin with the same response geometry. Divergent optima then
    # sort the same standing allelic effects differently, allowing G and the
    # resulting channels to diverge without mutation or drift.
    h2 = h1.copy()
    metric = np.linalg.inv(additive_covariance(base, effects) + 1e-8 * np.eye(3))
    rows: list[dict[str, float]] = []
    steps = int(round(generations / dt))
    for step in range(steps + 1):
        time = step * dt
        g1 = additive_covariance(p1, effects) + 1e-8 * np.eye(3)
        g2 = additive_covariance(p2, effects) + 1e-8 * np.eye(3)
        z1, z2 = trait_mean(p1, effects), trait_mean(p2, effects)
        beta1 = h1 @ (theta1 - z1)
        beta2 = h2 @ (theta2 - z2)
        angles = compare_channels(g1, h1, g2, h2, k=1, metric=metric)
        c1, c2 = channel_basis(g1, h1), channel_basis(g2, h2)
        rows.append(
            {
                "time": time,
                "channel_angle_deg": float(np.degrees(angles[0])),
                "gap_1": float(c1.eigenvalues[0] - c1.eigenvalues[1]),
                "gap_2": float(c2.eigenvalues[0] - c2.eigenvalues[1]),
                "mean_distance": float(np.linalg.norm(z2 - z1)),
            }
        )
        if step == steps:
            break
        dp1 = allele_frequency_velocity(p1, p2, effects, beta1, migration=migration)
        dp2 = allele_frequency_velocity(p2, p1, effects, beta2, migration=migration)
        p1 = np.clip(p1 + dt * dp1, 1e-5, 1.0 - 1e-5)
        p2 = np.clip(p2 + dt * dp2, 1e-5, 1.0 - 1e-5)
    return pd.DataFrame(rows)


def effective_migration(m0: float, ecological_log_load: float) -> float:
    """Phenomenological reduction m_eff=m0 exp(-L_E)."""
    if m0 < 0 or ecological_log_load < 0:
        raise ValueError("migration and load must be nonnegative")
    return float(m0 * math.exp(-ecological_log_load))


def coupled_isolation(ecological_load: float, dmi_count: float, *, eta: float = 1.0) -> float:
    """Combined multiplicative barrier 1-exp[-(L_E+eta D)]."""
    if min(ecological_load, dmi_count, eta) < 0:
        raise ValueError("arguments must be nonnegative")
    return float(1.0 - math.exp(-(ecological_load + eta * dmi_count)))


def simulate_dmi_handoff(
    *,
    duration: float = 10.0,
    dt: float = 0.01,
    m0: float = 0.2,
    geometric_strength: float = 0.055,
    substitution_rate: float = 0.34,
    selection: float = 0.12,
    incompatibility_probability: float = 0.025,
    eta: float = 0.75,
) -> pd.DataFrame:
    """Illustrative early ecological barrier coupled to later DMI accumulation.

    This is deliberately labelled phenomenological.  It demonstrates a possible
    handoff; it is not a derivation of DMI dynamics from S.
    """
    rows: list[dict[str, float | str]] = []
    for scenario, strength in [("no_geometric_precursor", 0.0), ("geometric_precursor", geometric_strength)]:
        k1 = k2 = 0.0
        steps = int(round(duration / dt))
        for step in range(steps + 1):
            time = step * dt
            le = strength * time**2
            meff = effective_migration(m0, le)
            establishment = selection / (selection + 2.0 * meff)
            dmi = incompatibility_probability * k1 * k2
            rows.append(
                {
                    "scenario": scenario,
                    "time": time,
                    "ecological_log_load": le,
                    "effective_migration": meff,
                    "substitutions_lineage_1": k1,
                    "substitutions_lineage_2": k2,
                    "dmi_count": dmi,
                    "intrinsic_isolation": 1.0 - math.exp(-eta * dmi),
                    "total_isolation": coupled_isolation(le, dmi, eta=eta),
                }
            )
            if step < steps:
                rate = substitution_rate * establishment
                k1 += dt * rate
                k2 += dt * rate
    return pd.DataFrame(rows)
