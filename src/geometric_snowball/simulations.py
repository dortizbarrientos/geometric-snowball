"""Small reproducible experiments illustrating the corrected theory."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .channels import channel_basis, local_chordal_coefficient
from .dynamics import simulate_dmi_handoff
from .hybrid_fitness import covariance_cross_load, ridge_covariance_load
from .linalg import chordal_distance_sq, principal_angles


def local_quadratic_experiment() -> pd.DataFrame:
    """Verify the deterministic t^2 channel-distance law."""
    s0 = np.diag([3.0, 1.65, 0.8, 0.35])
    a1 = np.array([[0.0, 0.12, 0.03, 0.0], [0.12, 0.0, 0.04, 0.02], [0.03, 0.04, 0.0, 0.01], [0.0, 0.02, 0.01, 0.0]])
    a2 = np.array([[0.0, -0.09, 0.07, 0.02], [-0.09, 0.0, -0.02, 0.04], [0.07, -0.02, 0.0, -0.03], [0.02, 0.04, -0.03, 0.0]])
    coefficient = local_chordal_coefficient(s0, a1, a2, k=2)
    rows = []
    for time in np.geomspace(1e-4, 0.12, 80):
        _, v1 = np.linalg.eigh(s0 + time * a1)
        _, v2 = np.linalg.eigh(s0 + time * a2)
        q1, q2 = v1[:, -2:], v2[:, -2:]
        exact = chordal_distance_sq(principal_angles(q1, q2))
        rows.append({"time": time, "exact_chordal_distance_sq": exact, "quadratic_prediction": coefficient * time**2})
    return pd.DataFrame(rows)


def ridge_load_experiment() -> pd.DataFrame:
    """Verify the exact ridge angle-to-load calibration."""
    n, k = 4, 2
    q1 = np.eye(n)[:, :k]
    delta_g, delta_h = 1.4, 1.1
    g_perp, h_parallel = 0.3, 0.45
    rows = []
    for degrees in np.linspace(0.0, 80.0, 41):
        a = np.deg2rad(degrees)
        r1 = np.array([[np.cos(a), 0.0], [0.0, np.cos(0.7 * a)], [np.sin(a), 0.0], [0.0, np.sin(0.7 * a)]])
        q2, _ = np.linalg.qr(r1)
        p1, p2 = q1 @ q1.T, q2[:, :k] @ q2[:, :k].T
        g1 = g_perp * np.eye(n) + delta_g * p1
        g2 = g_perp * np.eye(n) + delta_g * p2
        h_perp = h_parallel + delta_h
        h1 = h_perp * np.eye(n) - delta_h * p1
        h2 = h_perp * np.eye(n) - delta_h * p2
        angles = principal_angles(q1, q2[:, :k])
        rows.append(
            {
                "angle_deg": degrees,
                "exact_covariance_load": covariance_cross_load(g1, g2, h1, h2),
                "ridge_theorem_prediction": ridge_covariance_load(delta_g, delta_h, angles),
            }
        )
    return pd.DataFrame(rows)


def selection_vs_drift_experiment(*, replicates: int = 2500, seed: int = 614) -> pd.DataFrame:
    """Contrast deterministic O(t^2) with unbiased diffusive O(t) squared displacement."""
    rng = np.random.default_rng(seed)
    times = np.geomspace(1e-3, 1.0, 80)
    deterministic_rate = 0.21
    diffusion = 0.016
    rows = []
    for t in times:
        deterministic = (deterministic_rate * t) ** 2
        increments = rng.normal(scale=np.sqrt(2.0 * diffusion * t), size=(replicates, 2))
        stochastic = float(np.mean(np.sum(increments**2, axis=1)))
        rows.extend(
            [
                {"process": "deterministic selection", "time": t, "mean_squared_channel_distance": deterministic},
                {"process": "unbiased diffusion", "time": t, "mean_squared_channel_distance": stochastic},
            ]
        )
    return pd.DataFrame(rows)


def dmi_handoff_experiment() -> pd.DataFrame:
    return simulate_dmi_handoff()
