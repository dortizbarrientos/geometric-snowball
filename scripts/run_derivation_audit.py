#!/usr/bin/env python3
"""Transparent numerical audit of the main identities in the theory."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from geometric_snowball.channels import (  # noqa: E402
    channel_basis,
    channel_operator,
    channel_operator_derivative,
    local_chordal_coefficient,
    sqrt_derivative,
)
from geometric_snowball.dynamics import mean_response_identity  # noqa: E402
from geometric_snowball.hybrid_fitness import (  # noqa: E402
    covariance_cross_load,
    direct_covariance_benchmark_load,
    expected_negative_log_fitness,
    gaussian_mean_fitness,
    ridge_covariance_load,
)
from geometric_snowball.linalg import (  # noqa: E402
    chordal_distance_sq,
    metric_principal_angles,
    principal_angles,
    projector_distance_sq,
    spd_inverse_sqrt,
    spd_sqrt,
    volume_overlap_deficit,
)


@dataclass
class Check:
    name: str
    error: float
    tolerance: float
    explanation: str

    @property
    def passed(self) -> bool:
        return bool(abs(self.error) <= self.tolerance)


def random_spd(rng: np.random.Generator, n: int, floor: float = 0.4) -> np.ndarray:
    a = rng.normal(size=(n, n))
    return 0.5 * (a @ a.T + (a @ a.T).T) + floor * np.eye(n)


def run_checks() -> list[Check]:
    rng = np.random.default_rng(20260807)
    checks: list[Check] = []

    g = random_spd(rng, 5)
    h = random_spd(rng, 5)
    root = spd_sqrt(g)
    inverse = spd_inverse_sqrt(g)
    s = channel_operator(g, h)
    checks.append(Check("square_root_identity", np.linalg.norm(root @ root - g), 1e-10, "G^(1/2)G^(1/2)=G"))
    checks.append(Check("whitening_identity", np.linalg.norm(inverse @ g @ inverse - np.eye(5)), 1e-10, "G^(-1/2)GG^(-1/2)=I"))
    checks.append(Check("spectrum_GH_equals_S", np.max(np.abs(np.sort(np.linalg.eigvals(g @ h).real) - np.sort(np.linalg.eigvalsh(s)))), 1e-9, "GH and S are similar and share eigenvalues"))

    gmax = np.linalg.eigh(g)[1][:, -1]
    response = channel_basis(g, 1.7 * np.eye(5)).phenotype_basis[:, 0]
    cosine = abs(float(gmax @ response) / np.linalg.norm(response))
    checks.append(Check("Schluter_limit", 1.0 - cosine, 1e-11, "H proportional to I recovers g_max"))

    q1, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    q2, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    angles = principal_angles(q1, q2)
    checks.append(Check("projector_chordal_identity", chordal_distance_sq(angles) - projector_distance_sq(q1, q2), 1e-11, "sum sin^2 psi equals half projector distance"))
    volume_direct = 1.0 - np.linalg.det(q1.T @ q2) ** 2
    checks.append(Check("volume_overlap_identity", volume_overlap_deficit(angles) - volume_direct, 1e-11, "determinant and principal-angle expressions agree"))

    g1, g2 = random_spd(rng, 4), random_spd(rng, 4)
    h1, h2 = random_spd(rng, 4), random_spd(rng, 4)
    u1 = channel_basis(g1, h1, k=2).phenotype_basis
    u2 = channel_basis(g2, h2, k=2).phenotype_basis
    metric = np.linalg.inv(0.5 * (g1 + g2))
    baseline = metric_principal_angles(u1, u2, metric)
    transform = rng.normal(size=(4, 4))
    while abs(np.linalg.det(transform)) < 0.2:
        transform = rng.normal(size=(4, 4))
    inv = np.linalg.inv(transform)
    g1p, g2p = transform @ g1 @ transform.T, transform @ g2 @ transform.T
    h1p, h2p = inv.T @ h1 @ inv, inv.T @ h2 @ inv
    metricp = inv.T @ metric @ inv
    changed = metric_principal_angles(channel_basis(g1p, h1p, k=2).phenotype_basis, channel_basis(g2p, h2p, k=2).phenotype_basis, metricp)
    checks.append(Check("common_metric_coordinate_invariance", np.max(np.abs(baseline - changed)), 1e-9, "angles do not depend on nonsingular trait reparameterisation"))

    s0 = np.diag([3.2, 1.7, 0.8, 0.35])
    a1 = rng.normal(size=(4, 4)); a1 = 0.5 * (a1 + a1.T)
    a2 = rng.normal(size=(4, 4)); a2 = 0.5 * (a2 + a2.T)
    coefficient = local_chordal_coefficient(s0, a1, a2, k=2)
    eps = 1e-5
    q1e = np.linalg.eigh(s0 + eps * a1)[1][:, -2:]
    q2e = np.linalg.eigh(s0 + eps * a2)[1][:, -2:]
    finite = chordal_distance_sq(principal_angles(q1e, q2e)) / eps**2
    checks.append(Check("local_quadratic_coefficient", abs(finite - coefficient) / coefficient, 4e-4, "finite difference verifies ||K||_F^2"))
    checks.append(Check("commuting_rotation_zero", local_chordal_coefficient(s0, np.diag([0.1, 0.0, -0.1, 0.02]), np.diag([-0.1, 0.2, 0.03, -0.04]), k=2), 1e-14, "cross-channel forcing is required"))

    dg = 0.05 * (rng.normal(size=(5, 5))); dg = 0.5 * (dg + dg.T)
    dh = 0.05 * (rng.normal(size=(5, 5))); dh = 0.5 * (dh + dh.T)
    x = sqrt_derivative(g, dg)
    checks.append(Check("sqrt_derivative_Sylvester", np.linalg.norm(root @ x + x @ root - dg), 1e-10, "R X + X R = dG"))
    eps = 1e-6
    operator_fd = (channel_operator(g + eps * dg, h + eps * dh) - s) / eps
    operator_exact = channel_operator_derivative(g, h, dg, dh)
    checks.append(Check("operator_derivative", np.linalg.norm(operator_fd - operator_exact) / np.linalg.norm(operator_exact), 4e-6, "finite difference verifies dS"))

    mu, theta = rng.normal(scale=0.2, size=4), rng.normal(scale=0.2, size=4)
    p, hh = 0.18 * random_spd(rng, 4), 0.25 * random_spd(rng, 4)
    z = rng.multivariate_normal(mu, p, size=180_000)
    d = z - theta
    mc_log = 0.5 * np.mean(np.einsum("bi,ij,bj->b", d, hh, d))
    exact_log = expected_negative_log_fitness(mu, p, theta, hh)
    checks.append(Check("expected_quadratic_form", abs(mc_log - exact_log) / exact_log, 9e-3, "Monte Carlo verifies mean quadratic load"))
    mc_w = np.mean(np.exp(-0.5 * np.einsum("bi,ij,bj->b", d, hh, d)))
    exact_w = gaussian_mean_fitness(mu, p, theta, hh)
    checks.append(Check("Gaussian_mean_fitness", abs(mc_w - exact_w) / exact_w, 1e-2, "Monte Carlo verifies exact arithmetic mean fitness"))

    g1, g2 = random_spd(rng, 5), random_spd(rng, 5)
    h1, h2 = random_spd(rng, 5), random_spd(rng, 5)
    r = 0.02 * rng.normal(size=(5, 5)); r = 0.5 * (r + r.T)
    ph = 0.5 * (g1 + g2) + r
    direct = direct_covariance_benchmark_load(g1, g2, h1, h2, ph)
    compact = covariance_cross_load(g1, g2, h1, h2, segregation_excess=r)
    checks.append(Check("covariance_curvature_trace_identity", direct - compact, 1e-10, "direct and compact hybrid-load forms agree"))

    q1, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    q2, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    p1, p2 = q1 @ q1.T, q2 @ q2.T
    dg_ridge, dh_ridge = 1.5, 1.2
    rg1, rg2 = 0.3 * np.eye(6) + dg_ridge * p1, 0.3 * np.eye(6) + dg_ridge * p2
    rh1, rh2 = 1.7 * np.eye(6) - dh_ridge * p1, 1.7 * np.eye(6) - dh_ridge * p2
    exact = covariance_cross_load(rg1, rg2, rh1, rh2)
    predicted = ridge_covariance_load(dg_ridge, dh_ridge, principal_angles(q1, q2))
    checks.append(Check("ridge_theorem", exact - predicted, 1e-10, "curvature-weighted load equals angle formula under ridge assumptions"))

    effects = rng.normal(scale=0.05, size=(250, 4))
    pi = rng.uniform(0.08, 0.92, size=250)
    pj = rng.uniform(0.08, 0.92, size=250)
    beta = rng.normal(size=4)
    direct, matrix = mean_response_identity(pi, pj, effects, beta, migration=0.025)
    checks.append(Check("standing_variation_response", np.linalg.norm(direct - matrix), 1e-11, "allele-frequency and G beta plus migration forms agree"))

    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=ROOT / "build" / "audit")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    checks = run_checks()
    with (args.outdir / "derivation_audit.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["check", "error", "tolerance", "passed", "explanation"])
        for check in checks:
            writer.writerow([check.name, f"{check.error:.16e}", f"{check.tolerance:.4e}", check.passed, check.explanation])
    lines = ["Geometric snowball derivation audit", "=" * 37, ""]
    for check in checks:
        lines.append(f"{'PASS' if check.passed else 'FAIL'}  {check.name}: error={check.error:.4e}; tolerance={check.tolerance:.4e}")
        lines.append(f"      {check.explanation}")
    lines.extend(["", f"Overall: {sum(c.passed for c in checks)}/{len(checks)} checks passed."])
    report = "\n".join(lines) + "\n"
    (args.outdir / "derivation_audit.txt").write_text(report)
    print(report)
    if not all(check.passed for check in checks):
        raise SystemExit("one or more derivation checks failed")


if __name__ == "__main__":
    main()
