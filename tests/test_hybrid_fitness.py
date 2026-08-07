import numpy as np

from geometric_snowball.hybrid_fitness import (
    covariance_cross_load,
    direct_covariance_benchmark_load,
    expected_negative_log_fitness,
    gaussian_mean_fitness,
    ridge_covariance_load,
)
from geometric_snowball.linalg import principal_angles


def test_expected_quadratic_identity_by_monte_carlo() -> None:
    rng = np.random.default_rng(11)
    mu = np.array([0.2, -0.1])
    theta = np.array([-0.1, 0.05])
    p = np.array([[0.3, 0.05], [0.05, 0.2]])
    h = np.array([[0.8, -0.1], [-0.1, 1.1]])
    exact = expected_negative_log_fitness(mu, p, theta, h)
    z = rng.multivariate_normal(mu, p, size=250_000)
    d = z - theta
    monte_carlo = 0.5 * np.mean(np.einsum("bi,ij,bj->b", d, h, d))
    np.testing.assert_allclose(monte_carlo, exact, rtol=7e-3)


def test_gaussian_mean_fitness_by_monte_carlo() -> None:
    rng = np.random.default_rng(12)
    mu = np.array([0.2, -0.1])
    theta = np.array([-0.1, 0.05])
    p = np.array([[0.2, 0.03], [0.03, 0.15]])
    h = np.array([[0.7, -0.08], [-0.08, 0.9]])
    exact = gaussian_mean_fitness(mu, p, theta, h)
    z = rng.multivariate_normal(mu, p, size=300_000)
    d = z - theta
    monte_carlo = np.mean(np.exp(-0.5 * np.einsum("bi,ij,bj->b", d, h, d)))
    np.testing.assert_allclose(monte_carlo, exact, rtol=5e-3)


def test_covariance_curvature_trace_identity() -> None:
    rng = np.random.default_rng(13)
    a = rng.normal(size=(3, 3)); g1 = a @ a.T + np.eye(3)
    a = rng.normal(size=(3, 3)); g2 = a @ a.T + np.eye(3)
    a = rng.normal(size=(3, 3)); h1 = a @ a.T + np.eye(3)
    a = rng.normal(size=(3, 3)); h2 = a @ a.T + np.eye(3)
    r = 0.02 * rng.normal(size=(3, 3)); r = 0.5 * (r + r.T)
    hybrid = 0.5 * (g1 + g2) + r
    direct = direct_covariance_benchmark_load(g1, g2, h1, h2, hybrid)
    compact = covariance_cross_load(g1, g2, h1, h2, segregation_excess=r)
    np.testing.assert_allclose(direct, compact, atol=1e-11)


def test_ridge_theorem() -> None:
    rng = np.random.default_rng(14)
    q1, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    q2, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    p1, p2 = q1 @ q1.T, q2 @ q2.T
    delta_g, delta_h = 1.5, 1.1
    g1 = 0.3 * np.eye(5) + delta_g * p1
    g2 = 0.3 * np.eye(5) + delta_g * p2
    h1 = 1.6 * np.eye(5) - delta_h * p1
    h2 = 1.6 * np.eye(5) - delta_h * p2
    exact = covariance_cross_load(g1, g2, h1, h2)
    predicted = ridge_covariance_load(delta_g, delta_h, principal_angles(q1, q2))
    np.testing.assert_allclose(exact, predicted, atol=1e-11)


def test_angle_without_covariance_curvature_coupling_need_not_cost_fitness() -> None:
    g = np.diag([1.2, 0.5])
    h1 = np.diag([0.5, 1.3])
    h2 = np.diag([1.3, 0.5])
    np.testing.assert_allclose(covariance_cross_load(g, g, h1, h2), 0.0)
