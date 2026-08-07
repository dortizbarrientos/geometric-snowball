import numpy as np
import pytest

from geometric_snowball.channels import (
    channel_basis,
    channel_operator,
    channel_operator_derivative,
    compare_channels,
    davis_kahan_bound,
    local_chordal_coefficient,
    sqrt_derivative,
)
from geometric_snowball.linalg import chordal_distance_sq, principal_angles, spd_sqrt


def test_spectrum_matches_GH() -> None:
    g = np.array([[1.2, 0.3], [0.3, 0.7]])
    h = np.array([[0.8, -0.1], [-0.1, 1.4]])
    s = channel_operator(g, h)
    np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(s)), np.sort(np.linalg.eigvals(g @ h).real), atol=1e-12)


def test_schluter_limit_recovers_gmax() -> None:
    g = np.array([[1.2, 0.4], [0.4, 0.6]])
    values, vectors = np.linalg.eigh(g)
    gmax = vectors[:, np.argmax(values)]
    channel = channel_basis(g, 2.0 * np.eye(2)).phenotype_basis[:, 0]
    cosine = abs(float(gmax @ channel) / np.linalg.norm(channel))
    np.testing.assert_allclose(cosine, 1.0, atol=1e-12)


def test_unequal_G_requires_common_metric() -> None:
    g1 = np.diag([1.0, 0.5])
    g2 = np.diag([0.8, 0.7])
    with pytest.raises(ValueError, match="common"):
        compare_channels(g1, np.eye(2), g2, np.eye(2))
    angles = compare_channels(g1, np.eye(2), g2, np.eye(2), metric=np.eye(2))
    assert angles.shape == (1,)


def test_local_quadratic_coefficient() -> None:
    s0 = np.diag([3.0, 1.7, 0.8, 0.3])
    rng = np.random.default_rng(9)
    a1 = rng.normal(size=(4, 4)); a1 = 0.5 * (a1 + a1.T)
    a2 = rng.normal(size=(4, 4)); a2 = 0.5 * (a2 + a2.T)
    coefficient = local_chordal_coefficient(s0, a1, a2, k=2)
    eps = 1e-5
    _, v1 = np.linalg.eigh(s0 + eps * a1)
    _, v2 = np.linalg.eigh(s0 + eps * a2)
    exact = chordal_distance_sq(principal_angles(v1[:, -2:], v2[:, -2:])) / eps**2
    np.testing.assert_allclose(exact, coefficient, rtol=3e-4)


def test_commuting_perturbation_does_not_rotate() -> None:
    s0 = np.diag([3.0, 1.7, 0.8])
    a1 = np.diag([0.1, -0.2, 0.05])
    a2 = np.diag([-0.1, 0.15, -0.02])
    np.testing.assert_allclose(local_chordal_coefficient(s0, a1, a2, k=1), 0.0, atol=1e-15)
    assert davis_kahan_bound(s0, s0 + a2, k=1) > 0.0


def test_square_root_derivative() -> None:
    g = np.array([[1.3, 0.2], [0.2, 0.8]])
    dg = np.array([[0.1, -0.04], [-0.04, 0.06]])
    x = sqrt_derivative(g, dg)
    root = spd_sqrt(g)
    np.testing.assert_allclose(root @ x + x @ root, dg, atol=1e-12)
    eps = 1e-6
    np.testing.assert_allclose((spd_sqrt(g + eps * dg) - root) / eps, x, rtol=2e-6, atol=2e-7)


def test_channel_operator_derivative() -> None:
    g = np.array([[1.1, 0.1], [0.1, 0.7]])
    h = np.array([[0.8, -0.05], [-0.05, 1.2]])
    dg = np.array([[0.05, 0.01], [0.01, -0.02]])
    dh = np.array([[-0.02, 0.03], [0.03, 0.04]])
    analytic = channel_operator_derivative(g, h, dg, dh)
    eps = 1e-6
    finite = (channel_operator(g + eps * dg, h + eps * dh) - channel_operator(g, h)) / eps
    np.testing.assert_allclose(analytic, finite, rtol=3e-6, atol=3e-7)
