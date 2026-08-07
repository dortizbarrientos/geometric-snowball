import numpy as np

from geometric_snowball.linalg import (
    chordal_distance_sq,
    metric_principal_angles,
    principal_angles,
    projector_distance_sq,
    spd_inverse_sqrt,
    spd_sqrt,
    volume_overlap_deficit,
)


def random_spd(seed: int, n: int = 4) -> np.ndarray:
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, n))
    return a @ a.T + 0.5 * np.eye(n)


def test_spd_square_root_and_whitening() -> None:
    g = random_spd(1)
    root = spd_sqrt(g)
    inverse = spd_inverse_sqrt(g)
    np.testing.assert_allclose(root @ root, g, atol=1e-11)
    np.testing.assert_allclose(inverse @ g @ inverse, np.eye(4), atol=1e-11)


def test_principal_angle_projector_identity() -> None:
    rng = np.random.default_rng(2)
    q1, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    q2, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    angles = principal_angles(q1, q2)
    np.testing.assert_allclose(chordal_distance_sq(angles), projector_distance_sq(q1, q2), atol=1e-12)


def test_volume_overlap_endpoints() -> None:
    assert volume_overlap_deficit([0.0, 0.0]) == 0.0
    np.testing.assert_allclose(volume_overlap_deficit([np.pi / 2]), 1.0)


def test_common_metric_angles_are_coordinate_invariant() -> None:
    rng = np.random.default_rng(3)
    u1 = rng.normal(size=(4, 2))
    u2 = rng.normal(size=(4, 2))
    m = random_spd(4)
    baseline = metric_principal_angles(u1, u2, m)
    transform = rng.normal(size=(4, 4))
    while abs(np.linalg.det(transform)) < 0.2:
        transform = rng.normal(size=(4, 4))
    inv = np.linalg.inv(transform)
    changed = metric_principal_angles(transform @ u1, transform @ u2, inv.T @ m @ inv)
    np.testing.assert_allclose(baseline, changed, atol=1e-10)
