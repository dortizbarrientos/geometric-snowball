"""Linear-algebra foundations for the geometric snowball framework."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.float64]


def as_symmetric(a: ArrayLike) -> Array:
    """Return the symmetric part of a square matrix."""
    x = np.asarray(a, dtype=float)
    if x.ndim != 2 or x.shape[0] != x.shape[1]:
        raise ValueError("expected a square matrix")
    return 0.5 * (x + x.T)


def assert_spd(a: ArrayLike, *, tol: float = 1e-12) -> Array:
    """Validate and return a symmetric positive-definite matrix."""
    x = as_symmetric(a)
    smallest = float(np.linalg.eigvalsh(x).min())
    if smallest <= tol:
        raise ValueError(f"matrix is not positive definite; minimum eigenvalue={smallest:.3e}")
    return x


def eigh_descending(a: ArrayLike) -> tuple[Array, Array]:
    """Symmetric eigendecomposition ordered from largest to smallest eigenvalue."""
    values, vectors = np.linalg.eigh(as_symmetric(a))
    order = np.argsort(values)[::-1]
    return values[order], vectors[:, order]


def spd_sqrt(a: ArrayLike) -> Array:
    """Unique symmetric positive-definite square root."""
    x = assert_spd(a)
    values, vectors = np.linalg.eigh(x)
    return as_symmetric(vectors @ np.diag(np.sqrt(values)) @ vectors.T)


def spd_inverse_sqrt(a: ArrayLike) -> Array:
    """Inverse of the symmetric positive-definite square root."""
    x = assert_spd(a)
    values, vectors = np.linalg.eigh(x)
    return as_symmetric(vectors @ np.diag(values ** -0.5) @ vectors.T)


def orthonormalize(columns: ArrayLike) -> Array:
    """Euclidean-orthonormal basis for the supplied column space."""
    x = np.asarray(columns, dtype=float)
    if x.ndim != 2 or x.shape[0] < x.shape[1]:
        raise ValueError("columns must form a tall or square matrix")
    q, r = np.linalg.qr(x)
    signs = np.where(np.diag(r) < 0.0, -1.0, 1.0)
    return q * signs


def metric_orthonormalize(columns: ArrayLike, metric: ArrayLike) -> Array:
    """Return Q spanning columns with Q.T @ metric @ Q = I."""
    u = np.asarray(columns, dtype=float)
    m = assert_spd(metric)
    if u.ndim != 2 or u.shape[0] != m.shape[0]:
        raise ValueError("basis and metric dimensions do not agree")
    gram = as_symmetric(u.T @ m @ u)
    return u @ spd_inverse_sqrt(gram)


def principal_angles(q1: ArrayLike, q2: ArrayLike) -> Array:
    """Principal angles between Euclidean-orthonormal column spaces."""
    a = np.asarray(q1, dtype=float)
    b = np.asarray(q2, dtype=float)
    if a.shape[0] != b.shape[0] or a.shape[1] != b.shape[1]:
        raise ValueError("bases must have the same ambient and subspace dimensions")
    singular = np.linalg.svd(a.T @ b, compute_uv=False)
    return np.sort(np.arccos(np.clip(singular, 0.0, 1.0)))


def metric_principal_angles(u1: ArrayLike, u2: ArrayLike, metric: ArrayLike) -> Array:
    """Coordinate-invariant principal angles under one declared shared metric."""
    m = assert_spd(metric)
    q1 = metric_orthonormalize(u1, m)
    q2 = metric_orthonormalize(u2, m)
    singular = np.linalg.svd(q1.T @ m @ q2, compute_uv=False)
    return np.sort(np.arccos(np.clip(singular, 0.0, 1.0)))


def chordal_distance_sq(angles: ArrayLike) -> float:
    """Squared chordal distance, sum_j sin^2(psi_j)."""
    psi = np.asarray(angles, dtype=float)
    return float(np.sum(np.sin(psi) ** 2))


def geodesic_distance(angles: ArrayLike) -> float:
    """Grassmann geodesic distance, sqrt(sum_j psi_j^2)."""
    psi = np.asarray(angles, dtype=float)
    return float(np.linalg.norm(psi))


def volume_overlap_deficit(angles: ArrayLike) -> float:
    """One minus squared subspace-volume overlap; a geometry, not a fitness measure."""
    psi = np.asarray(angles, dtype=float)
    return float(1.0 - np.prod(np.cos(psi) ** 2))


def projector_distance_sq(q1: ArrayLike, q2: ArrayLike) -> float:
    """Squared chordal distance expressed as half the projector Frobenius distance."""
    a = np.asarray(q1, dtype=float)
    b = np.asarray(q2, dtype=float)
    p1, p2 = a @ a.T, b @ b.T
    return float(0.5 * np.linalg.norm(p1 - p2, ord="fro") ** 2)


def orthogonal_complement(q: ArrayLike) -> Array:
    """Euclidean-orthonormal complement of an orthonormal basis."""
    basis = np.asarray(q, dtype=float)
    _, _, vh = np.linalg.svd(basis.T, full_matrices=True)
    return vh[basis.shape[1] :, :].T
