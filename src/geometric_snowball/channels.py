"""Additive-response channels and their perturbation geometry."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike

from .linalg import (
    Array,
    as_symmetric,
    assert_spd,
    eigh_descending,
    metric_principal_angles,
    principal_angles,
    spd_sqrt,
)


@dataclass(frozen=True)
class ChannelBasis:
    """Leading response subspace in whitened and phenotype coordinates."""

    operator: Array
    eigenvalues: Array
    whitened_basis: Array
    phenotype_basis: Array


def channel_operator(g: ArrayLike, h: ArrayLike) -> Array:
    """Symmetric response operator S = G^(1/2) H G^(1/2)."""
    gg = assert_spd(g)
    hh = assert_spd(h)
    if gg.shape != hh.shape:
        raise ValueError("G and H must have the same shape")
    root = spd_sqrt(gg)
    return as_symmetric(root @ hh @ root)


def channel_basis(g: ArrayLike, h: ArrayLike, *, k: int = 1) -> ChannelBasis:
    """Leading k-dimensional channel of S, with phenotype map U=G^(1/2)V."""
    gg = assert_spd(g)
    n = gg.shape[0]
    if not 1 <= k < n:
        raise ValueError("k must satisfy 1 <= k < number of traits")
    operator = channel_operator(gg, h)
    values, vectors = eigh_descending(operator)
    v = vectors[:, :k]
    u = spd_sqrt(gg) @ v
    return ChannelBasis(operator, values, v, u)


def compare_channels(
    g1: ArrayLike,
    h1: ArrayLike,
    g2: ArrayLike,
    h2: ArrayLike,
    *,
    k: int = 1,
    metric: ArrayLike | None = None,
    equality_tol: float = 1e-10,
) -> Array:
    """Compare channels without mixing incompatible population-specific whitenings.

    If G1 and G2 are equal, their shared whitened coordinates are used.  If they
    differ, a common phenotype-space metric must be supplied explicitly.
    """
    gg1 = assert_spd(g1)
    gg2 = assert_spd(g2)
    c1 = channel_basis(gg1, h1, k=k)
    c2 = channel_basis(gg2, h2, k=k)
    if np.allclose(gg1, gg2, rtol=equality_tol, atol=equality_tol):
        return principal_angles(c1.whitened_basis, c2.whitened_basis)
    if metric is None:
        raise ValueError("unequal G matrices require a declared common phenotype-space metric")
    return metric_principal_angles(c1.phenotype_basis, c2.phenotype_basis, metric)


def local_rotation_matrix(s0: ArrayLike, a1: ArrayLike, a2: ArrayLike, *, k: int) -> Array:
    """First-order relative channel-rotation matrix K for S_i(t)=S0+t A_i+O(t^2)."""
    base = assert_spd(s0)
    values, vectors = eigh_descending(base)
    n = base.shape[0]
    if not 1 <= k < n:
        raise ValueError("k must satisfy 1 <= k < n")
    v, w = vectors[:, :k], vectors[:, k:]
    forcing = as_symmetric(np.asarray(a2, float) - np.asarray(a1, float))
    out = np.empty((n - k, k), dtype=float)
    for a in range(n - k):
        for b in range(k):
            gap = values[b] - values[k + a]
            if abs(gap) <= 1e-12:
                raise ValueError("selected channel is not spectrally separated")
            out[a, b] = w[:, a] @ forcing @ v[:, b] / gap
    return out


def local_chordal_coefficient(s0: ArrayLike, a1: ArrayLike, a2: ArrayLike, *, k: int) -> float:
    """Coefficient ||K||_F^2 in d_C^2(t)=t^2||K||_F^2+O(t^3)."""
    k_matrix = local_rotation_matrix(s0, a1, a2, k=k)
    return float(np.linalg.norm(k_matrix, ord="fro") ** 2)


def davis_kahan_bound(s1: ArrayLike, s2: ArrayLike, *, k: int) -> float:
    """Simple spectral-gap upper bound for the largest channel angle."""
    base = as_symmetric(s1)
    values, _ = eigh_descending(base)
    if not 1 <= k < len(values):
        raise ValueError("invalid k")
    gap = values[k - 1] - values[k]
    if gap <= 0.0:
        return float("inf")
    return float(np.linalg.norm(as_symmetric(s2) - base, ord=2) / gap)


def sqrt_derivative(g: ArrayLike, g_dot: ArrayLike) -> Array:
    """Fréchet derivative D sqrt(G)[G_dot] via R X + X R = G_dot."""
    gg = assert_spd(g)
    direction = as_symmetric(g_dot)
    root = spd_sqrt(gg)
    values, vectors = np.linalg.eigh(root)
    b = vectors.T @ direction @ vectors
    y = b / (values[:, None] + values[None, :])
    return as_symmetric(vectors @ y @ vectors.T)


def channel_operator_derivative(
    g: ArrayLike,
    h: ArrayLike,
    g_dot: ArrayLike,
    h_dot: ArrayLike,
) -> Array:
    """Derivative of S when both standing variation and curvature evolve."""
    gg = assert_spd(g)
    hh = assert_spd(h)
    root = spd_sqrt(gg)
    x = sqrt_derivative(gg, g_dot)
    return as_symmetric(x @ hh @ root + root @ as_symmetric(h_dot) @ root + root @ hh @ x)
