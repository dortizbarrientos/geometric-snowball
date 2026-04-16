"""
geometry.py  --  Channel-operator and spectral geometry helpers
===============================================================

This module contains the core mathematical operations that implement the
additive channels framework for one population:

    G                     additive genetic variance matrix (pxp, positive definite)
    H                     fitness curvature matrix near the optimum (pxp, PD)
    H_tilde               rank-1 surrogate for H, built from the selection gradient beta
                              H_tilde = (beta beta^T) / ||beta||^2 + eps * I
    S_tilde               = G^{1/2} H_tilde G^{1/2}     (symmetric, PD)
    g_max                 leading eigenvector of G                   ~ Schluter limit
    channel               leading eigenvector of S_tilde, pulled back to phenotype space
    rotation              acute angle between g_max and the channel  (degrees)
    lambda_1, lambda_2    two largest eigenvalues of G                (for spectral gap)
    ratio                 lambda_1 / lambda_2                         (channel anisotropy)
    delta                 lambda_1 - lambda_2                         (raw spectral gap)

The functions are deliberately small and pure: given numpy arrays in, they
return numpy arrays out. No file I/O, no global state.
"""

import numpy as np


# ----------------------------------------------------------------------
# SYMMETRIC MATRIX SQUARE ROOT VIA EIGENDECOMPOSITION
# ----------------------------------------------------------------------
def matrix_sqrt_sym(M, floor=1e-12):
    """
    Compute the principal (positive-definite) square root of a symmetric
    positive-definite matrix M.

    Parameters
    ----------
    M : array (p,p)
        Symmetric matrix.
    floor : float
        Small positive floor applied to eigenvalues to guard against
        numerical negatives in near-singular matrices.

    Returns
    -------
    M_half : array (p,p)
        Symmetric matrix such that M_half @ M_half approx= M.
    """
    # eigh is guaranteed to return real eigenvalues for symmetric input
    vals, vecs = np.linalg.eigh((M + M.T) / 2)
    vals = np.maximum(vals, floor)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.T


# ----------------------------------------------------------------------
# ACUTE ANGLE BETWEEN TWO VECTORS, IN DEGREES
# ----------------------------------------------------------------------
def acute_angle_deg(u, v):
    """
    Return the acute angle (0 to 90 degrees) between two non-zero vectors.

    We take the absolute value of the cosine because an eigenvector is
    defined only up to sign, so "parallel" and "anti-parallel" are the
    same direction for our purposes.
    """
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    c = np.clip(abs(u @ v), -1.0, 1.0)
    return float(np.degrees(np.arccos(c)))


# ----------------------------------------------------------------------
# GEOMETRY OF ONE POPULATION
# ----------------------------------------------------------------------
def compute_channel_geometry(G, beta, eps=1e-4):
    """
    Given one G-matrix draw and a selection gradient vector, compute the
    spectral geometry that enters the channel-operator framework.

    Parameters
    ----------
    G : array (p,p)
        Additive genetic covariance matrix for one posterior draw.
    beta : array (p,)
        Standardized phenotypic selection gradient (used to build the
        rank-1 surrogate curvature matrix H_tilde).
    eps : float
        Positive scalar added as eps * I to H_tilde so that the rank-1
        surrogate is positive definite (and therefore S_tilde is as well).

    Returns
    -------
    dict with the following keys:
        lambda_1, lambda_2 : float
            Two largest eigenvalues of G.
        ratio : float
            lambda_1 / lambda_2, the channel-anisotropy ratio of G.
        delta : float
            lambda_1 - lambda_2, the raw spectral gap of G.
        rotation_deg : float
            Acute angle between g_max (leading eigenvector of G) and the
            channel (leading eigenvector of S_tilde, pulled back to
            phenotype space via G^{1/2}).
        was_PD : bool
            False if the input G had a non-positive eigenvalue (these
            draws are excluded from summaries).
    """
    # Eigendecomposition of G (symmetric, expected PD)
    vals, vecs = np.linalg.eigh((G + G.T) / 2)
    idx = np.argsort(vals)[::-1]      # descending
    vals = vals[idx]
    vecs = vecs[:, idx]

    if vals[-1] <= 1e-10:
        return dict(lambda_1=np.nan, lambda_2=np.nan,
                    ratio=np.nan, delta=np.nan,
                    rotation_deg=np.nan, was_PD=False)

    lambda_1, lambda_2 = vals[0], vals[1]

    # G^{1/2} from the same decomposition
    G_half = vecs @ np.diag(np.sqrt(vals)) @ vecs.T

    # H_tilde = (beta beta^T) / ||beta||^2 + eps I
    beta_norm = beta / np.linalg.norm(beta)
    p = G.shape[0]
    H_tilde = np.outer(beta_norm, beta_norm) + eps * np.eye(p)

    # S_tilde = G^{1/2} H_tilde G^{1/2}
    S_tilde = G_half @ H_tilde @ G_half

    # Leading eigenvector of S_tilde (whitened); pull back with G^{1/2}
    _, S_vecs = np.linalg.eigh((S_tilde + S_tilde.T) / 2)
    channel_whitened = S_vecs[:, -1]            # eigh returns ascending
    channel = G_half @ channel_whitened

    # Rotation g_max (leading eigenvector of G, in phenotype space) -> channel
    g_max = vecs[:, 0]
    rotation = acute_angle_deg(g_max, channel)

    return dict(
        lambda_1=float(lambda_1),
        lambda_2=float(lambda_2),
        ratio=float(lambda_1 / lambda_2),
        delta=float(lambda_1 - lambda_2),
        rotation_deg=rotation,
        was_PD=True,
    )


# ----------------------------------------------------------------------
# POSTERIOR SUMMARIES
# ----------------------------------------------------------------------
def posterior_summary(x):
    """
    Posterior median and 95% credible interval for a 1-D array,
    ignoring NaNs.

    Returns
    -------
    dict(median, lo95, hi95, n_used)
    """
    x = np.asarray(x)
    x = x[~np.isnan(x)]
    if x.size == 0:
        return dict(median=np.nan, lo95=np.nan, hi95=np.nan, n_used=0)
    return dict(
        median=float(np.median(x)),
        lo95=float(np.quantile(x, 0.025)),
        hi95=float(np.quantile(x, 0.975)),
        n_used=int(x.size),
    )
