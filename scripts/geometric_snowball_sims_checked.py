"""
geometric_snowball_sims_checked.py
=============================================================================
"The Geometric Snowball of Speciation" — Ortiz-Barrientos, James, Burrage

PURPOSE
-------
Generates the theoretical figures (Figs 1, 2, S1, S2, S3) for the paper.
All figures are produced from simulation and mathematical computation;
no empirical data are used (see geometric_snowball_senecio.py for the
Senecio analysis).

FIGURES PRODUCED AND THEIR BIOLOGICAL MEANING
----------------------------------------------

Figure 1 — Two panels showing how channel misalignment grows and why it
  generates reproductive isolation:
  Panel a: Channel misalignment theta (degrees) vs divergence time, for three
           ecological contrasts (small, moderate, large Delta H). Shows that
           stronger ecological differences cause faster channel rotation.
  Panel b: Geometric incompatibility I_geom vs theta for several values of k
           (the effective number of channels). Shows that I_geom grows as
           theta^2 for small angles (the snowball), and accelerates as k grows.

Figure 2 — Two panels showing the classical recovery and the boundary signature:
  Panel a: Schluter recovery. When fitness curvature is isotropic (H = hI),
           the leading channel coincides with g_max (Schluter's line of least
           resistance, theta = 0). As curvature becomes anisotropic, the
           channel rotates away from g_max. The degree of rotation depends on
           how elliptical G is — more elliptical G resists rotation more.
  Panel b: Boundary signature. Populations near a "channel boundary" (spectral
           gap delta ≈ 0) accumulate channel misalignment fastest per unit of
           ecological perturbation. The Davis-Kahan bound sin(theta) <= |Delta S|/delta
           predicts a 1/delta sensitivity, which this panel verifies numerically
           across 120 random G-matrix configurations. Colour = additivity index Ag.

Supplementary Figure S1 — Two-trait worked example. Shows how the leading
  channel rotates away from g_max as curvature anisotropy grows, for
  different levels of genetic correlation r between the two traits.

Supplementary Figure S2 — Davis-Kahan bound verification in the two-trait case:
  Panel a: Contour surface of the bound over (g_perp, delta) parameter space.
  Panel b: Scatter plot of observed sin(theta) vs the predicted bound across
           120 random two-trait configurations. The bound is always satisfied
           (all points below the 45-degree line).

Supplementary Figure S3 — Boundary signature diagnostics across five-trait
  simulations:
  Panel a: Degeneracy D = lambda2/lambda1 vs observed sensitivity.
  Panel b: Variance toward the boundary (V_perp) vs sensitivity.
  Panel c: Additivity Ag vs D, with the exact 2D identity D=(1-Ag)/Ag.

THE CENTRAL MATHEMATICAL CORRECTION (why S, not G@H)
------------------------------------------------------
Lande's equation gives the mean phenotypic change per generation as:
    Delta z_bar = G * beta
where beta = -H*(z_bar - theta) is the selection gradient. The local adaptive
dynamics are therefore governed by the matrix G*H.

But G*H is NOT symmetric in general. That means it does not admit an
orthonormal eigendecomposition, and its "eigenvectors" are not a valid
orthonormal basis for computing angles between subspaces.

The correct object is the symmetric similarity transform S = G^{1/2} H G^{1/2}.
S shares eigenvalues with G*H (they are related by a congruence transform), but
S has a real orthonormal eigendecomposition. Its eigenvectors — the "genetic
channels" — are orthonormal in the genetically whitened metric, and standard
perturbation theory (the Davis-Kahan theorem) applies.

All angle calculations in this script therefore use the eigenspaces of S.
Where phenotype-space channel directions are needed for plotting or
interpretation, they are mapped back via u = G^{1/2} * v.

OUTPUT DIRECTORY
----------------
Figures are written to a directory resolved in this order:
  1. First command-line argument:  python3 script.py ~/Desktop/figures
  2. Environment variable:         SNOWBALL_OUTDIR=~/Desktop/figures python3 script.py
  3. Current working directory:    ./figures  (created automatically)

USAGE
-----
python geometric_snowball_sims_checked.py               # save to ./figures/
python geometric_snowball_sims_checked.py ~/Desktop     # save to ~/Desktop
SNOWBALL_OUTDIR=~/myfigures python geometric_snowball_sims_checked.py

DEPENDENCIES
------------
Python >= 3.8, numpy, matplotlib

SANITY CHECKS
-------------
Running this script also prints sanity checks that verify three fundamental
properties of the framework: (1) Schluter recovery when H is isotropic,
(2) I_geom endpoints at theta=0 and theta=90°, (3) small-angle approximation,
and (4) Davis-Kahan bound satisfaction.
=============================================================================

AUTHORS
-------
D. Ortiz-Barrientos, M.E. James, K. Burrage (University of Queensland / QUT)
Code annotations: Claude (Anthropic), April 2026
Repository: github.com/dortizbarrientos/geometric-snowball
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path
from typing import Tuple

import matplotlib
import numpy as np


# ---------------------------------------------------------------------------
# Output directory — works on any machine
# ---------------------------------------------------------------------------
def _resolve_outdir() -> Path:
    if len(sys.argv) > 1:
        d = Path(sys.argv[1]).expanduser()
    elif "SNOWBALL_OUTDIR" in os.environ:
        d = Path(os.environ["SNOWBALL_OUTDIR"])
    else:
        d = Path.cwd() / "figures"
    d.mkdir(parents=True, exist_ok=True)
    return d

OUTDIR = _resolve_outdir()
print(f"Figures will be saved to: {OUTDIR}")

# Backend: non-interactive (file output only).
# Remove this line if you want figures to pop up on screen.
matplotlib.use("Agg")

import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Plot style
# -----------------------------------------------------------------------------
plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 1.0,
        "figure.dpi": 300,
    }
)

# Okabe–Ito colour-blind-safe palette
C = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#7F7F7F",
}

MM_TO_IN = 1 / 25.4
COL1 = 89 * MM_TO_IN
COL2 = 183 * MM_TO_IN


# -----------------------------------------------------------------------------
# Linear algebra helpers for symmetric positive-definite matrices
# -----------------------------------------------------------------------------
def symmetrise(A: np.ndarray) -> np.ndarray:
    """Return the symmetric part of a square matrix."""
    return 0.5 * (A + A.T)



def spectral_decomp_spd(A: np.ndarray, tol: float = 1e-12) -> Tuple[np.ndarray, np.ndarray]:
    """
    Eigen-decomposition for a symmetric positive-definite matrix.

    Returns eigenvalues in descending order together with the corresponding
    eigenvectors as columns.
    """
    A = symmetrise(A)
    evals, evecs = np.linalg.eigh(A)
    if np.min(evals) <= tol:
        raise ValueError(
            f"Matrix is not numerically SPD. Smallest eigenvalue = {np.min(evals):.3e}"
        )
    order = np.argsort(evals)[::-1]
    return evals[order], evecs[:, order]



def spd_sqrt(A: np.ndarray) -> np.ndarray:
    """Symmetric square root of an SPD matrix."""
    evals, evecs = spectral_decomp_spd(A)
    return evecs @ np.diag(np.sqrt(evals)) @ evecs.T



def orthonormalise_columns(U: np.ndarray) -> np.ndarray:
    """
    Orthonormalise the columns of U with a QR decomposition.

    For symmetric eigensystems this is mostly a safety step against numerical
    drift.
    """
    Q, _ = np.linalg.qr(U)
    return Q


# -----------------------------------------------------------------------------
# Channel geometry
# -----------------------------------------------------------------------------
def channel_operator(G: np.ndarray, H: np.ndarray) -> np.ndarray:
    """
    Compute the symmetric channel operator S = G^{1/2} H G^{1/2}.

    BIOLOGICAL INTERPRETATION
    -------------------------
    S encodes both the available additive genetic variance (G) and the local
    curvature of the fitness landscape (H) into a single symmetric matrix.
    Its eigenvectors are the "additive genetic channels" — the directions in
    genetically whitened trait space along which adaptation proceeds most
    efficiently. The leading eigenvector corresponds to the direction where
    genetic variance and stabilising selection reinforce each other most.

    When H = h*I (isotropic fitness surface), S = h*G and the channels
    reduce to the eigenvectors of G, recovering Schluter's lines of least
    resistance.

    MATHEMATICAL DETAIL
    -------------------
    G^{1/2} is the symmetric matrix square root of G (see spd_sqrt).
    The output is symmetrised numerically to guard against floating-point drift.
    """
    G_half = spd_sqrt(G)
    S = G_half @ H @ G_half
    return symmetrise(S)



def leading_channel_basis(
    G: np.ndarray,
    H: np.ndarray,
    k: int = 1,
    space: str = "whitened",
) -> np.ndarray:
    """
    Return the leading k-dimensional channel basis.

    Parameters
    ----------
    G, H : ndarray
        SPD additive-genetic and curvature matrices.
    k : int
        Number of leading channel directions.
    space : {"whitened", "phenotype"}
        In whitened space, basis vectors are orthonormal eigenvectors of
        S = G^{1/2} H G^{1/2}.
        In phenotype space, vectors are mapped back by u = G^{1/2} v and then
        orthonormalised in the ordinary Euclidean metric for display.

    Notes
    -----
    WHICH SPACE TO USE?
    Use "whitened" space for all channel-angle computations — it is the space
    where S has an orthonormal eigendecomposition and the Davis-Kahan bound
    directly applies.

    Use "phenotype" space only for PLOTTING or BIOLOGICAL INTERPRETATION:
    to show which trait combinations the channel points along, or to compare
    with g_max from the original G-matrix. The mapping back to phenotype space
    is u = G^{1/2} * v, then orthonormalised.

    The Schluter limit (H = h*I) gives channels that coincide with the
    eigenvectors of G in phenotype space — directly comparable to the
    empirical PCVA divergence vectors used in the Senecio analysis.
    """
    S = channel_operator(G, H)
    _, V = spectral_decomp_spd(S)
    V_k = V[:, :k]

    if space == "whitened":
        return orthonormalise_columns(V_k)

    if space == "phenotype":
        G_half = spd_sqrt(G)
        U_k = G_half @ V_k
        return orthonormalise_columns(U_k)

    raise ValueError("space must be 'whitened' or 'phenotype'")



def principal_angles(U1: np.ndarray, U2: np.ndarray) -> np.ndarray:
    """
    Principal angles between two subspaces whose columns are orthonormal bases.

    Returns the angles in radians, sorted from smallest to largest.
    """
    M = U1.T @ U2
    singular_values = np.linalg.svd(M, compute_uv=False)
    singular_values = np.clip(singular_values, 0.0, 1.0)
    angles = np.arccos(singular_values)
    return np.sort(angles)



def largest_principal_angle(U1: np.ndarray, U2: np.ndarray) -> float:
    """Largest principal angle in radians."""
    return float(principal_angles(U1, U2)[-1])



def vector_angle(u: np.ndarray, v: np.ndarray) -> float:
    """Acute angle between two vectors in radians."""
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    dot = np.clip(abs(float(u @ v)), 0.0, 1.0)
    return float(np.arccos(dot))


# -----------------------------------------------------------------------------
# Incompatibility measures
# -----------------------------------------------------------------------------
def incompatibility_general(angles_rad: np.ndarray) -> float:
    """
    General geometric incompatibility:

        I_geom = 1 - prod_j cos^2(psi_j)
    """
    return float(1.0 - np.prod(np.cos(angles_rad) ** 2))



def incompatibility_equal_angle(theta_rad: np.ndarray | float, k: int) -> np.ndarray | float:
    """
    Equal-angle special case:

        I_geom = 1 - cos^(2k)(theta)
    """
    return 1.0 - np.cos(theta_rad) ** (2 * k)


# -----------------------------------------------------------------------------
# Random SPD matrix constructors
# -----------------------------------------------------------------------------
def random_orthogonal(n: int, seed: int | None = None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    return Q



def spd_from_eigensystem(eigenvalues: np.ndarray, seed: int | None = None) -> np.ndarray:
    """Construct a random SPD matrix with prescribed eigenvalues."""
    Q = random_orthogonal(len(eigenvalues), seed=seed)
    A = Q @ np.diag(eigenvalues) @ Q.T
    return symmetrise(A)


def make_curvature_matrix(n: int, anisotropy: float, seed: int | None = None) -> np.ndarray:
    """
    Build an SPD curvature matrix with mean eigenvalue one.

    anisotropy = 0 gives the isotropic case H = I.
    Larger anisotropy increases curvature differences across trait directions.
    """
    # Eigenvalues rise smoothly away from 1, then are normalised back to mean 1.
    base = np.linspace(0.0, 1.0, n)
    h_vals = 1.0 + anisotropy * base
    h_vals = h_vals / np.mean(h_vals)
    return spd_from_eigensystem(h_vals, seed=seed)



def make_curvature_matrix_at_45(G: np.ndarray, anisotropy: float) -> np.ndarray:
    """
    Curvature matrix whose strongest curvature direction lies at 45 degrees
    between the first two eigenvectors of G.

    The matrix is rescaled so that its mean eigenvalue is one. This lets us use
    \\|H - I\\|_2 as an honest anisotropy measure on the x-axis.
    """
    _, evecs_G = spectral_decomp_spd(G)
    e1 = evecs_G[:, 0]
    e2 = evecs_G[:, 1]
    v_45 = (e1 + e2) / math.sqrt(2.0)

    n = G.shape[0]
    H = np.eye(n) + anisotropy * np.outer(v_45, v_45)
    H = symmetrise(H)
    H = H / (np.trace(H) / n)
    return H


# -----------------------------------------------------------------------------
# Simulation pieces used in the figures
# -----------------------------------------------------------------------------
def simulate_channel_divergence(
    G: np.ndarray,
    H1: np.ndarray,
    H2_final: np.ndarray,
    n_steps: int = 200,
    k: int = 1,
) -> np.ndarray:
    """
    Track channel misalignment as one population's curvature matrix moves from
    H1 to H2_final along a linear path in SPD space.

    Because the SPD cone is convex, H(t) = (1-a)H1 + aH2 remains SPD.
    """
    basis_ref = leading_channel_basis(G, H1, k=k, space="whitened")
    thetas_deg = []

    for t in range(n_steps + 1):
        alpha = t / n_steps
        H_t = (1.0 - alpha) * H1 + alpha * H2_final
        basis_t = leading_channel_basis(G, H_t, k=k, space="whitened")
        theta = largest_principal_angle(basis_ref, basis_t)
        thetas_deg.append(np.degrees(theta))

    return np.asarray(thetas_deg)



def schluter_rotation_from_isotropy(G: np.ndarray, anisotropy: float) -> tuple[float, float]:
    """
    Return two quantities for the 45-degree curvature construction:

    1. the actual anisotropy \\|H - I\\|_2, and
    2. the angle between the phenotype-space leading channel and g_max.

    When anisotropy = 0 the result is exactly (0, 0), recovering Schluter's limit.
    """
    _, evecs_G = spectral_decomp_spd(G)
    gmax = evecs_G[:, 0]
    H = make_curvature_matrix_at_45(G, anisotropy)
    channel_vec = leading_channel_basis(G, H, k=1, space="phenotype")[:, 0]
    delta_H = np.linalg.norm(H - np.eye(G.shape[0]), ord=2)
    rotation = np.degrees(vector_angle(gmax, channel_vec))
    return float(delta_H), float(rotation)


# =============================================================================
# FIGURE 1: Channel misalignment growth and the incompatibility curve
# =============================================================================
# BIOLOGICAL STORY
# Panel a (left): How fast do channels rotate apart?
#   - Simulates one population's curvature matrix H moving from isotropic
#     (H = I, Schluter limit) to progressively anisotropic curvature.
#   - Three ecological contrasts (small/moderate/large Delta H) show that
#     stronger environmental differences cause faster channel rotation.
#   - The Schluter limit (dashed horizontal) is theta = 0 at all times —
#     no channel rotation when both populations share the same curvature.
#
# Panel b (right): How does channel rotation translate to isolation?
#   - Plots I_geom = 1 - cos^{2k}(theta) for k = 1,2,4,8 effective channels.
#   - The quadratic approximation I ~ k*theta^2 (dotted line) shows that the
#     growth is nonlinear — this is the SNOWBALL: slow start, then explosive.
#   - More dimensions (larger k) means faster accumulation of incompatibility.
# =============================================================================
def make_figure1() -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(COL2, COL2 * 0.42))
    ax1a, ax1b = axes

    # Moderately elliptical G with fixed total scale.
    n = 6
    G = spd_from_eigensystem(np.array([0.85, 0.62, 0.43, 0.30, 0.21, 0.15]), seed=42)
    H1 = np.eye(n)

    scenarios = [
        (0.5, C["sky"], "Weak ($\\Delta H$ small)"),
        (1.5, C["orange"], "Moderate ($\\Delta H$ medium)"),
        (3.0, C["red"], "Strong ($\\Delta H$ large)"),
    ]

    time = np.linspace(0.0, 1.0, 201)
    for idx, (anisotropy, colour, label) in enumerate(scenarios):
        H2 = make_curvature_matrix(n=n, anisotropy=anisotropy, seed=100 + idx)
        thetas_deg = simulate_channel_divergence(G, H1, H2, n_steps=200, k=1)
        ax1a.plot(time, thetas_deg, color=colour, lw=1.2, label=label)

    ax1a.axhline(
        0.0,
        color=C["black"],
        lw=0.6,
        ls="--",
        alpha=0.5,
        label="Schluter limit ($H \\propto I$)",
    )
    ax1a.set_xlabel("Divergence time (arbitrary units)")
    ax1a.set_ylabel("Channel misalignment $\\theta_{12}$ (degrees)")
    ax1a.set_xlim(0.0, 1.0)
    ax1a.set_ylim(-2.0, 80.0)
    ax1a.legend(loc="upper left", frameon=False)
    ax1a.spines["top"].set_visible(False)
    ax1a.spines["right"].set_visible(False)
    ax1a.text(-0.18, 1.04, "a", transform=ax1a.transAxes, fontsize=9, fontweight="bold", va="top")

    theta_deg = np.linspace(0.0, 90.0, 500)
    theta_rad = np.radians(theta_deg)
    k_vals = [1, 2, 4, 8]
    k_cols = [C["sky"], C["green"], C["orange"], C["red"]]
    k_labels = ["$k=1$", "$k=2$", "$k=4$", "$k=8$"]

    for k, colour, label in zip(k_vals, k_cols, k_labels):
        I = incompatibility_equal_angle(theta_rad, k)
        ax1b.plot(theta_deg, I, color=colour, lw=1.2, label=label)

    theta_ref_deg = np.linspace(0.0, 30.0, 100)
    theta_ref_rad = np.radians(theta_ref_deg)
    # Small-angle benchmark for k = 1: I ~ theta^2 when theta is in radians.
    I_quad = theta_ref_rad**2
    ax1b.plot(theta_ref_deg, I_quad, "k:", lw=0.8, alpha=0.6, label="$\\sim \\theta^2$")

    ax1b.set_xlabel("Channel misalignment $\\theta_{12}$ (degrees)")
    ax1b.set_ylabel("Geometric incompatibility $I_{\\mathrm{geom}}$")
    ax1b.set_xlim(0.0, 90.0)
    ax1b.set_ylim(0.0, 1.05)
    ax1b.set_xticks([0, 30, 60, 90])
    ax1b.legend(loc="upper left", frameon=False, ncol=2)
    ax1b.spines["top"].set_visible(False)
    ax1b.spines["right"].set_visible(False)
    ax1b.text(-0.18, 1.04, "b", transform=ax1b.transAxes, fontsize=9, fontweight="bold", va="top")
    ax1b.annotate(
        "Geometric\nsnowball",
        xy=(42, 0.50),
        xytext=(58, 0.30),
        fontsize=6,
        color=C["orange"],
        arrowprops=dict(arrowstyle="->", color=C["orange"], lw=0.7),
    )

    fig.tight_layout(pad=0.5)
    outfile = OUTDIR / "figure1_snowball_corrected.pdf"
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    return outfile


# =============================================================================
# FIGURE 2: Schluter recovery and boundary signature
# =============================================================================
# BIOLOGICAL STORY
# Panel a (left): Schluter recovery.
#   When the fitness surface is isotropic (the curvature is the same in every
#   trait direction), the leading channel exactly equals g_max — Schluter's
#   line of least resistance. As curvature becomes anisotropic (the fitness
#   surface has ridges or valleys in specific trait combinations), the channel
#   rotates away from g_max. How far it rotates depends on how elliptical
#   G is: a G-matrix with a strong dominant eigenvalue (strongly channelled
#   genetic variance) resists rotation; a near-spherical G-matrix rotates
#   easily. This panel recovers and extends the classical Schluter result.
#
# Panel b (right): Boundary signature.
#   The Davis-Kahan theorem gives a hard bound: sin(theta) <= |Delta S|/delta,
#   where delta = lambda_1 - lambda_2 is the spectral gap of S. Populations
#   near a "channel boundary" (delta ≈ 0, two eigenvalues nearly equal) can
#   accumulate large channel misalignment from a small ecological perturbation.
#   This panel verifies that across 120 random G-matrix configurations, the
#   observed sin(theta)/|Delta S| (y-axis) follows the predicted 1/delta law
#   (dashed curve). Points are coloured by the additivity index Ag = lambda_1
#   / sum(lambda), which is high (green) when one channel clearly dominates
#   and low (red) when channels are nearly degenerate.
# =============================================================================
def make_figure2() -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(COL2, COL2 * 0.42))
    ax2a, ax2b = axes

    # Panel (a): recovery of Schluter's limit.
    n = 4
    g_shapes = [
        (np.array([1.0, 0.6, 0.3, 0.1]), C["blue"], "Elliptical $\\mathbf{G}$"),
        (np.array([0.5, 0.45, 0.42, 0.40]), C["green"], "Near-spherical $\\mathbf{G}$"),
        (np.array([2.0, 0.4, 0.2, 0.1]), C["red"], "Highly elliptical $\\mathbf{G}$"),
    ]
    anisotropy_values = np.linspace(0.0, 4.0, 100)

    for seed, (evals_G, colour, label) in enumerate(g_shapes, start=7):
        G = spd_from_eigensystem(evals_G, seed=seed)
        xs, ys = zip(*(schluter_rotation_from_isotropy(G, a) for a in anisotropy_values))
        ax2a.plot(xs, ys, color=colour, lw=1.2, label=label)

    ax2a.axvline(0.0, color=C["gray"], lw=0.6, ls=":")
    ax2a.text(0.03, 2.0, "Schluter's\nlimit", fontsize=5.5, color=C["gray"], va="bottom")
    ax2a.set_xlabel("Fitness-landscape anisotropy $\\|\\mathbf{H}-\\mathbf{I}\\|_2$")
    ax2a.set_ylabel("Channel rotation from $\\mathbf{g}_{\\max}$ (degrees)")
    ax2a.set_xlim(0.0, 1.6)
    ax2a.set_ylim(-1.0, 50.0)
    ax2a.legend(loc="upper left", frameon=False)
    ax2a.spines["top"].set_visible(False)
    ax2a.spines["right"].set_visible(False)
    ax2a.text(-0.18, 1.04, "a", transform=ax2a.transAxes, fontsize=9, fontweight="bold", va="top")

    # Panel (b): exact Davis–Kahan sensitivity under the corrected geometry.
    rng = np.random.default_rng(13)
    n = 5
    spectral_gaps = []
    sensitivities = []
    A_g_values = []

    n_sims = 120
    for _ in range(n_sims):
        lam = np.sort(rng.exponential(scale=0.5, size=n))[::-1]
        lam = lam / lam.sum() * n  # standardise total additive scale
        G = spd_from_eigensystem(lam, seed=int(rng.integers(0, 1_000_000)))
        G_half = spd_sqrt(G)

        H1 = np.eye(n)
        raw = rng.normal(size=(n, n))
        dH = symmetrise(raw)
        dH_norm = np.linalg.norm(dH, ord=2)
        if dH_norm == 0.0:
            continue
        dH = 0.25 * dH / dH_norm
        H2 = H1 + dH

        # Ensure H2 remains SPD. If not, shift along the identity.
        min_eval_H2 = np.min(np.linalg.eigvalsh(symmetrise(H2)))
        if min_eval_H2 <= 1e-8:
            H2 = H2 + (abs(min_eval_H2) + 1e-3) * np.eye(n)

        S1 = channel_operator(G, H1)
        evals_S1, _ = spectral_decomp_spd(S1)
        delta = evals_S1[0] - evals_S1[1]
        if delta <= 1e-8:
            continue

        V1 = leading_channel_basis(G, H1, k=1, space="whitened")
        V2 = leading_channel_basis(G, H2, k=1, space="whitened")
        theta = largest_principal_angle(V1, V2)

        delta_S = G_half @ (H2 - H1) @ G_half
        delta_S_norm = np.linalg.norm(delta_S, ord=2)
        if delta_S_norm <= 0.0:
            continue

        sensitivity = np.sin(theta) / delta_S_norm
        A_g = evals_S1[0] / np.sum(evals_S1)

        spectral_gaps.append(delta)
        sensitivities.append(sensitivity)
        A_g_values.append(A_g)

    spectral_gaps = np.asarray(spectral_gaps)
    sensitivities = np.asarray(sensitivities)
    A_g_values = np.asarray(A_g_values)

    sc = ax2b.scatter(
        spectral_gaps,
        sensitivities,
        c=A_g_values,
        cmap="RdYlGn",
        vmin=float(np.min(A_g_values)),
        vmax=float(np.max(A_g_values)),
        s=18,
        alpha=0.8,
        linewidths=0.3,
        edgecolors="white",
    )

    delta_smooth = np.geomspace(
        max(1e-3, float(np.min(spectral_gaps))),
        float(np.max(spectral_gaps)), 200
    )
    ax2b.plot(
        delta_smooth,
        1.0 / delta_smooth,
        color=C["black"],
        lw=1.0,
        ls="--",
        label="Davis--Kahan bound $1/\\delta$",
    )

    cbar = fig.colorbar(sc, ax=ax2b, shrink=0.7, pad=0.02)
    cbar.set_label("Additivity index $\\mathcal{A}_g$", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    ax2b.set_xscale("log")
    ax2b.set_yscale("log")
    ax2b.set_xlabel("Spectral gap $\\delta$ of $\\mathbf{S}$ (log scale)")
    ax2b.set_ylabel(
        "Observed $\\sin\\theta_{12}/\\|\\Delta\\mathbf{S}\\|_2$ (log scale)"
    )
    ax2b.legend(loc="upper right", frameon=False)
    ax2b.spines["top"].set_visible(False)
    ax2b.spines["right"].set_visible(False)
    ax2b.text(-0.18, 1.04, "b", transform=ax2b.transAxes, fontsize=9, fontweight="bold", va="top")

    fig.tight_layout(pad=0.5)
    outfile = OUTDIR / "figure2_boundary_corrected.pdf"
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    return outfile


# -----------------------------------------------------------------------------
# Supplementary Figure S1 — 2-trait channel rotation vs H anisotropy
# Shows phi = angle(leading channel, g_max) as curvature ratio h22/h11 grows.
# We avoid r = 0 exactly because spherical G makes g_max undefined.
# -----------------------------------------------------------------------------
def make_figS1() -> Path:
    fig, ax = plt.subplots(figsize=(COL1, COL1 * 0.85))

    r_vals = [0.1, 0.4, 0.7, 0.9]
    r_cols = [C["sky"], C["green"], C["orange"], C["red"]]
    h22_range = np.linspace(1.0, 5.0, 200)

    for r, colour in zip(r_vals, r_cols):
        G2 = np.array([[1.0, r], [r, 1.0]])
        rotations = []
        for h22 in h22_range:
            H2 = np.diag([1.0, h22])
            evals_G, evecs_G = spectral_decomp_spd(G2)
            gmax = evecs_G[:, 0]
            ch = leading_channel_basis(G2, H2, k=1, space="phenotype")[:, 0]
            rotations.append(np.degrees(vector_angle(gmax, ch)))
        ax.plot(h22_range, rotations, color=colour, lw=1.2, label=f"$r={r}$")

    ax.axhline(0.0, color=C["gray"], lw=0.6, ls="--", label="Schluter limit")
    ax.set_xlabel("Curvature ratio $h_{22}/h_{11}$")
    ax.set_ylabel("Channel rotation from $\\mathbf{g}_{\\max}$ (degrees)")
    ax.set_xlim(1.0, 5.0)
    ax.set_ylim(-1.0, 50.0)
    ax.legend(loc="upper left", frameon=False, title="Genetic correlation $r$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("S1", loc="left", fontsize=8, fontweight="bold")

    fig.tight_layout(pad=0.5)
    outfile = OUTDIR / "figS1_channel_rotation.pdf"
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    return outfile


# -----------------------------------------------------------------------------
# Supplementary Figure S2 — Davis-Kahan bound: surface + scatter verification
# Panel a: bound surface sin(theta) <= g11 * dh / delta over (g11, delta).
# Panel b: observed sin(theta) vs bound across 120 two-trait simulations.
# -----------------------------------------------------------------------------
def make_figS2() -> Path:
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(COL2, COL2 * 0.45))

    # --- panel a: exact 2D diagonal worked example ---
    # G = diag(g_parallel, g_perp), H1 = I, H2 = I + diag(0, Delta h).
    # Then Delta S = diag(0, g_perp Delta h), so
    #   sin(theta_12) <= g_perp |Delta h| / delta.
    g_perp_vals = np.linspace(0.05, 1.5, 100)
    delta_vals = np.linspace(0.05, 2.0, 100)
    GPERP, DELTA = np.meshgrid(g_perp_vals, delta_vals)
    BOUND = np.clip(GPERP * 0.3 / DELTA, 0.0, 1.0)
    cf = ax_a.contourf(GPERP, DELTA, BOUND, levels=20, cmap="YlOrRd")
    ax_a.contour(
        GPERP, DELTA, BOUND,
        levels=[0.1, 0.3, 0.5, 0.7],
        colors="white", linewidths=0.5, alpha=0.7,
    )
    fig.colorbar(cf, ax=ax_a, shrink=0.85).set_label(
        "Bound on $\\sin\\theta_{12}$", fontsize=6
    )
    ax_a.set_xlabel(r"$g_{\perp}$ (variance in perturbed direction)")
    ax_a.set_ylabel("Spectral gap $\\delta$ of $S$")
    ax_a.set_title("2D diagonal baseline, $|\\Delta h|=0.3$", fontsize=6)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.text(-0.15, 1.06, "a", transform=ax_a.transAxes,
              fontsize=9, fontweight="bold", va="top")

    # --- panel b: scatter verification using the full matrix bound ---
    rng = np.random.default_rng(99)
    sins_obs, bounds_true = [], []
    for _ in range(120):
        g11 = rng.uniform(0.1, 1.5)
        g22 = rng.uniform(0.05, 0.8)
        r = rng.uniform(-0.5, 0.5)
        cov = r * math.sqrt(g11 * g22)
        G2 = np.array([[g11, cov], [cov, g22]])
        if np.linalg.det(G2) < 1e-6:
            continue
        H1 = np.eye(2)
        dh = rng.uniform(0.05, 0.5)
        dH = np.diag([0.0, dh])
        H2 = H1 + dH
        S1 = channel_operator(G2, H1)
        dS = channel_operator(G2, H2) - S1
        evals_S1 = spectral_decomp_spd(S1)[0]
        delta = evals_S1[0] - evals_S1[1]
        if delta < 1e-8:
            continue
        bnd = np.linalg.norm(dS, ord=2) / delta
        V1 = leading_channel_basis(G2, H1, k=1, space="whitened")
        V2 = leading_channel_basis(G2, H2, k=1, space="whitened")
        th = largest_principal_angle(V1, V2)
        sins_obs.append(float(np.sin(th)))
        bounds_true.append(float(bnd))

    sins_obs = np.array(sins_obs)
    bounds_true = np.array(bounds_true)
    lim = max(bounds_true.max(), sins_obs.max()) * 1.05

    ax_b.scatter(bounds_true, sins_obs, s=12, alpha=0.7,
                 color=C["blue"], edgecolors="white", linewidths=0.3)
    ax_b.plot([0, lim], [0, lim], "k--", lw=0.8, label="Bound (equality)")
    ax_b.set_xlabel("Davis--Kahan bound ($\\|\\Delta S\\|_2/\\delta$)")
    ax_b.set_ylabel("Observed $\\sin\\theta_{12}$")
    ax_b.set_xlim(0, lim)
    ax_b.set_ylim(0, lim)
    ax_b.legend(frameon=False, fontsize=5.5)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    pct = 100.0 * float(np.mean(sins_obs <= bounds_true))
    ax_b.text(0.55, 0.12, f"Bound satisfied: {pct:.0f}%",
              transform=ax_b.transAxes, fontsize=5.5, color=C["blue"])
    ax_b.text(-0.15, 1.06, "b", transform=ax_b.transAxes,
              fontsize=9, fontweight="bold", va="top")

    fig.tight_layout(pad=0.5)
    outfile = OUTDIR / "figS2_DK_bound.pdf"
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    return outfile


# -----------------------------------------------------------------------------
# Supplementary Figure S3 — boundary signature diagnostics (n=5, 100 sims)
# Panel a: degeneracy D = lambda2/lambda1 vs theta growth rate.
# Panel b: V_perp (variance toward boundary) vs theta growth rate.
# Panel c: additivity A_g vs D; dashed = exact 2D identity D=(1-A)/A.
# -----------------------------------------------------------------------------
def make_figS3() -> Path:
    rng = np.random.default_rng(21)
    n = 5
    D_list, Ag_list, Vp_list, sens_list = [], [], [], []

    for _ in range(100):
        lam = np.sort(rng.exponential(scale=0.4, size=n))[::-1] + 0.05
        G = spd_from_eigensystem(lam, seed=int(rng.integers(0, 1_000_000)))
        H1 = np.eye(n)
        S = channel_operator(G, H1)
        evals_S, evecs_S = spectral_decomp_spd(S)
        D = evals_S[1] / evals_S[0]
        Ag = evals_S[0] / np.sum(evals_S)

        # V_perp: phenotype-space variance in second channel direction.
        G_half = spd_sqrt(G)
        q2 = evecs_S[:, 1]
        u2 = G_half @ q2
        u2 /= np.linalg.norm(u2)
        Vp = float(u2 @ G @ u2)

        # Exact DK-style sensitivity for a small curvature perturbation.
        raw = rng.normal(size=(n, n))
        dH = 0.5 * (raw + raw.T)
        dH = 0.3 * dH / np.linalg.norm(dH, ord=2)
        H2 = H1 + dH
        min_eval_H2 = np.min(np.linalg.eigvalsh(symmetrise(H2)))
        if min_eval_H2 <= 1e-8:
            H2 = H2 + (abs(min_eval_H2) + 1e-3) * np.eye(n)
        V1 = leading_channel_basis(G, H1, k=1, space="whitened")
        V2 = leading_channel_basis(G, H2, k=1, space="whitened")
        theta = largest_principal_angle(V1, V2)
        delta_S_norm = np.linalg.norm(G_half @ (H2 - H1) @ G_half, ord=2)
        sensitivity = np.sin(theta) / delta_S_norm

        D_list.append(D)
        Ag_list.append(Ag)
        Vp_list.append(Vp)
        sens_list.append(sensitivity)

    D_arr = np.array(D_list)
    Ag_arr = np.array(Ag_list)
    Vp_arr = np.array(Vp_list)
    s_arr = np.array(sens_list)

    fig, axes = plt.subplots(1, 3, figsize=(COL2, COL2 * 0.38))
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # panel a
    axes[0].scatter(D_arr, s_arr, s=14, alpha=0.8,
                    c=Ag_arr, cmap="RdYlGn", vmin=0.1, vmax=0.6,
                    edgecolors="white", linewidths=0.3)
    axes[0].axvspan(0.7, 1.0, alpha=0.07, color=C["red"])
    axes[0].text(0.72, s_arr.max() * 0.90, "Near\nboundary",
                 fontsize=5, color=C["red"])
    axes[0].set_xlabel("Degeneracy $\\mathcal{D} = \\lambda_2/\\lambda_1$ of $S$")
    axes[0].set_ylabel("Observed $\\sin\\theta_{12}/\\|\\Delta S\\|_2$")
    axes[0].text(-0.22, 1.06, "a", transform=axes[0].transAxes,
                 fontsize=9, fontweight="bold", va="top")

    # panel b
    axes[1].scatter(Vp_arr, s_arr, s=14, alpha=0.8,
                    c=D_arr, cmap="RdYlBu_r", vmin=0.2, vmax=0.95,
                    edgecolors="white", linewidths=0.3)
    axes[1].set_xlabel(r"$V_{\perp}$ (variance toward boundary)")
    axes[1].set_ylabel("Observed $\\sin\\theta_{12}/\\|\\Delta S\\|_2$")
    axes[1].text(-0.22, 1.06, "b", transform=axes[1].transAxes,
                 fontsize=9, fontweight="bold", va="top")

    # panel c — A_g vs D with exact 2D identity
    Ag_sm = np.linspace(Ag_arr.min(), Ag_arr.max(), 100)
    axes[2].scatter(Ag_arr, D_arr, s=14, alpha=0.8,
                    color=C["purple"], edgecolors="white", linewidths=0.3)
    axes[2].plot(Ag_sm, (1.0 - Ag_sm) / Ag_sm, "k--", lw=1.0,
                 label="$\\mathcal{D}=(1{-}\\mathcal{A})/\\mathcal{A}$\n(exact in 2D)")
    axes[2].set_xlabel("Additivity $\\mathcal{A}_g$")
    axes[2].set_ylabel("Degeneracy $\\mathcal{D}$")
    axes[2].legend(frameon=False, fontsize=5.5)
    axes[2].text(-0.22, 1.06, "c", transform=axes[2].transAxes,
                 fontsize=9, fontweight="bold", va="top")

    fig.tight_layout(pad=0.5)
    outfile = OUTDIR / "figS3_boundary_signature.pdf"
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)
    return outfile


# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
def run_sanity_checks() -> None:
    """
    Run four algebraic checks that verify the core properties of the framework.

    These checks are NOT statistical tests — they are mathematical identities
    that MUST hold exactly (to numerical precision). If any check fails, there
    is a bug in the implementation.

    WHAT EACH CHECK VERIFIES (plain language)
    -----------------------------------------
    1. Schluter recovery: when the fitness surface is isotropic (H = hI, the
       same curvature in all directions), the genetic channel should coincide
       exactly with g_max — Schluter's line of least resistance. If this check
       fails, the channel operator or the phenotype-space mapping is broken.

    2. I_geom endpoints: the geometric overlap deficit must be exactly 0 when
       the two channel subspaces are perfectly aligned (theta = 0), and exactly
       1 when they are perfectly orthogonal (theta = 90°, k = 2 channels). This
       verifies the definition of I_geom = 1 - prod cos^2(psi_j).

    3. Small-angle approximation: for small theta, I_geom ≈ k*theta^2. The
       ratio exact/approx should be close to 1 for theta = 5°. This verifies
       the Taylor expansion used in the main text to derive the quadratic
       (snowball) scaling.

    4. Davis-Kahan: the observed sin(theta_12) must not exceed the bound
       ||Delta S||_2 / delta. This verifies that the perturbation theory is
       being applied correctly and that the spectral gap is computed from the
       right matrix (S, not G*H).
    """
    print("=== SANITY CHECKS ===")

    # 1. Schluter recovery: when H = h I, the phenotype-space leading channel
    # is parallel to g_max.
    G = np.diag([1.0, 0.5, 0.2])
    H = 2.0 * np.eye(3)
    gmax = spectral_decomp_spd(G)[1][:, 0]
    channel_vec = leading_channel_basis(G, H, k=1, space="phenotype")[:, 0]
    angle_deg = np.degrees(vector_angle(gmax, channel_vec))
    print(f"Schluter recovery angle = {angle_deg:.8f} deg [expected: 0]")

    # 2. Equal-angle incompatibility endpoints.
    I0 = incompatibility_equal_angle(0.0, k=2)
    I90 = incompatibility_equal_angle(np.pi / 2, k=2)
    print(f"I(0, k=2) = {I0:.8f} [expected: 0]")
    print(f"I(90 deg, k=2) = {I90:.8f} [expected: 1]")

    # 3. Small-angle expansion: I ~ k theta^2 for small theta.
    theta_small = np.radians(5.0)
    for k in [1, 2, 4]:
        exact = incompatibility_equal_angle(theta_small, k)
        approx = k * theta_small**2
        ratio = exact / approx
        print(
            f"small-angle check, k={k}: exact={exact:.8f}, approx={approx:.8f}, ratio={ratio:.5f}"
        )

    # 4. Davis–Kahan style behaviour for the symmetric operator.
    # Use an off-diagonal curvature perturbation so the angle is non-zero.
    G = np.diag([1.0, 0.6, 0.3])
    H1 = np.eye(3)
    dH = np.array([[0.00, 0.10, 0.00],
                   [0.10, -0.05, 0.03],
                   [0.00, 0.03, 0.02]])
    H2 = H1 + dH
    min_eval_H2 = np.min(np.linalg.eigvalsh(symmetrise(H2)))
    if min_eval_H2 <= 1e-8:
        H2 = H2 + (abs(min_eval_H2) + 1e-3) * np.eye(3)
    S1 = channel_operator(G, H1)
    S2 = channel_operator(G, H2)
    V1 = leading_channel_basis(G, H1, k=1, space="whitened")
    V2 = leading_channel_basis(G, H2, k=1, space="whitened")
    theta = largest_principal_angle(V1, V2)
    evals_S1, _ = spectral_decomp_spd(S1)
    delta = evals_S1[0] - evals_S1[1]
    bound = np.linalg.norm(S2 - S1, ord=2) / delta
    print(
        f"Davis–Kahan check: sin(theta)={np.sin(theta):.8f}, bound={bound:.8f}, holds={np.sin(theta) <= bound}"
    )


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    fig1   = make_figure1()
    fig2   = make_figure2()
    fig_s1 = make_figS1()
    fig_s2 = make_figS2()
    fig_s3 = make_figS3()
    run_sanity_checks()
    print()
    print(f"Saved: {fig1}")
    print(f"Saved: {fig2}")
    print(f"Saved: {fig_s1}")
    print(f"Saved: {fig_s2}")
    print(f"Saved: {fig_s3}")


if __name__ == "__main__":
    main()
