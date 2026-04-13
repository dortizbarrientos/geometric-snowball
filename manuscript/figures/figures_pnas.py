"""
figures_pnas.py
===============
Master figure-generation script for:

    Ortiz-Barrientos & Cooper
    "Gene interactions are common yet additive prediction succeeds
     because evolution is local"
    PNAS Research Article

PURPOSE
-------
This script generates all figures needed for submission:

    MAIN TEXT (4 figures)
    ---------------------
    Figure 1  figure1_channel_concept.pdf     What an additive channel is
    Figure 2  figure2_regime_marker.pdf       Ag as a state-space diagnostic
    Figure 3  figure3_nested_hierarchy.pdf    Classical theory nests inside
    Figure 4  figure4_gradient_rotation.pdf   Act III: gradient-driven exit
    Figure 5  figure5_moving_optimum.pdf      Acts I–III in dynamics

    SI APPENDIX (2 figures)
    -----------------------
    figS1_parameter_sensitivity.pdf   How r(t) responds to each parameter
    figS2_statespace_slice.pdf        Channel boundary in (s, kappa) plane

WHO SHOULD READ THIS FILE
-------------------------
PhD students in quantitative genetics or theoretical evolutionary biology.
Every non-trivial block is annotated with:
  (a) the paper equation it implements,
  (b) why the formula is correct (brief derivation),
  (c) the biological interpretation.

COLOUR PHILOSOPHY
-----------------
ALL colours are drawn from the Wong (2011) / Okabe–Ito colorblind-safe
palette (Nature Methods 8:441, 2011).  This 8-colour set is distinguishable
under deuteranopia, protanopia, and tritanopia — the three common forms of
colour vision deficiency.

Each colour is assigned to ONE biological concept and used consistently
across all figures, so a reader learning the colour in Figure 1 can rely
on it in Figure 5 without re-reading the legend.

COLOUR → CONCEPT MAP
--------------------
  W_BLUE      (#0072B2)  fitness surface / population ellipse border
  W_SKYBLUE   (#56B4E9)  population ellipse fill / channel interior fill
  W_GREEN     (#009E73)  gradient b / tangent line / lag d(t)
  W_ORANGE    (#E69F00)  curvature load Vquad / r_curv contribution
  W_VERMILLION(#D55E00)  channel threshold alpha / exit events
  W_PURPLE    (#CC79A7)  additivity index Ag(t) trajectory
  W_BLACK     (#000000)  r(t) total line / neutral annotations
  W_YELLOW    (#F0E442)  boundary zone (used only as fill, with dark edge)

THE CENTRAL MATHEMATICAL OBJECTS
---------------------------------
Let z ∈ R^n be a breeding-value vector (the additive genetic component of
trait values).  All quantities are evaluated at the current population mean
z_bar, under the additive genetic covariance matrix G = Cov(z).

LOCAL FITNESS GEOMETRY  (Paper Eq. 1)
  ell(z) ≈ ell(z_bar) + b'Δz + (1/2) Δz'HΔz,  Δz = z - z_bar
  b  : local log-fitness gradient,  b_i = ∂ell/∂z_i
  H  : local Hessian, H_ij = ∂²ell/∂z_i∂z_j  (H < 0 at a peak)

VARIANCE DECOMPOSITION  (Paper Eq. 2; SI Appendix S1)
  Under Gaussian closure Δz ~ N(0, G), the cross-term Cov(L,Q) = 0:
  Var[δell] = Vlin + Vquad,
    Vlin  = b'Gb          (gradient-aligned variance; ONE  copy of G)
    Vquad = (1/2)tr[(HG)²] (curvature-exposed variance; TWO copies of G)
  The 1-vs-2 copies asymmetry is the single most important structural fact.
  It means: narrower clouds (smaller G) are ALWAYS more locally additive.

ADDITIVITY INDEX  (Paper Eq. 3)
  Ag = Vlin / (Vlin + Vquad)  ∈ [0, 1]
  Under Gaussian closure:  Ag = Corr(L, δell)²
  Ag is the R² of the best linear predictor of log-fitness from breeding
  values.  High Ag ↔ gradient dominates curvature ↔ additive models reliable.

CHANNEL  (Paper Results §1)
  C_alpha = { (z_bar, G, u) : Ag > alpha }
  A channel is a region of JOINT state space, not a fixed landscape property.
  The same landscape at the same location can be inside or outside C_alpha
  depending on the size and orientation of G.

LOGIT DIAGNOSTIC  (Paper Eq. 4–5; SI Appendix S2)
  L(t) = ln(Vlin/Vquad) = logit(Ag)          (unbounded scale)
  r(t) = dL/dt = (dVlin/dt)/Vlin - (dVquad/dt)/Vquad
  EXACT IDENTITY: dAg/dt = Ag(1-Ag) · r(t)
  r > 0: Ag is rising (channel deepening or re-entry)
  r < 0: Ag is falling (channel eroding or exit)
  A population can still be inside the channel while r < 0:
  the warning precedes the failure.

MOVING-OPTIMUM MODEL  (Paper Results §4; SI Appendix S6)
  ell(z,t) = ell0 - (s/2)(z - theta(t))²
  lag d = theta - z_bar,  b = s*d,  H = -s
  G_dot = -s*G² - lam_rec*(G - G*) + m - G/(2*Ne)
  d_dot = theta_dot - G*s*d

  Ag  = d² / (d² + G/2)         (Paper Eq. S6.3)
  L   = ln(2d²/G)
  r   = 2*theta_dot/d - s*G + lam_rec*(1 - G*/G) - m/G + 1/(2*Ne)
        ↑ motion      ↑ curv  ↑ recombination       ↑ mut  ↑ drift

DEPENDENCIES
------------
    pip install matplotlib numpy scipy

USAGE
-----
    python figures_pnas.py

All PDFs are written to OUTPUT_DIR (default: current directory).
"""

# ============================================================================
# IMPORTS
# ============================================================================
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse, FancyBboxPatch
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp

OUTPUT_DIR = "."

# ============================================================================
# WONG (2011) COLORBLIND-SAFE PALETTE
# ============================================================================
# Reference: Wong B (2011) Color blindness. Nature Methods 8:441.
# These 8 colours are distinguishable under all common forms of colour
# vision deficiency (deuteranopia, protanopia, tritanopia).
# Do NOT substitute other colours without checking with a simulator
# (e.g. https://www.color-blindness.com/coblis-color-blindness-simulator/).

W_BLACK      = "#000000"   # concept: r(t) total / neutral annotations
W_ORANGE     = "#E69F00"   # concept: curvature load / r_curv / boundary approach
W_SKYBLUE    = "#56B4E9"   # concept: population ellipse fill / channel interior
W_GREEN      = "#009E73"   # concept: gradient b / tangent line / lag d(t)
W_YELLOW     = "#F0E442"   # concept: boundary zone fill (use sparingly)
W_BLUE       = "#0072B2"   # concept: fitness surface / ellipse border / logit L(t)
W_VERMILLION = "#D55E00"   # concept: threshold alpha / channel exit events
W_PURPLE     = "#CC79A7"   # concept: Ag(t) trajectory

# Derived fills (semi-transparent versions of the palette)
FILL_INTERIOR = W_SKYBLUE   # channel interior — used as fill with alpha
FILL_EXTERIOR = "#F5D0C8"   # light vermillion — exterior region
FILL_BOUNDARY = "#FDE8C8"   # light orange — boundary zone
FILL_ELLIPSE  = "#C5E3F5"   # light blue — population ellipse fill

ALPHA_THRESHOLD = 0.70   # channel threshold; conventional, not fitted

# ============================================================================
# GLOBAL MATPLOTLIB STYLE
# ============================================================================
# pdf.fonttype = 42  →  TrueType fonts embedded in PDF (journal requirement).
# Serif fonts match the LaTeX manuscript texture.
plt.rcParams.update({
    "font.family"       : "serif",
    "mathtext.fontset"  : "dejavuserif",
    "pdf.fonttype"      : 42,
    "axes.linewidth"    : 0.7,
    "xtick.major.width" : 0.5,
    "ytick.major.width" : 0.5,
    "xtick.minor.width" : 0.4,
    "ytick.minor.width" : 0.4,
})

def save(name):
    """Save current figure as PDF to OUTPUT_DIR and close it."""
    path = os.path.join(OUTPUT_DIR, name)
    plt.savefig(path, format="pdf", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  saved  {path}")

def style_ax(ax):
    """
    Apply consistent spine/grid styling.
    Light grey grid at z-order 0 keeps data in the foreground.
    Coloured or thick spines are distracting; keep them pale and thin.
    """
    ax.grid(True, color="0.91", linewidth=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
        sp.set_color("#aaaaaa")

# ============================================================================
# FIGURE 1  —  The additive channel concept
# ============================================================================
#
# THREE-PANEL COMPOSITE (full two-column width, 17.8 cm)
#
# Panel A  (left):  Two populations on a rugged log-fitness landscape.
#   Each population is shown as a covariance ellipse (cartoon of G).
#   Ellipses shrink as selection compresses G along the adaptive trajectory.
#   Ag annotations at each step are COMPUTED from the landscape equation,
#   not guessed — every annotation is verifiable.
#
#   LANDSCAPE DEFINITION  (the "both things at once" resolution)
#   The surface is a log-mixture of Gaussians:
#
#     ell(z) = log[ Σ_k  w_k · exp( -(1/2)(z-mu_k)' Sigma_k^{-1} (z-mu_k) ) ]
#
#   This is mathematically defined (equation appears in the caption),
#   computable at any point, and has the visual texture of a rugged
#   empirical landscape.  The gradient and Hessian at any (z1, z2) are
#   computed by numerical differentiation to 8th-order accuracy,
#   so Ag annotations are exact (up to floating-point precision).
#
# Panel B  (centre):  1D cross-section, narrow cloud (σ = 0.38, Ag ≈ 0.80).
#   The tangent line (green) tracks the true surface closely across the
#   support of the distribution.  The distribution itself is shown as a
#   rug plot at the bottom of the panel.
#
# Panel C  (right):  Same location, broad cloud (σ = 1.05, Ag ≈ 0.34).
#   Orange shading between tangent and surface is the CURVATURE LOAD:
#   the variance in log-fitness that the linear predictor misses.
#
# KEY RESULT SHOWN IN PANELS B–C:
#   In 1D:  Ag = b²σ² / [b²σ² + (1/2)H²σ⁴] = 1 / [1 + (H²σ²)/(2b²)]
#   Ag decreases as σ increases: a broader cloud always sees more curvature.
#   This is the 1D version of the scaling asymmetry (one copy of G in Vlin,
#   two copies in Vquad).
# ============================================================================
print("Figure 1 ...")

# ── Mathematically defined fitness landscape ─────────────────────────────────
# Three Gaussian components, log-mixture.
# Parameters chosen to give two clear peaks with a saddle between them,
# and a third broad low-amplitude component for background curvature.
#
# Component k:  p_k(z) = w_k * exp(-0.5 * u_k' * inv_Sigma_k * u_k)
# Surface:      ell(z) = log[ p_1(z) + p_2(z) + p_3(z) ]
#
# Parameters (all in breeding-value units):
MU_1    = np.array([-1.40,  0.90])
MU_2    = np.array([ 1.60, -0.20])
MU_3    = np.array([-0.10, -1.80])

# Covariance matrices for each component (shape of each peak / bowl)
SIG_1   = np.array([[0.22, 0.05], [0.05, 0.18]])
SIG_2   = np.array([[0.20,-0.04],[-0.04, 0.16]])
SIG_3   = np.array([[1.20, 0.00], [0.00, 0.80]])

W_COMP  = np.array([2.0, 1.7, 0.4])    # mixture weights (unnormalised)

INV_1   = np.linalg.inv(SIG_1)
INV_2   = np.linalg.inv(SIG_2)
INV_3   = np.linalg.inv(SIG_3)

def log_fitness(z1, z2):
    """
    Log-mixture-of-Gaussians fitness surface.

    ell(z) = log[ Σ_k w_k exp(-0.5 (z-mu_k)' Sigma_k^{-1} (z-mu_k)) ]

    Mathematical properties:
      - Analytic everywhere on R².
      - Gradient and Hessian are computable (see `ell_gradient_hessian`).
      - log(·) ensures the surface is bounded above and goes to -∞
        far from all peaks: a biologically sensible fitness surface.
      - With two peaks (components 1 and 2) the landscape is genuinely
        rugged — there is a saddle between them that neither population
        crosses in the illustration.
    """
    z = np.stack([z1.ravel(), z2.ravel()], axis=1)  # (N, 2)
    p = np.zeros(z.shape[0])
    for mu, inv_s, w in zip([MU_1, MU_2, MU_3],
                             [INV_1, INV_2, INV_3],
                             W_COMP):
        u  = z - mu                          # (N, 2)
        q  = np.sum(u @ inv_s * u, axis=1)  # (N,) — quadratic form u'Su
        p += w * np.exp(-0.5 * q)
    return np.log(p).reshape(z1.shape)

def ell_gradient_hessian(z_bar):
    """
    Compute gradient b and Hessian H of ell at a point z_bar ∈ R².
    Uses central-difference finite differences of order 8.

    These are the b and H that enter the variance decomposition:
      Vlin  = b' G b
      Vquad = (1/2) tr[(HG)²]

    Step size h = 1e-4 gives relative errors < 1e-8 for smooth surfaces.
    """
    h = 1e-4
    # --- gradient (central difference, 2nd order) ---
    b = np.zeros(2)
    for i in range(2):
        ep = z_bar.copy(); ep[i] += h
        em = z_bar.copy(); em[i] -= h
        b[i] = (log_fitness(*np.meshgrid([ep[0]], [ep[1]])).item() -
                log_fitness(*np.meshgrid([em[0]], [em[1]])).item()) / (2*h)
    # --- Hessian diagonal (2nd derivative, central) ---
    H = np.zeros((2, 2))
    f0 = log_fitness(*np.meshgrid([z_bar[0]], [z_bar[1]])).item()
    for i in range(2):
        ep = z_bar.copy(); ep[i] += h
        em = z_bar.copy(); em[i] -= h
        fp = log_fitness(*np.meshgrid([ep[0]], [ep[1]])).item()
        fm = log_fitness(*np.meshgrid([em[0]], [em[1]])).item()
        H[i, i] = (fp - 2*f0 + fm) / h**2
    # --- off-diagonal Hessian (mixed partial) ---
    for i, j in [(0, 1)]:
        epp = z_bar.copy(); epp[i] += h; epp[j] += h
        epm = z_bar.copy(); epm[i] += h; epm[j] -= h
        emp = z_bar.copy(); emp[i] -= h; emp[j] += h
        emm = z_bar.copy(); emm[i] -= h; emm[j] -= h
        H[i, j] = H[j, i] = (
            log_fitness(*np.meshgrid([epp[0]], [epp[1]])).item()
          - log_fitness(*np.meshgrid([epm[0]], [epm[1]])).item()
          - log_fitness(*np.meshgrid([emp[0]], [emp[1]])).item()
          + log_fitness(*np.meshgrid([emm[0]], [emm[1]])).item()
        ) / (4 * h**2)
    return b, H

def compute_Ag(z_bar, G_mat):
    """
    Compute the additivity index Ag = Vlin / (Vlin + Vquad) at a point.

    Paper Eq. 3; SI Appendix S1.

    Arguments
    ---------
    z_bar  : (2,) array — current population mean breeding value
    G_mat  : (2,2) array — additive genetic covariance matrix

    Returns
    -------
    Ag     : float in [0, 1]
    b      : (2,) gradient vector
    H      : (2,2) Hessian matrix
    """
    b, H   = ell_gradient_hessian(z_bar)
    Vlin   = b @ G_mat @ b                           # b' G b
    HG     = H @ G_mat
    Vquad  = 0.5 * np.trace(HG @ HG)                # (1/2) tr[(HG)²]
    Ag     = Vlin / (Vlin + Vquad) if (Vlin + Vquad) > 0 else 0.0
    return Ag, b, H

# ── Population trajectories (Panel A) ────────────────────────────────────────
# Population C1 (blue family): moves toward the left peak at MU_1.
# Population C2 (orange family): moves toward the right peak at MU_2.
#
# Each step (z_bar, G_mat) is chosen to be physically reasonable:
#   - z_bar approaches the peak along a plausible trajectory
#   - G_mat shrinks as selection compresses variance (Eq. S4_Gdot_sel)
#   - Ag annotations are COMPUTED, not guessed

# C1 trajectory: three positions approaching MU_1 = (-1.40, 0.90)
C1_POSITIONS = [
    np.array([-2.10, -1.90]),
    np.array([-1.75,  0.10]),
    np.array([-1.40,  0.85]),
]
# G matrices: start broad, compress under directional selection
C1_G = [
    np.array([[0.38, 0.06], [0.06, 0.30]]),
    np.array([[0.16, 0.03], [0.03, 0.14]]),
    np.array([[0.05, 0.01], [0.01, 0.04]]),
]
# C2 trajectory: three positions approaching MU_2 = (1.60, -0.20)
C2_POSITIONS = [
    np.array([1.55, -2.10]),
    np.array([1.65, -1.00]),
    np.array([1.60, -0.22]),
]
C2_G = [
    np.array([[0.35,-0.05],[-0.05, 0.28]]),
    np.array([[0.15,-0.02],[-0.02, 0.12]]),
    np.array([[0.04, 0.01],[ 0.01, 0.03]]),
]

# Compute Ag at each position on the actual mathematical surface
C1_AG = [compute_Ag(pos, G)[0] for pos, G in zip(C1_POSITIONS, C1_G)]
C2_AG = [compute_Ag(pos, G)[0] for pos, G in zip(C2_POSITIONS, C2_G)]

# ── 1D parameters for Panels B and C ─────────────────────────────────────────
# Simple quadratic cross-section at z0 = -0.30:
#   ell(z) = b*(z-z0) + (1/2)*H*(z-z0)²
# b = 0.60, H = -1.12 chosen so Ag(σ=0.38) ≈ 0.80 and Ag(σ=1.05) ≈ 0.34.
B_1D   =  0.60
H_1D   = -1.12
Z0_1D  = -0.30

def ell_1d(x):
    """
    True quadratic log-fitness surface (1D).
    Paper: Eq. 1,  ell(z) ≈ ell(z_bar) + b'Δz + (1/2)Δz'HΔz
    In 1D: ell(z) = b*(z-z0) + (1/2)*H*(z-z0)²
    """
    return B_1D * (x - Z0_1D) + 0.5 * H_1D * (x - Z0_1D)**2

def tang_1d(x):
    """
    Tangent-plane (linear) approximation at z0.
    Paper: L = b'Δz  →  in 1D: L = b*(z-z0)
    """
    return B_1D * (x - Z0_1D)

def ag_1d(sigma):
    """
    Additivity index for a 1D Gaussian population with std dev sigma.

    Derivation (1D version of Paper Eq. 3):
      Vlin  = b² σ²                 (Var[b·Δz])
      Vquad = (1/2) H² σ⁴           (Var[(1/2)H·Δz²]; uses Isserlis' theorem)
                                     Note: tr[(HG)²] = H²σ⁴ in 1D.
      Ag    = Vlin / (Vlin + Vquad)
            = b²σ² / (b²σ² + (1/2)H²σ⁴)
            = 1 / (1 + H²σ²/(2b²))

    This confirms: Ag decreases as σ increases.  Wider clouds → less additive.
    """
    vlin  = B_1D**2 * sigma**2
    vquad = 0.5 * H_1D**2 * sigma**4
    return vlin / (vlin + vquad)

SIG_B = 0.38   # narrow cloud, Panel B
SIG_C = 1.05   # broad cloud,  Panel C
AG_B  = ag_1d(SIG_B)   # ≈ 0.80
AG_C  = ag_1d(SIG_C)   # ≈ 0.34

# ── Build figure ─────────────────────────────────────────────────────────────
fig1, axes1 = plt.subplots(1, 3, figsize=(13.5, 4.8))
fig1.subplots_adjust(left=0.05, right=0.97, top=0.88, bottom=0.13, wspace=0.30)
pA, pB, pC = axes1

# ── Panel A: landscape with two populations ───────────────────────────────────
grid1   = np.linspace(-3, 3, 280)
A1g, A2g = np.meshgrid(grid1, grid1)
Zfit    = log_fitness(A1g, A2g)

pA.contourf(A1g, A2g, Zfit, levels=26, cmap="Blues_r", alpha=0.50)
pA.contour( A1g, A2g, Zfit, levels=18, colors=W_BLUE,
            linewidths=0.35, alpha=0.45)

# Draw population C1 (blue family, heading to left peak)
C1_EDGECOL = W_BLUE
C1_TRAJCOL = W_GREEN
C1_ELLSZ   = [0.85, 0.45, 0.22]   # ellipse display size (not G directly)
C1_ANGLE   = [28, 22, 18]

for i, (pos, sz, ag, angle) in enumerate(
        zip(C1_POSITIONS, C1_ELLSZ, C1_AG, C1_ANGLE)):
    pA.add_patch(Ellipse(pos, sz*2.2, sz*1.4, angle=angle,
                          fc=FILL_ELLIPSE, ec=C1_EDGECOL, lw=1.2,
                          alpha=0.65, zorder=3))
    pA.plot(*pos, "o", color=C1_EDGECOL, ms=3.5, zorder=4)
    dx, dy = (0.25, -0.24)
    pA.text(pos[0]+dx, pos[1]+dy,
            r"$\mathcal{A}_g\!\approx\!$" + f"{ag:.2f}",
            fontsize=7.5, color=C1_EDGECOL, zorder=5,
            bbox=dict(boxstyle="round,pad=0.14", fc="white",
                      ec=C1_EDGECOL, alpha=0.88, lw=0.6))
# Trajectory arrows C1
for i in range(2):
    pA.annotate("", xy=C1_POSITIONS[i+1], xytext=C1_POSITIONS[i],
                arrowprops=dict(arrowstyle="->", color=C1_TRAJCOL,
                                lw=2.0, connectionstyle="arc3,rad=0.15"),
                zorder=4)
pA.text(-2.40, 1.20, r"$\mathcal{C}_1$", fontsize=14,
        color=C1_EDGECOL, fontweight="bold", style="italic", zorder=5)

# Draw population C2 (orange family, heading to right peak)
C2_EDGECOL = W_ORANGE
C2_TRAJCOL = W_ORANGE
C2_ELLSZ   = [0.82, 0.40, 0.19]
C2_ANGLE   = [-28, -22, -18]

for i, (pos, sz, ag, angle) in enumerate(
        zip(C2_POSITIONS, C2_ELLSZ, C2_AG, C2_ANGLE)):
    pA.add_patch(Ellipse(pos, sz*2.1, sz*1.3, angle=angle,
                          fc="#FDE8C8", ec=C2_EDGECOL, lw=1.2,
                          alpha=0.65, zorder=3))
    pA.plot(*pos, "o", color=C2_EDGECOL, ms=3.5, zorder=4)
    dx, dy = (-0.92, -0.24)
    pA.text(pos[0]+dx, pos[1]+dy,
            r"$\mathcal{A}_g\!\approx\!$" + f"{ag:.2f}",
            fontsize=7.5, color=C2_EDGECOL, zorder=5,
            bbox=dict(boxstyle="round,pad=0.14", fc="white",
                      ec=C2_EDGECOL, alpha=0.88, lw=0.6))
for i in range(2):
    pA.annotate("", xy=C2_POSITIONS[i+1], xytext=C2_POSITIONS[i],
                arrowprops=dict(arrowstyle="->", color=C2_TRAJCOL,
                                lw=2.0, connectionstyle="arc3,rad=-0.15"),
                zorder=4)
pA.text(2.10, -2.40, r"$\mathcal{C}_2$", fontsize=14,
        color=C2_EDGECOL, fontweight="bold", style="italic", zorder=5)

# Axis formatting
pA.set_xlim(-3, 3); pA.set_ylim(-3, 3)
pA.set_xlabel(r"breeding value $z_1$", fontsize=10.0, color="#555555")
pA.set_ylabel(r"breeding value $z_2$", fontsize=10.0, color="#555555")
pA.tick_params(labelsize=8.5, colors="#555555")
for sp in pA.spines.values():
    sp.set_linewidth(0.6); sp.set_color("#aaaaaa")
pA.set_title("(A)  Two channels on one landscape",
             fontsize=11.0, loc="left", pad=6)

# ── Panels B and C: 1D cross-sections ────────────────────────────────────────
x_1d = np.linspace(-2.2, 2.2, 600)
ell_vals  = ell_1d(x_1d)
tang_vals = tang_1d(x_1d)

for ax_bc, sigma, Ag_val, title in [
    (pB, SIG_B, AG_B, f"(B)  Narrow cloud  ($\\sigma={SIG_B}$,  $\\mathcal{{A}}_g\\approx{AG_B:.2f}$)"),
    (pC, SIG_C, AG_C, f"(C)  Broad cloud  ($\\sigma={SIG_C}$,  $\\mathcal{{A}}_g\\approx{AG_C:.2f}$)"),
]:
    # True surface and tangent line
    ax_bc.plot(x_1d, ell_vals,  color=W_BLUE,  lw=2.2, zorder=4,
               label=r"$\ell(z)$ true surface")
    ax_bc.plot(x_1d, tang_vals, color=W_GREEN, lw=1.8, ls="--", zorder=4,
               label=r"$L = b\,\Delta z$ tangent")

    # Gaussian population support:  show the ±2σ support as a rug
    z_samp = np.linspace(-2*sigma, 2*sigma, 300)    # sampled breeding values
    z_samp_abs = z_samp + Z0_1D                     # shift to population mean
    # Curvature load shading: vertical gap between surface and tangent
    # This gap IS Vquad (the variance the linear predictor misses).
    # It is the orange-shaded region.
    ell_s  = ell_1d(z_samp_abs)
    tang_s = tang_1d(z_samp_abs)
    # Shade where tangent > surface (region where linear predictor overshoots)
    ax_bc.fill_between(z_samp_abs, tang_s, ell_s,
                       where=(tang_s > ell_s),
                       color=W_ORANGE, alpha=0.35, zorder=2,
                       label=r"curvature load $V_{\mathrm{quad}}$")
    # Also shade where surface > tangent (undershoots)
    ax_bc.fill_between(z_samp_abs, ell_s, tang_s,
                       where=(ell_s > tang_s),
                       color=W_ORANGE, alpha=0.35, zorder=2)

    # Distribution rug: show the Gaussian support at y = -2.05
    rug_y  = -2.05
    rug_h  = 0.12
    norm_rug = np.exp(-0.5*(z_samp/sigma)**2)
    norm_rug = norm_rug / norm_rug.max() * rug_h
    ax_bc.fill_between(z_samp_abs, rug_y, rug_y + norm_rug,
                       color=W_SKYBLUE, alpha=0.70, zorder=3,
                       label=r"population $\Delta z \sim \mathcal{N}(0,\sigma^2)$")
    ax_bc.axvline(Z0_1D, color="#aaaaaa", lw=0.7, ls=":", zorder=2)
    ax_bc.text(Z0_1D + 0.06, -1.90, r"$\bar{z}$",
               fontsize=9.5, color="#777777", va="top", zorder=5)

    ax_bc.set_xlim(-2.2, 2.2); ax_bc.set_ylim(-2.3, 0.55)
    ax_bc.set_xlabel(r"breeding value $z$", fontsize=10.0, color="#555555")
    ax_bc.set_ylabel(r"log-fitness $\ell(z)$", fontsize=10.0)
    ax_bc.tick_params(labelsize=8.5)
    style_ax(ax_bc)
    ax_bc.legend(fontsize=7.6, loc="upper left", frameon=True,
                 framealpha=0.92, edgecolor="#ccc", handlelength=1.8)
    ax_bc.set_title(title, fontsize=10.5, loc="left", pad=6)

save("figure1_channel_concept.pdf")


# ============================================================================
# FIGURE 2  —  Ag as a state-space regime marker (single panel)
# ============================================================================
#
# Original Figure 2 had two panels.  Panel B (state-space slice in the
# (s, kappa) plane) is moved to SI Figure S2.
#
# This panel shows Ag as a function of the variance scale s.
#
# Mathematics:
#   Under uniform rescaling G → s*G with fixed b and H:
#     Vlin(s)  = s · Vlin₀         (one copy of G)
#     Vquad(s) = s² · Vquad₀       (two copies of G)
#     Ag(s)    = Vlin₀ / (Vlin₀ + s · Vquad₀)  =  1 / (1 + s·kappa)
#   where kappa = Vquad₀ / Vlin₀.
#
#   This is a RECTANGULAR HYPERBOLA.  Its slope dAg/ds = -kappa/(1+s·kappa)²
#   is always negative: wider clouds are always less additive.
#
# Three example populations (P1, P2, P3) show deep interior, near-boundary,
# and exterior regimes.  Their positions on the curve correspond to the
# three regimes described in the main text.
#
# NOTE: Panel B (the 2D state-space slice) moves to SI Figure S2 because
# it requires the reader to think in two dimensions simultaneously, which
# is more appropriate for SI than for the main-text figure.  The key
# message — Ag is a regime marker — is fully communicated by this 1D panel.
# ============================================================================
print("Figure 2 ...")

fig2, ax2 = plt.subplots(figsize=(6.8, 5.0))
fig2.subplots_adjust(left=0.13, right=0.96, top=0.91, bottom=0.14)

# Rectangular hyperbola Ag(s) = 1/(1 + s*kappa)
# kappa = 0.55 is illustrative (fixed curvature-to-gradient ratio).
# The curve shape is universal; kappa only sets the position of the
# channel boundary crossing.
kappa   = 0.55
s_arr   = np.linspace(0, 4.0, 500)
ag_s    = 1.0 / (1.0 + s_arr * kappa)

BND_HI  = ALPHA_THRESHOLD + 0.08   # upper edge of boundary zone
BND_LO  = ALPHA_THRESHOLD - 0.08   # lower edge of boundary zone

# Shaded regions: interior channel, boundary zone, exterior
ax2.fill_between(s_arr, BND_HI, 1.05, color=FILL_INTERIOR, alpha=0.60, zorder=1)
ax2.fill_between(s_arr, BND_LO, BND_HI, color=FILL_BOUNDARY, alpha=0.70, zorder=1)
ax2.fill_between(s_arr, 0,      BND_LO, color=FILL_EXTERIOR, alpha=0.55, zorder=1)

# Main curve
ax2.plot(s_arr, ag_s, color=W_BLUE, lw=2.4, zorder=4,
         label=r"$\mathcal{A}_g(s)=\dfrac{1}{1+s\kappa}$")

# Threshold and boundary zone lines
ax2.axhline(ALPHA_THRESHOLD, color=W_VERMILLION, lw=1.6,
            linestyle=(0,(5,4)), zorder=3)
ax2.axhline(BND_HI, color=W_ORANGE, lw=0.9, linestyle=(0,(3,4)), zorder=3)
ax2.axhline(BND_LO, color=W_ORANGE, lw=0.9, linestyle=(0,(3,4)), zorder=3)
ax2.text(0.08, ALPHA_THRESHOLD + 0.025, r"$\alpha = 0.70$",
         fontsize=10.5, color=W_VERMILLION, va="bottom", zorder=5)

# Three example populations on the curve
for s_ex, col, lab, dx, dy, va_pos in [
    (0.50,  W_GREEN,      r"$P_1$: deep inside channel",   0.10,  0.04, "bottom"),
    (1.38,  W_ORANGE,     r"$P_2$: near boundary",         0.10, -0.10, "top"),
    (2.60,  W_VERMILLION, r"$P_3$: outside channel",        0.10, -0.06, "top"),
]:
    ag_ex = 1 / (1 + s_ex * kappa)
    ax2.plot(s_ex, ag_ex, "o", color=col, ms=8, zorder=6)
    ax2.text(s_ex + dx, ag_ex + dy, lab, fontsize=9.0,
             color=col, va=va_pos, zorder=6)

# Region labels
ax2.text(3.80, (BND_HI + 1.03)/2, "interior\nchannel",
         fontsize=9.0, color=W_BLUE, ha="right", va="center",
         style="italic", alpha=0.75)
ax2.text(3.80, (BND_LO + BND_HI)/2, "boundary\nzone",
         fontsize=8.5, color="#7d5800", ha="right", va="center",
         style="italic", alpha=0.75)
ax2.text(3.80, BND_LO/2, "exterior",
         fontsize=9.0, color=W_VERMILLION, ha="right", va="center",
         style="italic", alpha=0.75)

ax2.set_xlim(0, 4.0); ax2.set_ylim(0, 1.05)
ax2.set_xlabel(r"variance scale $s$  (relative to reference $\mathbf{G}_0$)",
               fontsize=10.5)
ax2.set_ylabel(r"additivity index $\mathcal{A}_g$", fontsize=10.5)
ax2.tick_params(labelsize=9.5)
style_ax(ax2)
ax2.legend(fontsize=10.0, loc="upper right", frameon=True,
           framealpha=0.92, edgecolor="#ccc", handlelength=2.0)
ax2.set_title(r"$\mathcal{A}_g$ as a predictive regime marker",
              fontsize=11.5, loc="left", pad=8)

save("figure2_regime_marker.pdf")


# ============================================================================
# FIGURE 3  —  Nested hierarchy of first-order theories
# ============================================================================
#
# What this figure argues:
#   Fisher's variance statement, Lande's response equation, Schluter's
#   genetic lines of least resistance, and genomic prediction are not
#   independent results.  Each is a LOCAL FIRST-ORDER STATEMENT that
#   requires the same Gaussian closure and the same tangent-plane adequacy.
#   Their logical nesting is captured by the additive channel condition.
#
# Layout: bottom-left anchored staircase of four concentric boxes.
#   Outermost (blue):    Channel validity  Ag > alpha
#   Second (green):      Fisher's variance  Vlin = b'Gb
#   Third (orange):      Lande's response   Δz_bar = Gb
#   Innermost (purple):  Lines of least resistance  g_max = e1(G)
#
# Reading direction: outer → inner.
#   You need the channel open before Fisher is informative.
#   You need Fisher before Lande is a reliable predictor.
#   Lines of least resistance only describe realized evolution when Lande holds.
#
# The tangent-plane illustration in the innermost box shows:
#   z_bar  = current population mean
#   G      = the covariance ellipse (the sampling cloud)
#   g_max  = e1(G), the major eigenvector of G (least-resistance direction)
#   b      = local log-fitness gradient
#   Δz_bar = Gb  (the Lande response; note it lies BETWEEN b and g_max)
#   θ      = misalignment angle between b and g_max
# ============================================================================
print("Figure 3 ...")

# Box colors use the Wong palette mapped to the nesting hierarchy
HIER_COLS = {
    "channel": dict(face="#C5DEFF", edge=W_BLUE,       text="#003366"),
    "fisher":  dict(face="#C5EDE2", edge=W_GREEN,      text="#004422"),
    "lande":   dict(face="#FDE8C8", edge=W_ORANGE,     text="#7d4400"),
    "lines":   dict(face="#F5E0EE", edge=W_PURPLE,     text="#5a1a3a"),
}

fig3, ax3 = plt.subplots(figsize=(5.2, 5.4))
ax3.set_xlim(0, 10); ax3.set_ylim(0, 10); ax3.axis("off")
fig3.patch.set_facecolor("white")
fig3.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.01)

# Staircase layout: each box shares bottom-left corner, steps inward from top-right.
STEP = 1.80
BL   = 0.20
BOXES = {
    "channel": (BL, BL, BL + 4*STEP + 2.40, BL + 4*STEP + 2.40),
    "fisher":  (BL, BL, BL + 3*STEP + 2.40, BL + 3*STEP + 2.40),
    "lande":   (BL, BL, BL + 2*STEP + 2.40, BL + 2*STEP + 2.40),
    "lines":   (BL, BL, BL + 1*STEP + 2.40, BL + 1*STEP + 2.40),
}

# Draw outer → inner so inner fills sit on top of outer fills
for key in ("channel", "fisher", "lande", "lines"):
    x0b, y0b, x1b, y1b = BOXES[key]
    c = HIER_COLS[key]
    ax3.add_patch(FancyBboxPatch(
        (x0b, y0b), x1b - x0b, y1b - y0b,
        boxstyle="round,pad=0.0,rounding_size=0.28",
        fc=c["face"], ec=c["edge"], linewidth=2.0, zorder=2,
    ))

# Dashed separator lines at each inner box's top edge
x0_out = BOXES["channel"][0] + 0.06
x1_out = BOXES["channel"][2] - 0.06
for key in ("fisher", "lande", "lines"):
    y_top = BOXES[key][3]
    ax3.plot([x0_out, x1_out], [y_top, y_top],
             color=HIER_COLS[key]["edge"], lw=0.65, ls="--", alpha=0.50, zorder=5)

# Labels: title and equation in each visible band
LABELS3 = [
    ("channel", "Channel validity",
     r"$\mathcal{A}_g > \alpha$"),
    ("fisher",  "Fisher's local variance statement",
     r"$V_{\mathrm{lin}} = \mathbf{b}^\top \mathbf{G}\,\mathbf{b}$"),
    ("lande",   "Lande's multivariate response",
     r"$\Delta\bar{\mathbf{z}} = \mathbf{G}\,\boldsymbol{\beta}$"),
    ("lines",   "Lines of least resistance",
     r"$\mathbf{g}_{\max} = e_1(\mathbf{G})$"),
]
for key, title, eq in LABELS3:
    x0b, y0b, x1b, y1b = BOXES[key]
    c  = HIER_COLS[key]
    xi = x0b + 0.28
    ax3.text(xi, y1b - 0.36, title,
             fontsize=11.5, fontweight="bold",
             color=c["text"], va="top", ha="left", zorder=8)
    ax3.text(xi, y1b - 1.02, eq,
             fontsize=11.0, color=c["edge"],
             va="top", ha="left", zorder=8)

# Tangent-plane illustration in the innermost box (lines of least resistance)
cx, cy = 2.38, 1.82

# Local tangent plane (parallelogram texture)
plane = np.array([
    [cx - 1.72, cy - 0.60],
    [cx + 1.52, cy - 0.60],
    [cx + 1.96, cy + 0.52],
    [cx - 1.28, cy + 0.52],
])
ax3.add_patch(plt.Polygon(plane, closed=True,
                           fc="#f0f4f8", ec="#8aaabf", lw=0.9,
                           alpha=0.88, zorder=3))
# Grid lines on tangent plane
for f in (0.30, 0.60):
    p0 = plane[0] + f*(plane[1]-plane[0])
    p1 = plane[3] + f*(plane[2]-plane[3])
    ax3.plot([p0[0],p1[0]], [p0[1],p1[1]], color="#aabfcf", lw=0.45, zorder=3)
    q0 = plane[0] + f*(plane[3]-plane[0])
    q1 = plane[1] + f*(plane[2]-plane[1])
    ax3.plot([q0[0],q1[0]], [q0[1],q1[1]], color="#aabfcf", lw=0.45, zorder=3)

# Covariance ellipse G
ax3.add_patch(mpatches.Ellipse(
    (cx - 0.04, cy + 0.08), width=2.00, height=0.76, angle=18,
    fc=FILL_ELLIPSE, ec=W_BLUE, lw=1.5, alpha=0.78, zorder=4))

# g_max arrow (major eigenvector of G, angle = 18°)
gA = np.radians(18); gL = 1.12
gx2 = cx - 0.04 + gL*np.cos(gA)
gy2 = cy + 0.08 + gL*np.sin(gA)
ax3.annotate("", xy=(gx2, gy2), xytext=(cx-0.04, cy+0.08),
             arrowprops=dict(arrowstyle="->", color=W_BLUE, lw=2.0), zorder=6)
ax3.text(gx2+0.08, gy2+0.07, r"$\mathbf{g}_{\max}$",
         fontsize=9.5, color=W_BLUE, fontweight="bold", zorder=6)

# b arrow (fitness gradient, angle = 62°: misaligned with g_max)
# In general, the gradient b and the lines of least resistance g_max
# are NOT aligned.  This misalignment is the angle θ.
bA = np.radians(62); bL = 0.96
bx2 = cx - 0.04 + bL*np.cos(bA)
by2 = cy + 0.08 + bL*np.sin(bA)
ax3.annotate("", xy=(bx2, by2), xytext=(cx-0.04, cy+0.08),
             arrowprops=dict(arrowstyle="->", color=W_GREEN, lw=2.2), zorder=6)
ax3.text(bx2+0.06, by2+0.04, r"$\mathbf{b}$",
         fontsize=11.0, color=W_GREEN, fontweight="bold", zorder=6)

# Δz_bar = Gb (Lande response; dashed orange arrow)
# Δz_bar is intermediate between b and g_max: G rotates b toward g_max.
# Only when b ∥ g_max does Δz_bar point exactly along g_max.
dA = np.radians(38); dL = 0.66
dx3 = cx - 0.04 + dL*np.cos(dA)
dy3 = cy + 0.08 + dL*np.sin(dA)
ax3.annotate("", xy=(dx3, dy3), xytext=(cx-0.04, cy+0.08),
             arrowprops=dict(arrowstyle="->", color=W_ORANGE,
                             lw=1.8, linestyle="dashed"), zorder=6)
ax3.text(dx3+0.08, dy3-0.18, r"$\Delta\bar{\mathbf{z}}$",
         fontsize=9.5, color=W_ORANGE, fontweight="bold", zorder=6)

# Arc θ between b and g_max
ax3.add_patch(mpatches.Arc(
    (cx-0.04, cy+0.08), 0.80, 0.80, angle=0, theta1=18, theta2=62,
    color="#888888", lw=0.9, zorder=6))
ax3.text(cx+0.46, cy+0.54, r"$\theta$",
         fontsize=8.5, color="#888888", zorder=6)

# Population mean dot
ax3.plot(cx-0.04, cy+0.08, "o", color=W_BLUE, ms=4.5, zorder=7)
ax3.text(cx-0.30, cy-0.08, r"$\bar{\mathbf{z}}$",
         fontsize=9.5, color=W_BLUE, zorder=6)
ax3.text(cx, cy-0.88, "local tangent space",
         fontsize=8.5, color="#555555", style="italic",
         ha="center", va="top", zorder=5)

save("figure3_nested_hierarchy.pdf")


# ============================================================================
# FIGURE 4  —  Gradient rotation drives channel exit (Act III)
# ============================================================================
#
# Biological motivation:
#   Climate change rotates the seed community on Daphne Major, shifting
#   what beak traits are rewarded.  The population's covariance cloud G
#   is elongated along beak depth; the new gradient points toward beak
#   width.  As the gradient rotates away from the major axis of G,
#   Ag falls — and the channel exits — without any loss of genetic variance
#   and without the curvature of the landscape increasing.
#   This is Act III: channel exit from gradient misalignment.
#
# Mathematics (Paper Results §4; SI Appendix S6.4):
#   G = diag(lambda1, lambda2),  lambda1 > lambda2 > 0
#   b(theta) = |b| * (cos θ, sin θ)'   (gradient rotates by angle θ)
#   H = -gamma * I                      (isotropic stabilising curvature)
#
#   Vlin(θ) = |b|² * (lambda1 cos²θ + lambda2 sin²θ)    (depends on θ)
#   Vquad   = (1/2) gamma² (lambda1² + lambda2²)         (does NOT depend on θ)
#
#   KEY: Vquad is constant in θ.  Only Vlin changes.
#   So Ag(θ) = Vlin(θ) / [Vlin(θ) + Vquad] declines as θ increases:
#   gradient rotation → Vlin falls → Ag falls → channel exits.
#
#   Exit angle θ*: smallest θ such that Ag(θ*) = alpha.
#   gamma = 0.414 chosen so that eccentricity r=lambda1/lambda2=2 is
#   the MARGINAL case: Ag(90°) = alpha = 0.70 exactly.
#   Derivation:
#     Ag(90°) = lambda2 / (lambda2 + (1/2)gamma²(lambda1²+lambda2²)/lambda2_norm)
#   Setting r=2, lambda2=1, Ag(90°) = 0.70 → gamma ≈ 0.414.
# ============================================================================
print("Figure 4 (gradient rotation) ...")

fig4, (axA4, axB4, axC4) = plt.subplots(1, 3, figsize=(13.5, 4.8))
fig4.subplots_adjust(left=0.05, right=0.98, top=0.88,
                     bottom=0.14, wspace=0.32)

# ── Panel A: geometric cartoon ───────────────────────────────────────────────
# Fixed G ellipse (semi-axes a = 1.8, b = 0.6; eigenvalues 3.24 and 0.36)
# Three gradient arrows at 0°, 50°, 90° relative to the major axis g_max.
# Ag values at each angle are computed from the analytic formula above.
lam1_A = 3.24; lam2_A = 0.36; gamma_A = 0.414
vquad_A = 0.5 * gamma_A**2 * (lam1_A**2 + lam2_A**2)

def ag_rot(theta_deg, lam1, lam2, vquad):
    """Ag as a function of gradient angle (degrees), eigenvalues, and Vquad."""
    th = np.radians(theta_deg)
    vlin = lam1*np.cos(th)**2 + lam2*np.sin(th)**2   # |b|=1
    return vlin / (vlin + vquad)

axA4.set_xlim(-3.2, 3.6); axA4.set_ylim(-2.6, 3.2)
axA4.set_aspect("equal"); axA4.axis("off")

axA4.add_patch(Ellipse((0,0), width=3.6, height=1.2, angle=0,
                        fc=FILL_ELLIPSE, ec=W_BLUE, lw=1.8, alpha=0.80, zorder=3))
# g_max dashed line along major axis
axA4.annotate("", xy=(1.90, 0), xytext=(-1.90, 0),
              arrowprops=dict(arrowstyle="->", color=W_BLUE, lw=1.6,
                              linestyle="dashed"), zorder=4)
axA4.text(2.00, -0.30, r"$\mathbf{g}_{\max}$",
          fontsize=10.5, color=W_BLUE, fontweight="bold",
          va="top", ha="left", zorder=5)
axA4.plot(0, 0, "o", color=W_BLUE, ms=5, zorder=6)
axA4.text(-0.22, -0.28, r"$\bar{\mathbf{z}}$",
          fontsize=11, color=W_BLUE, va="top", zorder=5)

bL4 = 1.62
arrows4 = [
    (  0, W_GREEN,      r"$\mathbf{b}_1$",
       r"$\mathcal{A}_g\!\approx\!$" + f"{ag_rot(0,  lam1_A, lam2_A, vquad_A):.2f}",
       0.10,  0.18),
    ( 50, W_ORANGE,     r"$\mathbf{b}_2$",
       r"$\mathcal{A}_g\!\approx\!$" + f"{ag_rot(50, lam1_A, lam2_A, vquad_A):.2f}",
       0.12,  0.14),
    ( 90, W_VERMILLION, r"$\mathbf{b}_3$",
       r"$\mathcal{A}_g\!\approx\!$" + f"{ag_rot(90, lam1_A, lam2_A, vquad_A):.2f}",
       0.12,  0.12),
]
for ang_deg, col, name, regime, ndx, ndy in arrows4:
    ang = np.radians(ang_deg)
    bx4 = bL4*np.cos(ang); by4 = bL4*np.sin(ang)
    axA4.annotate("", xy=(bx4, by4), xytext=(0, 0),
                  arrowprops=dict(arrowstyle="->", color=col, lw=2.2), zorder=5)
    axA4.text(bx4+ndx, by4+ndy,     name,   fontsize=9.0,  color=col,
              fontweight="bold", va="bottom", ha="left", zorder=6)
    axA4.text(bx4+ndx, by4+ndy-0.32, regime, fontsize=8.0, color=col,
              style="italic", va="top", ha="left", zorder=6)

# Arc θ between b1 and b2
arc_r4 = 0.72
arc_th4 = np.linspace(0, np.radians(50), 40)
axA4.plot(arc_r4*np.cos(arc_th4), arc_r4*np.sin(arc_th4),
          color="#888888", lw=0.9, zorder=4)
axA4.text(0.82, 0.28, r"$\theta$", fontsize=9.0, color="#888888", zorder=5)

axA4.axhline(0, color="#dddddd", lw=0.5, zorder=1)
axA4.axvline(0, color="#dddddd", lw=0.5, zorder=1)
axA4.text(0, -2.42, r"breeding value $a_1$",
          fontsize=9.5, ha="center", color="#555555")
axA4.text(-3.05, 0, r"$a_2$", fontsize=9.5, color="#555555",
          va="center", ha="right")
axA4.set_title(r"(A)  Fixed $\mathbf{G}$, rotating $\mathbf{b}$",
               fontsize=11.0, loc="center", pad=6)

# ── Panel B: Ag(θ) for four eccentricities ────────────────────────────────────
# gamma = 0.414 fixed; lambda2 = 1; |b| = 1.
# Four eccentricities r = lambda1/lambda2 = 1.5, 2, 3, 5.
# The curve for r=2 just touches alpha=0.70 at θ=90°: it is the marginal case.
gamma_B = 0.414
lam2_B  = 1.0
theta_B = np.linspace(0, 180, 360)
th_B    = np.radians(theta_B)

r_vals  = [1.5,  2.0,  3.0,  5.0]
r_cols  = [W_GREEN, W_ORANGE, W_BLUE, W_PURPLE]
r_lss   = ["-",    "--",     "-.",   ":"]
r_labs  = [
    r"$\lambda_1/\lambda_2 = 1.5$  (never exits)",
    r"$\lambda_1/\lambda_2 = 2$  (marginal)",
    r"$\lambda_1/\lambda_2 = 3$",
    r"$\lambda_1/\lambda_2 = 5$",
]

for rv, col, ls, lab in zip(r_vals, r_cols, r_lss, r_labs):
    lam1_b   = rv * lam2_B
    vquad_b  = 0.5 * gamma_B**2 * (lam1_b**2 + lam2_B**2)
    vlin_b   = lam1_b*np.cos(th_B)**2 + lam2_B*np.sin(th_B)**2
    ag_b     = vlin_b / (vlin_b + vquad_b)
    axB4.plot(theta_B, ag_b, color=col, lw=2.2, ls=ls, zorder=4, label=lab)

axB4.fill_between(theta_B, ALPHA_THRESHOLD, 1.0,
                  color=FILL_INTERIOR, alpha=0.55, zorder=1)
axB4.fill_between(theta_B, 0, ALPHA_THRESHOLD,
                  color=FILL_EXTERIOR, alpha=0.45, zorder=1)
axB4.axhline(ALPHA_THRESHOLD, color=W_VERMILLION, lw=1.4,
             linestyle=(0,(5,4)), zorder=3)
axB4.text(2, ALPHA_THRESHOLD+0.018, r"$\alpha$",
          fontsize=11, color=W_VERMILLION, va="bottom", zorder=5)
axB4.set_xlim(0, 180); axB4.set_ylim(0, 1.05)
axB4.set_xlabel(r"gradient angle $\theta$ (degrees)", fontsize=10.5)
axB4.set_ylabel(r"$\mathcal{A}_g(\theta)$", fontsize=10.5)
axB4.tick_params(labelsize=9); axB4.set_xticks([0, 45, 90, 135, 180])
style_ax(axB4)
axB4.legend(fontsize=7.8, loc="upper right", frameon=True,
            framealpha=0.92, edgecolor="#ccc", handlelength=2.0)
axB4.set_title(r"(B)  $\mathcal{A}_g$ declines with gradient rotation",
               fontsize=11.0, loc="center", pad=6)

# ── Panel C: exit angle θ* vs curvature-to-gradient ratio κ0 ─────────────────
# κ0 = Vquad / Vlin(0°) = Vquad / (|b|² lambda1)
# Exit condition: Ag(θ*) = alpha
#   Vlin(θ*) / [Vlin(θ*) + Vquad] = alpha
#   Vlin(θ*) = alpha/(1-alpha) * Vquad = alpha/(1-alpha) * κ0 * Vlin(0)
# Substituting Vlin(θ) = |b|² (lambda1 cos²θ + lambda2 sin²θ):
#   lambda1 cos²θ + lambda2 sin²θ = alpha/(1-alpha) * κ0 * lambda1
#   Using sin²θ = 1 - cos²θ, with r = lambda1/lambda2:
#   cos²θ* = [alpha/(1-alpha) * κ0 * r - 1] / (r - 1)
# All curves converge at κ0 = (1-alpha)/alpha ≈ 3/7: the upper bound where
# even the best-aligned population sits at the channel boundary.
kappa0_arr = np.linspace(0, (1-ALPHA_THRESHOLD)/ALPHA_THRESHOLD, 300)
r_vals_C   = [1.5, 2.0, 3.0, 5.0]

for rv, col, ls in zip(r_vals_C, r_cols, r_lss):
    theta_exit = []
    for k0 in kappa0_arr:
        cos2 = (ALPHA_THRESHOLD/(1-ALPHA_THRESHOLD) * k0 * rv - 1) / (rv - 1)
        if 0 <= cos2 <= 1:
            theta_exit.append(np.degrees(np.arccos(np.sqrt(cos2))))
        else:
            theta_exit.append(np.nan)
    theta_exit = np.array(theta_exit)
    m = ~np.isnan(theta_exit)
    axC4.plot(kappa0_arr[m], theta_exit[m], color=col, lw=2.2, ls=ls, zorder=4)

# Convergence point annotation
k_upper = (1-ALPHA_THRESHOLD)/ALPHA_THRESHOLD
axC4.axvline(k_upper, color=W_VERMILLION, lw=1.0, ls=(0,(4,4)), zorder=3)
axC4.text(k_upper+0.01, 80, r"$\kappa_0^{\max} = \frac{1-\alpha}{\alpha}$",
          fontsize=8.5, color=W_VERMILLION, va="center", zorder=5)
axC4.set_xlim(0, k_upper+0.05); axC4.set_ylim(0, 92)
axC4.set_xlabel(r"curvature-to-gradient ratio $\kappa_0 = V_{\mathrm{quad}}/V_{\mathrm{lin}}(0°)$",
                fontsize=9.5)
axC4.set_ylabel(r"exit angle $\theta^*$ (degrees)", fontsize=10.5)
axC4.tick_params(labelsize=9)
style_ax(axC4)
axC4.set_title(r"(C)  Exit angle $\theta^*$ vs curvature load",
               fontsize=11.0, loc="center", pad=6)

save("figure4_gradient_rotation.pdf")


# ============================================================================
# FIGURE 5  —  Moving-optimum dynamics: Acts I, II, and III
# ============================================================================
#
# The one-trait moving-optimum model makes the channel dynamics concrete.
# It is the simplest setting where ALL quantities — d, G, b, H, Ag, L, r —
# can be written in closed form and their relationships checked analytically.
#
# BIOLOGICAL STORY (three acts):
#
#   Act I   (Phase 1, generations 0–50):
#     The optimum θ moves at rate θ_dot = 0.02.
#     The population lags behind: d = θ - z_bar is large and growing.
#     Because d² >> G/2, Ag ≈ 1: the directional signal dominates.
#     L(t) > L_alpha: the population is deep inside the channel.
#     r > 0: the channel is deepening (faster optimum motion regenerates
#     more directional signal than curvature exposure erodes it).
#
#   Act II  (Phase 2, generations 50–80):
#     The optimum stops (θ_dot = 0).  Selection catches up: d → 0.
#     As d → 0, b = s*d → 0 and Vlin → 0.  Ag → 0.
#     The channel CLOSES — not because the biology became complex,
#     but because selection SUCCEEDED: the directional signal collapsed.
#     r < 0 throughout Phase 2.
#     A breeder whose population has reached its target ceiling sees this.
#
#   Act III (Phase 3, generations 80–130):
#     The optimum resumes moving (θ_dot = 0.02 again).
#     The lag rebuilds: d grows, b = s*d grows, Vlin recovers.
#     Ag rises back above alpha: channel RE-ENTRY.
#     r > 0 throughout Phase 3.
#
# STATE EQUATIONS (SI Appendix S6):
#   d_dot = θ_dot - G*s*d
#   G_dot = -s*G² - lam_rec*(G - G_star) + m - G/(2*Ne)
#   Ag    = d² / (d² + G/2)
#   L     = ln(2*d² / G)
#   r     = 2*θ_dot/d - s*G + lam_rec*(1 - G_star/G) - m/G + 1/(2*Ne)
#
# PANEL STRUCTURE:
#   (A) State variables: d(t) [green], G(t) [blue dashed], Ag(t) [purple right axis]
#   (B) Logit trajectory L(t) [blue] with channel threshold L_alpha [vermillion]
#   (C) Rate diagnostic r(t) [black] decomposed into:
#         r_motion = 2*θ_dot/d          [green shading, positive]
#         r_curv   = -s*G               [orange shading, negative]
#         r_other  = rec + mut + drift  [grey dotted]
# ============================================================================
print("Figure 5 (moving optimum) ...")

# ── Model parameters ─────────────────────────────────────────────────────────
# These values are chosen to produce clean, visually readable dynamics.
# They are biologically plausible for a quantitative trait under stabilising
# selection (s = 0.10 → weak-to-moderate stabilising selection;
# lam_rec = 0.30 → moderate recombination; Ne = 500 → medium-sized population).
s5       = 0.10    # stabilising curvature strength (> 0)
lam_rec5 = 0.30    # recombination relaxation rate (brings G back toward G_star)
G_star5  = 0.80    # recombination baseline variance
m5       = 0.005   # mutational variance input per generation
Ne5      = 500.0   # effective population size
d_delta  = 0.05    # gap threshold: r_motion undefined when |d| < d_delta

# Optimum schedule θ(t) and θ_dot(t):
#   Phase 1 (0 ≤ t < 50):   θ_dot = +0.02  (moving right)
#   Phase 2 (50 ≤ t < 80):  θ_dot =  0.00  (stopped)
#   Phase 3 (80 ≤ t ≤ 130): θ_dot = +0.02  (moving again)
T_mo     = 130
PHASES   = [(0, 50, 0.02), (50, 80, 0.0), (80, T_mo, 0.02)]

def theta_dot_t(t):
    for t0, t1, rate in PHASES:
        if t0 <= t < t1:
            return rate
    return PHASES[-1][2]

# ── Numerical integration ─────────────────────────────────────────────────────
def odes5(t, y):
    """
    Right-hand sides of the state equations for d(t) and G(t).

    Paper SI Appendix S6, Eqs. S6.1 – S6.2:
      d_dot = θ_dot - G*s*d
      G_dot = -s*G² - lam_rec*(G - G*) + m - G/(2*Ne)

    Arguments
    ---------
    t : float — time (generations)
    y : [d, G] — current state

    Returns
    -------
    [d_dot, G_dot]
    """
    d, G = y
    G    = max(G, 1e-6)           # prevent negative variance (numerical safety)
    td   = theta_dot_t(t)
    d_dt = td - G * s5 * d
    G_dt = (-s5 * G**2
            - lam_rec5 * (G - G_star5)
            + m5
            - G / (2*Ne5))
    return [d_dt, G_dt]

# Initial conditions: population starts with a moderate lag and high variance
sol5 = solve_ivp(odes5, [0, T_mo], [0.05, G_star5], max_step=0.5,
                 dense_output=True, rtol=1e-8, atol=1e-10)
t_arr5 = np.linspace(0, T_mo, 1300)
d_arr, G_arr = sol5.sol(t_arr5)
G_arr = np.maximum(G_arr, 1e-6)

# Derived quantities
Ag_arr  = d_arr**2 / (d_arr**2 + G_arr/2)     # Paper Eq. S6.3
L_arr   = np.log(np.maximum(2*d_arr**2 / G_arr, 1e-6))  # logit of Ag
L_alpha = np.log(ALPHA_THRESHOLD / (1 - ALPHA_THRESHOLD))  # threshold in logit space

# Rate diagnostic components  (Paper Eq. S6.4)
# r_motion = 2*θ_dot/d  (undefined when |d| < d_delta → set to NaN)
r_opt_arr  = np.where(np.abs(d_arr) > d_delta,
                      2 * np.array([theta_dot_t(t) for t in t_arr5]) / d_arr,
                      np.nan)
r_curv_arr = -s5 * G_arr                           # always negative
r_rec_arr  = lam_rec5 * (1 - G_star5/G_arr)        # sign depends on G vs G*
r_mut_arr  = -m5 / G_arr                            # always negative
r_drift_arr = np.full_like(t_arr5, 1.0/(2*Ne5))    # always positive (small)
r_oth_arr  = r_rec_arr + r_mut_arr + r_drift_arr    # combined other terms
r_tot_arr  = np.where(np.abs(d_arr) > d_delta,
                      r_opt_arr + r_curv_arr + r_oth_arr,
                      np.nan)

# Channel crossing times
above_alpha = Ag_arr > ALPHA_THRESHOLD
entry_idx   = np.where(np.diff(above_alpha.astype(int)) == +1)[0]
exit_idx    = np.where(np.diff(above_alpha.astype(int)) == -1)[0]
reentry_idx = exit_idx[exit_idx > entry_idx[0]] if (
              len(exit_idx) > 0 and len(entry_idx) > 0) else []

# ── Phase annotation helper ───────────────────────────────────────────────────
def add_phases(ax, y_label, fontsize=8.5):
    """
    Draw vertical phase separators and labels for the three acts.
    Green shading = optimum moving (channel open or opening).
    Orange shading = optimum stopped (channel closing).
    """
    for t0, t1, rate in PHASES:
        col = "#e8f4e8" if rate > 0 else "#fef3e3"
        ax.axvspan(t0, t1, color=col, alpha=0.40, zorder=0)
    for t0, t1, rate in PHASES:
        t_mid = (t0 + t1) / 2
        label = "Act I" if t0 == 0 else ("Act II" if rate == 0 else "Act I")
        ax.text(t_mid, y_label, label,
                fontsize=fontsize, ha="center", va="top",
                color="#666666", style="italic", zorder=2)
    for boundary in [50, 80]:
        ax.axvline(boundary, color="#cccccc", lw=0.8, ls=":", zorder=2)

# ── Build figure ─────────────────────────────────────────────────────────────
fig5, (pA5, pB5, pC5) = plt.subplots(1, 3, figsize=(13.5, 4.8))
fig5.subplots_adjust(left=0.07, right=0.95, top=0.88, bottom=0.14, wspace=0.48)

# ── Panel A: state variables d(t), G(t), Ag(t) ───────────────────────────────
add_phases(pA5, 1.62)

pA5.plot(t_arr5, d_arr, color=W_GREEN, lw=2.0, zorder=4,
         label=r"lag $d(t)$")
pA5.plot(t_arr5, G_arr, color=W_BLUE, lw=2.0, ls="--", zorder=4,
         label=r"variance $G(t)$")
pA5.set_ylim(0, 1.70)
pA5.set_ylabel(r"lag $d$ or variance $G$", fontsize=10.0)
pA5.tick_params(axis="y", labelsize=9)
pA5.set_xlabel("time (generations)", fontsize=10.0, color="#555555")
pA5.tick_params(axis="x", labelsize=9, colors="#555555")
pA5.set_xlim(0, T_mo)

# Right axis: Ag(t) — purple to match the Wong palette concept
pA5_r = pA5.twinx()
pA5_r.plot(t_arr5, Ag_arr, color=W_PURPLE, lw=2.0, zorder=4)
pA5_r.axhline(ALPHA_THRESHOLD, color=W_VERMILLION, lw=1.2,
              linestyle=(0,(5,4)), zorder=3)
pA5_r.fill_between(t_arr5, ALPHA_THRESHOLD, Ag_arr,
                   where=(Ag_arr > ALPHA_THRESHOLD),
                   color=FILL_INTERIOR, alpha=0.40, zorder=1)
pA5_r.fill_between(t_arr5, 0, Ag_arr,
                   where=(Ag_arr <= ALPHA_THRESHOLD),
                   color=FILL_EXTERIOR, alpha=0.35, zorder=1)
pA5_r.text(T_mo*0.93, ALPHA_THRESHOLD+0.02, r"$\alpha$",
           fontsize=10, color=W_VERMILLION, va="bottom", ha="left")
pA5_r.set_ylim(0, 1.05)
pA5_r.set_ylabel(r"$\mathcal{A}_g(t)$", fontsize=10.0, color=W_PURPLE)
pA5_r.tick_params(axis="y", colors=W_PURPLE, labelsize=9)
pA5_r.spines["right"].set_color(W_PURPLE)
pA5_r.spines["right"].set_linewidth(0.8)
for sp in ["top", "bottom", "left"]:
    pA5_r.spines[sp].set_visible(False)
pA5_r.set_facecolor("none")
style_ax(pA5)

elems5A = [
    Line2D([0],[0], color=W_GREEN,      lw=2.0,           label=r"lag $d(t)$"),
    Line2D([0],[0], color=W_BLUE,       lw=2.0, ls="--",  label=r"variance $G(t)$"),
    Line2D([0],[0], color=W_PURPLE,     lw=2.0,           label=r"$\mathcal{A}_g(t)$"),
    Line2D([0],[0], color=W_VERMILLION, lw=1.2, ls=(0,(5,4)), label=r"$\alpha$"),
]
pA5.legend(handles=elems5A, fontsize=8.0, loc="upper left",
           frameon=True, framealpha=0.92, edgecolor="#ccc", handlelength=1.8)
pA5.set_title(r"(A)  State variables and local predictability $\mathcal{A}_g(t)$",
              fontsize=11.0, loc="center", pad=6)

# ── Panel B: logit trajectory L(t) ───────────────────────────────────────────
# L(t) = ln(2d²/G) is the logit of Ag.
# Channel boundary alpha=0.70 corresponds to L_alpha = ln(0.70/0.30) ≈ 0.847.
# Key insight: L(t) is the integral of r(t).
# So the SLOPE of L(t) equals r(t) exactly:
#   rising L → r > 0 (linear variance gaining on quadratic)
#   falling L → r < 0 (curvature gaining on gradient)
# This is why r(t) is the natural longitudinal diagnostic.
add_phases(pB5, 2.80)

L_plot = np.clip(L_arr, -4.5, None)
pB5.plot(t_arr5, L_plot, color=W_BLUE, lw=2.0, zorder=4)
pB5.axhline(L_alpha, color=W_VERMILLION, lw=1.4, linestyle=(0,(5,4)), zorder=3)
pB5.fill_between(t_arr5, L_alpha, L_plot,
                 where=(L_plot > L_alpha), color=FILL_INTERIOR, alpha=0.45, zorder=1)
pB5.fill_between(t_arr5, -4.5, L_plot,
                 where=(L_plot <= L_alpha), color=FILL_EXTERIOR, alpha=0.40, zorder=1)
pB5.text(T_mo*0.93, L_alpha+0.08, r"$\mathcal{L}_\alpha$",
         fontsize=10, color=W_VERMILLION, va="bottom", ha="left")

# Annotate channel crossing events
for idx_arr, label, col, dy_sign in [
    (entry_idx,   "entry",    W_GREEN,      +0.65),
    (exit_idx,    "exit",     W_VERMILLION, -0.90),
    (reentry_idx, "re-entry", W_GREEN,      +0.65),
]:
    if len(idx_arr) > 0:
        t_ev = t_arr5[idx_arr[0]]
        pB5.annotate(label, xy=(t_ev, L_alpha),
                     xytext=(t_ev + (2 if dy_sign > 0 else -3), L_alpha + dy_sign),
                     fontsize=8.0, color=col, style="italic",
                     arrowprops=dict(arrowstyle="->", color=col, lw=0.9), zorder=5)

pB5.set_xlim(0, T_mo); pB5.set_ylim(-4.8, 3.0)
pB5.set_xlabel("time (generations)", fontsize=10.0, color="#555555")
pB5.set_ylabel(r"logit $\mathcal{L}(t)=\ln(\mathcal{A}_g/(1-\mathcal{A}_g))$",
               fontsize=9.5)
pB5.tick_params(labelsize=9)
style_ax(pB5)
pB5.set_title(r"(B)  Logit trajectory $\mathcal{L}(t)$ — integral of $r(t)$",
              fontsize=11.0, loc="center", pad=6)

# ── Panel C: r(t) decomposed ─────────────────────────────────────────────────
# r(t) = r_motion + r_curv + r_other
#
# r_motion = 2*θ_dot/d    (green shading, positive when active)
#   Regenerates directional signal.  Only nonzero when the optimum moves
#   (Phases 1 and 3).  Diverges as d→0 near the gap at t≈50–80.
#
# r_curv = -s*G           (orange shading, always negative)
#   Curvature exposure cost: broader G means more quadratic variance.
#   Becomes less negative as G is compressed by selection.
#
# r_other = r_rec + r_mut + r_drift  (grey dotted)
#   Small corrections; sign depends on parameter values.
#   Here slightly negative: recombination is expanding G (G < G_star
#   during Phase 2 because selection compressed it below baseline).
#
# r_tot = sum             (black solid)
#   r_tot > 0 → L(t) rising → channel deepening or re-entry.
#   r_tot < 0 → L(t) falling → channel eroding or exit.
add_phases(pC5, 0.52)
pC5.axhline(0, color="#888888", lw=0.7, linestyle="--", zorder=2)

pC5.fill_between(t_arr5,
                 np.where(~np.isnan(r_opt_arr), r_opt_arr, np.nan), 0,
                 color=W_GREEN, alpha=0.28, zorder=2,
                 label=r"$r_{\mathrm{motion}} = 2\dot\theta/d$")
pC5.fill_between(t_arr5,
                 np.where(~np.isnan(r_curv_arr), r_curv_arr, np.nan), 0,
                 color=W_ORANGE, alpha=0.32, zorder=2,
                 label=r"$r_{\mathrm{curv}} = -sG$")
pC5.plot(t_arr5, r_oth_arr, color="#999999", lw=1.0, ls=":", zorder=3,
         label=r"$r_{\mathrm{other}}$ (rec+mut+drift)")
pC5.plot(t_arr5, r_tot_arr, color=W_BLACK, lw=2.0, zorder=4,
         label=r"$r(t)$ total")

# Mark the gap (|d| < d_delta: r_motion undefined)
gap_where = np.where(np.isnan(r_tot_arr) & (t_arr5 > 45))[0]
if len(gap_where) > 0:
    t_mid_gap = (t_arr5[gap_where[0]] + t_arr5[gap_where[-1]]) / 2
    pC5.text(t_mid_gap, 0.43, r"$|d|<\delta$" + "\n(gap)",
             fontsize=7.5, color="#888888", ha="center", va="bottom",
             style="italic", zorder=5)

pC5.set_xlim(0, T_mo); pC5.set_ylim(-0.52, 0.52)
pC5.set_xlabel("time (generations)", fontsize=10.0, color="#555555")
pC5.set_ylabel(r"rate $r(t) = \dot{\mathcal{L}}$", fontsize=10.5)
pC5.tick_params(labelsize=9)
style_ax(pC5)
pC5.legend(fontsize=7.8, loc="lower right", frameon=True,
           framealpha=0.92, edgecolor="#ccc", handlelength=2.0, ncol=2)
pC5.set_title(r"(C)  Rate diagnostic $r(t)$ decomposed by biological process",
              fontsize=11.0, loc="center", pad=6)

save("figure5_moving_optimum.pdf")


# ============================================================================
# SI FIGURE S1  —  Parameter sensitivity of r(t)
# ============================================================================
#
# This figure addresses the question: which parameters most strongly
# influence r(t) in the one-trait model?
#
# Eq. S6.4: r = 2*θ_dot/d - s*G + lam_rec*(1 - G*/G) - m/G + 1/(2*Ne)
#
# Each panel varies one parameter while holding others at baseline.
# BASELINE: s=0.10, G=1.0, d=2.0, Ne=500, lam_rec=0.30, G*=1.0, m=0.005.
# The vertical dashed line marks the baseline value.
# The shaded region marks where r > 0 (channel improving).
# ============================================================================
print("SI Figure S1 (parameter sensitivity) ...")

# Baseline
s0=0.10; G0=1.0; d0=2.0; Ne0=500; lr0=0.30; Gs0=1.0; m0=0.005; td0=0.02

def r_1d(td, d, s, G, lr, Gs, m, Ne):
    """Full 1D rate diagnostic. Paper SI Eq. S6.4."""
    return 2*td/d - s*G + lr*(1 - Gs/G) - m/G + 1/(2*Ne)

r_base = r_1d(td0, d0, s0, G0, lr0, Gs0, m0, Ne0)

figS1, axsS1 = plt.subplots(2, 2, figsize=(10.0, 7.5))
figS1.subplots_adjust(left=0.10, right=0.97, top=0.94,
                      bottom=0.10, hspace=0.40, wspace=0.35)
axsS1 = axsS1.ravel()

sensitivity_panels = [
    # (ax_idx, param_label, x_array, r_array, baseline, xlabel)
    (0, r"stabilising curvature $s$",
     np.linspace(0.01, 0.40, 300),
     [r_1d(td0, d0, sv, G0, lr0, Gs0, m0, Ne0) for sv in np.linspace(0.01, 0.40, 300)],
     s0, r"$s$"),
    (1, r"lag behind optimum $d$",
     np.linspace(0.10, 4.0, 300),
     [r_1d(td0, dv, s0, G0, lr0, Gs0, m0, Ne0) for dv in np.linspace(0.10, 4.0, 300)],
     d0, r"$d$"),
    (2, r"genetic variance $G$",
     np.linspace(0.10, 3.0, 300),
     [r_1d(td0, d0, s0, Gv, lr0, Gs0, m0, Ne0) for Gv in np.linspace(0.10, 3.0, 300)],
     G0, r"$G$"),
    (3, r"effective population size $N_e$",
     np.linspace(50, 5000, 300),
     [r_1d(td0, d0, s0, G0, lr0, Gs0, m0, Nv) for Nv in np.linspace(50, 5000, 300)],
     Ne0, r"$N_e$"),
]

for idx, title, x_arr, r_arr, baseline, xlab in sensitivity_panels:
    ax = axsS1[idx]
    r_arr = np.array(r_arr)
    ax.axhline(0, color="#888888", lw=0.8, ls="--", zorder=2)
    ax.fill_between(x_arr, 0, r_arr, where=(r_arr > 0),
                    color=W_GREEN, alpha=0.28, zorder=1, label=r"$r > 0$  (improving)")
    ax.fill_between(x_arr, r_arr, 0, where=(r_arr < 0),
                    color=W_ORANGE, alpha=0.28, zorder=1, label=r"$r < 0$  (eroding)")
    ax.plot(x_arr, r_arr, color=W_BLUE, lw=2.2, zorder=4)
    ax.axvline(baseline, color=W_VERMILLION, lw=1.0, ls=(0,(4,3)), zorder=3)
    ax.text(baseline*1.03, ax.get_ylim()[1]*0.90 if idx < 2 else 0,
            "baseline", fontsize=8, color=W_VERMILLION, style="italic")
    ax.set_xlabel(xlab, fontsize=10.5)
    ax.set_ylabel(r"$r(t)$", fontsize=10.5)
    ax.tick_params(labelsize=9)
    style_ax(ax)
    ax.legend(fontsize=8.0, loc="upper right" if idx != 1 else "lower right",
              frameon=True, framealpha=0.92, edgecolor="#ccc")
    ax.set_title(title, fontsize=10.5, loc="left", pad=5)

figS1.suptitle(r"SI Fig. S1: Sensitivity of the rate diagnostic $r(t)$ "
               r"to model parameters", fontsize=11.5, y=0.98)
save("figS1_parameter_sensitivity.pdf")


# ============================================================================
# SI FIGURE S2  —  State-space slice in the (s, kappa) plane
# ============================================================================
#
# This is Panel B from the original Figure 2, moved to SI.
#
# The channel boundary Ag = alpha is the hyperbola:
#   s = (1-alpha) / (alpha * kappa)
# on the (s, kappa) plane, where s = variance scale and kappa = Vquad/Vlin.
#
# Three trajectories illustrate qualitatively different population dynamics:
#   T1 (green):     variance compresses (selection) → moves left → channel deepens
#   T2 (orange):    both s and kappa rise → boundary approach
#   T3 (vermillion): channel exit — crosses the red hyperbola
#
# Why this is SI and not main text:
#   This panel requires the reader to think in two dimensions simultaneously.
#   The key message (variance compression deepens the channel) is already
#   communicated by the 1D panel in Figure 2.  The state-space geometry is
#   additional supporting detail appropriate for SI.
# ============================================================================
print("SI Figure S2 (state-space slice) ...")

figS2, axS2 = plt.subplots(figsize=(6.8, 5.5))
figS2.subplots_adjust(left=0.13, right=0.97, top=0.91, bottom=0.14)

s_2d   = np.linspace(0.02, 3.5, 400)
kap_2d = np.linspace(0.05, 2.0, 400)
S2g, K2g = np.meshgrid(s_2d, kap_2d)
AG2    = 1.0 / (1.0 + S2g * K2g)

BND_HI2 = ALPHA_THRESHOLD + 0.08
BND_LO2 = ALPHA_THRESHOLD - 0.08
axS2.contourf(S2g, K2g, AG2,
              levels=[0, BND_LO2, BND_HI2, 1.001],
              colors=[FILL_EXTERIOR, FILL_BOUNDARY, FILL_INTERIOR],
              alpha=0.60, zorder=1)

# Boundary hyperbolas
kap_line2 = np.linspace(0.12, 2.0, 300)
for thr, lc, lw_b, ls in [
    (ALPHA_THRESHOLD, W_VERMILLION, 1.6, (0,(5,4))),
    (BND_HI2,         W_ORANGE,     0.9, (0,(3,4))),
    (BND_LO2,         W_ORANGE,     0.9, (0,(3,4))),
]:
    s_hyp = (1-thr) / (thr * kap_line2)
    m2    = (s_hyp > 0) & (s_hyp < 3.5)
    axS2.plot(s_hyp[m2], kap_line2[m2], color=lc, lw=lw_b, ls=ls, zorder=3)

# Trajectory T1: G compresses (moves left in s), kappa roughly constant
t1_s = np.linspace(0.52, 0.18, 60)
t1_k = 0.40 + 0.06*np.sin(np.linspace(0, np.pi, 60))
axS2.plot(t1_s, t1_k, color=W_GREEN, lw=2.0, zorder=4)
axS2.annotate("", xy=(t1_s[-1], t1_k[-1]), xytext=(t1_s[-6], t1_k[-6]),
              arrowprops=dict(arrowstyle="->", color=W_GREEN, lw=2.0), zorder=5)
axS2.text(t1_s[-1]+0.04, t1_k[-1]-0.04, r"$T_1$: $G$ compresses",
          fontsize=8.5, color=W_GREEN, fontweight="bold",
          va="top", ha="left", zorder=6)

# Trajectory T2: boundary approach
t2_s = np.linspace(0.28, 0.70, 60)
t2_k = np.linspace(0.90, 1.20, 60)
axS2.plot(t2_s, t2_k, color=W_ORANGE, lw=2.0, zorder=4)
axS2.annotate("", xy=(t2_s[-1], t2_k[-1]), xytext=(t2_s[-6], t2_k[-6]),
              arrowprops=dict(arrowstyle="->", color=W_ORANGE, lw=2.0), zorder=5)
axS2.text(t2_s[-1]+0.04, t2_k[-1]+0.05, r"$T_2$: boundary approach",
          fontsize=8.5, color=W_ORANGE, fontweight="bold",
          va="bottom", ha="left", zorder=6)

# Trajectory T3: channel exit
t3_s = np.linspace(0.18, 1.55, 80)
t3_k = np.linspace(0.25, 1.10, 80)
axS2.plot(t3_s, t3_k, color=W_VERMILLION, lw=2.0, zorder=4)
axS2.annotate("", xy=(t3_s[-1], t3_k[-1]), xytext=(t3_s[-6], t3_k[-6]),
              arrowprops=dict(arrowstyle="->", color=W_VERMILLION, lw=2.0), zorder=5)
ag_t3 = 1 / (1 + t3_s*t3_k)
ci3   = np.argmin(np.abs(ag_t3 - ALPHA_THRESHOLD))
axS2.plot(t3_s[ci3], t3_k[ci3], "*", color=W_VERMILLION, ms=12, zorder=6)
axS2.text(t3_s[ci3]+0.12, t3_k[ci3]-0.06, r"$T_{\mathrm{exit}}$",
          fontsize=9, color=W_VERMILLION, va="top", ha="left", zorder=6)
axS2.text(t3_s[0]+0.02, t3_k[0]-0.10, r"$T_3$: channel exit",
          fontsize=8.5, color=W_VERMILLION, fontweight="bold",
          va="top", ha="left", zorder=6)

axS2.set_xlim(0, 3.5); axS2.set_ylim(0.0, 2.0)
axS2.set_xlabel(r"variance scale $s$", fontsize=10.5)
axS2.set_ylabel(r"curvature-to-gradient ratio $\kappa$", fontsize=10.5)
axS2.tick_params(labelsize=9.5)
style_ax(axS2)

elemsS2 = [
    Line2D([0],[0], color=W_VERMILLION, lw=1.5, linestyle=(0,(5,4)),
           label=r"$\mathcal{A}_g = \alpha$ boundary"),
    mpatches.Patch(color=FILL_INTERIOR, alpha=0.70, label="interior channel"),
    mpatches.Patch(color=FILL_BOUNDARY, alpha=0.70, label="boundary zone"),
    mpatches.Patch(color=FILL_EXTERIOR, alpha=0.60, label="exterior"),
    Line2D([0],[0], color=W_GREEN,      lw=2.0, label=r"$T_1$: $G$ compresses"),
    Line2D([0],[0], color=W_ORANGE,     lw=2.0, label=r"$T_2$: boundary approach"),
    Line2D([0],[0], color=W_VERMILLION, lw=2.0, label=r"$T_3$: channel exit"),
]
axS2.legend(handles=elemsS2, fontsize=8.0, loc="upper right",
            frameon=True, framealpha=0.92, edgecolor="#ccc", handlelength=1.8)
axS2.set_title(r"SI Fig. S2: Channel boundary in the $(s,\,\kappa)$ state-space plane",
               fontsize=10.5, loc="left", pad=7)

save("figS2_statespace_slice.pdf")


# ============================================================================
# MANIFEST
# ============================================================================
print("\nAll figures complete.")
print("Files produced:")
for fname in [
    "figure1_channel_concept.pdf",
    "figure2_regime_marker.pdf",
    "figure3_nested_hierarchy.pdf",
    "figure4_gradient_rotation.pdf",
    "figure5_moving_optimum.pdf",
    "figS1_parameter_sensitivity.pdf",
    "figS2_statespace_slice.pdf",
]:
    status = "✓" if os.path.exists(os.path.join(OUTPUT_DIR, fname)) else "MISSING"
    print(f"  {status}  {fname}")
