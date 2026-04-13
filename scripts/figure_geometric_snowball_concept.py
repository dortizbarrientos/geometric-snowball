#!/usr/bin/env python3
"""
figure_geometric_snowball_concept.py
=============================================================================
"The Geometric Snowball of Speciation" — Ortiz-Barrientos, James, Burrage

PURPOSE
-------
Generates Figure 1 of the paper: the conceptual overview of the geometric
snowball. The figure has ONE main panel — the snowball growth curve — with
THREE small inset diagrams placed on the curve at early, intermediate, and
late stages of divergence.

WHAT THE FIGURE IS SHOWING (for the biologist)
-----------------------------------------------
When two populations adapt to different environments, each population evolves
a preferred "direction" of additive genetic response — its GENETIC CHANNEL.
This channel is the direction in multivariate trait space where genetic
variance and selective pressure reinforce each other most strongly.

As the two populations diverge, their channels rotate apart. When you cross
them, the hybrid inherits a mix of additive variance from both channels.
Some of that variance points toward one parental optimum; some toward the
other. None of it is perfectly aligned with either environment. The result
is reduced hybrid fitness.

The key quantity is I_geom — the geometric overlap deficit — which measures
how much channel overlap has been lost. When channels point in the same
direction (theta = 0): I_geom = 0, no cost, hybrids can be vigorous.
As channels rotate apart: I_geom grows as theta^2 for small angles,
then accelerates — the GEOMETRIC SNOWBALL.

THE THREE INSETS
----------------
Each inset is a "V-tree" diagram: ancestor at top, two diverged populations
in the middle (shown as coloured ellipses representing their G-matrices),
and the hybrid at the bottom (ellipse coloured by the fitness cost I_geom).

The insets sit directly ON the snowball curve, so their y-position reads
off the incompatibility at that stage. The reader can see simultaneously:
  - The shape of the hybrid G-matrix (how misaligned the additive variance is)
  - The value of I_geom (where on the curve that stage falls)
  - The channel arrows (the direction that has rotated)

MATHEMATICAL ANNOTATIONS (for the programmer)
----------------------------------------------
The G-matrix is shown as an ellipse. The semi-axes are:
  - Major axis length = g_dom (dominant eigenvalue sqrt ~ variance along g_max)
  - Minor axis length = g_min (minor eigenvalue sqrt ~ variance orthogonal to g_max)

The channel direction u is a unit vector in the direction of g_max, rotated
by theta/2 symmetrically about the vertical for the two populations.

The hybrid G-matrix is V_H = 0.5*(G1 + G2). Its shape is computed from the
eigendecomposition of V_H. Its colour interpolates from grey (I_geom = 0,
no cost) to vermillion (I_geom = 1, complete isolation).

OUTPUT
------
Saves figure_geometric_snowball_concept.pdf to the output directory, resolved
in this order:
  1. First command-line argument
  2. Environment variable SNOWBALL_OUTDIR
  3. ./figures/ (created if needed)

USAGE
-----
python figure_geometric_snowball_concept.py              # save to ./figures/
python figure_geometric_snowball_concept.py ~/Desktop    # save to ~/Desktop
SNOWBALL_OUTDIR=~/myfigures python figure_geometric_snowball_concept.py

DEPENDENCIES
------------
Python >= 3.8, numpy, matplotlib
=============================================================================

AUTHORS
-------
D. Ortiz-Barrientos, M.E. James, K. Burrage (University of Queensland / QUT)
Code annotations: Claude (Anthropic), April 2026
Repository: github.com/dortizbarrientos/geometric-snowball
"""

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse

# ── Non-interactive backend: saves PDFs without opening a window ──────────────
# Remove this line (or set matplotlib.use("TkAgg")) to see the figure on screen.
matplotlib.use("Agg")


# =============================================================================
# SECTION 0: OUTPUT DIRECTORY
# =============================================================================
# Portable resolution — works on any machine, no hardcoded paths.

def _resolve_outdir() -> Path:
    """Return the output directory, creating it if necessary."""
    if len(sys.argv) > 1:
        d = Path(sys.argv[1]).expanduser()
    elif "SNOWBALL_OUTDIR" in os.environ:
        d = Path(os.environ["SNOWBALL_OUTDIR"]).expanduser()
    else:
        d = Path.cwd() / "figures"
    d.mkdir(parents=True, exist_ok=True)
    return d

OUTDIR = _resolve_outdir()
print(f"Figures will be saved to: {OUTDIR}")


# =============================================================================
# SECTION 1: PLOT STYLE
# =============================================================================
# Nature Ecology & Evolution: Arial/Helvetica, 7pt body, 6pt labels, PDF type 42.

plt.rcParams.update({
    "font.family"      : "sans-serif",
    "font.sans-serif"  : ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset" : "dejavusans",
    "pdf.fonttype"     : 42,        # editable text in Illustrator/Inkscape
    "axes.linewidth"   : 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size" : 3,
    "ytick.major.size" : 3,
    "font.size"        : 7,
    "axes.labelsize"   : 7,
    "xtick.labelsize"  : 6,
    "ytick.labelsize"  : 6,
})

# Okabe-Ito colour-blind-safe palette
C_BLUE    = "#0072B2"   # Population 1 (Dune-like)
C_SKY     = "#56B4E9"   # Population 1 fill
C_ORANGE  = "#E69F00"   # Population 2 (Headland-like)
C_ORANGE2 = "#C06000"   # Population 2 edge
C_VERM    = "#D55E00"   # High-incompatibility colour (deep vermillion)


# =============================================================================
# SECTION 2: V-TREE INSET FUNCTION
# =============================================================================

def draw_v_snapshot(ax, theta_deg: float, stage_name: str,
                    g_dom: float = 0.24, g_min: float = 0.042) -> None:
    """
    Draw one V-tree inset diagram showing the state of adaptation and
    hybridisation at a given channel misalignment angle.

    BIOLOGICAL INTERPRETATION
    -------------------------
    A V-tree represents the evolutionary history of a cross:
      - Ancestor (grey dot, top): the common ancestor before divergence
      - Two populations (coloured ellipses, middle): the parental populations
        that have adapted to different environments
      - Hybrid (bottom ellipse): the F2 cross between the parents

    The G-matrix ellipse of each parent has:
      - Major axis = direction of maximum additive variance (g_max = the
        leading genetic channel in the Schluter limit)
      - Minor axis = direction of minimum additive variance

    KEY DESIGN PRINCIPLE — populations are placed ALONG their channel
    directions from the ancestor:
        p_i = ancestor - d * u_i

    This is biologically accurate: populations diverge by adapting along
    their lines of least resistance (the channels). The V-tree branches
    therefore point in the channel directions. At theta = 0 both channels
    coincide and the populations are at the same position (complete ellipse
    overlap). As theta grows the populations splay apart and their ellipses
    lose overlap — exactly mirroring the loss of channel alignment.

    The hybrid G-matrix V_H = ½(G₁ + G₂) is the additive-genetic
    variance expected in an F2 hybrid population under Mendelian segregation
    without dominance. Its ellipse shape shows how the hybrid variance is
    spread across directions:
      - theta = 0:  V_H = G₁ = G₂ — hybrid looks like either parent
      - theta = 90°: V_H is circular — variance equally spread everywhere
    The total variance (trace) is conserved regardless of theta.

    The hybrid ellipse COLOUR encodes I_geom: grey = no fitness cost (theta
    near 0), vermillion = high fitness cost (theta near 90°).

    PARAMETERS
    ----------
    ax : matplotlib Axes
    theta_deg : float  — channel misalignment angle in degrees
    stage_name : str   — label above the inset (e.g. "Early", "Late")
    g_dom : float      — sqrt of dominant eigenvalue of G (major axis half-length)
    g_min : float      — sqrt of minor eigenvalue of G (minor axis half-length)
    """
    ax.set_xlim(-0.88, 0.88)
    ax.set_ylim(-0.88, 0.88)
    ax.set_aspect("equal")
    ax.axis("off")

    theta = np.radians(theta_deg)

    # I_geom for k=2 channels (equal-angle formula: I = 1 - cos^{2k}(theta))
    I_geom = 1.0 - np.cos(theta) ** 4

    # ── Channel vectors ───────────────────────────────────────────────────────
    # Symmetric about the vertical: total misalignment = theta.
    u1 = np.array([ np.sin(theta / 2),  np.cos(theta / 2)])   # pop 1 channel
    u2 = np.array([-np.sin(theta / 2),  np.cos(theta / 2)])   # pop 2 channel
    w1 = np.array([-u1[1],  u1[0]])    # minor axis of pop 1
    w2 = np.array([-u2[1],  u2[0]])    # minor axis of pop 2

    def angle_of(u):
        return np.degrees(np.arctan2(u[1], u[0]))

    # ── Population positions — placed ALONG their channels from the ancestor ──
    # This is the key fix: p_i = ancestor − d·u_i
    #   • At theta=0: u1=u2=[0,1], so p1=p2 — ellipses fully overlap
    #   • At theta=20°: tiny separation — slight overlap (early divergence) ✓
    #   • At theta=50°: moderate gap (intermediate divergence) ✓
    #   • At theta=80°: large gap (late divergence, little overlap) ✓
    # The V-tree branches now point ALONG the channel directions — biologically
    # consistent with adaptation proceeding along lines of least resistance.
    anc = np.array([0.00, 0.58])
    d   = 0.30                         # displacement distance along channel
    p1  = anc - d * u1
    p2  = anc - d * u2
    hyb = np.array([0.00, -0.56])      # hybrid below both parents

    # ── V-tree lines ──────────────────────────────────────────────────────────
    # Ancestor → branching node (vertical drop), then diagonal arms to each
    # population. The arms point along the channel directions.
    ax.plot(*anc, "o", color="0.38", ms=3.8, zorder=5)
    node = anc - np.array([0, 0.12])   # branching node below ancestor
    ax.plot([anc[0], node[0]], [anc[1], node[1]], "-", color="0.52", lw=0.9, zorder=2)
    ax.plot([node[0], p1[0]+g_dom*0.5*u1[0]],
            [node[1], p1[1]+g_dom*0.5*u1[1]], "-", color="0.52", lw=0.9, zorder=2)
    ax.plot([node[0], p2[0]+g_dom*0.5*u2[0]],
            [node[1], p2[1]+g_dom*0.5*u2[1]], "-", color="0.52", lw=0.9, zorder=2)

    # ── Parental G-matrix ellipses and channel arrows ─────────────────────────
    # ORIENTATION: width=2·g_dom (major), height=2·g_min (minor), angle=angle_of(u).
    # matplotlib Ellipse places 'width' along the body x-axis before rotation,
    # so after rotating by angle_of(u) the major axis points along u. ✓
    for pos, u_ch, fill, edge in [
        (p1, u1, C_SKY,    C_BLUE),
        (p2, u2, C_ORANGE, C_ORANGE2),
    ]:
        ax.add_patch(Ellipse(
            pos, 2 * g_dom, 2 * g_min, angle=angle_of(u_ch),
            fc=fill, ec=edge, lw=1.1, alpha=0.8, zorder=4
        ))
        s = g_dom * 0.82
        ax.annotate(
            "", xy=(pos[0] + s * u_ch[0], pos[1] + s * u_ch[1]),
            xytext=(pos[0] - s * u_ch[0], pos[1] - s * u_ch[1]),
            arrowprops=dict(arrowstyle="->", color=edge, lw=1.3), zorder=6
        )

    # ── Hybrid G-matrix ellipse ───────────────────────────────────────────────
    # V_H = ½(G₁+G₂): additive variance in F2 hybrids under Mendelian
    # segregation, no dominance. Eigenvalues computed numerically.
    # Verified properties: V_H = G at theta=0; V_H circular at theta=90°;
    # tr(V_H) = tr(G) always (total variance conserved).
    G1 = g_dom**2 * np.outer(u1, u1) + g_min**2 * np.outer(w1, w1)
    G2 = g_dom**2 * np.outer(u2, u2) + g_min**2 * np.outer(w2, w2)
    V_hybrid = 0.5 * (G1 + G2)

    evals, evecs = np.linalg.eigh(V_hybrid)
    evals = np.abs(evals)[::-1]
    evecs = evecs[:, ::-1]
    ha = np.sqrt(evals[0])
    hb = np.sqrt(evals[1])

    grey       = np.array([0.70, 0.70, 0.70])
    vermillion = np.array([0.835, 0.369, 0.0])
    hybrid_colour = tuple((1 - I_geom) * grey + I_geom * vermillion)

    # Hybrid ellipse: width=2·ha (larger eigenvalue along x before rotation),
    # height=2·hb, angle=direction of first eigenvector. ✓
    ax.add_patch(Ellipse(
        hyb, 2 * ha, 2 * hb,
        angle=np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0])),
        fc=hybrid_colour, ec="0.30", lw=1.2, alpha=0.65, zorder=4
    ))

    # ── Labels ────────────────────────────────────────────────────────────────
    ax.text(0.0, 0.85, stage_name, ha="center", va="bottom",
            fontsize=6.5, color="0.30", style="italic")
    ax.text(0.0, -0.82,
            r"$\theta_{12} = " + f"{theta_deg}^\\circ$",
            ha="center", va="top", fontsize=6.5, color="0.20", fontweight="bold")


# =============================================================================
# SECTION 3: MAIN FIGURE
# =============================================================================

def make_figure(outfile: str = "figure_geometric_snowball_concept.pdf") -> None:
    """
    Build the full conceptual figure and save it as a PDF.

    LAYOUT
    ------
    One main axes (the snowball growth curve) plus three inset axes
    placed at the positions on the curve corresponding to the three
    snapshot stages.

    SNOWBALL CURVES
    ---------------
    The equal-angle formula I_geom = 1 - cos^{2k}(theta) is plotted for
    k=1 (one effective channel) and k=2 (two effective channels), where
    theta grows linearly with divergence time (omega*t = pi*t/2 maps
    t in [0,1] to theta in [0, 90 degrees]).

    For small angles, I_geom ≈ k*theta^2 ≈ k*(omega*t)^2, so isolaton
    grows QUADRATICALLY — the snowball. A linear reference line shows
    what linear accumulation would look like. The region where I_geom
    lies below the linear reference is shaded grey ("slow start"); the
    region where it exceeds the linear reference is shaded vermillion
    ("accelerating isolation").

    Biology: the slow start mirrors Orr's snowball — both mechanisms have
    a slow build-up when divergence is young. The acceleration phase is
    when the geometry of adaptation begins compounding.
    """

    # ── Snapshot positions ────────────────────────────────────────────────────
    # These are the three theta values (in degrees) at which the insets sit.
    theta_snaps = [20, 50, 80]
    # Map theta to "time" axis: t = theta / 90
    t_snaps = [th / 90 for th in theta_snaps]
    # I_geom values at each snapshot (k=2 formula)
    I_snaps = [1 - np.cos(np.radians(th)) ** 4 for th in theta_snaps]

    # ── Figure dimensions ─────────────────────────────────────────────────────
    # 183mm = full double-column width in Nature journals
    mm = 1.0 / 25.4
    fig = plt.figure(figsize=(183 * mm, 118 * mm))

    # ── Main axes ─────────────────────────────────────────────────────────────
    AX_L, AX_B, AX_W, AX_H = 0.09, 0.09, 0.875, 0.850
    ax = fig.add_axes([AX_L, AX_B, AX_W, AX_H])

    # ── Snowball curves ───────────────────────────────────────────────────────
    t       = np.linspace(0, 1, 600)
    omega_t = t * np.pi / 2        # theta in radians, grows from 0 to pi/2

    I_linear = t                                   # linear reference (no acceleration)
    I_k1     = 1 - np.cos(omega_t) ** 2           # k=1 channel (one trait dimension)
    I_k2     = 1 - np.cos(omega_t) ** 4           # k=2 channels (two trait dimensions)

    # Shading: where I_k2 < linear → slow start (grey)
    #          where I_k2 > linear → acceleration (light vermillion)
    ax.grid(True, color="0.93", lw=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
        sp.set_color("#c0c0c0")

    ax.fill_between(t, I_k2, I_linear,
                    where=(I_k2 < I_linear), color="0.83", alpha=0.32, zorder=1)
    ax.fill_between(t, I_k2, I_linear,
                    where=(I_k2 >= I_linear), color=C_VERM, alpha=0.13, zorder=1)

    ax.plot(t, I_linear, "--", color="0.44", lw=1.9, zorder=2, label="Linear reference")
    ax.plot(t, I_k1, "-", color=C_SKY, lw=1.5, alpha=0.72, zorder=3,
            label=r"Geometric snowball $k=1$")
    ax.plot(t, I_k2, "-", color=C_BLUE, lw=2.4, zorder=4,
            label=r"Geometric snowball $k=2$")

    # Dots on the curve at each snapshot position
    for t_s, I_s in zip(t_snaps, I_snaps):
        ax.plot(t_s, I_s, "o", color=C_BLUE, ms=4.5, zorder=7, clip_on=False)

    # ── Annotations on the main curve ────────────────────────────────────────
    # NOTE: "Accelerating isolation" annotation removed.
    # The curve has an inflection at t=1/3 (I≈0.44); it decelerates above that.
    # Any annotation position either sits in the wrong zone or looks like the
    # slow-start region. The formula box I≈k(ωt)² already conveys the
    # quadratic scaling argument without ambiguity.

    # Formula box — the key mathematical relationship
    ax.text(0.54, 0.95,
            r"$I_{\rm geom} \approx k\,(\omega t)^2$",
            ha="center", va="top", fontsize=7.5, color=C_BLUE,
            bbox=dict(fc="white", ec=C_BLUE, lw=0.6, pad=2.2,
                      boxstyle="round,pad=0.28"))

    ax.set_xlabel("Ecological divergence time", fontsize=8.5, labelpad=3)
    ax.set_ylabel(r"Hybrid fitness cost $I_{\rm geom}$", fontsize=8.5, labelpad=3)
    ax.set_xlim(-0.015, 1.015)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.legend(fontsize=6.5, loc="lower right",
              bbox_to_anchor=(1.00, 0.52), framealpha=0.93)

    # ── Place the three inset V-tree diagrams ─────────────────────────────────
    # Each inset is a small axes placed in figure-level coordinates.
    # INSET_W and INSET_H are in figure-fraction units.
    INSET_W, INSET_H = 0.155, 0.280

    # (left, bottom) coordinates for each inset in figure-fraction units
    inset_positions = [
        (AX_L + 0.015,          AX_B + AX_H * 0.580),   # Early: top-left
        (AX_L + AX_W * 0.365,   AX_B + 0.015),           # Intermediate: bottom-middle
        (AX_L + AX_W * 0.685,   AX_B + 0.015),           # Late: bottom-right
    ]
    stage_labels = ["Early", "Intermediate", "Late"]

    for (il, ib), theta, label in zip(inset_positions, theta_snaps, stage_labels):
        ax_inset = fig.add_axes([il, ib, INSET_W, INSET_H])
        draw_v_snapshot(ax_inset, theta, label)
        ax_inset.patch.set_alpha(0.0)   # transparent background

    # ── Save ─────────────────────────────────────────────────────────────────
    outpath = OUTDIR / outfile
    plt.savefig(outpath, format="pdf", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {outpath}")
    plt.close(fig)


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    print("Building integrated snowball concept figure ...")
    make_figure("figure_geometric_snowball_concept.pdf")
