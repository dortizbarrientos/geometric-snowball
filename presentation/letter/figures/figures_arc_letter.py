# =============================================================================
# figures_arc_letter.py
# Python 3 — Figures for "A Letter to the ARC Panellist"
# Daniel & Mark Cooper
#
# Concept
# -------
# Six figures covering every major conceptual stage of the letter:
#   Fig 1  (Section I)   — G matrix: family resemblance & genetic cloud
#   Fig 2  (Section I)   — Fitness surface shapes: what γ means
#   Fig 3  (Section I)   — Additive channel vs curvature load
#   Fig 4  (Section III) — Darwin's finches: normal year vs drought
#   Fig 5  (Section IV)  — Sorghum root depth & the yield cliff
#   Fig 6  (Section V)   — Channel switching: position and heading
#
# Dependencies
# ------------
#   pip install numpy matplotlib scipy
#
# Usage
# -----
#   python figures_arc_letter.py
#
# Output
# ------
#   Six PDF files (300 dpi), one per figure, in the working directory.
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.special import expit          # logistic sigmoid σ(x) = 1/(1+e^{-x})

# ─────────────────────────────────────────────────────────────────────────────
# 0.  Global style
#     All figures share one style block so they look like a single publication.
# ─────────────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    'font.family'       : 'DejaVu Sans',
    'font.size'         : 9,
    'axes.spines.top'   : False,
    'axes.spines.right' : False,
    'axes.linewidth'    : 0.6,
    'xtick.major.width' : 0.6,
    'ytick.major.width' : 0.6,
    'xtick.major.size'  : 3,
    'ytick.major.size'  : 3,
    'figure.dpi'        : 150,
    'savefig.dpi'       : 300,
    'savefig.bbox'      : 'tight',
    'savefig.facecolor' : 'white',
})

# ─── Semantic colour palette ─────────────────────────────────────────────────
#   The same colour encodes the same concept across all six figures.

CG     = '#1D9E75'   # teal-green  : G matrix, genetic cloud, surviving individuals
CGMAX  = '#085041'   # dark teal   : g_max arrow
CFIT   = '#D85A30'   # coral       : fitness surface, γ matrix, severe environment
CAG    = '#7F77DD'   # soft purple : Ag index, channel-switching dynamics
CMILD  = '#378ADD'   # sky blue    : mild-stress environment
CCLIFF = '#A32D2D'   # dark red    : fitness threshold / cliff
CGRAY  = '#888780'   # mid-gray    : axis text, annotations
CLGRAY = '#D3D1C7'   # light-gray  : background scatter, gridlines


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Shared drawing helpers
#     These are called in every figure to guarantee visual consistency.
# ─────────────────────────────────────────────────────────────────────────────

def g_ellipse(ax, cx, cy, rx, ry, angle_deg,
              color=CG, alpha_fill=0.10, lw=1.4):
    """
    Draw the G-matrix ellipse as a filled region + coloured edge.

    Parameters
    ----------
    ax         : matplotlib Axes
    cx, cy     : centre of the ellipse (data coordinates)
    rx, ry     : semi-axes: rx along the long axis, ry along the short axis
    angle_deg  : rotation of the long axis from the x-axis (degrees)
    color      : edge and fill colour (default: teal-green for G)
    alpha_fill : transparency of the filled interior
    lw         : edge line width

    Notes for the biologist
    -----------------------
    The G-matrix ellipse summarises where genetic variation is concentrated
    in two-trait space.  The long axis points in the direction families
    differ most (g_max).  The ellipse is the 2-sigma contour of the
    bivariate genetic distribution.
    """
    kw = dict(angle=angle_deg, linewidth=lw)
    ax.add_patch(Ellipse((cx, cy), 2*rx, 2*ry,
                         color=color, alpha=alpha_fill, **kw))
    ax.add_patch(Ellipse((cx, cy), 2*rx, 2*ry,
                         fill=False, edgecolor=color, **kw))


def gmax_arrow(ax, cx, cy, angle_deg, half_len,
               color=CGMAX, lw=1.8, fontsize=9):
    """
    Draw a centred arrow labelled g_max pointing in the long-axis direction.

    The arrow runs from (cx - 0.4·half_len) in the angle_deg direction
    to (cx + half_len) in that direction, with an arrowhead at the far end.
    """
    rad      = np.radians(angle_deg)
    dx, dy   = np.cos(rad) * half_len, np.sin(rad) * half_len
    ax.annotate('',
                xy=(cx + dx, cy + dy),
                xytext=(cx - dx * 0.4, cy - dy * 0.4),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                mutation_scale=10))
    ax.text(cx + dx * 1.22, cy + dy * 1.22,
            r'$\mathbf{g}_{\mathrm{max}}$',
            color=color, fontsize=fontsize, ha='center', va='center',
            fontweight='bold')


def panel_label(ax, letter, x=-0.10, y=1.05):
    """Bold panel label (A, B …) placed outside the upper-left corner."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top', ha='left',
            color='#2C2C2A')


def style_ax(ax, xlabel='', ylabel='', title=''):
    """Apply minimal, clean axis styling."""
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9, color=CGRAY)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9, color=CGRAY)
    if title:
        ax.set_title(title, fontsize=9.5, pad=5,
                     color='#2C2C2A', fontweight='bold')
    ax.tick_params(labelsize=8, colors=CGRAY, width=0.6)
    for sp in ['bottom', 'left']:
        ax.spines[sp].set_color(CLGRAY)


def gmax_from_cov(cov):
    """
    Return (angle_degrees, long_semi_axis, short_semi_axis, unit_vector)
    for the leading eigenvector of a 2x2 covariance matrix.

    The eigenvector is normalised so its angle lies in (-90°, 90°],
    ensuring the arrow points into the upper half of trait space.

    Notes for the biologist
    -----------------------
    The leading eigenvector of G is g_max: the direction of maximum
    genetic variation.  The eigenvalue is the variance in that direction.
    """
    eigvals, eigvecs = np.linalg.eigh(np.array(cov, dtype=float))
    idx     = np.argmax(eigvals)
    vec     = eigvecs[:, idx]
    angle   = np.degrees(np.arctan2(vec[1], vec[0]))
    # Normalise to (-90, 90] so the arrow points 'forward'
    if angle > 90:
        angle -= 180
    elif angle <= -90:
        angle += 180
    rx  = 2.0 * np.sqrt(eigvals[idx])        # 2-sigma semi-axis
    ry  = 2.0 * np.sqrt(eigvals[1 - idx])
    u1  = np.array([np.cos(np.radians(angle)),
                    np.sin(np.radians(angle))])
    return angle, rx, ry, u1


# =============================================================================
# FIGURE 1 — The G Matrix
# =============================================================================
# Placement in letter: after "Visualise G as an ellipse drawn in the space
#                       of all traits…"
#
# Panel A  single-trait parent-offspring scatter → heritability
# Panel B  two-trait scatter with G ellipse + g_max arrow
#
# Key message: family resemblance is measurable and the G matrix captures
# both the heritability of each trait (Panel A) and their genetic correlation
# (Panel B, ellipse orientation).
# =============================================================================

def make_fig1():
    rng  = np.random.default_rng(42)
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.5))

    # ── Panel A : parent-offspring scatter for a single trait ─────────────
    ax   = axes[0]
    h2   = 0.65     # heritability = slope of parent-offspring regression
    n    = 80

    # Generate parent values, then offspring = h² × parent + residual
    parent    = rng.normal(0, 1, n)
    offspring = h2 * parent + np.sqrt(1 - h2**2) * rng.normal(0, 1, n)

    ax.scatter(parent, offspring, s=18, color=CLGRAY, alpha=0.75,
               linewidths=0, zorder=2)

    # Regression line (slope = heritability)
    x_line = np.linspace(-2.8, 2.8, 200)
    ax.plot(x_line, h2 * x_line, color=CG, lw=1.8, zorder=3)
    ax.axhline(0, color=CLGRAY, lw=0.5)
    ax.axvline(0, color=CLGRAY, lw=0.5)

    # Annotate slope
    ax.annotate(rf'slope $= h^2 = {h2}$',
                xy=(1.8, h2 * 1.8), xytext=(0.3, -1.5),
                fontsize=8, color=CG,
                arrowprops=dict(arrowstyle='->', color=CG, lw=0.8))

    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    ax.set_aspect('equal')
    style_ax(ax,
             xlabel='Parent trait value',
             ylabel='Offspring trait value',
             title='A  Family resemblance')
    ax.text(0.04, 0.97,
            'Each point = one family\n'
            'Slope = heritability $h^2$',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    panel_label(ax, 'A')

    # ── Panel B : two-trait G ellipse ─────────────────────────────────────
    ax   = axes[1]

    # Genetic covariance matrix (two traits genetically correlated)
    cov_G = [[1.6, 1.0],
             [1.0, 0.9]]

    pts = rng.multivariate_normal([0, 0], cov_G, n)
    ax.scatter(pts[:, 0], pts[:, 1], s=18, color=CLGRAY, alpha=0.75,
               linewidths=0, zorder=2)

    angle, rx, ry, _ = gmax_from_cov(cov_G)
    g_ellipse(ax, 0, 0, rx, ry, angle)
    gmax_arrow(ax, 0, 0, angle, rx * 0.70)

    # Label the ellipse
    ax.text(-2.6, 2.6, r'$\mathbf{G}$ ellipse', color=CG,
            fontsize=8.5, ha='left', va='top')

    ax.axhline(0, color=CLGRAY, lw=0.5)
    ax.axvline(0, color=CLGRAY, lw=0.5)
    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    ax.set_aspect('equal')
    style_ax(ax,
             xlabel='Trait 1 (e.g. beak depth)',
             ylabel='Trait 2 (e.g. beak width)',
             title='B  The G matrix ellipse')
    ax.text(0.04, 0.97,
            r'Long axis = $\mathbf{g}_{\max}$' '\n'
            'Direction families differ most',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    panel_label(ax, 'B')

    fig.tight_layout(w_pad=3.0)
    fig.savefig('fig1_G_matrix.pdf')
    print('Saved  fig1_G_matrix.pdf')
    plt.close(fig)


# =============================================================================
# FIGURE 2 — Fitness Surface Shapes: What γ Means
# =============================================================================
# Placement in letter: after "…that sudden steepening is what we mean when
#                       we say the fitness surface bends."
#
# Four panels showing how the curvature term γ changes the shape of the
# fitness–trait relationship: flat, stabilising, disruptive, cliff.
#
# Key message: γ is not a vague notion of 'selection' — it is the second
# derivative of the fitness surface.  Its sign and magnitude tell you
# exactly how the surface bends.
# =============================================================================

def make_fig2():
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.2),
                             gridspec_kw={'wspace': 0.35})

    x = np.linspace(-2.5, 2.5, 300)

    # Each entry: (title, fitness_function, gamma_label, colour, footnote)
    panels = [
        (
            'Linear\n(γ ≈ 0)',
            lambda x: 0.55 * x,
            r'$\gamma \approx 0$',
            CMILD,
            'Proportional: prediction\nreliable throughout'
        ),
        (
            'Stabilising\n(γ < 0)',
            lambda x: -0.55 * x**2 + 0.6,
            r'$\gamma < 0$',
            CG,
            'Intermediate trait values\nfavoured (peak fitness)'
        ),
        (
            'Disruptive\n(γ > 0)',
            lambda x: 0.50 * x**2 - 0.5,
            r'$\gamma > 0$',
            CFIT,
            'Extremes favoured\n(valley fitness)'
        ),
        (
            'Cliff\n(|γ| large)',
            lambda x: expit(5.0 * (x - 0.6)) - 0.5,
            r'$|\gamma|$ large',
            CCLIFF,
            'Abrupt threshold: small\ntrait difference → large\nfitness collapse'
        ),
    ]

    for ax, (title, fn, glbl, col, note) in zip(axes, panels):
        y = fn(x)
        ax.plot(x, y, color=col, lw=2.2, zorder=3)
        ax.fill_between(x, y, y.min() - 0.4,
                        color=col, alpha=0.09, zorder=2)
        ax.axhline(0, color=CLGRAY, lw=0.5)

        # γ label in upper-right corner
        ax.text(0.97, 0.96, glbl,
                transform=ax.transAxes, fontsize=9.5,
                ha='right', va='top', color=col, fontweight='bold')

        # Explanatory note in lower-left
        ax.text(0.04, 0.04, note,
                transform=ax.transAxes, fontsize=7.5,
                color=CGRAY, va='bottom', linespacing=1.4)

        ax.set_xlim(-2.6, 2.6)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['bottom'].set_color(CLGRAY)
        ax.spines['left'].set_color(CLGRAY)
        ax.set_title(title, fontsize=9.5, pad=5,
                     color='#2C2C2A', fontweight='bold')

        # Add threshold marker on cliff panel
        if 'Cliff' in title:
            ax.axvline(0.6, color=CCLIFF, lw=1.2, ls='--', alpha=0.7)
            ax.text(0.64, 0.54, 'threshold',
                    color=CCLIFF, fontsize=7.5,
                    transform=ax.transAxes)

    axes[0].set_ylabel('Fitness (reproductive success)',
                       fontsize=9, color=CGRAY)

    # Shared x-axis label via figure text
    fig.text(0.52, -0.04, 'Trait value →',
             ha='center', fontsize=9, color=CGRAY)

    fig.suptitle('Four shapes of the fitness surface — what γ captures',
                 y=1.02, fontsize=10, fontweight='bold', color='#2C2C2A')
    fig.savefig('fig2_gamma_shapes.pdf')
    print('Saved  fig2_gamma_shapes.pdf')
    plt.close(fig)


# =============================================================================
# FIGURE 3 — Additive Channel vs Curvature Load
# =============================================================================
# Placement in letter: after the curvature-load paragraph, before the table.
#
# Both panels show the same two-trait space.  The G ellipse and g_max are
# identical.  What changes is the fitness surface underneath.
#
# Panel A  g_max crosses gently varying contours → additive channel
# Panel B  g_max crosses a sharp fitness cliff   → curvature load
#
# Key message: the channel/load distinction is not about G or γ alone —
# it is about whether g_max and the high-γ region coincide.
# =============================================================================

def make_fig3():
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.0))

    lim  = 2.8
    res  = 300
    x1, x2 = np.meshgrid(np.linspace(-lim, lim, res),
                         np.linspace(-lim, lim, res))

    # g_max at 45°: u1 = (1,1)/√2, u2 = (-1,1)/√2
    angle_deg = 45.0
    rad = np.radians(angle_deg)
    u1  = np.array([np.cos(rad), np.sin(rad)])    # g_max direction
    u2  = np.array([-np.sin(rad), np.cos(rad)])   # perpendicular

    # Projections of every grid point onto g_max (proj1) and perp (proj2)
    proj1 = x1 * u1[0] + x2 * u1[1]
    proj2 = x1 * u2[0] + x2 * u2[1]

    # G ellipse parameters (long axis at 45°)
    rx_G, ry_G = 1.7, 0.6

    # ── Panel A : additive channel ────────────────────────────────────────
    # Fitness is nearly linear in the g_max direction (no curvature along it)
    # and varies smoothly perpendicular to it.
    # γ along g_max ≈ 0 → additive channel open.
    w_channel = 0.5 * proj2 + 0.3 * proj1   # linear; curvature = 0

    ax = axes[0]
    ax.contourf(x1, x2, w_channel,
                levels=12, cmap='YlGn', alpha=0.45, zorder=1)
    ax.contour(x1, x2, w_channel,
               levels=8, colors='white', linewidths=0.4, alpha=0.6, zorder=2)

    g_ellipse(ax, 0, 0, rx_G, ry_G, angle_deg)
    gmax_arrow(ax, 0, 0, angle_deg, rx_G * 0.75)

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    style_ax(ax,
             xlabel='Trait 1', ylabel='Trait 2',
             title='A  Additive channel open')
    ax.text(0.04, 0.97,
            r'$\mathbf{g}_{\max}$ runs along gentle contours' '\n'
            r'$\gamma$ small in $\mathbf{g}_{\max}$ direction' '\n'
            'Proportional prediction reliable',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    ax.text(0.5, -0.08, 'colour = fitness (dark green: high)',
            transform=ax.transAxes, ha='center', fontsize=7.5, color=CGRAY)
    panel_label(ax, 'A')

    # ── Panel B : curvature load ──────────────────────────────────────────
    # Fitness is a sigmoid in the g_max direction: near-zero below threshold,
    # moderate above.  This creates a sharp cliff that crosses g_max.
    # γ large in g_max direction → curvature load.
    threshold = 0.5
    w_load = expit(4.5 * (proj1 - threshold))

    ax = axes[1]
    ax.contourf(x1, x2, w_load,
                levels=12, cmap='RdYlGn', alpha=0.45, zorder=1)
    ax.contour(x1, x2, w_load,
               levels=8, colors='white', linewidths=0.4, alpha=0.5, zorder=2)

    # Cliff line: points where proj1 = threshold → x2 = threshold·√2 - x1
    x_cl = np.linspace(-lim, lim, 300)
    y_cl = threshold * np.sqrt(2) - x_cl
    mask = (y_cl > -lim) & (y_cl < lim)
    ax.plot(x_cl[mask], y_cl[mask],
            color=CCLIFF, lw=2.0, ls='--', alpha=0.9, zorder=5)
    ax.text(lim - 0.1, threshold * np.sqrt(2) - lim + 0.2,
            'fitness cliff', color=CCLIFF, fontsize=7.5,
            ha='right', va='bottom', style='italic')

    g_ellipse(ax, 0, 0, rx_G, ry_G, angle_deg)
    gmax_arrow(ax, 0, 0, angle_deg, rx_G * 0.75)

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    style_ax(ax,
             xlabel='Trait 1', ylabel='Trait 2',
             title='B  Curvature load')
    ax.text(0.04, 0.97,
            r'$\mathbf{g}_{\max}$ runs into the cliff' '\n'
            r'$\gamma$ large in $\mathbf{g}_{\max}$ direction' '\n'
            'Proportional prediction fails',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    ax.text(0.5, -0.08, 'colour = fitness (green: high, red: low)',
            transform=ax.transAxes, ha='center', fontsize=7.5, color=CGRAY)
    panel_label(ax, 'B')

    fig.tight_layout(w_pad=3.0)
    fig.savefig('fig3_channels.pdf')
    print('Saved  fig3_channels.pdf')
    plt.close(fig)


# =============================================================================
# FIGURE 4 — Darwin's Finches: Normal Year vs Drought Year
# =============================================================================
# Placement in letter: after the finch section's curvature-load conclusion,
#                      before "The lesson."
#
# Both panels share the same G ellipse and g_max.  The fitness surface
# underneath changes between years — this is the core of the finch example.
#
# Panel A  normal year: small soft seeds; fitness surface nearly flat; γ small
# Panel B  drought year: hard seeds only; cliff in g_max direction; γ large
#          Individuals are coloured: green = survived, red cross = did not
#
# Key message: nothing in the inheritance system changed between years.
# What changed was the shape of the fitness surface relative to g_max.
# =============================================================================

def make_fig4():
    rng  = np.random.default_rng(7)
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.0))

    lim  = 2.8
    res  = 300
    x1, x2 = np.meshgrid(np.linspace(-lim, lim, res),
                         np.linspace(-lim, lim, res))

    # G matrix for beak depth × beak width (two positively correlated traits)
    cov_G = [[1.4, 0.9],
             [0.9, 0.9]]
    angle, rx, ry, u1 = gmax_from_cov(cov_G)
    u2 = np.array([-u1[1], u1[0]])     # perpendicular to g_max

    # Project grid onto g_max (proj1) and its perpendicular (proj2)
    proj1 = x1 * u1[0] + x2 * u1[1]
    proj2 = x1 * u2[0] + x2 * u2[1]

    # Simulate 90 individual birds as a bivariate normal cloud
    birds  = rng.multivariate_normal([0, 0], cov_G, 90)

    # ── Panel A : normal year (small, soft seeds) ─────────────────────────
    # Fitness surface is nearly linear — small γ.
    w_normal = 0.4 * proj1 + 0.15 * proj2

    ax = axes[0]
    ax.contourf(x1, x2, w_normal,
                levels=10, cmap='YlGn', alpha=0.40, zorder=1)
    ax.contour(x1, x2, w_normal,
               levels=8, colors='white', linewidths=0.4, alpha=0.5, zorder=2)
    ax.scatter(birds[:, 0], birds[:, 1], s=18, color=CG, alpha=0.55,
               linewidths=0, zorder=3)
    g_ellipse(ax, 0, 0, rx, ry, angle)
    gmax_arrow(ax, 0, 0, angle, rx * 0.75)

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    style_ax(ax,
             xlabel='Beak depth', ylabel='Beak width',
             title='A  Normal year — small soft seeds')
    ax.text(0.04, 0.97,
            r'$\gamma$ small — surface nearly flat' '\n'
            r'along $\mathbf{g}_{\max}$' '\n'
            r'$\mathcal{A}_g$ small — additive channel open',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    panel_label(ax, 'A')

    # ── Panel B : drought year (hard seeds only — mechanical threshold) ────
    # Fitness collapses below a critical beak-force threshold.
    # The cliff runs perpendicular to g_max and crosses the main cloud.
    threshold = 0.3
    w_drought = expit(5.0 * (proj1 - threshold))

    ax = axes[1]
    ax.contourf(x1, x2, w_drought,
                levels=10, cmap='RdYlGn', alpha=0.40, zorder=1)
    ax.contour(x1, x2, w_drought,
               levels=8, colors='white', linewidths=0.4, alpha=0.5, zorder=2)

    # Cliff line: proj1 = threshold → x2 = (threshold - x1·u1[0]) / u1[1]
    x_cl = np.linspace(-lim, lim, 300)
    y_cl = (threshold - x_cl * u1[0]) / u1[1]
    mask = (y_cl > -lim) & (y_cl < lim)
    ax.plot(x_cl[mask], y_cl[mask],
            color=CCLIFF, lw=2.0, ls='--', alpha=0.9, zorder=5)
    ax.text(lim - 0.1, (threshold - (lim - 0.1) * u1[0]) / u1[1] + 0.1,
            'seed-cracking\nthreshold',
            color=CCLIFF, fontsize=7.0, ha='right', va='bottom', style='italic')

    # Colour birds by whether they clear the threshold
    bird_proj1 = birds[:, 0] * u1[0] + birds[:, 1] * u1[1]
    survived   = bird_proj1 >= threshold
    ax.scatter(birds[~survived, 0], birds[~survived, 1], s=22,
               color=CCLIFF, alpha=0.55, linewidths=0.5, zorder=4,
               marker='x', label='Did not survive')
    ax.scatter(birds[survived, 0], birds[survived, 1], s=18,
               color=CG, alpha=0.60, linewidths=0, zorder=4,
               label='Survived')

    g_ellipse(ax, 0, 0, rx, ry, angle)
    gmax_arrow(ax, 0, 0, angle, rx * 0.75)

    ax.legend(fontsize=7, loc='lower right', frameon=False,
              labelcolor=CGRAY, markerscale=1.2)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    style_ax(ax,
             xlabel='Beak depth', ylabel='Beak width',
             title='B  Drought year — hard seeds only')
    ax.text(0.04, 0.97,
            r'$\gamma$ large in $\mathbf{g}_{\max}$ direction' '\n'
            'Cliff crosses main genetic cloud\n'
            r'$\mathcal{A}_g$ large — additive prediction fails',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    panel_label(ax, 'B')

    fig.tight_layout(w_pad=3.0)
    fig.savefig('fig4_finches.pdf')
    print('Saved  fig4_finches.pdf')
    plt.close(fig)


# =============================================================================
# FIGURE 5 — Sorghum: Root Depth and the Yield Cliff
# =============================================================================
# Placement in letter: after "…the surface has bent: sharply, suddenly,
#                       and in exactly the place the population's variation
#                       is concentrated."
#
# Panel A  mild stress: yield rises proportionally with root depth;
#          no threshold; γ small; additive channel open
# Panel B  severe terminal drought: abrupt cliff at 60 cm (subsoil moisture);
#          near-zero yield below; γ large; same lines, very different outcome
#
# Key message: the same breeding programme, the same genomic model, and the
# same lines — but the failure is purely geometric.
# =============================================================================

def make_fig5():
    rng    = np.random.default_rng(99)
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8))

    root_depth  = np.linspace(20, 100, 300)   # cm
    n_lines     = 65
    line_depths = rng.uniform(28, 95, n_lines)  # root depths of breeding lines

    # ── Panel A : mild stress ─────────────────────────────────────────────
    def yield_mild(d, noise=0):
        """Proportional yield–root-depth curve under mild stress."""
        return 2.0 + 0.065 * (d - 20) + noise

    ax   = axes[0]
    y_m  = yield_mild(root_depth)

    ax.plot(root_depth, y_m, color=CMILD, lw=2.2, zorder=3)
    ax.fill_between(root_depth, y_m - 0.35, y_m + 0.35,
                    color=CMILD, alpha=0.12, zorder=2)

    # Scatter of individual breeding lines
    noise_m = rng.normal(0, 0.28, n_lines)
    ax.scatter(line_depths, yield_mild(line_depths, noise_m),
               s=22, color=CMILD, alpha=0.55, linewidths=0, zorder=4)

    style_ax(ax,
             xlabel='Root depth (cm)',
             ylabel='Grain yield (t/ha)',
             title='A  Mild stress — additive channel open')
    ax.set_xlim(18, 102); ax.set_ylim(1.2, 8.0)
    ax.text(0.04, 0.97,
            'Yield rises proportionally\nwith root depth\n'
            r'$\gamma$ small — genomic selection reliable',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    panel_label(ax, 'A')

    # ── Panel B : severe terminal drought ─────────────────────────────────
    threshold = 60.0   # cm — depth of subsoil moisture layer in these soils

    def yield_severe(d, noise=0):
        """
        Yield as a logistic function of root depth relative to the subsoil
        moisture threshold.  Below threshold: near-zero.  Above: moderate.
        The sigmoid mimics the sharp mechanical threshold described in the text.
        """
        return 3.8 * expit(0.50 * (d - threshold)) + 0.25 + noise

    ax   = axes[1]
    y_s  = yield_severe(root_depth)

    ax.plot(root_depth, y_s, color=CFIT, lw=2.2, zorder=3)
    ax.fill_between(root_depth, y_s - 0.30, y_s + 0.30,
                    color=CFIT, alpha=0.10, zorder=2)

    # Threshold line
    ax.axvline(threshold, color=CCLIFF, lw=1.5, ls='--', alpha=0.85, zorder=5)
    ax.text(threshold + 1.0, 7.5,
            f'{threshold:.0f} cm\nsubsoil\nmoisture',
            color=CCLIFF, fontsize=7.5, va='top', ha='left')

    # Same breeding lines — now coloured by which side of the threshold they fall
    noise_s = rng.normal(0, 0.24, n_lines)
    y_lines = yield_severe(line_depths, noise_s)
    above   = line_depths >= threshold
    ax.scatter(line_depths[above], y_lines[above], s=22,
               color=CG, alpha=0.60, linewidths=0, zorder=4,
               label='Roots reached moisture layer')
    ax.scatter(line_depths[~above], y_lines[~above], s=22,
               color=CCLIFF, alpha=0.55, linewidths=0.5, zorder=4,
               marker='x', label='Roots fell short')

    # Annotate the cliff — show two lines only 16 cm apart but very different
    d_lo, d_hi = threshold - 8, threshold + 8
    y_lo = yield_severe(d_lo)
    y_hi = yield_severe(d_hi)
    ax.annotate('',
                xy=(d_hi, y_hi), xytext=(d_lo, y_lo),
                arrowprops=dict(arrowstyle='<->', color='#2C2C2A',
                                lw=1.0, mutation_scale=8))
    ax.text(d_lo - 1.5, (y_lo + y_hi) / 2,
            '16 cm apart\n→ large yield\ngap',
            fontsize=7.0, color='#2C2C2A', va='center', ha='right')

    ax.legend(fontsize=7, loc='upper left', frameon=False, labelcolor=CGRAY)
    style_ax(ax,
             xlabel='Root depth (cm)',
             ylabel='Grain yield (t/ha)',
             title='B  Severe drought — curvature load')
    ax.set_xlim(18, 102); ax.set_ylim(0.1, 8.0)
    ax.text(0.55, 0.12,
            r'$\gamma$ large — cliff at threshold' '\n'
            r'$\mathcal{A}_g$ large — same model misfires',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='bottom', linespacing=1.5)
    panel_label(ax, 'B')

    fig.tight_layout(w_pad=3.0)
    fig.savefig('fig5_sorghum.pdf')
    print('Saved  fig5_sorghum.pdf')
    plt.close(fig)


# =============================================================================
# FIGURE 6 — Channel Switching: Position and Heading
# =============================================================================
# Placement in letter: after "The key conceptual advance is that you can tell
#                       which direction you are moving before you cross the
#                       boundary."
#
# Panel A  Ag over time — the "dashboard" metaphor.
#          Ag is small (channel open), rises (channel closing), then falls
#          (new channel opens after disturbance).  The heading arrow marks
#          the moment of warning — Ag rising before the cliff is reached.
#
# Panel B  Population cloud migrating through trait space toward a fixed cliff.
#          Three snapshots: far from cliff (generation 1), approaching
#          (generation 2), and post-disturbance (generation 3 — narrowed,
#          re-aligned, new channel open).
#
# Key message: Ag tells you where you are.  Its rate of change tells you
# where you are going.  The second instrument gives you warning in advance.
# =============================================================================

def make_fig6():
    rng  = np.random.default_rng(55)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))

    # ── Panel A : Ag over time ────────────────────────────────────────────
    ax   = axes[0]
    t    = np.linspace(0, 1, 300)

    # Ag profile: starts low, rises as population drifts toward cliff,
    # peaks at disturbance, then falls as new channel opens.
    ag = (
        0.08
        + 0.80 * np.exp(-18.0 * (t - 0.50)**2)
    )
    ag = np.clip(ag, 0.02, 0.95)

    # Shaded background regions
    threshold_ag = 0.35   # notional threshold between open and loaded
    ax.axhspan(0.0,          threshold_ag, color=CG,   alpha=0.07, zorder=1)
    ax.axhspan(threshold_ag, 1.0,          color=CFIT, alpha=0.07, zorder=1)
    ax.axhline(threshold_ag, color=CGRAY,
               lw=0.7, ls=':', alpha=0.6, zorder=2)

    ax.plot(t, ag, color=CAG, lw=2.2, zorder=3)
    ax.fill_between(t, ag, 0,
                    color=CAG, alpha=0.12, zorder=2)

    # Phase labels
    ax.text(0.12, 0.26, 'Channel open',
            ha='center', fontsize=7.5, color=CG,
            transform=ax.transAxes, va='center', style='italic')
    ax.text(0.50, 0.92, 'Disturbance\n(channel closed)',
            ha='center', fontsize=7.5, color=CCLIFF,
            transform=ax.transAxes, va='top')
    ax.text(0.88, 0.26, 'New channel\nopen',
            ha='center', fontsize=7.5, color=CG,
            transform=ax.transAxes, va='center', style='italic')

    # "Heading" arrow: rising Ag is the warning signal
    i0, i1 = int(0.33 * 299), int(0.48 * 299)
    ax.annotate('',
                xy=(t[i1], ag[i1]),
                xytext=(t[i0], ag[i0]),
                arrowprops=dict(arrowstyle='->', color=CAG, lw=1.4))
    ax.text(0.385, ag[int(0.405*299)] + 0.04,
            'heading:\nchannel closing',
            ha='center', fontsize=7.0, color=CAG, style='italic')

    # Vertical marker at disturbance peak
    t_peak = t[np.argmax(ag)]
    ax.axvline(t_peak, color=CCLIFF, lw=0.8, ls='--', alpha=0.6)

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks([0, t_peak, 1.0])
    ax.set_xticklabels(['Generation 0', 'Disturbance', ''],
                       fontsize=8, color=CGRAY)
    ax.set_yticks([0, threshold_ag, 1.0])
    ax.set_yticklabels(['0', 'threshold', '1'],
                       fontsize=8, color=CGRAY)
    style_ax(ax,
             ylabel=r'Additivity index $\mathcal{A}_g$',
             title=r'A  $\mathcal{A}_g$ over time')
    ax.text(0.96, threshold_ag + 0.02, r'$\mathcal{A}_g$ threshold',
            color=CGRAY, fontsize=7, va='bottom', ha='right',
            transform=ax.get_yaxis_transform())
    panel_label(ax, 'A')

    # ── Panel B : population cloud migrating toward fitness cliff ─────────
    ax  = axes[1]
    lim = 3.0

    # Fixed fitness cliff: points where x1 + x2 = cliff_level
    cliff_level = 1.5
    xx, yy = np.meshgrid(np.linspace(-lim, lim, 200),
                         np.linspace(-lim, lim, 200))
    proj_bg = (xx + yy) / np.sqrt(2)
    w_bg    = expit(3.5 * (proj_bg - cliff_level / np.sqrt(2)))
    ax.contourf(xx, yy, w_bg, levels=10, cmap='RdYlGn', alpha=0.28, zorder=1)

    # Cliff line
    x_cl = np.linspace(-lim, lim, 300)
    y_cl = cliff_level - x_cl
    mask = (y_cl > -lim) & (y_cl < lim)
    ax.plot(x_cl[mask], y_cl[mask],
            color=CCLIFF, lw=2.0, ls='--', alpha=0.9, zorder=5)
    ax.text(lim - 0.1, cliff_level - lim + 0.15,
            'fitness cliff', color=CCLIFF, fontsize=7.5,
            ha='right', va='bottom', style='italic')

    # Three time-point G ellipses — same g_max direction, migrating cloud
    # Generation 1: far from cliff, full cloud
    # Generation 2: cloud centre near threshold — approaching
    # Generation 3: post-selection, narrowed cloud, past cliff, new channel
    snapshots = [
        # (centre_x, centre_y, rx, ry, colour, label_offset_x, label_offset_y, label)
        (-1.6, -1.6, 0.90, 0.30, CMILD,  0.0, -0.55, 'Generation 1\n(channel open)'),
        (-0.2, -0.2, 0.90, 0.30, CAG,    0.0, -0.55, 'Generation 2\n(Ag rising — warning)'),
        ( 1.4,  1.4, 0.55, 0.20, CG,     0.0, -0.45, 'Generation 3\n(new channel open)'),
    ]

    centres_x, centres_y = [], []
    for (cx, cy, rxv, ryv, col, lox, loy, lbl) in snapshots:
        g_ellipse(ax, cx, cy, rxv, ryv, angle_deg=45,
                  color=col, alpha_fill=0.14, lw=1.5)
        # g_max mini-arrow
        r = np.radians(45)
        dx, dy = np.cos(r) * rxv * 0.85, np.sin(r) * ryv * 0.85 * (ryv / rxv) * (rxv / 0.30)
        # simpler: fixed direction, scaled by rxv
        dx, dy = np.cos(r) * rxv * 0.8, np.sin(r) * rxv * 0.8
        ax.annotate('', xy=(cx+dx, cy+dy), xytext=(cx-dx*0.4, cy-dy*0.4),
                    arrowprops=dict(arrowstyle='->', color=col, lw=1.3,
                                    mutation_scale=8))
        ax.text(cx + lox, cy + loy, lbl,
                color=col, fontsize=7.0, ha='center', va='top')
        centres_x.append(cx); centres_y.append(cy)

    # Trajectory arrow connecting the three centres
    ax.plot(centres_x, centres_y, color=CGRAY, lw=0.7, ls=':', zorder=3, alpha=0.7)
    for i in range(len(centres_x) - 1):
        ax.annotate('',
                    xy=(centres_x[i+1], centres_y[i+1]),
                    xytext=(centres_x[i], centres_y[i]),
                    arrowprops=dict(arrowstyle='->', color=CGRAY,
                                    lw=0.9, alpha=0.6, mutation_scale=8))

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    style_ax(ax,
             xlabel='Trait 1', ylabel='Trait 2',
             title='B  Population cloud migrating through trait space')
    ax.text(0.04, 0.97,
            'Same g_max direction throughout\n'
            'Cliff is fixed — cloud moves\n'
            'Post-selection cloud: narrowed,\nre-aligned, new channel open',
            transform=ax.transAxes, fontsize=7.5, color=CGRAY,
            va='top', linespacing=1.5)
    panel_label(ax, 'B')

    fig.tight_layout(w_pad=3.5)
    fig.savefig('fig6_switching.pdf')
    print('Saved  fig6_switching.pdf')
    plt.close(fig)


# =============================================================================
# Entry point
# =============================================================================

if __name__ == '__main__':
    print('Generating figures for "A Letter to the ARC Panellist"…\n')
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    make_fig5()
    make_fig6()
    print('\nDone.  Six PDF files written to the working directory.')
