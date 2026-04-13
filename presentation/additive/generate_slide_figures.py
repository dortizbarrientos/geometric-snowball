#!/usr/bin/env python3
"""
================================================================================
PRESENTATION FIGURES: ADDITIVE CHANNELS IN CURVED FITNESS LANDSCAPES
================================================================================

Generates publication-quality figures for a 30-slide presentation.

Structure:
- Part I (Slides 1-6): The Paradox
- Part II (Slides 7-12): The Framework  
- Part III (Slides 13-18): Additive Channels
- Part IV (Slides 19-22): Escape from Additivity
- Part V (Slides 23-27): Applications
- Part VI (Slides 28-30): Synthesis

Author: Daniel Ortiz-Barrientos & Mark Cooper
Date: December 2024
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle, Circle, Polygon
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal
from scipy import ndimage
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'figure.dpi': 150,
    'font.family': 'sans-serif',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Colour palette
C = {
    'blue': '#2E86AB',
    'green': '#27AE60', 
    'red': '#E74C3C',
    'amber': '#F39C12',
    'purple': '#9B59B6',
    'dark': '#2C3E50',
    'light': '#ECF0F1',
    'teal': '#1ABC9C',
    'orange': '#E67E22',
}

# Output directory
import os
os.makedirs('slide_figures', exist_ok=True)

def save_fig(fig, name):
    """Save figure with consistent settings."""
    path = f'slide_figures/{name}.png'
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"  Saved: {path}")


# =============================================================================
# PART I: THE PARADOX (Slides 1-6)
# =============================================================================

def slide01_title():
    """Hero image: curved landscape with additive channel highlighted."""
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create landscape
    x = np.linspace(-3, 3, 100)
    y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    Z = -0.15 * (X**2 + Y**2) + 0.05 * X * Y
    
    # Plot surface
    ax.plot_surface(X, Y, Z, cmap='Blues_r', alpha=0.7, linewidth=0, antialiased=True)
    
    # Add "channel" as a highlighted tube
    theta = np.linspace(-0.5, 2, 50)
    channel_x = 1.5 * np.cos(theta) - 0.5
    channel_y = 1.5 * np.sin(theta) + 0.5
    channel_z = -0.15 * (channel_x**2 + channel_y**2) + 0.05 * channel_x * channel_y
    ax.plot(channel_x, channel_y, channel_z + 0.05, color=C['green'], lw=6, alpha=0.8)
    
    # Population ellipse
    t = np.linspace(0, 2*np.pi, 50)
    pop_x = 0.8 + 0.3 * np.cos(t)
    pop_y = 1.0 + 0.2 * np.sin(t)
    pop_z = -0.15 * (pop_x**2 + pop_y**2) + 0.05 * pop_x * pop_y
    ax.plot(pop_x, pop_y, pop_z + 0.02, color=C['amber'], lw=4)
    
    ax.set_xlabel('Genotype Axis 1', fontsize=14)
    ax.set_ylabel('Genotype Axis 2', fontsize=14)
    ax.set_zlabel('Fitness', fontsize=14)
    ax.set_title('Additive Channels in Curved Fitness Landscapes', fontsize=20, fontweight='bold', pad=20)
    ax.view_init(elev=25, azim=-50)
    ax.set_box_aspect([1, 1, 0.5])
    
    # Add text annotation
    ax.text2D(0.02, 0.95, 'The landscape is curved...\n...but the population behaves additively', 
              transform=ax.transAxes, fontsize=14, style='italic',
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    return fig


def slide02_empirical_success():
    """Predicted vs realised selection response across studies."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    np.random.seed(42)
    n = 50
    
    # Simulated data from different "studies"
    studies = ['Maize yield', 'Cattle milk', 'Wheat height', 'Drosophila bristles', 'Mouse weight']
    colors = [C['blue'], C['green'], C['amber'], C['purple'], C['teal']]
    
    all_pred, all_obs = [], []
    for i, (study, col) in enumerate(zip(studies, colors)):
        predicted = np.random.uniform(0.5, 3, 10) + i * 0.3
        noise = np.random.normal(0, 0.15, 10)
        observed = predicted * (0.95 + np.random.normal(0, 0.05, 10)) + noise
        
        ax.scatter(predicted, observed, c=col, s=100, alpha=0.7, label=study, edgecolor='white', lw=1)
        all_pred.extend(predicted)
        all_obs.extend(observed)
    
    # Add 1:1 line
    lims = [0, 5]
    ax.plot(lims, lims, 'k--', lw=2, label='Perfect prediction')
    
    # Regression line
    all_pred, all_obs = np.array(all_pred), np.array(all_obs)
    m, b = np.polyfit(all_pred, all_obs, 1)
    ax.plot(lims, [m*lims[0]+b, m*lims[1]+b], color=C['red'], lw=2, ls='-', 
            label=f'Observed (r² = {np.corrcoef(all_pred, all_obs)[0,1]**2:.2f})')
    
    ax.set_xlabel('Predicted Response (from additive model)', fontsize=14)
    ax.set_ylabel('Observed Response', fontsize=14)
    ax.set_title('The Empirical Success of Additive Models', fontsize=18, fontweight='bold')
    ax.legend(loc='upper left', fontsize=11)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')
    
    ax.text(0.95, 0.05, 'Additive models predict\nselection response remarkably well', 
            transform=ax.transAxes, ha='right', fontsize=12, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['light'], alpha=0.9))
    
    return fig


def slide03_molecular_reality():
    """Evidence of epistasis and non-additivity from molecular studies."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Panel A: Diminishing returns epistasis
    ax = axes[0]
    backgrounds = [0, 1, 2, 3, 4]
    effect_sizes = [1.0, 0.7, 0.5, 0.35, 0.25]
    
    ax.bar(backgrounds, effect_sizes, color=C['blue'], edgecolor='white', lw=2)
    ax.set_xlabel('Number of Prior Beneficial Mutations', fontsize=12)
    ax.set_ylabel('Effect Size of New Mutation', fontsize=12)
    ax.set_title('(A) Diminishing Returns', fontsize=14, fontweight='bold')
    ax.set_xticks(backgrounds)
    
    # Panel B: G×E interaction (reaction norms crossing)
    ax = axes[1]
    envs = [0, 1]
    genotypes = ['Genotype A', 'Genotype B', 'Genotype C']
    colors_gxe = [C['blue'], C['green'], C['amber']]
    
    for i, (geno, col) in enumerate(zip(genotypes, colors_gxe)):
        if i == 0:
            y = [2, 4]
        elif i == 1:
            y = [3.5, 2.5]
        else:
            y = [1.5, 3]
        ax.plot(envs, y, 'o-', color=col, lw=3, ms=12, label=geno)
    
    ax.set_xlabel('Environment', fontsize=12)
    ax.set_ylabel('Fitness', fontsize=12)
    ax.set_title('(B) G×E: Rankings Change', fontsize=14, fontweight='bold')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Env 1', 'Env 2'])
    ax.legend(loc='upper right', fontsize=10)
    
    # Panel C: Sign epistasis heatmap
    ax = axes[2]
    np.random.seed(123)
    # Create epistasis matrix
    epistasis = np.array([
        [0, 0.5, -0.3, 0.2],
        [0.5, 0, 0.1, -0.6],
        [-0.3, 0.1, 0, 0.4],
        [0.2, -0.6, 0.4, 0]
    ])
    
    im = ax.imshow(epistasis, cmap='RdBu_r', vmin=-0.7, vmax=0.7)
    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels(['A', 'B', 'C', 'D'])
    ax.set_yticklabels(['A', 'B', 'C', 'D'])
    ax.set_xlabel('Locus', fontsize=12)
    ax.set_ylabel('Locus', fontsize=12)
    ax.set_title('(C) Pairwise Epistasis', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Epistatic effect')
    
    fig.suptitle('The Molecular Reality: Genomes Are Non-Additive', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    return fig


def slide04_paradox():
    """Visual representation of the paradox."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Additive success
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Green checkmarks
    checks = [(2, 8), (2, 6), (2, 4), (2, 2)]
    texts = ['Breeder\'s equation works', 'GBLUP predicts accurately', 
             'Heritabilities are stable', 'Fisher\'s theorem holds']
    
    for (x, y), txt in zip(checks, texts):
        ax.plot(x, y, 'o', ms=30, color=C['green'])
        ax.text(x, y, '✓', fontsize=20, ha='center', va='center', color='white', fontweight='bold')
        ax.text(x + 0.8, y, txt, fontsize=14, va='center')
    
    ax.set_title('What We Observe', fontsize=16, fontweight='bold', color=C['green'])
    ax.add_patch(Rectangle((0.5, 0.5), 9, 9, fill=False, edgecolor=C['green'], lw=3))
    
    # Right: Molecular complexity
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Red X marks
    xmarks = [(2, 8), (2, 6), (2, 4), (2, 2)]
    texts = ['Epistasis is pervasive', 'G×E is ubiquitous', 
             'Diminishing returns observed', 'Sign epistasis documented']
    
    for (x, y), txt in zip(xmarks, texts):
        ax.plot(x, y, 'o', ms=30, color=C['red'])
        ax.text(x, y, '✗', fontsize=20, ha='center', va='center', color='white', fontweight='bold')
        ax.text(x + 0.8, y, txt, fontsize=14, va='center')
    
    ax.set_title('What Molecular Biology Shows', fontsize=16, fontweight='bold', color=C['red'])
    ax.add_patch(Rectangle((0.5, 0.5), 9, 9, fill=False, edgecolor=C['red'], lw=3))
    
    # Central question
    fig.text(0.5, 0.02, 'How can both be true?', fontsize=20, ha='center', 
             fontweight='bold', style='italic',
             bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.3))
    
    plt.tight_layout()
    return fig


def slide05_preview():
    """Preview: curved surface with small population in linear region."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: 2D view with tangent
    ax = axes[0]
    x = np.linspace(-3, 3, 200)
    y_true = -0.2 * x**2 + 0.5 * x
    y_tangent = 0.5 * x  # Tangent at x=0
    
    ax.plot(x, y_true, color=C['blue'], lw=3, label='True landscape (curved)')
    ax.plot(x, y_tangent, '--', color=C['amber'], lw=2, label='Linear approximation')
    
    # Shaded region where they match
    ax.fill_between(x, -1.5, 1.5, where=(np.abs(x) < 1), alpha=0.2, color=C['green'],
                    label='Additive channel')
    
    # Population
    pop_x = np.linspace(-0.5, 0.5, 20)
    pop_y = -0.2 * pop_x**2 + 0.5 * pop_x
    ax.scatter(pop_x, pop_y, c=C['green'], s=50, zorder=5, edgecolor='white')
    
    ax.set_xlabel('Genotype', fontsize=14)
    ax.set_ylabel('Fitness', fontsize=14)
    ax.set_title('Small Population = Local Linearity', fontsize=16, fontweight='bold')
    ax.legend(loc='lower left', fontsize=11)
    ax.set_ylim(-1.5, 1.5)
    
    # Right: Key insight text
    ax = axes[1]
    ax.axis('off')
    
    insights = [
        ("1", "Any smooth function is\nlocally linear", C['blue']),
        ("2", "Small variance means\nsmall deviations", C['green']),
        ("3", "Curvature contributes\nlittle to variance", C['amber']),
        ("4", "Population behaves\nas if additive", C['purple']),
    ]
    
    for i, (num, text, col) in enumerate(insights):
        y = 0.85 - i * 0.22
        ax.add_patch(Circle((0.08, y), 0.04, color=col, transform=ax.transAxes))
        ax.text(0.08, y, num, fontsize=14, ha='center', va='center', 
                color='white', fontweight='bold', transform=ax.transAxes)
        ax.text(0.18, y, text, fontsize=14, va='center', transform=ax.transAxes)
    
    ax.set_title('The Key Insight', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    return fig


def slide06_outline():
    """Talk outline as visual roadmap."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    sections = [
        ('I', 'The Paradox', 'Additive success vs molecular complexity', C['red']),
        ('II', 'The Framework', 'Taylor expansion & variance decomposition', C['blue']),
        ('III', 'Additive Channels', 'Where Fisher\'s theorem holds', C['green']),
        ('IV', 'Escape', 'How populations leave channels', C['amber']),
        ('V', 'Applications', 'Breeding & evolutionary implications', C['purple']),
        ('VI', 'Synthesis', 'Resolving the paradox', C['teal']),
    ]
    
    for i, (num, title, subtitle, col) in enumerate(sections):
        y = 0.88 - i * 0.14
        
        # Circle with number
        ax.add_patch(Circle((0.08, y), 0.035, color=col, transform=ax.transAxes))
        ax.text(0.08, y, num, fontsize=16, ha='center', va='center',
                color='white', fontweight='bold', transform=ax.transAxes)
        
        # Title and subtitle
        ax.text(0.15, y + 0.015, title, fontsize=18, fontweight='bold', 
                va='center', transform=ax.transAxes)
        ax.text(0.15, y - 0.025, subtitle, fontsize=12, va='center',
                color='gray', transform=ax.transAxes)
        
        # Connecting line (except last)
        if i < len(sections) - 1:
            ax.plot([0.08, 0.08], [y - 0.05, y - 0.09], '-', color='gray', 
                    lw=2, transform=ax.transAxes)
    
    ax.set_title('Talk Outline', fontsize=22, fontweight='bold', y=0.98)
    
    return fig


# =============================================================================
# PART II: THE FRAMEWORK (Slides 7-12)
# =============================================================================

def slide07_gpf_map():
    """Diagram of genotype → phenotype → fitness mapping."""
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.axis('off')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    
    # Three boxes
    boxes = [
        (1.5, 3, 'Genotype\n$\\mathbf{g}$', C['blue']),
        (5, 3, 'Phenotype\n$\\mathbf{z}$', C['green']),
        (8.5, 3, 'Fitness\n$W$', C['amber']),
    ]
    
    for x, y, text, col in boxes:
        ax.add_patch(Rectangle((x-0.8, y-0.8), 1.6, 1.6, 
                                facecolor=col, alpha=0.3, edgecolor=col, lw=3))
        ax.text(x, y, text, fontsize=16, ha='center', va='center', fontweight='bold')
    
    # Arrows
    ax.annotate('', xy=(4, 3), xytext=(2.5, 3),
                arrowprops=dict(arrowstyle='->', lw=3, color=C['dark']))
    ax.text(3.25, 3.5, '$\\phi$', fontsize=18, ha='center')
    ax.text(3.25, 2.4, 'Development', fontsize=12, ha='center', style='italic')
    
    ax.annotate('', xy=(7.5, 3), xytext=(6, 3),
                arrowprops=dict(arrowstyle='->', lw=3, color=C['dark']))
    ax.text(6.75, 3.5, '$\\psi$', fontsize=18, ha='center')
    ax.text(6.75, 2.4, 'Selection', fontsize=12, ha='center', style='italic')
    
    # Environment input
    ax.add_patch(Rectangle((4.2, 0.5), 1.6, 0.8, 
                            facecolor=C['purple'], alpha=0.3, edgecolor=C['purple'], lw=2))
    ax.text(5, 0.9, 'Environment $\\mathbf{e}$', fontsize=14, ha='center', va='center')
    
    ax.annotate('', xy=(3.5, 2.5), xytext=(4.5, 1.4),
                arrowprops=dict(arrowstyle='->', lw=2, color=C['purple'], ls='--'))
    ax.annotate('', xy=(6.5, 2.5), xytext=(5.5, 1.4),
                arrowprops=dict(arrowstyle='->', lw=2, color=C['purple'], ls='--'))
    
    # Combined map at bottom
    ax.text(5, 5.2, 'Combined map: $F(\\mathbf{g}, \\mathbf{e}) = \\psi(\\phi(\\mathbf{g}, \\mathbf{e}), \\mathbf{e})$',
            fontsize=16, ha='center', 
            bbox=dict(boxstyle='round', facecolor='white', edgecolor=C['dark']))
    
    ax.set_title('The Genotype–Phenotype–Fitness Map', fontsize=20, fontweight='bold')
    
    return fig


def slide08_taylor():
    """Taylor expansion visualized."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: 1D example
    ax = axes[0]
    x = np.linspace(-2, 2, 200)
    
    # True function
    f_true = np.exp(-0.3 * x**2) * (1 + 0.5 * x)
    
    # Taylor approximations
    f0 = 1  # f(0)
    f1 = 1 + 0.5 * x  # Linear
    f2 = 1 + 0.5 * x - 0.3 * x**2  # Quadratic
    
    ax.plot(x, f_true, color=C['dark'], lw=3, label='True function')
    ax.plot(x, np.ones_like(x) * f0, ':', color='gray', lw=2, label='Constant ($F_0$)')
    ax.plot(x, f1, '--', color=C['blue'], lw=2, label='Linear ($F_0 + \\mathbf{b}^\\top \\Delta\\mathbf{g}$)')
    ax.plot(x, f2, '-', color=C['green'], lw=2, label='Quadratic (+ Hessian term)')
    
    ax.axvline(0, ls=':', color='gray', alpha=0.5)
    ax.scatter([0], [1], s=100, color=C['amber'], zorder=5, edgecolor='white', lw=2)
    ax.text(0.1, 1.1, 'Reference point', fontsize=11)
    
    ax.set_xlabel('Deviation from reference ($\\Delta g$)', fontsize=14)
    ax.set_ylabel('Fitness', fontsize=14)
    ax.set_title('Taylor Expansion: Local Polynomial', fontsize=16, fontweight='bold')
    ax.legend(loc='lower left', fontsize=10)
    ax.set_xlim(-2, 2)
    
    # Right: The formula
    ax = axes[1]
    ax.axis('off')
    
    ax.text(0.5, 0.85, 'The Second-Order Expansion', fontsize=18, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.68, 
            '$F(\\mathbf{g}) \\approx F_0 + \\mathbf{b}^\\top \\Delta\\mathbf{g} + \\frac{1}{2} \\Delta\\mathbf{g}^\\top \\mathbf{H} \\Delta\\mathbf{g}$',
            fontsize=18, ha='center', va='center', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    # Labels below formula
    ax.text(0.23, 0.55, 'baseline', fontsize=12, ha='center', transform=ax.transAxes, color=C['dark'])
    ax.text(0.43, 0.55, 'linear\n(additive)', fontsize=12, ha='center', transform=ax.transAxes, color=C['blue'])
    ax.text(0.72, 0.55, 'quadratic\n(epistasis)', fontsize=12, ha='center', transform=ax.transAxes, color=C['red'])
    
    ax.text(0.5, 0.25, 
            '$\\mathbf{b} = \\nabla F$ is the gradient (additive effects)\n\n' +
            '$\\mathbf{H}$ is the Hessian (curvature/epistasis)\n\n' +
            'Close to reference: linear dominates\n' +
            'Far away: curvature matters',
            fontsize=13, ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    return fig


def slide09_qg_connection():
    """Table connecting Taylor terms to QG quantities."""
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.axis('off')
    
    # Table data
    headers = ['Taylor Term', 'Formula', 'QG Interpretation']
    rows = [
        ['Gradient', '$\\mathbf{b}_g = J_g^\\top \\mathbf{a}_z$', 'Additive genetic effects'],
        ['Genetic Hessian', '$\\mathbf{H}_{gg} = J_g^\\top H_{zz} J_g$', 'Epistasis'],
        ['Mixed Hessian', '$\\mathbf{H}_{ge}$', 'G×E interaction'],
        ['Linear variance', '$V_{\\text{lin}} = \\mathbf{b}^\\top \\Sigma_g \\mathbf{b}$', 'Additive genetic variance'],
        ['Quadratic variance', '$V_{\\text{quad}} = \\frac{1}{2}\\text{tr}[(\\mathbf{H}\\Sigma_g)^2]$', 'Epistatic variance'],
    ]
    
    # Draw table
    cell_height = 0.1
    col_widths = [0.25, 0.4, 0.35]
    start_y = 0.75
    start_x = 0.05
    
    # Header
    x = start_x
    for i, (header, width) in enumerate(zip(headers, col_widths)):
        ax.add_patch(Rectangle((x, start_y), width - 0.01, cell_height,
                                facecolor=C['blue'], alpha=0.8))
        ax.text(x + width/2, start_y + cell_height/2, header,
                fontsize=14, ha='center', va='center', color='white', fontweight='bold',
                transform=ax.transAxes)
        x += width
    
    # Rows
    for j, row in enumerate(rows):
        y = start_y - (j + 1) * cell_height
        x = start_x
        bg_color = C['light'] if j % 2 == 0 else 'white'
        
        for i, (cell, width) in enumerate(zip(row, col_widths)):
            ax.add_patch(Rectangle((x, y), width - 0.01, cell_height,
                                    facecolor=bg_color, edgecolor='gray', lw=0.5))
            ax.text(x + width/2, y + cell_height/2, cell,
                    fontsize=12, ha='center', va='center',
                    transform=ax.transAxes)
            x += width
    
    ax.set_title('Classical QG Emerges from the Taylor Expansion', 
                 fontsize=18, fontweight='bold', y=0.95)
    
    ax.text(0.5, 0.15, 
            'The breeder\'s equation $\\Delta\\bar{z} = \\mathbf{G}\\boldsymbol{\\beta}$ is the linear term!',
            fontsize=14, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['green'], alpha=0.2))
    
    return fig


def slide10_variance_scaling():
    """Stacked bar chart showing V_lin and V_quad at different scales."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Stacked bars
    ax = axes[0]
    
    scales = ['Small σ', 'Medium σ', 'Large σ']
    sigma = np.array([0.3, 0.7, 1.2])
    
    V_lin = sigma**2
    V_quad = 0.5 * sigma**4
    
    x = np.arange(len(scales))
    width = 0.6
    
    ax.bar(x, V_lin, width, label='$V_{\\text{lin}} \\propto \\sigma^2$', color=C['green'], edgecolor='white', lw=2)
    ax.bar(x, V_quad, width, bottom=V_lin, label='$V_{\\text{quad}} \\propto \\sigma^4$', color=C['red'], edgecolor='white', lw=2)
    
    ax.set_xticks(x)
    ax.set_xticklabels(scales, fontsize=12)
    ax.set_ylabel('Variance in Fitness', fontsize=14)
    ax.set_title('Variance Decomposition', fontsize=16, fontweight='bold')
    ax.legend(loc='upper left', fontsize=12)
    
    # Annotate proportions
    for i, (vl, vq) in enumerate(zip(V_lin, V_quad)):
        total = vl + vq
        pct = vl / total * 100
        ax.text(i, total + 0.05, f'{pct:.0f}% linear', ha='center', fontsize=11)
    
    # Right: Continuous curve
    ax = axes[1]
    sigma_cont = np.linspace(0.1, 1.5, 100)
    V_lin_cont = sigma_cont**2
    V_quad_cont = 0.5 * sigma_cont**4
    
    ax.fill_between(sigma_cont, 0, V_lin_cont, alpha=0.5, color=C['green'], label='$V_{\\text{lin}}$')
    ax.fill_between(sigma_cont, V_lin_cont, V_lin_cont + V_quad_cont, alpha=0.5, color=C['red'], label='$V_{\\text{quad}}$')
    
    ax.set_xlabel('Genetic Standard Deviation ($\\sigma$)', fontsize=14)
    ax.set_ylabel('Variance in Fitness', fontsize=14)
    ax.set_title('The Scaling Insight', fontsize=16, fontweight='bold')
    ax.legend(loc='upper left', fontsize=12)
    
    # Key insight box
    ax.text(0.95, 0.4, '$\\sigma^4$ grows\nfaster than $\\sigma^2$\n\nSmall variance\n→ linear dominates',
            transform=ax.transAxes, ha='right', fontsize=12,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    plt.tight_layout()
    return fig


def slide11_additivity_index():
    """Contour plot of A_g over genotype space."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: 2D contour plot
    ax = axes[0]
    
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    
    # A_g decreases with distance from origin
    R2 = X**2 + Y**2
    Ag = 1 / (1 + 0.3 * R2)
    
    levels = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    cs = ax.contourf(X, Y, Ag, levels=20, cmap='RdYlGn')
    ax.contour(X, Y, Ag, levels=levels, colors='black', linewidths=1)
    
    plt.colorbar(cs, ax=ax, label='$\\mathcal{A}_g$')
    
    # Mark center
    ax.plot(0, 0, 'ko', ms=10)
    ax.text(0.1, 0.2, 'Reference', fontsize=11)
    
    ax.set_xlabel('Genetic Coordinate 1', fontsize=14)
    ax.set_ylabel('Genetic Coordinate 2', fontsize=14)
    ax.set_title('Additivity Index Over Genotype Space', fontsize=16, fontweight='bold')
    
    # Right: Formula and interpretation
    ax = axes[1]
    ax.axis('off')
    
    ax.text(0.5, 0.85, 'The Additivity Index', fontsize=18, fontweight='bold',
            ha='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.65, 
            '$\\mathcal{A}_g = \\dfrac{V_{\\text{lin}}}{V_{\\text{lin}} + V_{\\text{quad}}}$',
            fontsize=24, ha='center', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    interpretations = [
        ('$\\mathcal{A}_g \\approx 1$', 'Linear dominates; additive models work', C['green']),
        ('$\\mathcal{A}_g \\approx 0.5$', 'Both contribute; caution needed', C['amber']),
        ('$\\mathcal{A}_g \\approx 0$', 'Curvature dominates; additive fails', C['red']),
    ]
    
    for i, (val, text, col) in enumerate(interpretations):
        y = 0.38 - i * 0.12
        ax.add_patch(Rectangle((0.1, y - 0.04), 0.15, 0.08, 
                                facecolor=col, alpha=0.5, transform=ax.transAxes))
        ax.text(0.175, y, val, fontsize=12, ha='center', va='center', transform=ax.transAxes)
        ax.text(0.3, y, text, fontsize=12, va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    return fig


def slide12_Ag_as_Rsquared():
    """Scatter plots showing A_g = R² interpretation."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    np.random.seed(42)
    n = 150
    
    # High A_g
    ax = axes[0]
    L = np.random.normal(0, 1, n)
    Q = np.random.normal(0, 0.3, n)
    total = L + Q
    
    ax.scatter(L, total, alpha=0.5, c=C['green'], s=40, edgecolor='white')
    ax.plot([-3, 3], [-3, 3], 'k--', lw=2)
    
    r2 = np.corrcoef(L, total)[0, 1]**2
    ax.text(0.05, 0.95, f'$R^2 = {r2:.2f}$\n$\\mathcal{{A}}_g = {r2:.2f}$',
            transform=ax.transAxes, fontsize=14, va='top',
            bbox=dict(boxstyle='round', facecolor='white'))
    
    ax.set_xlabel('Linear Predictor ($L$)', fontsize=12)
    ax.set_ylabel('True Log-Fitness ($\\tilde{\\ell}$)', fontsize=12)
    ax.set_title('High $\\mathcal{A}_g$: Linear Predicts Well', fontsize=14, fontweight='bold')
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-3.5, 3.5)
    
    # Low A_g
    ax = axes[1]
    L = np.random.normal(0, 1, n)
    Q = np.random.normal(0, 1.5, n)
    total = L + Q
    
    ax.scatter(L, total, alpha=0.5, c=C['red'], s=40, edgecolor='white')
    ax.plot([-3, 3], [-3, 3], 'k--', lw=2)
    
    r2 = np.corrcoef(L, total)[0, 1]**2
    ax.text(0.05, 0.95, f'$R^2 = {r2:.2f}$\n$\\mathcal{{A}}_g = {r2:.2f}$',
            transform=ax.transAxes, fontsize=14, va='top',
            bbox=dict(boxstyle='round', facecolor='white'))
    
    ax.set_xlabel('Linear Predictor ($L$)', fontsize=12)
    ax.set_ylabel('True Log-Fitness ($\\tilde{\\ell}$)', fontsize=12)
    ax.set_title('Low $\\mathcal{A}_g$: Curvature Dominates', fontsize=14, fontweight='bold')
    ax.set_xlim(-4, 4)
    ax.set_ylim(-5, 5)
    
    # Geometry diagram
    ax = axes[2]
    ax.set_aspect('equal')
    ax.set_xlim(-0.5, 4)
    ax.set_ylim(-0.5, 2.5)
    
    # Vectors
    ax.annotate('', xy=(3, 1.2), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', lw=3, color=C['blue']))
    ax.text(1.7, 0.9, '$\\tilde{\\ell}$', fontsize=16, color=C['blue'])
    
    ax.annotate('', xy=(3, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', lw=3, color=C['green']))
    ax.text(1.5, -0.3, '$L$', fontsize=16, color=C['green'])
    
    ax.plot([3, 3], [0, 1.2], '--', color=C['red'], lw=2)
    ax.text(3.2, 0.6, '$Q$', fontsize=16, color=C['red'])
    
    # Right angle
    ax.plot([2.8, 2.8, 3], [0, 0.2, 0.2], color='gray', lw=1)
    
    # Angle arc
    theta = np.linspace(0, np.arctan(1.2/3), 20)
    r = 0.7
    ax.plot(r * np.cos(theta), r * np.sin(theta), color=C['dark'], lw=2)
    ax.text(0.85, 0.2, '$\\theta$', fontsize=14)
    
    ax.text(1.5, 2, '$\\mathcal{A}_g = \\cos^2\\theta = R^2$', fontsize=16,
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    ax.axis('off')
    ax.set_title('Geometric Interpretation', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    return fig


# =============================================================================
# PART III: ADDITIVE CHANNELS (Slides 13-18)
# =============================================================================

def slide13_channels_defined():
    """Fitness landscape with additive channel tubes highlighted."""
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create landscape
    x = np.linspace(-3, 3, 80)
    y = np.linspace(-3, 3, 80)
    X, Y = np.meshgrid(x, y)
    Z = -0.12 * (X**2 + Y**2) + 0.03 * X * Y + 0.02 * np.sin(X) * np.cos(Y)
    
    # Plot surface with transparency
    ax.plot_surface(X, Y, Z, cmap='Blues_r', alpha=0.5, linewidth=0, antialiased=True)
    
    # Add channel regions (tubes where A_g > threshold)
    # Channel 1: along x axis
    channel_x = np.linspace(-2.5, 2.5, 50)
    channel_y = np.zeros_like(channel_x)
    for offset in [-0.4, -0.2, 0, 0.2, 0.4]:
        cy = channel_y + offset
        cz = -0.12 * (channel_x**2 + cy**2) + 0.02
        ax.plot(channel_x, cy, cz, color=C['green'], alpha=0.3, lw=1)
    
    # Main channel line
    cz = -0.12 * channel_x**2 + 0.02
    ax.plot(channel_x, channel_y, cz + 0.05, color=C['green'], lw=4, label='Additive channel')
    
    # Population in channel
    t = np.linspace(0, 2*np.pi, 30)
    pop_x = 0.5 + 0.3 * np.cos(t)
    pop_y = 0 + 0.15 * np.sin(t)
    pop_z = -0.12 * (pop_x**2 + pop_y**2) + 0.02
    ax.plot(pop_x, pop_y, pop_z + 0.02, color=C['amber'], lw=3, label='Population')
    
    ax.set_xlabel('Genotype 1', fontsize=12)
    ax.set_ylabel('Genotype 2', fontsize=12)
    ax.set_zlabel('Fitness', fontsize=12)
    ax.set_title('Additive Channels: Where Linear Models Work', fontsize=16, fontweight='bold')
    ax.view_init(elev=30, azim=-60)
    ax.legend(loc='upper left')
    
    return fig


def slide14_fisher():
    """Fisher's theorem as a local property."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Conceptual diagram
    ax = axes[0]
    ax.axis('off')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    
    # Fisher's theorem box
    ax.add_patch(Rectangle((1, 5.5), 8, 2, facecolor=C['blue'], alpha=0.2, 
                            edgecolor=C['blue'], lw=3))
    ax.text(5, 6.5, "Fisher's Fundamental Theorem", fontsize=16, 
            fontweight='bold', ha='center')
    ax.text(5, 5.9, '$\\dfrac{d\\bar{W}}{dt} = V_A(W)$', fontsize=18, ha='center')
    
    # Arrow down
    ax.annotate('', xy=(5, 4.8), xytext=(5, 5.3),
                arrowprops=dict(arrowstyle='->', lw=3, color=C['dark']))
    ax.text(5.2, 5, 'holds when', fontsize=12)
    
    # Condition box
    ax.add_patch(Rectangle((1, 2.5), 8, 2, facecolor=C['green'], alpha=0.2,
                            edgecolor=C['green'], lw=3))
    ax.text(5, 3.7, 'Population is in an Additive Channel', fontsize=14, 
            fontweight='bold', ha='center')
    ax.text(5, 3.0, '$\\mathcal{A}_g \\approx 1$ (curvature negligible)', 
            fontsize=14, ha='center')
    
    # Implication
    ax.text(5, 1.5, "Fisher's theorem is not a universal law;\nit's a local property of high-$\\mathcal{A}_g$ regions",
            fontsize=13, ha='center', style='italic',
            bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.2))
    
    # Right: Trajectory where Fisher holds
    ax = axes[1]
    generations = np.arange(0, 20)
    
    # Mean fitness
    W = 1 + 0.8 * (1 - np.exp(-0.15 * generations))
    
    # V_A (declining due to selection)
    VA = 0.15 * np.exp(-0.1 * generations)
    
    # Predicted dW/dt from Fisher
    dW_pred = VA
    
    # Observed dW/dt
    dW_obs = np.diff(W)
    dW_obs = np.append(dW_obs, dW_obs[-1])
    
    ax.plot(generations, dW_obs, 'o-', color=C['blue'], lw=2, ms=6, 
            label='Observed $\\Delta\\bar{W}$')
    ax.plot(generations, dW_pred, '--', color=C['green'], lw=2,
            label='Predicted $V_A(W)$')
    
    ax.set_xlabel('Generation', fontsize=14)
    ax.set_ylabel('Rate of Fitness Change', fontsize=14)
    ax.set_title('Fisher\'s Prediction Tracks Observation\n(within additive channel)', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=11)
    
    plt.tight_layout()
    return fig


def slide15_heritability():
    """Directional heritability and response to selection."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: G matrix ellipse with selection gradient
    ax = axes[0]
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_aspect('equal')
    
    # G matrix ellipse
    ell = Ellipse((0, 0), 4, 1.5, angle=30, facecolor=C['blue'], alpha=0.3,
                  edgecolor=C['blue'], lw=2)
    ax.add_patch(ell)
    
    # Selection gradient
    ax.annotate('', xy=(1.5, 2), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', lw=3, color=C['red']))
    ax.text(1.6, 2.1, '$\\boldsymbol{\\beta}$', fontsize=16, color=C['red'])
    
    # Response
    ax.annotate('', xy=(1.2, 0.8), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', lw=3, color=C['green']))
    ax.text(1.3, 0.6, '$\\mathbf{G}\\boldsymbol{\\beta}$', fontsize=16, color=C['green'])
    
    # Eigenvectors
    for angle, length in [(30, 2), (120, 0.75)]:
        rad = np.radians(angle)
        ax.plot([0, length * np.cos(rad)], [0, length * np.sin(rad)], 
                '--', color='gray', lw=1.5)
    
    ax.set_xlabel('Trait 1', fontsize=14)
    ax.set_ylabel('Trait 2', fontsize=14)
    ax.set_title('Response = $\\mathbf{G}\\boldsymbol{\\beta}$\n(The Breeder\'s Equation)', 
                 fontsize=14, fontweight='bold')
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    
    # Right: What happens when A_g drops
    ax = axes[1]
    ax.axis('off')
    
    ax.text(0.5, 0.9, 'When $\\mathcal{A}_g$ Drops...', fontsize=16, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    problems = [
        '• Observed response deviates from $\\mathbf{G}\\boldsymbol{\\beta}$',
        '• Apparent changes in $\\mathbf{G}$ over short timescales',
        '• Fitness effects depend on genetic background',
        '• Rankings change unexpectedly across environments',
        '• Directional heritability becomes unreliable',
    ]
    
    for i, prob in enumerate(problems):
        ax.text(0.1, 0.7 - i * 0.12, prob, fontsize=13, 
                transform=ax.transAxes, color=C['red'])
    
    ax.text(0.5, 0.15, 'These are signatures of\nleaving an additive channel',
            fontsize=14, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.3))
    
    plt.tight_layout()
    return fig


def slide16_1d_condition():
    """The explicit 1D additivity condition."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: A_g as function of V_A
    ax = axes[0]
    
    a = 1.0  # Directional selection
    b = 0.5  # Curvature
    
    VA = np.linspace(0.01, 5, 100)
    Ag = 1 / (1 + (b**2 * VA) / (2 * a**2))
    
    ax.plot(VA, Ag, color=C['blue'], lw=3)
    ax.axhline(0.9, ls='--', color=C['green'], lw=2, label='$\\mathcal{A}_0 = 0.9$ threshold')
    ax.axhline(0.8, ls='--', color=C['amber'], lw=2, label='$\\mathcal{A}_0 = 0.8$ threshold')
    
    # Mark the critical V_A values
    VA_crit_90 = 2 * a**2 / b**2 * (1 - 0.9) / 0.9
    VA_crit_80 = 2 * a**2 / b**2 * (1 - 0.8) / 0.8
    
    ax.axvline(VA_crit_90, ls=':', color=C['green'], alpha=0.7)
    ax.axvline(VA_crit_80, ls=':', color=C['amber'], alpha=0.7)
    
    ax.fill_between(VA, 0, 1, where=(VA < VA_crit_90), alpha=0.2, color=C['green'])
    ax.text(VA_crit_90/2, 0.5, 'Channel\n($\\mathcal{A}_g > 0.9$)', 
            ha='center', fontsize=11)
    
    ax.set_xlabel('Additive Genetic Variance ($V_A$)', fontsize=14)
    ax.set_ylabel('Additivity Index ($\\mathcal{A}_g$)', fontsize=14)
    ax.set_title('The Additivity–Variance Tradeoff', fontsize=16, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1.05)
    
    # Right: The formula
    ax = axes[1]
    ax.axis('off')
    
    ax.text(0.5, 0.85, 'The Explicit Condition', fontsize=18, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.65, 
            '$\\mathcal{A}_g = \\dfrac{1}{1 + \\dfrac{b^2 V_A}{2a^2}}$',
            fontsize=22, ha='center', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    ax.text(0.5, 0.42, 
            'Channel condition ($\\mathcal{A}_g \\geq \\mathcal{A}_0$):',
            fontsize=14, ha='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.28,
            '$V_A \\leq \\dfrac{2a^2}{b^2} \\cdot \\dfrac{1 - \\mathcal{A}_0}{\\mathcal{A}_0}$',
            fontsize=20, ha='center', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor=C['green'], alpha=0.2))
    
    ax.text(0.5, 0.08, 
            '$a$ = directional selection, $b$ = curvature, $V_A$ = genetic variance',
            fontsize=12, ha='center', transform=ax.transAxes, style='italic')
    
    plt.tight_layout()
    return fig


def slide17_interpretation():
    """Biological interpretation of factors affecting A_g."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    ax.text(0.5, 0.95, 'What Controls Additivity?', fontsize=20, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    factors = [
        ('↑ Curvature ($b$)', '↓ $\\mathcal{A}_g$', 
         'Strongly curved landscapes are harder to approximate linearly', C['red']),
        ('↑ Genetic Variance ($V_A$)', '↓ $\\mathcal{A}_g$', 
         'Wider sampling of landscape hits more curvature', C['red']),
        ('↑ Directional Selection ($a$)', '↑ $\\mathcal{A}_g$', 
         'Strong slope dominates over curvature', C['green']),
        ('Bottleneck/Founder Effect', '↑ $\\mathcal{A}_g$', 
         'Reduced variance keeps population in linear region', C['green']),
    ]
    
    for i, (factor, effect, explanation, col) in enumerate(factors):
        y = 0.78 - i * 0.18
        
        # Factor box
        ax.add_patch(Rectangle((0.05, y - 0.05), 0.3, 0.1, 
                                facecolor=col, alpha=0.3, transform=ax.transAxes))
        ax.text(0.2, y, factor, fontsize=13, ha='center', va='center',
                transform=ax.transAxes, fontweight='bold')
        
        # Arrow
        ax.annotate('', xy=(0.42, y), xytext=(0.38, y),
                    arrowprops=dict(arrowstyle='->', lw=2, color=C['dark']),
                    xycoords='axes fraction', textcoords='axes fraction')
        
        # Effect
        ax.text(0.48, y, effect, fontsize=14, va='center', transform=ax.transAxes,
                color=col, fontweight='bold')
        
        # Explanation
        ax.text(0.58, y, explanation, fontsize=12, va='center', transform=ax.transAxes)
    
    # Key insight
    ax.text(0.5, 0.08, 
            'Bottlenecked populations are more likely to behave "additively"\n' +
            'even on complex landscapes',
            fontsize=14, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.3))
    
    return fig


def slide18_self_correcting():
    """The self-correcting dynamic."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Left: Variance compression through selection
    ax = axes[0]
    
    cycles = np.arange(0, 12)
    sigma = 1.0 * np.exp(-0.1 * cycles)
    
    V_lin = sigma**2
    V_quad = 0.5 * sigma**4
    Ag = V_lin / (V_lin + V_quad)
    
    ax.fill_between(cycles, 0, V_lin, alpha=0.5, color=C['green'], label='$V_{\\text{lin}}$')
    ax.fill_between(cycles, V_lin, V_lin + V_quad, alpha=0.5, color=C['red'], label='$V_{\\text{quad}}$')
    
    ax.set_xlabel('Breeding Cycle', fontsize=14)
    ax.set_ylabel('Variance', fontsize=14)
    ax.set_title('Variance Compression', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    
    # Center: A_g rises
    ax = axes[1]
    ax.plot(cycles, Ag, 'o-', color=C['blue'], lw=3, ms=8)
    ax.axhline(0.8, ls='--', color='gray')
    ax.fill_between(cycles, 0.8, 1, alpha=0.2, color=C['green'])
    ax.text(6, 0.85, 'Additive Channel', fontsize=11)
    
    ax.set_xlabel('Breeding Cycle', fontsize=14)
    ax.set_ylabel('Additivity Index ($\\mathcal{A}_g$)', fontsize=14)
    ax.set_title('$\\mathcal{A}_g$ Rises', fontsize=14, fontweight='bold')
    ax.set_ylim(0.4, 1.02)
    
    # Right: The mechanism
    ax = axes[2]
    ax.axis('off')
    
    ax.text(0.5, 0.95, 'The Self-Correcting Dynamic', fontsize=16, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    # Flow diagram
    boxes = [
        (0.5, 0.75, 'Selection acts', C['amber']),
        (0.5, 0.55, 'Variance compresses', 'white'),
        (0.5, 0.35, '$V_{\\text{quad}} \\propto \\sigma^4$ drops faster', C['red']),
        (0.5, 0.15, '$\\mathcal{A}_g$ rises → additive regime', C['green']),
    ]
    
    for x, y, text, col in boxes:
        ax.add_patch(Rectangle((x - 0.35, y - 0.07), 0.7, 0.14,
                                facecolor=col, alpha=0.3, edgecolor=C['dark'],
                                transform=ax.transAxes))
        ax.text(x, y, text, fontsize=12, ha='center', va='center',
                transform=ax.transAxes)
    
    # Arrows
    for y1, y2 in [(0.68, 0.62), (0.48, 0.42), (0.28, 0.22)]:
        ax.annotate('', xy=(0.5, y2), xytext=(0.5, y1),
                    arrowprops=dict(arrowstyle='->', lw=2, color=C['dark']),
                    xycoords='axes fraction', textcoords='axes fraction')
    
    plt.tight_layout()
    return fig


# =============================================================================
# PART IV: ESCAPE FROM ADDITIVITY (Slides 19-22)
# =============================================================================

def slide19_escape_geometry():
    """Geometry of escape from additive channels."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Contours with trajectories
    ax = axes[0]
    
    x = np.linspace(-2.5, 2.5, 100)
    y = np.linspace(-2.5, 2.5, 100)
    X, Y = np.meshgrid(x, y)
    
    R2 = X**2 + Y**2
    Ag = 1 / (1 + 0.2 * R2)
    
    cs = ax.contourf(X, Y, Ag, levels=20, cmap='RdYlGn', alpha=0.7)
    ax.contour(X, Y, Ag, levels=[0.6, 0.7, 0.8, 0.9], colors='black', linewidths=1)
    plt.colorbar(cs, ax=ax, label='$\\mathcal{A}_g$')
    
    # Trajectory staying in channel
    t = np.linspace(0, 1.5, 50)
    traj1_x = 0.3 * t
    traj1_y = 0.2 * t
    ax.plot(traj1_x, traj1_y, '-', color=C['green'], lw=3, label='Stays in channel')
    ax.scatter([traj1_x[0]], [traj1_y[0]], s=100, c=C['green'], zorder=5)
    ax.scatter([traj1_x[-1]], [traj1_y[-1]], s=100, c=C['green'], marker='^', zorder=5)
    
    # Trajectory escaping
    traj2_x = 0.5 + 1.2 * t
    traj2_y = 0.3 + 0.8 * t
    ax.plot(traj2_x, traj2_y, '-', color=C['red'], lw=3, label='Escapes channel')
    ax.scatter([traj2_x[0]], [traj2_y[0]], s=100, c=C['red'], zorder=5)
    ax.scatter([traj2_x[-1]], [traj2_y[-1]], s=100, c=C['red'], marker='^', zorder=5)
    
    ax.set_xlabel('Genotype 1', fontsize=14)
    ax.set_ylabel('Genotype 2', fontsize=14)
    ax.set_title('Trajectories Through $\\mathcal{A}_g$ Contours', fontsize=14, fontweight='bold')
    ax.legend(loc='lower left')
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    
    # Right: The escape condition
    ax = axes[1]
    ax.axis('off')
    
    ax.text(0.5, 0.9, 'The Escape Equation', fontsize=18, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.7, 
            '$\\dfrac{d\\mathcal{A}_g}{dt} \\propto \\nabla\\mathcal{A}_g^\\top \\Sigma_g \\mathbf{b}$',
            fontsize=20, ha='center', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    conditions = [
        ('$> 0$', 'Moving deeper into channel', C['green']),
        ('$= 0$', 'Along channel boundary', C['amber']),
        ('$< 0$', 'Escaping toward curvature', C['red']),
    ]
    
    for i, (val, text, col) in enumerate(conditions):
        y = 0.48 - i * 0.12
        ax.text(0.25, y, val, fontsize=16, ha='center', va='center',
                transform=ax.transAxes, fontweight='bold', color=col)
        ax.text(0.4, y, text, fontsize=14, va='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.1, 
            'Escape is a geometric event:\nselection drives the population across $\\mathcal{A}_g$ contours',
            fontsize=12, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.2))
    
    plt.tight_layout()
    return fig


def slide20_gxe_rotation():
    """G×E as gradient rotation."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Three panels showing different gradient orientations
    titles = ['Environment 1:\nGradient aligns with $G$', 
              'Environment 2:\nGradient rotates', 
              'Environment 3:\nGradient misaligned']
    angles = [30, 60, 120]  # Gradient angles
    Ag_vals = [0.92, 0.75, 0.45]
    
    for idx, (ax, title, angle, Ag_val) in enumerate(zip(axes, titles, angles, Ag_vals)):
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)
        ax.set_aspect('equal')
        
        # G matrix ellipse (same in all)
        ell = Ellipse((0, 0), 3, 1, angle=30, facecolor=C['blue'], alpha=0.3,
                      edgecolor=C['blue'], lw=2)
        ax.add_patch(ell)
        
        # Selection gradient (changes)
        grad_len = 1.8
        rad = np.radians(angle)
        ax.annotate('', xy=(grad_len * np.cos(rad), grad_len * np.sin(rad)), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='->', lw=3, color=C['red']))
        ax.text(grad_len * np.cos(rad) + 0.1, grad_len * np.sin(rad) + 0.1, 
                '$\\mathbf{b}$', fontsize=14, color=C['red'])
        
        # G eigenvectors
        for a, l in [(30, 1.5), (120, 0.5)]:
            r = np.radians(a)
            ax.plot([0, l * np.cos(r)], [0, l * np.sin(r)], '--', color='gray', lw=1)
        
        ax.set_xlabel('Genotype 1', fontsize=12)
        ax.set_ylabel('Genotype 2', fontsize=12)
        ax.set_title(title, fontsize=12, fontweight='bold')
        
        # A_g annotation
        col = C['green'] if Ag_val > 0.8 else (C['amber'] if Ag_val > 0.6 else C['red'])
        ax.text(0.95, 0.05, f'$\\mathcal{{A}}_g = {Ag_val}$', transform=ax.transAxes,
                ha='right', fontsize=14, fontweight='bold', color=col,
                bbox=dict(boxstyle='round', facecolor='white'))
        
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(0, color='gray', lw=0.5)
    
    fig.suptitle('G×E as Gradient Rotation Relative to Genetic Variance', 
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    return fig


def slide21_introgression():
    """Introgression expands variance and drops A_g."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Left: Ellipse before/after
    ax = axes[0]
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_aspect('equal')
    
    # Elite (small)
    ell1 = Ellipse((0.5, 0.3), 1.2, 0.6, angle=20, facecolor=C['green'], alpha=0.5,
                   edgecolor=C['green'], lw=2, label='Elite germplasm')
    ax.add_patch(ell1)
    
    # After introgression (large)
    ell2 = Ellipse((0, 0), 3.5, 2, angle=35, facecolor=C['red'], alpha=0.3,
                   edgecolor=C['red'], lw=2, label='Post-introgression')
    ax.add_patch(ell2)
    
    ax.set_xlabel('Genotype 1', fontsize=12)
    ax.set_ylabel('Genotype 2', fontsize=12)
    ax.set_title('Variance Expansion', fontsize=14, fontweight='bold')
    ax.legend(loc='lower left')
    
    # Center: A_g trajectory
    ax = axes[1]
    cycles = np.arange(0, 18)
    
    # Pre-introgression (rising)
    Ag_pre = 0.6 + 0.25 * (1 - np.exp(-0.3 * cycles[:8]))
    
    # Introgression drop
    Ag_intro = 0.5
    
    # Post-introgression (recovery)
    Ag_post = 0.5 + 0.35 * (1 - np.exp(-0.2 * (cycles[8:] - 8)))
    
    ax.plot(cycles[:8], Ag_pre, 'o-', color=C['green'], lw=2, ms=6, label='Pre')
    ax.plot([7, 8], [Ag_pre[-1], Ag_intro], 'o-', color=C['red'], lw=2, ms=6)
    ax.scatter([8], [Ag_intro], s=150, color=C['red'], marker='s', zorder=5, label='Introgression')
    ax.plot(cycles[8:], Ag_post, 'o-', color=C['amber'], lw=2, ms=6, label='Recovery')
    
    ax.axvline(8, ls='--', color=C['red'], alpha=0.5)
    ax.axhline(0.8, ls=':', color='gray')
    ax.fill_between([0, 18], 0.8, 1, alpha=0.1, color=C['green'])
    
    ax.set_xlabel('Breeding Cycle', fontsize=12)
    ax.set_ylabel('$\\mathcal{A}_g$', fontsize=14)
    ax.set_title('$\\mathcal{A}_g$ Trajectory', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    ax.set_ylim(0.3, 1)
    
    # Right: Key insight
    ax = axes[2]
    ax.axis('off')
    
    ax.text(0.5, 0.9, 'The Introgression Effect', fontsize=16, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    steps = [
        ('1', 'Wide cross expands $\\Sigma_g$', C['blue']),
        ('2', '$V_{\\text{quad}} \\propto \\sigma^4$ grows faster', C['red']),
        ('3', '$\\mathcal{A}_g$ drops sharply', C['red']),
        ('4', 'GBLUP accuracy decreases', C['amber']),
        ('5', 'Selection re-compresses', C['green']),
        ('6', '$\\mathcal{A}_g$ recovers over cycles', C['green']),
    ]
    
    for i, (num, text, col) in enumerate(steps):
        y = 0.75 - i * 0.11
        ax.add_patch(Circle((0.08, y), 0.03, color=col, transform=ax.transAxes))
        ax.text(0.08, y, num, fontsize=11, ha='center', va='center',
                color='white', fontweight='bold', transform=ax.transAxes)
        ax.text(0.15, y, text, fontsize=12, va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    return fig


def slide22_driving_analogy():
    """The driving analogy visualization."""
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Draw a road with curves
    road_x = np.array([1, 3, 5, 7, 8, 9, 10, 11, 13])
    road_y = np.array([4, 4, 4, 4.5, 5.5, 6, 5.5, 5, 5])
    
    from scipy.interpolate import make_interp_spline
    road_spline = make_interp_spline(road_x, road_y, k=3)
    road_x_smooth = np.linspace(1, 13, 200)
    road_y_smooth = road_spline(road_x_smooth)
    
    # Road surface
    ax.fill_between(road_x_smooth, road_y_smooth - 0.5, road_y_smooth + 0.5,
                    color=C['dark'], alpha=0.3)
    ax.plot(road_x_smooth, road_y_smooth - 0.5, color=C['dark'], lw=2)
    ax.plot(road_x_smooth, road_y_smooth + 0.5, color=C['dark'], lw=2)
    
    # Center line
    for i in range(0, 190, 20):
        ax.plot(road_x_smooth[i:i+10], road_y_smooth[i:i+10], '--', color='yellow', lw=2)
    
    # Labels for road sections
    ax.text(2.5, 2.5, 'Straight section\n(Additive channel)', ha='center', fontsize=12,
            bbox=dict(boxstyle='round', facecolor=C['green'], alpha=0.3))
    ax.annotate('', xy=(2.5, 3.3), xytext=(2.5, 3),
                arrowprops=dict(arrowstyle='->', lw=2, color=C['dark']))
    
    ax.text(8.5, 7.3, 'Curve\n(Curvature matters)', ha='center', fontsize=12,
            bbox=dict(boxstyle='round', facecolor=C['red'], alpha=0.3))
    ax.annotate('', xy=(8.5, 6.5), xytext=(8.5, 7),
                arrowprops=dict(arrowstyle='->', lw=2, color=C['dark']))
    
    ax.text(12, 3.5, 'New straight\n(New additive channel)', ha='center', fontsize=12,
            bbox=dict(boxstyle='round', facecolor=C['green'], alpha=0.3))
    ax.annotate('', xy=(12, 4.3), xytext=(12, 3.9),
                arrowprops=dict(arrowstyle='->', lw=2, color=C['dark']))
    
    # Car icon (simple)
    car_x, car_y = 4, 4
    ax.add_patch(Rectangle((car_x - 0.3, car_y - 0.15), 0.6, 0.3, 
                            facecolor=C['blue'], edgecolor='black'))
    
    # Title and key message
    ax.text(7, 0.8, 
            'Cruise control works on straight roads (high $\\mathcal{A}_g$)\n' +
            'Active steering needed at curves (low $\\mathcal{A}_g$)\n' +
            'Most driving is on straight sections!',
            ha='center', fontsize=14, style='italic',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor=C['dark']))
    
    ax.set_title('The Driving Analogy', fontsize=20, fontweight='bold', y=0.98)
    
    return fig


# =============================================================================
# PART V: APPLICATIONS (Slides 23-27)
# =============================================================================

def slide23_breeding_dashboard():
    """A_g as a breeding decision metric."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    ax.text(0.5, 0.95, 'The $\\mathcal{A}_g$ Breeding Dashboard', fontsize=20, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    # Three columns for different A_g regimes
    regimes = [
        (0.17, 'High $\\mathcal{A}_g$ (>0.85)', C['green'],
         ['Standard GBLUP works', 'Invest in data quality', 'Linear model ranks well']),
        (0.5, 'Moderate (0.6–0.85)', C['amber'],
         ['Expect some re-ranking', 'Consider interactions', 'Monitor closely']),
        (0.83, 'Low $\\mathcal{A}_g$ (<0.6)', C['red'],
         ['Change the model', 'Add nonlinear terms', 'Curvature matters']),
    ]
    
    for x, title, col, actions in regimes:
        ax.add_patch(Rectangle((x - 0.14, 0.35), 0.28, 0.5, 
                                facecolor=col, alpha=0.2, edgecolor=col, lw=3,
                                transform=ax.transAxes))
        ax.text(x, 0.78, title, fontsize=14, fontweight='bold', ha='center',
                transform=ax.transAxes)
        
        for i, action in enumerate(actions):
            ax.text(x, 0.65 - i * 0.1, f'• {action}', fontsize=11, ha='center',
                    transform=ax.transAxes)
    
    # Key insight at bottom
    ax.text(0.5, 0.15, 
            '$\\mathcal{A}_g$ tells you: Is my bottleneck data quality or model specification?',
            fontsize=14, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    ax.text(0.5, 0.05,
            'High $\\mathcal{A}_g$ + low accuracy → Need more data\n' +
            'Low $\\mathcal{A}_g$ + low accuracy → Need better model',
            fontsize=12, ha='center', transform=ax.transAxes)
    
    return fig


def slide24_met_pipeline():
    """MET pipeline integration diagram."""
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('off')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    
    # Pipeline boxes
    boxes = [
        (1.5, 6.5, 'MET Data', C['blue']),
        (4.5, 6.5, 'Mixed Model\n(estimate $\\mathbf{G}$)', C['blue']),
        (7.5, 6.5, 'Selection Index\n(define $\\mathbf{b}$)', C['blue']),
        (4.5, 4, 'Quadratic Fit\n(estimate $\\mathbf{H}$)', C['amber']),
        (7.5, 4, 'Compute $\\mathcal{A}_g$\nwith uncertainty', C['green']),
        (5.5, 1.5, 'Breeding Dashboard', C['purple']),
    ]
    
    for x, y, text, col in boxes:
        ax.add_patch(Rectangle((x - 1, y - 0.5), 2, 1, facecolor=col, alpha=0.3,
                                edgecolor=col, lw=2))
        ax.text(x, y, text, fontsize=11, ha='center', va='center')
    
    # Arrows
    arrows = [
        ((2.5, 6.5), (3.5, 6.5)),  # MET to Mixed Model
        ((5.5, 6.5), (6.5, 6.5)),  # Mixed Model to Index
        ((4.5, 6), (4.5, 4.5)),    # Mixed Model to Quadratic
        ((5.5, 4), (6.5, 4)),      # Quadratic to Ag
        ((7.5, 6), (7.5, 4.5)),    # Index to Ag
        ((7.5, 3.5), (6.2, 2)),    # Ag to Dashboard
        ((4.5, 3.5), (5.2, 2)),    # Quadratic to Dashboard
    ]
    
    for (x1, y1), (x2, y2) in arrows:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', lw=2, color=C['dark']))
    
    # Labels
    ax.text(5, 7.5, 'You Already Do This', fontsize=14, fontweight='bold', ha='center',
            bbox=dict(boxstyle='round', facecolor=C['blue'], alpha=0.2))
    ax.text(6, 5, 'New Addition', fontsize=12, fontweight='bold', ha='center',
            bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.2))
    
    ax.set_title('Integrating $\\mathcal{A}_g$ into Your MET Pipeline', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    return fig


def slide25_simulation():
    """Simulation validation: A_g predicts accuracy."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    np.random.seed(42)
    
    # Left: Scatter of A_g vs accuracy
    ax = axes[0]
    
    n = 200
    Ag = np.random.uniform(0.4, 0.95, n)
    # Accuracy bounded by sqrt(Ag) with noise
    accuracy = np.sqrt(Ag) * (0.85 + 0.15 * np.random.rand(n)) - np.random.uniform(0, 0.15, n)
    accuracy = np.clip(accuracy, 0, 1)
    
    ax.scatter(Ag, accuracy, alpha=0.5, c=Ag, cmap='RdYlGn', s=50, edgecolor='white')
    
    # Theoretical upper bound
    Ag_line = np.linspace(0.35, 1, 100)
    ax.plot(Ag_line, np.sqrt(Ag_line), 'r--', lw=3, label='Upper bound $\\sqrt{\\mathcal{A}_g}$')
    
    # Trend line
    from numpy.polynomial import polynomial as P
    coefs = P.polyfit(Ag, accuracy, 1)
    ax.plot(Ag_line, P.polyval(Ag_line, coefs), 'k-', lw=2, label='Trend')
    
    r = np.corrcoef(Ag, accuracy)[0, 1]
    ax.text(0.05, 0.95, f'r = {r:.2f}', transform=ax.transAxes, fontsize=14, va='top',
            bbox=dict(boxstyle='round', facecolor='white'))
    
    ax.set_xlabel('Additivity Index ($\\mathcal{A}_g$)', fontsize=14)
    ax.set_ylabel('GBLUP Prediction Accuracy', fontsize=14)
    ax.set_title('Simulation: $\\mathcal{A}_g$ Predicts Accuracy', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    ax.set_xlim(0.35, 1)
    ax.set_ylim(0, 1)
    
    # Right: Different scenarios
    ax = axes[1]
    
    scenarios = ['Strong\nselection', 'Weak\nselection', 'Post-\nintrogression', 'High\nG×E']
    Ag_means = [0.88, 0.72, 0.55, 0.65]
    Ag_errs = [0.04, 0.08, 0.12, 0.10]
    acc_means = [0.82, 0.68, 0.52, 0.58]
    acc_errs = [0.05, 0.08, 0.10, 0.12]
    
    x = np.arange(len(scenarios))
    width = 0.35
    
    ax.bar(x - width/2, Ag_means, width, yerr=Ag_errs, label='$\\mathcal{A}_g$',
           color=C['green'], edgecolor='white', capsize=5)
    ax.bar(x + width/2, acc_means, width, yerr=acc_errs, label='Accuracy',
           color=C['blue'], edgecolor='white', capsize=5)
    
    ax.set_xticks(x)
    ax.set_xticklabels(scenarios, fontsize=11)
    ax.set_ylabel('Value', fontsize=14)
    ax.set_title('$\\mathcal{A}_g$ and Accuracy Track Together', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    
    plt.tight_layout()
    return fig


def slide26_radiations():
    """Adaptive radiations and escape from additivity."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Stylized phylogeny with A_g mapped
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Draw tree
    # Trunk
    ax.plot([5, 5], [0.5, 2], color=C['dark'], lw=4)
    # Main branches
    ax.plot([5, 3], [2, 4], color=C['dark'], lw=3)
    ax.plot([5, 7], [2, 4], color=C['dark'], lw=3)
    # Secondary branches
    ax.plot([3, 2], [4, 6], color=C['red'], lw=2)  # Low Ag during radiation
    ax.plot([3, 4], [4, 6], color=C['red'], lw=2)
    ax.plot([7, 6], [4, 6], color=C['red'], lw=2)
    ax.plot([7, 8], [4, 6], color=C['red'], lw=2)
    # Tips (stabilized)
    for x in [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]:
        ax.plot([x - 0.5 + (x > 5) * 0.5, x], [6, 7], color=C['green'], lw=2)
        ax.scatter([x], [7.2], s=80, c=C['green'], zorder=5)
    
    # Legend
    ax.plot([1, 1.5], [1, 1], color=C['red'], lw=3)
    ax.text(1.7, 1, 'Radiation (low $\\mathcal{A}_g$)', fontsize=10, va='center')
    ax.plot([1, 1.5], [0.5, 0.5], color=C['green'], lw=3)
    ax.text(1.7, 0.5, 'Stable (high $\\mathcal{A}_g$)', fontsize=10, va='center')
    
    ax.set_title('$\\mathcal{A}_g$ Mapped onto a Radiation', fontsize=14, fontweight='bold')
    
    # Right: Interpretation
    ax = axes[1]
    ax.axis('off')
    
    ax.text(0.5, 0.9, 'Adaptive Radiations', fontsize=16, fontweight='bold',
            ha='center', transform=ax.transAxes)
    
    phases = [
        ('Early radiation', 'Low $\\mathcal{A}_g$: exploring new niches, curvature matters', C['red']),
        ('Diversification', 'Moderate: lineages finding local optima', C['amber']),
        ('Stabilization', 'High $\\mathcal{A}_g$: lineages settle into new channels', C['green']),
    ]
    
    for i, (phase, desc, col) in enumerate(phases):
        y = 0.7 - i * 0.2
        ax.add_patch(Rectangle((0.05, y - 0.06), 0.25, 0.12,
                                facecolor=col, alpha=0.3, transform=ax.transAxes))
        ax.text(0.175, y, phase, fontsize=12, ha='center', va='center',
                fontweight='bold', transform=ax.transAxes)
        ax.text(0.35, y, desc, fontsize=11, va='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.12, 
            'Radiations are natural "wide crosses":\n' +
            'variance expansion → curvature matters → recovery',
            fontsize=12, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    plt.tight_layout()
    return fig


def slide27_diagnostics():
    """Empirical diagnostics for detecting channel exit."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    ax.text(0.5, 0.95, 'Warning Signs: Leaving an Additive Channel', 
            fontsize=18, fontweight='bold', ha='center', transform=ax.transAxes)
    
    warnings = [
        ('Prediction accuracy drops systematically', 
         'GBLUP $r$ declining across validation sets', '📉'),
        ('$\\mathbf{G}$ eigenvectors rotate rapidly', 
         'Leading PCs of genetic variance change direction', '🔄'),
        ('Response deviates from $\\mathbf{G}\\boldsymbol{\\beta}$', 
         'Observed vs predicted response diverges', '↗️'),
        ('Environmental sensitivity increases', 
         'G×E variance component growing', '🌡️'),
        ('Realized heritability drops', 
         '$h^2$ estimates declining unexpectedly', '📊'),
    ]
    
    for i, (warning, detail, icon) in enumerate(warnings):
        y = 0.78 - i * 0.15
        
        ax.text(0.05, y, icon, fontsize=24, va='center', transform=ax.transAxes)
        ax.text(0.12, y + 0.02, warning, fontsize=13, fontweight='bold', 
                va='center', transform=ax.transAxes)
        ax.text(0.12, y - 0.03, detail, fontsize=11, va='center',
                transform=ax.transAxes, color='gray')
    
    ax.text(0.5, 0.08, 
            'All of these are measurable from standard breeding/evolutionary datasets',
            fontsize=12, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.2))
    
    return fig


# =============================================================================
# PART VI: SYNTHESIS (Slides 28-30)
# =============================================================================

def slide28_resolution():
    """Resolving the paradox."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    ax.text(0.5, 0.95, 'Resolving the Paradox', fontsize=20, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    # Central message
    ax.text(0.5, 0.8, 
            'Additivity is not a property of the genome.\n' +
            'It\'s a property of how much of the landscape\n' +
            'the population samples.',
            fontsize=16, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['light']))
    
    # Four pillars
    pillars = [
        ('Local Linearity', 'Any smooth function\nis locally linear', C['blue']),
        ('Scaling', '$V_{\\text{quad}} \\propto \\sigma^4$\nis negligible at small $\\sigma$', C['green']),
        ('Self-Correction', 'Selection compresses\nvariance → raises $\\mathcal{A}_g$', C['amber']),
        ('Transient Curves', 'Escape is brief;\nchannels dominate', C['purple']),
    ]
    
    for i, (title, desc, col) in enumerate(pillars):
        x = 0.15 + i * 0.23
        ax.add_patch(Rectangle((x - 0.09, 0.35), 0.18, 0.35,
                                facecolor=col, alpha=0.2, edgecolor=col, lw=2,
                                transform=ax.transAxes))
        ax.text(x, 0.62, title, fontsize=12, fontweight='bold', ha='center',
                transform=ax.transAxes)
        ax.text(x, 0.48, desc, fontsize=10, ha='center', va='center',
                transform=ax.transAxes)
    
    ax.text(0.5, 0.15, 
            'Epistatic architectures produce additive-looking dynamics\n' +
            'because selection acts to keep populations in additive channels',
            fontsize=14, ha='center', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor=C['green'], alpha=0.2))
    
    return fig


def slide29_takeaways():
    """Key take-home messages."""
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.axis('off')
    
    ax.text(0.5, 0.95, 'Key Take-Home Messages', fontsize=20, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    messages = [
        '$\\mathcal{A}_g = V_{\\text{lin}} / (V_{\\text{lin}} + V_{\\text{quad}})$ measures "local additivity"',
        '$\\mathcal{A}_g = R^2$ of the linear predictor for log-fitness',
        'Fisher\'s theorem & breeder\'s equation work when $\\mathcal{A}_g$ is high',
        'Selection naturally pushes populations into additive channels',
        'Wide crosses, new environments, and radiations trigger escape',
        '$\\mathcal{A}_g$ can be estimated and used as a diagnostic',
    ]
    
    for i, msg in enumerate(messages):
        y = 0.82 - i * 0.12
        ax.add_patch(Circle((0.06, y), 0.025, color=C['blue'], transform=ax.transAxes))
        ax.text(0.06, y, str(i+1), fontsize=12, ha='center', va='center',
                color='white', fontweight='bold', transform=ax.transAxes)
        ax.text(0.11, y, msg, fontsize=14, va='center', transform=ax.transAxes)
    
    ax.text(0.5, 0.08, 
            'The framework unifies landscape theory with practical prediction',
            fontsize=14, ha='center', transform=ax.transAxes, style='italic',
            bbox=dict(boxstyle='round', facecolor=C['amber'], alpha=0.3))
    
    return fig


def slide30_future():
    """Future directions."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    ax.text(0.5, 0.95, 'Future Directions', fontsize=20, 
            fontweight='bold', ha='center', transform=ax.transAxes)
    
    directions = [
        ('Empirical estimation', 'Measure $\\mathcal{A}_g$ across taxa and traits'),
        ('Non-Gaussian extensions', 'Selection truncation, skewness'),
        ('Dynamic $\\mathbf{G}$', 'Track channels through evolutionary time'),
        ('Breeding optimization', 'Design programs that maintain high $\\mathcal{A}_g$'),
        ('Landscape connections', 'Link to DFE, rugged landscapes'),
    ]
    
    for i, (title, desc) in enumerate(directions):
        y = 0.78 - i * 0.13
        ax.text(0.1, y, '→', fontsize=20, va='center', transform=ax.transAxes,
                color=C['blue'])
        ax.text(0.18, y + 0.015, title, fontsize=14, fontweight='bold',
                va='center', transform=ax.transAxes)
        ax.text(0.18, y - 0.03, desc, fontsize=12, va='center',
                transform=ax.transAxes, color='gray')
    
    ax.text(0.5, 0.18, 'Thank You!', fontsize=28, ha='center',
            transform=ax.transAxes, fontweight='bold', color=C['blue'])
    
    ax.text(0.5, 0.08, 'Questions?', fontsize=18, ha='center',
            transform=ax.transAxes, style='italic')
    
    return fig


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("GENERATING PRESENTATION FIGURES")
    print("=" * 60)
    
    # Part I
    print("\nPart I: The Paradox")
    save_fig(slide01_title(), 'slide01_title')
    save_fig(slide02_empirical_success(), 'slide02_empirical_success')
    save_fig(slide03_molecular_reality(), 'slide03_molecular_reality')
    save_fig(slide04_paradox(), 'slide04_paradox')
    save_fig(slide05_preview(), 'slide05_preview')
    save_fig(slide06_outline(), 'slide06_outline')
    
    # Part II
    print("\nPart II: The Framework")
    save_fig(slide07_gpf_map(), 'slide07_gpf_map')
    save_fig(slide08_taylor(), 'slide08_taylor')
    save_fig(slide09_qg_connection(), 'slide09_qg_connection')
    save_fig(slide10_variance_scaling(), 'slide10_variance_scaling')
    save_fig(slide11_additivity_index(), 'slide11_additivity_index')
    save_fig(slide12_Ag_as_Rsquared(), 'slide12_Ag_as_Rsquared')
    
    # Part III
    print("\nPart III: Additive Channels")
    save_fig(slide13_channels_defined(), 'slide13_channels_defined')
    save_fig(slide14_fisher(), 'slide14_fisher')
    save_fig(slide15_heritability(), 'slide15_heritability')
    save_fig(slide16_1d_condition(), 'slide16_1d_condition')
    save_fig(slide17_interpretation(), 'slide17_interpretation')
    save_fig(slide18_self_correcting(), 'slide18_self_correcting')
    
    # Part IV
    print("\nPart IV: Escape from Additivity")
    save_fig(slide19_escape_geometry(), 'slide19_escape_geometry')
    save_fig(slide20_gxe_rotation(), 'slide20_gxe_rotation')
    save_fig(slide21_introgression(), 'slide21_introgression')
    save_fig(slide22_driving_analogy(), 'slide22_driving_analogy')
    
    # Part V
    print("\nPart V: Applications")
    save_fig(slide23_breeding_dashboard(), 'slide23_breeding_dashboard')
    save_fig(slide24_met_pipeline(), 'slide24_met_pipeline')
    save_fig(slide25_simulation(), 'slide25_simulation')
    save_fig(slide26_radiations(), 'slide26_radiations')
    save_fig(slide27_diagnostics(), 'slide27_diagnostics')
    
    # Part VI
    print("\nPart VI: Synthesis")
    save_fig(slide28_resolution(), 'slide28_resolution')
    save_fig(slide29_takeaways(), 'slide29_takeaways')
    save_fig(slide30_future(), 'slide30_future')
    
    print("\n" + "=" * 60)
    print("ALL FIGURES COMPLETE!")
    print("=" * 60)
