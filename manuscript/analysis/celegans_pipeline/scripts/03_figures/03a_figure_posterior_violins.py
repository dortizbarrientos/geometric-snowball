#!/usr/bin/env python3
"""
03a_figure_posterior_violins.py
================================

PURPOSE
-------
Build the publication-quality figure showing the full posterior
distributions of (a) channel anisotropy (lambda_1/lambda_2 of G) and
(b) channel rotation (g_max to S_tilde-channel) across the four
C. elegans populations in the high-salt environment.

INPUT
-----
    data/processed/celegans_posterior_perdraw.csv

OUTPUT
------
    results/figures/celegans_posterior_violins.pdf
    results/figures/celegans_posterior_violins.png

DESIGN CHOICES
--------------
- Violins show full posterior density (kernel density estimate).
- Point = posterior median; vertical bar = 95% credible interval.
- Log10 axis for the ratio panel because ratios span a wide range
  (ancestor near 10, boundary replicates near 1).
- Okabe-Ito colorblind-safe palette.
- Fonts sized for a Nature Brief Communication double-column figure
  (183 mm wide, ~80 mm tall).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PROCESSED = ROOT / "data" / "processed"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# LOAD PER-DRAW DATA
# ----------------------------------------------------------------------
df = pd.read_csv(PROCESSED / "celegans_posterior_perdraw.csv")

# This figure focuses on the high-salt environment (primary selection env)
df = df[(df.environment == "NaCl") & (df.was_PD == True)].copy()

# ----------------------------------------------------------------------
# VISUAL SETUP
# ----------------------------------------------------------------------
COLORS = {
    "A6140": "#0072B2",   # dark blue = ancestor
    "GA150": "#D55E00",   # vermillion
    "GA250": "#CC79A7",   # purple/pink
    "GA450": "#E69F00",   # orange
}
POP_ORDER = ["A6140", "GA150", "GA250", "GA450"]
LABELS = {
    "A6140": "A6140\n(ancestor)",
    "GA150": "GA150",
    "GA250": "GA250",
    "GA450": "GA450",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6,
})
mm = 1.0 / 25.4


# ----------------------------------------------------------------------
# HELPER: a single violin panel
# ----------------------------------------------------------------------
def add_violin_panel(ax, data_by_pop, log10=False, boundary=None,
                     boundary_label=None, ylabel="", title="", label=""):
    positions = list(range(len(POP_ORDER)))
    for i, pop in enumerate(POP_ORDER):
        vals = data_by_pop[pop]
        if log10:
            vals = np.log10(vals)
        parts = ax.violinplot([vals], positions=[i], widths=0.7,
                              showmeans=False, showmedians=False,
                              showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(COLORS[pop])
            pc.set_alpha(0.45)
            pc.set_edgecolor(COLORS[pop])
            pc.set_linewidth(0.8)
        # Median + 95 % CI
        med = np.median(vals)
        lo, hi = np.quantile(vals, [0.025, 0.975])
        ax.plot([i, i], [lo, hi], color=COLORS[pop], lw=1.6)
        ax.scatter([i], [med], s=30, color="white",
                   edgecolor=COLORS[pop], linewidth=1.5, zorder=10)

    if boundary is not None:
        y = np.log10(boundary) if log10 else boundary
        ax.axhline(y, color="0.35", lw=0.6, ls="--", alpha=0.7)
        if boundary_label:
            ax.text(len(POP_ORDER) - 0.5, y, boundary_label,
                    fontsize=6.5, color="0.35",
                    va="bottom", ha="right")

    ax.set_xticks(positions)
    ax.set_xticklabels([LABELS[p] for p in POP_ORDER], fontsize=7.5)
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.set_title(title, fontsize=9, pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(0.02, 0.98, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top")


# ----------------------------------------------------------------------
# BUILD FIGURE
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(183 * mm, 80 * mm))

# Panel a: ratio
ratio_by_pop = {pop: df[df.population == pop]["ratio"].to_numpy()
                for pop in POP_ORDER}
add_violin_panel(
    axes[0], ratio_by_pop, log10=True,
    boundary=1.0, boundary_label="spectral-gap\nboundary",
    ylabel=r"$\log_{10}\left(\lambda_1 / \lambda_2\right)$",
    title=r"Channel anisotropy of $\mathbf{G}$",
    label="a",
)

# Panel b: rotation
rot_by_pop = {pop: df[df.population == pop]["rotation_deg"].to_numpy()
              for pop in POP_ORDER}
add_violin_panel(
    axes[1], rot_by_pop, log10=False,
    boundary=90.0, boundary_label="orthogonal limit",
    ylabel=r"rotation $\mathbf{g}_{\max} \to \tilde{\mathbf{S}}$-channel (degrees)",
    title="Channel rotation under fitness curvature",
    label="b",
)
axes[1].set_ylim(-3, 95)

plt.tight_layout()

out_pdf = FIG_DIR / "celegans_posterior_violins.pdf"
out_png = FIG_DIR / "celegans_posterior_violins.png"
plt.savefig(out_pdf, format="pdf", dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(out_png, format="png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close()

print(f"Wrote {out_pdf}")
print(f"Wrote {out_png}")
