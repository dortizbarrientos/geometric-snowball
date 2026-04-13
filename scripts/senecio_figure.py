"""
senecio_figure.py
=================
Publication-quality figure summarising S. lautus Dune–Headland pair analyses.

Four panels:
  a  PCA ordination of all S. lautus ecotypes (population means)
     → Shows where Dune and Headland sit in the full radiation
  b  Phenotypic divergence per pair (standardised Euclidean Δμ)
     → Proxy for channel rotation θ₁₂
  c  Environmental divergence per pair (standardised Euclidean ΔE)
     → Proxy for fitness landscape contrast ΔH
  d  Env. divergence vs. Pheno. divergence scatter
     → Tests the prediction that ecological contrast ∝ phenotypic divergence
     → Pair 4 and 3 expected to be most divergent
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr
import shutil

# ── Palette ───────────────────────────────────────────────────────────────────
BLUE    = "#0072B2"
SKY     = "#56B4E9"
ORANGE  = "#E69F00"
VERM    = "#D55E00"
GREEN   = "#009E73"
PURPLE  = "#CC79A7"
YELLOW  = "#F0E442"
GREY    = "#999999"

# Ecotype colour map
ECO_COLS = {
    'D.Dune.N'     : BLUE,
    'G.Head_N'     : ORANGE,
    'H.Head_S'     : VERM,
    'L.SeaBird'    : PURPLE,
    'A.Alpine_E'   : GREEN,
    'B.Alpine_P'   : SKY,
    'C.Brigalow'   : '#8B4513',
    'F.Lacustrinus': '#4B0082',
    'I.Depressicola': '#006400',
    'J.Lanceolatus': '#8B008B',
    'K.Mallee'     : '#FF8C00',
    'M.Tableland'  : '#2F4F4F',
    'N.Woodland'   : GREY,
}

ECO_LABELS = {
    'D.Dune.N'     : 'Dune (N)',
    'G.Head_N'     : 'Headland (N)',
    'H.Head_S'     : 'Headland (S)',
    'L.SeaBird'    : 'Sea Bird Colony',
    'A.Alpine_E'   : 'Alpine (E)',
    'B.Alpine_P'   : 'Alpine (P)',
    'C.Brigalow'   : 'Brigalow',
    'F.Lacustrinus': 'Lacustrinus',
    'I.Depressicola': 'Depressicola',
    'J.Lanceolatus': 'Lanceolatus',
    'K.Mallee'     : 'Mallee',
    'M.Tableland'  : 'Tableland',
    'N.Woodland'   : 'Woodland',
}

plt.rcParams.update({
    "font.family"      : "serif",
    "mathtext.fontset" : "dejavuserif",
    "pdf.fonttype"     : 42,
    "axes.linewidth"   : 0.7,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

# ── Load pre-computed CSVs ─────────────────────────────────────────────────────
pca_df   = pd.read_csv('/home/claude/pca_pops.csv')
merged   = pd.read_csv('/home/claude/pair_summary.csv', index_col='Pair')
loadings = pd.read_csv('/home/claude/pca_loadings.csv', index_col=0)

# ── Clade colours for pairs ───────────────────────────────────────────────────
CLADE_COL = {'eastern': BLUE, 'southeastern': ORANGE}
CLADE_MRK = {'eastern': 'o',  'southeastern': 's'}


def make_figure(outfile="senecio_pair_analysis.pdf"):

    mm = 1/25.4
    fig = plt.figure(figsize=(183*mm, 190*mm))
    gs  = GridSpec(3, 2, figure=fig,
                   hspace=0.48, wspace=0.38,
                   left=0.10, right=0.96,
                   top=0.95,  bottom=0.07)

    ax_pca  = fig.add_subplot(gs[0, :])   # full-width top
    ax_phen = fig.add_subplot(gs[1, 0])
    ax_env  = fig.add_subplot(gs[1, 1])
    ax_scat = fig.add_subplot(gs[2, :])   # full-width bottom

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL a — PCA ordination
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_pca

    # Plot non-focal ecotypes first (smaller, greyed)
    focal = {'D.Dune.N', 'G.Head_N', 'H.Head_S'}
    for eco in pca_df['Ecotype'].unique():
        sub = pca_df[pca_df['Ecotype'] == eco]
        col = ECO_COLS.get(eco, GREY)
        if eco not in focal:
            ax.scatter(sub['PC1'], sub['PC2'],
                       color=col, s=28, alpha=0.50, zorder=2,
                       marker='o', label=ECO_LABELS.get(eco, eco))
        else:
            # Focal ecotypes: larger, solid
            mkr = 'o' if 'Dune' in eco else ('^' if 'Head_N' in eco else 's')
            ax.scatter(sub['PC1'], sub['PC2'],
                       color=col, s=60, alpha=0.88, zorder=4,
                       edgecolors='white', linewidths=0.4,
                       marker=mkr, label=ECO_LABELS.get(eco, eco))

    # Draw connecting lines between paired Dune and Headland populations
    PAIR_POPS = {
        1: ('DN0','HN0'), 2: ('DN1','HN1'), 3: ('DN3','HN2'),
        4: ('DN4','HN5'), 5: ('DN12','HN14'), 6: ('DN23','HS21'),
        7: ('DN32','HN12'), 8: ('DN14','HN15'),
    }
    for pid, (dp, hp) in PAIR_POPS.items():
        rd = pca_df[pca_df['Population'] == dp]
        rh = pca_df[pca_df['Population'] == hp]
        if rd.empty or rh.empty:
            continue
        x = [float(rd['PC1'].iloc[0]), float(rh['PC1'].iloc[0])]
        y = [float(rd['PC2'].iloc[0]), float(rh['PC2'].iloc[0])]
        clade = merged.loc[pid, 'Clade'] if pid in merged.index else 'eastern'
        ax.plot(x, y, '-', color=CLADE_COL.get(clade, GREY),
                lw=0.9, alpha=0.55, zorder=3)
        # Pair number at midpoint
        mx, my = np.mean(x), np.mean(y)
        ax.text(mx, my + 0.08, str(pid),
                ha='center', va='bottom', fontsize=5.5,
                color=CLADE_COL.get(clade, GREY))

    # Biplot arrows for top-loading traits
    scale = 2.4
    arrow_traits = ['TotalHeight', 'Circularity', 'Nindents',
                    'Area', 'IndentWidth']
    arrow_labels = {'TotalHeight':'Height', 'Circularity':'Circularity',
                    'Nindents':'N indents', 'Area':'Leaf area',
                    'IndentWidth':'Indent width'}
    for tr in arrow_traits:
        dx = loadings.loc[tr, 'PC1_load'] * scale
        dy = loadings.loc[tr, 'PC2_load'] * scale
        ax.annotate('', xy=(dx, dy), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='->', color='0.50',
                                   lw=0.8), zorder=5)
        ax.text(dx*1.12, dy*1.12, arrow_labels.get(tr, tr),
                ha='center', va='center', fontsize=5.2, color='0.40')

    ax.axhline(0, color='0.85', lw=0.5, ls='--', zorder=1)
    ax.axvline(0, color='0.85', lw=0.5, ls='--', zorder=1)
    ax.set_xlabel(f'PC1  (36.5% variance)', fontsize=7.5)
    ax.set_ylabel(f'PC2  (25.6% variance)', fontsize=7.5)
    ax.set_title(r'$\mathbf{a}$  Morphological ordination of $S.\,lautus$ ecotypes',
                 fontsize=8, loc='left', pad=3)
    ax.tick_params(labelsize=6.5)

    # Legend (two columns: focal highlighted)
    handles, labels = ax.get_legend_handles_labels()
    # Reorder: focal first
    focal_labels = ['Dune (N)', 'Headland (N)', 'Headland (S)']
    order = ([labels.index(l) for l in focal_labels if l in labels] +
             [i for i, l in enumerate(labels) if l not in focal_labels])
    ax.legend([handles[i] for i in order],
              [labels[i] for i in order],
              fontsize=5.0, ncol=4, loc='lower right',
              framealpha=0.90, edgecolor='0.82',
              handlelength=0.8, handletextpad=0.35,
              labelspacing=0.28, columnspacing=0.6,
              markerscale=0.8)

    for sp in ax.spines.values():
        sp.set_linewidth(0.6); sp.set_color('#c0c0c0')
    ax.grid(True, color='0.93', lw=0.5, zorder=0)

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL b — Phenotypic divergence per pair
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_phen
    pairs_sorted = merged.sort_values('Pair').index
    x     = np.arange(len(pairs_sorted))
    cols  = [CLADE_COL[merged.loc[p, 'Clade']] for p in pairs_sorted]
    bars  = ax.bar(x, merged.loc[pairs_sorted, 'Pheno_dist'],
                   color=cols, alpha=0.82, width=0.65,
                   edgecolor='white', linewidth=0.5, zorder=3)

    # Annotate pair numbers above bars
    for xi, pid in zip(x, pairs_sorted):
        ax.text(xi, merged.loc[pid, 'Pheno_dist'] + 0.05,
                str(pid), ha='center', va='bottom', fontsize=6.0,
                color='0.25')

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"P{p}\n{merged.loc[p,'Dune_pop']}\n{merged.loc[p,'Head_pop']}"
         for p in pairs_sorted],
        fontsize=5.2, rotation=0, ha='center')
    ax.set_ylabel(r'Phenotypic divergence $\Delta\mu$' + '\n(standardised Euclidean)',
                  fontsize=7.0)
    ax.set_title(r'$\mathbf{b}$  Phenotypic divergence  (proxy: $\theta_{12}$)',
                 fontsize=8, loc='left', pad=3)
    ax.tick_params(labelsize=6.5)
    ax.set_ylim(0, merged['Pheno_dist'].max() * 1.25)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6); sp.set_color('#c0c0c0')
    ax.grid(True, axis='y', color='0.93', lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    # Clade legend
    ax.legend(handles=[
        mpatches.Patch(fc=BLUE,   ec='white', label='Eastern clade'),
        mpatches.Patch(fc=ORANGE, ec='white', label='Southeastern clade'),
    ], fontsize=6.0, loc='upper left', framealpha=0.90, edgecolor='0.82')

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL c — Environmental divergence per pair
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_env

    # Stacked horizontal decomposition: show ΔpH and ΔTOC as component bars
    # (signed: Dune−Headland; Dune typically has higher pH, lower TOC)
    env_s = merged.loc[pairs_sorted]

    # Normalise soil variables to max abs value for comparability
    max_env = merged['Env_dist'].max()
    bars2 = ax.bar(x, env_s['Env_dist'],
                   color=cols, alpha=0.82, width=0.65,
                   edgecolor='white', linewidth=0.5, zorder=3)

    # Small overlay: ΔpH magnitude as hatched portion
    delta_pH_norm = np.abs(env_s['Delta_pH']) / (soil_pH_range := 3.0)
    delta_TOC_norm = np.abs(env_s['Delta_TOC']) / (soil_TOC_range := 7.0)

    for xi, pid in zip(x, pairs_sorted):
        ax.text(xi, env_s.loc[pid, 'Env_dist'] + 0.05,
                str(pid), ha='center', va='bottom', fontsize=6.0,
                color='0.25')

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"P{p}\n{merged.loc[p,'Dune_pop']}\n{merged.loc[p,'Head_pop']}"
         for p in pairs_sorted],
        fontsize=5.2, rotation=0, ha='center')
    ax.set_ylabel(r'Environmental divergence $\Delta E$' + '\n(standardised Euclidean)',
                  fontsize=7.0)
    ax.set_title(r'$\mathbf{c}$  Environmental divergence  (proxy: $\Delta H$)',
                 fontsize=8, loc='left', pad=3)
    ax.tick_params(labelsize=6.5)
    ax.set_ylim(0, merged['Env_dist'].max() * 1.25)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6); sp.set_color('#c0c0c0')
    ax.grid(True, axis='y', color='0.93', lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    # ── ΔpH and ΔTOC inset annotation ──────────────────────────────────────
    ax.legend(handles=[
        mpatches.Patch(fc=BLUE,   ec='white', label='Eastern clade'),
        mpatches.Patch(fc=ORANGE, ec='white', label='Southeastern clade'),
    ], fontsize=6.0, loc='upper right', framealpha=0.90, edgecolor='0.82')

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL d — Scatter: Env divergence vs Phenotypic divergence
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_scat

    for pid in merged.index:
        row   = merged.loc[pid]
        col   = CLADE_COL[row['Clade']]
        mkr   = CLADE_MRK[row['Clade']]
        ax.scatter(row['Env_dist'], row['Pheno_dist'],
                   color=col, marker=mkr, s=70, zorder=5,
                   edgecolors='white', linewidths=0.5, alpha=0.90)
        ax.text(row['Env_dist'] + 0.08, row['Pheno_dist'] + 0.05,
                f"Pair {pid}\n({row['Dune_pop']}/{row['Head_pop']})",
                ha='left', va='bottom', fontsize=5.5, color='0.30')

    # Pearson r
    r, pval = pearsonr(merged['Env_dist'], merged['Pheno_dist'])
    # Regression line
    x_fit = np.linspace(merged['Env_dist'].min()*0.9,
                        merged['Env_dist'].max()*1.05, 100)
    b1, b0 = np.polyfit(merged['Env_dist'], merged['Pheno_dist'], 1)
    ax.plot(x_fit, b0 + b1*x_fit, '--', color='0.55', lw=1.3, zorder=2)

    # Pearson r annotation — bottom-left, away from legend
    pstr = f"p = {pval:.3f}" if pval >= 0.001 else "p < 0.001"
    ax.text(0.04, 0.08,
            rf'$r = {r:.2f}$,  {pstr}  ($n = 8$ pairs)',
            transform=ax.transAxes,
            ha='left', va='bottom', fontsize=7.5, color='0.20',
            bbox=dict(fc='white', ec='0.80', lw=0.6, pad=2.5,
                      boxstyle='round,pad=0.3'))

    # Prediction zone annotation — bottom-right
    ax.text(0.97, 0.04,
            'Large ΔE, large Δμ\n→ predicted high ' + r'$I_{\rm geom}$',
            transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6.5, color=BLUE,
            style='italic')

    ax.set_xlabel(r'Environmental divergence $\Delta E$  (proxy for $\Delta H$)',
                  fontsize=8)
    ax.set_ylabel(r'Phenotypic divergence $\Delta\mu$  (proxy for $\theta_{12}$)',
                  fontsize=8)
    ax.set_title(
        r'$\mathbf{d}$  Environmental vs. phenotypic divergence per pair'
        '\n' + r'(prediction: $\Delta E \propto \Delta\mu \propto \theta_{12}$'
        r'$\Rightarrow I_{\rm geom} \approx k\,\theta_{12}^2$)',
        fontsize=8, loc='left', pad=3)
    ax.tick_params(labelsize=7)
    ax.legend(handles=[
        mpatches.Patch(fc=BLUE,   ec='white', label='Eastern clade'),
        mpatches.Patch(fc=ORANGE, ec='white', label='Southeastern clade'),
    ], fontsize=6.5, loc='upper right',
       framealpha=0.90, edgecolor='0.82')
    for sp in ax.spines.values():
        sp.set_linewidth(0.6); sp.set_color('#c0c0c0')
    ax.grid(True, color='0.93', lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    # ── Save ──────────────────────────────────────────────────────────────────
    plt.savefig(outfile, format='pdf', dpi=300,
                bbox_inches='tight', facecolor='white')
    print(f"  Saved: {outfile}")


if __name__ == "__main__":
    # Re-run analysis to get merged DataFrame in scope
    import subprocess, sys
    subprocess.run([sys.executable, '/home/claude/senecio_analysis.py'], check=True)
    make_figure("senecio_pair_analysis.pdf")
    shutil.copy("senecio_pair_analysis.pdf",
                "/mnt/user-data/outputs/senecio_pair_analysis.pdf")
    shutil.copy(__file__,
                "/mnt/user-data/outputs/senecio_figure.py")
    shutil.copy("/home/claude/senecio_analysis.py",
                "/mnt/user-data/outputs/senecio_analysis.py")
    print("Done.")
