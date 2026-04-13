"""
senecio_pair_analysis.py
========================
Analyses of Senecio lautus Dune–Headland pairs using available
soil chemistry and leaf morphology data.

CONTEXT
-------
These analyses are preparation for computing the geometric snowball
of speciation (I_geom = 1 − cos^(2k)(θ₁₂)) once G and H matrices
become available from the DH breeding design.

What we can compute NOW:
  1. Phenotypic divergence (Δμ) between paired populations — proxies
     for the degree of adaptive divergence; expected to correlate with θ₁₂.
  2. Environmental divergence (ΔE) from soil chemistry — proxies for
     the difference in fitness landscapes (ΔH); pairs with large ΔE
     are expected to have large I_geom.
  3. A multivariate ordination of all populations to visualise the
     structure of the radiation.
  4. Pair-level summary table ready to receive I_geom once G/H data arrive.

POPULATION MAPPING (soil code → leaf code)
-------------------------------------------
Soil pair 1 : D00→DN0,  H00→HN0   [eastern]
Soil pair 2 : D01→DN1,  H01→HN1   [eastern]
Soil pair 3 : D03→DN3,  H02→HN2   [eastern]
Soil pair 4 : D04→DN4,  H05→HN5   [eastern]
Soil pair 5 : D12→DN12, H14→HN14  [eastern]
Soil pair 6 : D23→DN23, H21→HS21  [southeastern]
Soil pair 7 : D32→DN32, H12→HN12  [southeastern]
Soil pair 8 : D14→DN14, H15→HN15  [southeastern]
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.spatial.distance import mahalanobis
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import shutil, warnings
warnings.filterwarnings("ignore")

# ── Colour palette (Okabe-Ito) ────────────────────────────────────────────────
BLUE    = "#0072B2"
ORANGE  = "#E69F00"
VERM    = "#D55E00"
GREEN   = "#009E73"
SKY     = "#56B4E9"
YELLOW  = "#F0E442"
PURPLE  = "#CC79A7"

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
    "pdf.fonttype": 42,
    "axes.linewidth": 0.7,
})

# ── Load data ─────────────────────────────────────────────────────────────────
leaf = pd.read_excel('/mnt/user-data/uploads/data_AL_Arch_Leaf_ARVersion.xlsx')
soil_all = pd.read_excel(
    '/mnt/user-data/uploads/All_S__lautus_soil_data.xlsx',
    sheet_name='All soil data')

# ── Population-level leaf trait means ─────────────────────────────────────────
# Traits chosen for:
#   (a) biological relevance to dune-headland adaptation
#   (b) no missing data in either ecotype group
TRAITS = ['TotalHeight', 'MSD', 'MSL_W', 'Area',
          'Circularity', 'Nindents', 'IndentWidth', 'IndentDepth', 'SB']

TRAIT_LABELS = {
    'TotalHeight' : 'Plant height (mm)',
    'MSD'         : 'Max. stem diam. (mm)',
    'MSL_W'       : 'Leaf width ratio',
    'Area'        : 'Leaf area (mm²)',
    'Circularity' : 'Leaf circularity',
    'Nindents'    : 'No. leaf indents',
    'IndentWidth' : 'Indent width (mm)',
    'IndentDepth' : 'Indent depth (mm)',
    'SB'          : 'Stem branches',
}

# Compute population means
pop_means = (leaf.dropna(subset=TRAITS)
             .groupby(['Ecotype','Population'])[TRAITS]
             .mean())

# ── Pair mapping ──────────────────────────────────────────────────────────────
# Each tuple: (pair_id, clade, dune_pop, head_pop)
PAIRS = [
    (1, 'eastern',      'DN0',  'HN0'),
    (2, 'eastern',      'DN1',  'HN1'),
    (3, 'eastern',      'DN3',  'HN2'),
    (4, 'eastern',      'DN4',  'HN5'),
    (5, 'eastern',      'DN12', 'HN14'),
    (6, 'southeastern', 'DN23', 'HS21'),
    (7, 'southeastern', 'DN32', 'HN12'),
    (8, 'southeastern', 'DN14', 'HN15'),
]

# ── Soil chemistry variables ──────────────────────────────────────────────────
# Selected for ecological relevance and minimal missingness
SOIL_VARS = [
    'pH Value', 'Electrical Conductivity @ 25_C',
    'Total Organic Carbon', 'Cation Exchange Capacity',
    'Phosphorus', 'Total Nitrogen as N',
    'Exchangeable Calcium', 'Exchangeable Magnesium',
]

SOIL_LABELS = {
    'pH Value'                        : 'pH',
    'Electrical Conductivity @ 25_C'  : 'Elec. conduct.',
    'Total Organic Carbon'            : 'Organic C',
    'Cation Exchange Capacity'        : 'CEC',
    'Phosphorus'                      : 'Phosphorus',
    'Total Nitrogen as N'             : 'Total N',
    'Exchangeable Calcium'            : 'Exch. Ca',
    'Exchangeable Magnesium'          : 'Exch. Mg',
}

# Soil pair mapping (soil pop code → pair)
SOIL_PAIR_MAP = {
    ('D00','H00'): 1, ('D01','H01'): 2, ('D03','H02'): 3,
    ('D04','H05'): 4, ('D12','H14'): 5, ('D23','H21'): 6,
    ('D32','H12'): 7, ('D14','H15'): 8,
}

# ── Step 1: Phenotypic divergence per pair ────────────────────────────────────
# Δμ = Euclidean distance in standardised trait space.
# Standardisation uses the grand SD across all Dune+Headland individuals,
# so each trait contributes on the same scale regardless of units.

dh_leaf = leaf[leaf['Ecotype'].isin(
    ['D.Dune.N', 'G.Head_N', 'H.Head_S'])].dropna(subset=TRAITS)

grand_mean = dh_leaf[TRAITS].mean()
grand_std  = dh_leaf[TRAITS].std()

def std_means(pop):
    """Standardised trait vector for one population (numpy array)."""
    mask = pop_means.index.get_level_values('Population') == pop
    if mask.sum() == 0:
        return None
    raw = pop_means[mask][TRAITS].values[0]   # numpy array
    return (raw - grand_mean.values) / grand_std.values

pheno_rows = []
for pid, clade, dp, hp in PAIRS:
    mv_d = std_means(dp)
    mv_h = std_means(hp)
    if mv_d is None or mv_h is None:
        print(f"  Pair {pid}: missing pop ({dp} or {hp})")
        continue
    delta = mv_d - mv_h
    dist  = np.linalg.norm(delta)
    # delta is a numpy array; index by integer position
    pheno_rows.append({
        'Pair'            : pid,
        'Clade'           : clade,
        'Dune_pop'        : dp,
        'Head_pop'        : hp,
        'Pheno_dist'      : dist,
        'Delta_Height'    : delta[TRAITS.index('TotalHeight')],
        'Delta_MSD'       : delta[TRAITS.index('MSD')],
        'Delta_Area'      : delta[TRAITS.index('Area')],
        'Delta_Circ'      : delta[TRAITS.index('Circularity')],
        'Delta_Nindents'  : delta[TRAITS.index('Nindents')],
    })

pheno_df = pd.DataFrame(pheno_rows).set_index('Pair')
print("Phenotypic divergence per pair:\n", pheno_df[['Clade','Pheno_dist']].round(3))

# ── Step 2: Environmental divergence per pair ─────────────────────────────────
# ΔE = Euclidean distance in standardised soil-chemistry space.
# This proxies ΔH (the difference in fitness landscapes) in the absence of H matrices.

soil_dune = soil_all[soil_all['Ecotype'] == 'Dune'].set_index('Population')
soil_head = soil_all[soil_all['Ecotype'] == 'Headland'].set_index('Population')

# Standardise across all soil samples
soil_matrix = soil_all[soil_all['Ecotype'].isin(
    ['Dune','Headland'])][SOIL_VARS].apply(pd.to_numeric, errors='coerce')
soil_std = soil_matrix.std()
soil_mean = soil_matrix.mean()

env_rows = []
for pid, clade, dp, hp in PAIRS:
    # Retrieve soil population codes from the inverse of the pair map
    soil_d_code = f"D{dp[2:]}"   # DN0→D00, DN1→D01 etc
    # Fix: leaf DN0→soil D00 mapping needs zero-padding
    # DN0→D00, DN1→D01, DN3→D03, DN4→D04, DN12→D12
    # HN0→H00, HN1→H01, HN2→H02, HN5→H05, HN14→H14
    # HS21→H21, HN12→H12, HN15→H15, DN23→D23, DN32→D32, DN14→D14
    lp_to_soil = {
        'DN0':'D00','DN1':'D01','DN3':'D03','DN4':'D04','DN12':'D12',
        'DN23':'D23','DN32':'D32','DN14':'D14',
        'HN0':'H00','HN1':'H01','HN2':'H02','HN5':'H05','HN14':'H14',
        'HS21':'H21','HN12':'H12','HN15':'H15',
    }
    sc_d = lp_to_soil.get(dp)
    sc_h = lp_to_soil.get(hp)
    if sc_d not in soil_dune.index or sc_h not in soil_head.index:
        print(f"  Pair {pid}: soil data missing for {sc_d} or {sc_h}")
        continue
    sv_d = pd.to_numeric(soil_dune.loc[sc_d, SOIL_VARS], errors='coerce')
    sv_h = pd.to_numeric(soil_head.loc[sc_h, SOIL_VARS], errors='coerce')
    delta_e = (sv_d - sv_h) / soil_std
    env_dist = np.linalg.norm(delta_e.dropna())
    delta_pH = float(sv_d['pH Value']) - float(sv_h['pH Value'])
    delta_TOC = float(sv_d['Total Organic Carbon']) - float(sv_h['Total Organic Carbon'])
    delta_N   = float(sv_d['Total Nitrogen as N'])  - float(sv_h['Total Nitrogen as N'])
    env_rows.append({
        'Pair'       : pid,
        'Clade'      : clade,
        'Env_dist'   : env_dist,
        'Delta_pH'   : delta_pH,
        'Delta_TOC'  : delta_TOC,
        'Delta_N'    : delta_N,
    })

env_df = pd.DataFrame(env_rows).set_index('Pair')
print("\nEnvironmental divergence per pair:\n", env_df[['Clade','Env_dist','Delta_pH','Delta_TOC']].round(3))

# ── Step 3: Merge and compute proxy I_geom ranking ───────────────────────────
# Without G and H matrices we cannot compute I_geom exactly.
# However, under the geometric snowball model:
#   θ₁₂ ∝ Δμ (phenotypic divergence drives channel rotation)
#   A ∝ ΔE  (environmental contrast drives the amplitude of RF)
#   I_geom ≈ k·θ₁₂²  →  proxy I_geom ∝ Pheno_dist²
# A combined proxy: Proxy_I = (Pheno_dist / max(Pheno_dist))²
# This ranks pairs by predicted isolation under the geometric snowball.

merged = pheno_df.join(env_df[['Env_dist','Delta_pH','Delta_TOC','Delta_N']])
merged['Proxy_I'] = (merged['Pheno_dist'] / merged['Pheno_dist'].max())**2
merged = merged.sort_values('Pheno_dist', ascending=False)

print("\nMerged pair summary (ranked by phenotypic divergence):")
print(merged[['Clade','Dune_pop','Head_pop','Pheno_dist',
              'Env_dist','Delta_pH','Delta_TOC','Proxy_I']].round(3))

# ── Step 4: PCA of leaf morphology across ALL ecotypes ───────────────────────
# All populations, all ecotypes — to show where Dune and Headland sit
# within the full S. lautus radiation.
all_pops = leaf.dropna(subset=TRAITS).groupby(['Ecotype','Population'])[TRAITS].mean()
scaler   = StandardScaler()
X        = scaler.fit_transform(all_pops[TRAITS])
pca      = PCA(n_components=2)
pcs      = pca.fit_transform(X)
pca_df   = pd.DataFrame(pcs, columns=['PC1','PC2'],
                         index=all_pops.index)
pca_df['Ecotype']    = pca_df.index.get_level_values('Ecotype')
pca_df['Population'] = pca_df.index.get_level_values('Population')

var_expl = pca.explained_variance_ratio_ * 100

# Return merged and pca_df for the figure
print(f"\nPCA: PC1={var_expl[0]:.1f}%, PC2={var_expl[1]:.1f}%")

# Save intermediates for figure script
merged.to_csv('/home/claude/pair_summary.csv')
pca_df.to_csv('/home/claude/pca_pops.csv')
all_pops.to_csv('/home/claude/pop_means.csv')

# Store pca loadings
loadings = pd.DataFrame(pca.components_.T, index=TRAITS,
                        columns=['PC1_load','PC2_load'])
loadings.to_csv('/home/claude/pca_loadings.csv')
print("Loadings:\n", loadings.round(3))
print("\nVariance explained:", var_expl.round(2))
print("Analysis complete.")
