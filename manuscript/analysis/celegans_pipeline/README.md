# Caenorhabditis elegans channel-operator posterior analysis

**A reproducible pipeline for the worm section of**
*"The Geometric Snowball of Speciation"* (Ortiz-Barrientos, James & Burrage, in prep.)

This repository re-analyses the multivariate quantitative-genetic data of
Mallard, Afonso & Teotônio (*eLife* **12**, e80993, 2023) using the
additive-channels framework of Ortiz-Barrientos & Cooper
(in prep.). The goal is to test, with full posterior uncertainty, whether
sustained directional selection erodes the spectral gap of the **G**
matrix and destabilises the leading response channel defined by
**S̃** = **G**^1/2^ **H̃ G**^1/2^.

---

## What the analysis does, in plain language

Four populations of *Caenorhabditis elegans* are available: a single
**ancestor** (A6140) and three **evolved replicates** that were adapted
to high salt for 50 generations (GA150, GA250, GA450). For every
population, Mallard et al. measured 10,000 MCMC posterior draws of the
7×7 additive genetic variance matrix **G** covering six locomotion
transitions and body size. They also reported a standardised selection
gradient **β** — the direction in trait space that fitness pulls on.

The additive-channels framework says that adaptation is steered not
by **G** alone but by the symmetric operator

> **S̃** = **G**^1/2^ **H̃ G**^1/2^

where **H̃** is the local fitness curvature. Because **H̃** is rarely
measured directly, we build a rank-1 surrogate from **β**:

> **H̃** = **ββ**^⊤^ / ‖**β**‖^2^ + ε **I**

and track two quantities across every posterior draw:

1. **Channel anisotropy of G**: the eigenvalue ratio λ₁ / λ₂.
   A large ratio means genetic variation concentrates along a single
   direction; a ratio close to 1 means two directions are nearly equal
   and the leading eigenvector is unstable.

2. **Channel rotation**: the angle between **g**_max (leading eigenvector
   of **G**) and the leading eigenvector of **S̃** pulled back to
   phenotype space. A small angle means fitness curvature reinforces
   genetic structure (Schluter's line of least resistance); a large
   angle means selection pulls the response direction away from
   **g**_max.

The Davis–Kahan theorem predicts that rotation becomes unstable as the
spectral gap shrinks. So the prediction is: evolved replicates that
have burned through their spectral anisotropy should show poorly
identified rotations — wide credible intervals — not because the data
are noisy, but because the leading eigenvector itself is no longer well
defined.

The pipeline runs every MCMC draw through this calculation, extracts
posterior median + 95% credible interval for each quantity in each
population, and produces a publication-ready figure and table.

---

## Repository layout

```
celegans_pipeline/
├── run_pipeline.sh              ← orchestrator: runs the full pipeline end-to-end
├── requirements.txt             ← Python dependencies
├── requirements.R               ← R dependencies (none beyond base R 4.x)
├── README.md                    ← this file
├── LICENSE
├── .gitignore
│
├── data/
│   ├── raw/                     ← PRIMARY DATA from Mallard et al. (2023)
│   │   ├── G_matrices_high_salt.RData        ← MCMC chains, high salt  (44 MB)
│   │   ├── G_matrices_low_salt.RData         ← MCMC chains, low salt   (44 MB)
│   │   ├── Combined_Selection_Tables.txt     ← β tables and selection differentials
│   │   ├── *_posterior_mode.txt               ← point-estimate matrices (sanity check)
│   │   └── Phenotypic_data.RData             ← individual-level phenotypes
│   │
│   └── processed/               ← derived, regenerable files
│       ├── beta_ancestor_highsalt_standardized.csv
│       ├── beta_G_evolved_populations.csv
│       ├── selection_differentials.csv
│       ├── G_matrices_posterior_mode_long.csv
│       ├── G_chains_high_salt_long.csv       (~70 MB)
│       ├── G_chains_low_salt_long.csv        (~70 MB)
│       └── celegans_posterior_perdraw.csv    (per-draw geometry for every pop)
│
├── scripts/
│   ├── 01_parse/
│   │   ├── 01a_parse_selection_tables.py     ← split the stacked β tables
│   │   ├── 01b_parse_posterior_mode_matrices.py  ← read the 9 point-estimate matrices
│   │   └── 01c_extract_posterior_chains.R    ← read MCMC chains from RData to CSV
│   │
│   ├── 02_analysis/
│   │   ├── geometry.py          ← core math: matsqrt, channel geometry, summaries
│   │   └── 02a_posterior_analysis.py          ← propagate posterior through S̃
│   │
│   └── 03_figures/
│       ├── 03a_figure_posterior_violins.py    ← the two-panel violin figure
│       └── 03b_table_celegans_summary.py      ← Word + LaTeX summary table
│
├── results/
│   ├── figures/
│   │   ├── celegans_posterior_violins.pdf
│   │   └── celegans_posterior_violins.png
│   │
│   └── tables/
│       ├── celegans_posterior_summary.csv     ← 40-row tidy summary (all pops × envs)
│       ├── table_celegans_summary.docx        ← publication-ready Word table
│       └── table_celegans_summary.tex         ← LaTeX version for Overleaf
│
└── docs/
    └── (space for supplementary notes)
```

---

## Quick start

```bash
# 1. Clone the repo and enter it
git clone https://github.com/dortizbarrientos/celegans-channel-operator.git
cd celegans-channel-operator

# 2. Install dependencies
pip install -r requirements.txt             # Python side
# R side: only base R is needed. No CRAN packages required.

# 3. Run the full pipeline (takes ~2 minutes on a laptop)
bash run_pipeline.sh
```

The outputs appear in `results/figures/` and `results/tables/`.

### A note about Git LFS

The two posterior-chain files in `data/raw/` are about 44 MB each —
below GitHub's 100 MB hard limit but above the 50 MB warning
threshold. For a clean push experience, initialise Git LFS **before**
your first commit:

```bash
git lfs install
git lfs track "data/raw/*.RData"
git add .gitattributes
```

This keeps the working repository small while still versioning the
raw data. If you choose not to use LFS, the repo will still work; you
will just see GitHub warnings on push.

---

## Software requirements

| Software | Version tested | Purpose |
|---|---|---|
| Python | ≥ 3.10 | Main analysis, figures, tables |
| R      | ≥ 4.0 | Only for reading the MCMCglmm S4 objects (one script) |
| numpy  | ≥ 1.24 | Linear algebra |
| pandas | ≥ 2.0  | Data wrangling |
| matplotlib | ≥ 3.7 | Figure rendering |
| python-docx | ≥ 1.0 | Word table generation |

No R packages beyond base R are required. No MCMCglmm install is
needed; we only read the already-fitted posterior chains from the
RData files supplied by Mallard et al.

---

## Pipeline stages in detail

### Stage 01: parsing

**01a — `scripts/01_parse/01a_parse_selection_tables.py`**

Parses the file `Combined_Selection_Tables.txt`, which contains three
stacked tables with different column structures:

| Table | Columns | Content |
|---|---|---|
| 1 | 6 | Standardised β (ancestor, high salt) — posterior mode and CIs |
| 2 | 17 | Evolved-population CIs and "Beta G" (genetic-regression β) |
| 3 | 7 | Selection differentials by environment |

The script auto-detects table boundaries by column count and splits
them into three tidy CSVs under `data/processed/`. The β used in every
downstream calculation is the Post.mode column of Table 1 — the
standardised phenotypic selection gradient — because that is the β that
enters Lande's equation Δ**z̄** = **Gβ**.

**01b — `scripts/01_parse/01b_parse_posterior_mode_matrices.py`**

Parses the nine 7×7 text files (four populations × two environments,
plus MA lines) each containing "value  [lower;upper]" cells. Produces
a long-format CSV with 441 rows (9 matrices × 49 cells). These are
point estimates only — used as a sanity check against the posterior
modes we compute from the MCMC chains later. The chains are the
primary data.

**01c — `scripts/01_parse/01c_extract_posterior_chains.R`**

The only R script in the pipeline. Reads the two `.RData` files:

- `G_matrices_high_salt.RData` → object `VCV_mat_NaCl` (list of 4 populations)
- `G_matrices_low_salt.RData`  → object `VCV_mat_NGM`  (same structure)

Each population entry contains a `VCV_Mat` slot that is an `mcmc`
object with **10,000 rows × 147 columns**. The 147 columns decompose as:

| Columns | Variance component | Used here? |
|---|---|---|
| 1–49   | `pop_label` — additive genetic G | **Yes** (this is what we analyse) |
| 50–98  | `date_str` — date-block random effect | No |
| 99–147 | `units` — residual covariance | No |

The script exports columns 1–49 for every (population, draw) as a
long-format CSV. After this step, the rest of the pipeline has no R
dependency.

### Stage 02: analysis

**02a — `scripts/02_analysis/02a_posterior_analysis.py`**

Main analysis. For each environment × population × draw:

1. Reshape the 49-element row into a 7×7 symmetric matrix **G**.
2. Eigendecompose **G**; record λ₁, λ₂, ratio λ₁/λ₂, gap δ = λ₁ − λ₂.
3. Compute **G**^1/2^ from the same decomposition.
4. Build the rank-1 surrogate **H̃** from **β** (held fixed across all
   populations — a conservative choice that forces any detected rotation
   to arise from divergence in **G** alone).
5. Compute **S̃** = **G**^1/2^ **H̃ G**^1/2^.
6. Eigendecompose **S̃**, pull the leading eigenvector back through
   **G**^1/2^, and measure the acute angle to **g**_max.

Outputs:

- `data/processed/celegans_posterior_perdraw.csv` (80,000 rows:
  4 pops × 2 envs × 10,000 draws)
- `results/tables/celegans_posterior_summary.csv` (40 rows: posterior
  median + 95% CI per pop × env × quantity)

The core math sits in `scripts/02_analysis/geometry.py`:

- `matrix_sqrt_sym(M)` — symmetric positive-definite square root
- `acute_angle_deg(u, v)` — acute angle in degrees (invariant to sign)
- `compute_channel_geometry(G, beta, eps)` — all per-draw quantities
- `posterior_summary(x)` — median + 95% CI with NaN handling

### Stage 03: outputs

**03a — `scripts/03_figures/03a_figure_posterior_violins.py`**

Two-panel violin plot for the high-salt environment:

- Panel (a): log₁₀(λ₁/λ₂) per population, with the spectral-gap
  boundary marked at 1.
- Panel (b): channel rotation in degrees, with the orthogonal limit
  marked at 90°.

Each violin shows the full posterior; the overlaid point is the median,
the bar is the 95% credible interval. Okabe-Ito palette, Arial fonts,
183 mm wide to fit a Nature Brief Communication double-column figure.
Both PDF and PNG are emitted.

**03b — `scripts/03_figures/03b_table_celegans_summary.py`**

Produces two versions of a manuscript-ready summary table:

- `table_celegans_summary.docx` (python-docx, styled header row, alternate
  row shading)
- `table_celegans_summary.tex` (booktabs format, for Overleaf)

Both include the percent gap erosion relative to the ancestor — a
direct quantification of the prediction.

---

## What the results show

Running the pipeline on the published Mallard et al. data produces the
following high-salt summary:

| Population | *n* | λ₁/λ₂ | Spectral gap δ | Erosion (%) | Rotation (deg) |
|---|---|---|---|---|---|
| A6140 (ancestor) | 362 | 7.22 [3.98, 12.22] | 0.819 [0.441, 1.306] | — | 6.5 [3.7, 12.1] |
| GA150 | 118 | 1.76 [1.12, 3.57] | 0.134 [0.026, 0.364] | 84 | 27.8 [9.1, 84.6] |
| GA250 | 121 | 1.79 [1.12, 3.45] | 0.118 [0.022, 0.301] | 86 | 32.8 [9.7, 85.5] |
| GA450 | 83  | 3.20 [1.44, 7.42] | 0.418 [0.119, 0.936] | 49 | 15.2 [5.8, 47.7] |

The ancestor sits well above the spectral-gap boundary. Two of three
evolved replicates (GA150, GA250) collapsed toward the boundary, with
their 95% credible intervals touching 1.12 — very close to the
channel-degeneracy threshold. The third replicate (GA450) remained at
an intermediate ratio. In the two boundary-approaching replicates, the
rotation 95% CIs span [9°, 85°] — almost the entire possible range.
This wide interval is not a bug: it is the Davis–Kahan bound made
flesh. When the spectral gap collapses, the leading eigenvector itself
becomes ill-defined, and the posterior correctly reports that it
cannot be pinned down from data. **The unidentifiability is the
prediction.**

The low-salt environment (NGM) shows less dramatic collapse — expected,
because these populations were selected in high salt, not low salt.
Both environments are analysed; the manuscript focuses on high salt.

---

## Reproducibility notes

- **Determinism.** The pipeline is deterministic given the raw data.
  There is no random seed, no bootstrap. Every number in the results
  is computed from the posterior chains analytically.
- **Incremental runs.** Each stage is idempotent: re-running overwrites
  the same files. The R chain extraction (stage 01c) is the slowest —
  about 30 s per environment. Comment out a stage in
  `run_pipeline.sh` to skip it if already complete.
- **Memory.** The per-draw analysis peaks at about 1 GB when the
  high-salt chain CSV is loaded into pandas.
- **File sizes.** `data/processed/G_chains_*.csv` are ~70 MB each.
  They are regenerable from the raw `.RData` files and are excluded by
  `.gitignore` — if you clone this repo you regenerate them by running
  `run_pipeline.sh`.

---

## Data provenance and citation

The raw data in `data/raw/` are republished here for reproducibility,
with thanks to Mallard, Afonso & Teotônio. They are available in their
original form at the Dryad repository associated with the *eLife* paper.

If you use this pipeline or its results, please cite:

> Mallard, F., Afonso, B. & Teotônio, H. (2023). Selection and the
> direction of phenotypic evolution. *eLife* **12**, e80993.

and, for the additive-channels framework:

> Ortiz-Barrientos, D. & Cooper, M. (in prep.). Additive channels in
> curved fitness landscapes.
>
> Ortiz-Barrientos, D., James, M. E. & Burrage, K. (in prep.). The
> geometric snowball of speciation.

---

## License

This pipeline is released under the MIT License. See `LICENSE`.
The raw data retain the license under which Mallard et al. originally
released them.

---

## Contact

Daniel Ortiz-Barrientos — d.ortizbarrientos@uq.edu.au
School of the Environment, The University of Queensland
