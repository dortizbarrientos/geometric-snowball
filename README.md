# The Geometric Snowball of Speciation

**Daniel Ortiz-Barrientos, Maddie James, and Kevin Burrage**  
The University of Queensland, Queensland University of Technology, and ARC CoE for Plant Success

---

Code and manuscript files for the paper:

> **The Geometric Snowball of Speciation**  
> Daniel Ortiz-Barrientos, Maddie E. James, and Kevin Burrage, in preparation.

When two populations adapt to different environments, their preferred directions of genetic response rotate apart, and hybrids inherit a mixture that suits neither environment. The resulting fitness cost grows quadratically with the angle between parental directions — a geometric snowball that starts at the first moment of ecological divergence and requires no new mutations.

---

## Repository structure

```
geometric-snowball-main/
│
├── manuscript/
│   ├── word/                               Manuscript drafts (current: snowball_ms_v9.docx)
│   ├── analysis/
│   │   ├── celegans_pipeline/              C. elegans reproducible pipeline
│   │   ├── senecio_pipeline/               Senecio reproducible pipeline
│   │   ├── otherScripts/                   Exploratory analysis scripts
│   │   └── README.md                       Analysis-level documentation
│   ├── figures/
│   │   ├── forKeynote/                     Individual PDF panels for compositing
│   │   ├── legends/                        Figure legend drafts
│   │   └── retired/                        Previous figure versions
│   ├── tables/
│   │   └── all_tables.docx                 Formatted manuscript tables
│   ├── abstracts/
│   │   └── abstract_options.docx           Abstract drafts
│   ├── supporting/
│   │   ├── supportingMethods.pdf           Supplementary methods
│   │   ├── Math_Boot_Camp.pdf              Mathematical appendix (61 pages)
│   │   └── channel_nonlinearity_tutorial.html   Interactive tutorial (open in browser)
│   └── playground/
│       ├── coupled_ratchet_sweep.py        Coupled-ratchet simulations
│       └── README_coupled_ratchet_sweep.md
│
├── data/                               Published datasets organised by source paper
│   ├── amNat/                              Walter et al. 2018 Am Nat (G-matrices, D)
│   ├── evolLetters/                        Walter et al. 2020 Evol Letters (field transplant)
│   ├── evol/                               Walter et al. 2016 Evolution (seed set, F4)
│   ├── evol2/                              Richards et al. (multi-generation transplant)
│   ├── elifeCelegans/                      Mallard et al. 2023 eLife (C. elegans G + β)
│   ├── ecol/                               Walter et al. 2018 Ecology
│   ├── newPhyt/                            Melo et al. 2014 New Phytologist (extrinsic RI)
│   ├── phyloandGenomics/                   James et al. 2021, 2025 (RI crosses, genomics)
│   ├── evolLetters-Ken/                    Thompson et al. (stickleback)
│   ├── sticklebacks/                       Schluter lab supplementary data
│   ├── refsDatasets/                       Reference datasets (Silene, rotation)
│   ├── All S. lautus soil data.xlsx        Soil chemistry at population sites
│   └── data_AL_Arch&Leaf_ARVersion.xlsx    Architecture and leaf data
│
├── extension/                          Extended theoretical materials
│   ├── math_bootcamp_appendix.tex/pdf      Full mathematical appendix (source + PDF)
│   ├── gph_tutorial_for_biologists.tex/pdf         Geometric tutorial (short version)
│   ├── gph_tutorial_for_biologists_detailed.tex/pdf  Geometric tutorial (long version)
│   ├── gph_landscape_manuscript_draft.tex/pdf      Extended landscape theory
│   ├── gph_math_figures.py                 Python code for math figures
│   ├── gph_math_figure_pack.pdf            Compiled math figures
│   └── gph_math_figure_notes.md            Notes on figure construction
│
├── playground/                         Exploratory work (not part of submission)
│   ├── sympatric_coupled_ratchet_sweep.py
│   └── section_revision_hybrid_cost.docx
│
├── LICENSE
└── README.md                           You are here
```

---

## Reproducing the analyses

The paper's empirical results come from two reproducible pipelines, one per study system. Each pipeline is self-contained: raw data, scripts, a shell orchestrator, and documentation. A single command runs everything end-to-end.

### Requirements

**R** (≥ 4.0) with packages:

```bash
Rscript -e 'install.packages(c("MCMCglmm", "coda", "readxl", "dplyr", "tidyr", "ggplot2", "expm", "Matrix"))'
```

**Python** (≥ 3.9) with packages:

```bash
pip install numpy pandas matplotlib xlrd openpyxl
```

**Git LFS** (for large RData files in the C. elegans pipeline only):

```bash
brew install git-lfs    # macOS; see https://git-lfs.com for other platforms
git lfs install
```

### Senecio pipeline

Located at `manuscript/analysis/senecio_pipeline/`. Estimates population-specific **G**-matrices for eight *S. lautus* populations, constructs the channel operator **S** = **G**^{1/2} **H** **G**^{1/2}, computes pairwise channel angles, and analyses the multigenerational transplant data.

```bash
cd manuscript/analysis/senecio_pipeline/

# Install R packages (one-time)
Rscript requirements.R

# Quick test — short MCMC chains, validates everything works (~5 minutes)
bash run_pipeline.sh --quick-test

# Full run — production MCMC chains (several hours for Stage 1)
bash run_pipeline.sh
```

Four stages, runnable individually with `--stage N`:

| Stage | What it does | Runtime |
|-------|--------------|---------|
| **1** | Fits population-specific **G** via MCMCglmm (8 populations × 2 trait blocks: 4×4 vegetative + 2×2 leaf PCs). Computes pairwise **D** from common-garden means and pairwise **G** geometry. | Hours (full) / ~5 min (quick) |
| **2** | Fits Lande–Arnold selection gradients (β, γ) per population × environment. Constructs **H** from quadratic gradients. Computes **S** = **G**^{1/2} **H** **G**^{1/2} and extracts channel angles. Cross-validates with soil-based rank-1 **H** proxy. | ~2 min |
| **3** | Produces summary figure: **G** anisotropy, pairwise g_max angles, and **S** channel angles across all populations and environments. | Seconds |
| **4** | Analyses the multigenerational transplant data (Richards et al.) and seed-set trajectory (Walter et al. 2020). Computes F1 heterosis, F2 breakdown, extrinsic RI per generation, and backcross recovery. | Seconds |

Outputs land in `data/processed/` (Stages 1–2) and `results/` (Stages 3–4). See the pipeline's own `README.md` and `docs/expected_outputs.md` for validated reference numbers.

**Tracing paper claims to output files:**

| Paper claim | Output file | Column |
|---|---|---|
| G-matrix anisotropy per population (Fig. 5a) | `leaf_PCs_G_summary.csv` | `lambda1_over_lambda2` |
| Pairwise g_max angles, Schluter limit (Fig. 5b) | `pairwise_G_geometry_leafPC.csv` | `gmax_angle_deg` |
| S channel angles, full operator (Fig. 5c) | `S_channel_angles_LA.csv` | `S_channel_angle_deg` |
| θ explains 87% of RI variance | Computed from `CrossesDH-ReciprocalsAveraged.txt` | `RIseedSet-SobelChen` vs θ |
| F1 heterosis (+39% vs midparent) | `transplant_multigenerational_summary.csv` | `L10_mean` |
| F2 breakdown (−29% vs F1) | same | same |
| RI sign reversal F1→F2 in Dune env | `transplant_extrinsic_RI.csv` | `RI_extrinsic` |
| Seed-set F2 crash (0.61 → 0.12) | `seedset_trajectory.csv` | `mean_PropFilled` |
| Backcross recovery (83% of parental) | `transplant_fitness_by_generation.csv` | `L10` for BCHD × HEAD |

### C. elegans pipeline

Located at `manuscript/analysis/celegans_pipeline/`. Tests the spectral-gap prediction using Mallard et al. (2023) data: 50 generations of adaptation to high salt in one ancestral and three evolved replicate populations.

```bash
cd manuscript/analysis/celegans_pipeline/

# Full run — deterministic, ~2 minutes
bash run_pipeline.sh
```

The pipeline parses posterior **G**-matrices from RData files, constructs a rank-1 surrogate **S** from the genetic selection gradient **β**, propagates uncertainty through the full MCMC posterior (10,000 draws × 4 populations × 2 environments), and produces a publication-ready violin figure with posterior credible intervals.

Key output numbers (high-salt environment):

| Population | λ₁/λ₂ | gmax → S rotation (°) | Spectral gap erosion |
|---|---|---|---|
| A6140 (ancestor) | 7.22 [3.98, 12.22] | 6.5 [3.7, 12.1] | — |
| GA150 | 1.76 [1.12, 3.57] | 27.8 [9.1, 84.6] | 84% |
| GA250 | 1.79 [1.12, 3.45] | 32.8 [9.7, 85.5] | 86% |
| GA450 | 3.20 [1.44, 7.42] | 15.2 [5.8, 47.7] | 49% |

### Figures

Manuscript figures are composited in Keynote from individual PDF panels. To regenerate the analysis panels:

```bash
# Senecio summary figure (Fig. 5 panels) — run pipeline Stage 3
cd manuscript/analysis/senecio_pipeline/
bash run_pipeline.sh --stage 3

# C. elegans figure — run the C. elegans pipeline
cd manuscript/analysis/celegans_pipeline/
bash run_pipeline.sh
```

### Mathematical derivations

The Mathematica supplement (`manuscript/analysis/GeometricSnowballSupplement.nb`) contains all derivations: the channel operator construction, the geometric overlap deficit, the Davis–Kahan perturbation bound, and the connection to Schluter's line of least resistance as a limiting case. The companion file `GeometricSnowballTools.wl` provides reusable functions for computing **S**, extracting channels, and measuring principal angles.

The extended mathematical appendix (`extension/math_bootcamp_appendix.tex`, 61 pages) provides a pedagogical treatment of the linear algebra and matrix perturbation theory underlying the framework, written for quantitative biologists. A compiled PDF is available at `manuscript/supporting/Math_Boot_Camp.pdf`.

---

## Companion papers

- **Paper 1:** Additive channels in curved fitness landscapes — *GENETICS*, under review (2026).  
  Ortiz-Barrientos D & Cooper M. [bioRxiv 2026.03.21.713332](https://doi.org/10.1101/2026.03.21.713332)

- **Paper 2:** This manuscript — in preparation.

- **Paper 3:** The Hidden Geometry of Breeding Constraints — *Nature Plants*, accepted.

---

## Interactive tutorials

A pedagogical HTML app explaining why hybrid cost scales nonlinearly with channel misalignment is included at `manuscript/supporting/channel_nonlinearity_tutorial.html`. Open it directly in any browser — no installation required. It teaches the cosine product structure, the sin²θ cost curve, and the spectral-gap heterogeneity prediction through interactive visualisations with sliders and annotated equations rendered with KaTeX.

Additional tutorial materials are in `extension/`, including a detailed geometric tutorial for biologists (`gph_tutorial_for_biologists_detailed.pdf`).

These tutorials were created with the assistance of Claude (Anthropic).

---

## Data sources

No new empirical data were generated. All data come from published sources:

| Dataset | Source | Location |
|---|---|---|
| *S. lautus* common-garden morphology | Walter et al. 2018 *Am Nat* | `data/amNat/` |
| *S. lautus* breeding-design offspring | Walter et al. 2018 *Am Nat* | `data/amNat/` |
| *S. lautus* field transplant (4 environments) | Walter et al. 2016, 2020 | `data/evol/`, `data/evolLetters/` |
| *S. lautus* multigenerational transplant | Richards et al. | `data/evol2/` |
| *S. lautus* soil chemistry | Roda et al. 2013 *Evolution* | `data/` |
| *S. lautus* reproductive isolation crosses | James et al. 2025 *Mol Ecol* | `data/phyloandGenomics/` |
| *S. lautus* extrinsic RI (Cabarita Beach) | Melo et al. 2014 *New Phyt* | `data/newPhyt/` |
| *C. elegans* G-matrices and selection gradients | Mallard et al. 2023 *eLife* | `data/elifeCelegans/` |

---

## Software

| Tool | Version | Purpose |
|---|---|---|
| R | ≥ 4.0 | G-matrix estimation (MCMCglmm), Lande–Arnold regressions |
| Python | ≥ 3.9 | Figure generation, transplant analysis, pipeline orchestration |
| Mathematica | ≥ 13 | Symbolic derivations, supplement notebook |
| LaTeX | TeX Live 2024 | Extended appendix, tutorial documents |
| Git LFS | ≥ 3.0 | Large file storage (C. elegans RData files) |

---

## License

Code: MIT — see `LICENSE`.

Manuscript text: © Ortiz-Barrientos, James & Burrage, all rights reserved until published.
