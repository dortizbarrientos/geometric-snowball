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
├── manuscript/                         Paper and analysis
│   ├── word/
│   │   └── snowball_ms_v9.docx             Current manuscript draft
│   ├── analysis/
│   │   ├── celegans_pipeline/              C. elegans channel-operator pipeline
│   │   ├── senecio_pipeline/               Senecio channel-operator pipeline
│   │   ├── GeometricSnowballSupplement.nb  Mathematica supplement (derivations)
│   │   └── GeometricSnowballTools.wl       Mathematica function library
│   ├── figures/                            Publication figures (Keynote composites)
│   ├── tables/                             Formatted tables
│   └── supporting/                         Supplementary methods PDF
│
├── data/                               Published datasets (by source paper)
│   ├── amNat/                              Walter et al. 2018 Am Nat (G-matrices)
│   ├── evolLetters/                        Walter et al. 2020 Evol Letters (transplant)
│   ├── evol/                               Walter et al. 2016 Evolution (seed set)
│   ├── evol2/                              Richards et al. (transplant, multi-generation)
│   ├── elifeCelegans/                      Mallard et al. 2023 eLife (C. elegans G + β)
│   ├── ecol/                               Walter et al. 2018 Ecology
│   ├── newPhyt/                            Melo et al. 2014 New Phytologist (extrinsic RI)
│   ├── phyloandGenomics/                   James et al. 2021 Evolution, 2021 MBE, 2025 Mol Ecol
│   ├── evolLetters-Ken/                    Thompson et al. (stickleback)
│   ├── sticklebacks/                       Schluter lab supplementary data
│   └── refsDatasets/                       Reference datasets (Silene, rotation)
│
├── scripts/                            Standalone analysis scripts
│   ├── figure_geometric_snowball_concept.py   Conceptual figure (Fig. 1)
│   ├── geometric_snowball_sims_checked.py     Simulation sweeps (Fig. 2)
│   └── gmatricesDH/                           Original validated G-matrix R scripts
│       ├── population_G_D_snowball_analysis.R
│       ├── senecio_channel_operator_analysis.R
│       └── population_G_D_results/            Pre-computed outputs
│
├── extension/                          Extended materials
│   ├── math_bootcamp_appendix.tex/pdf      Mathematical appendix (61 pages)
│   ├── gph_tutorial_for_biologists.tex/pdf Geometric tutorial for biologists
│   └── gph_landscape_manuscript_draft.tex  Extended landscape theory
│
├── LICENSE
└── README.md                           You are here
```

---

## Reproducing the analyses

The paper's empirical results come from two reproducible pipelines, one per study system. Each pipeline is self-contained: it includes all raw data, all scripts, and a single shell command that runs everything end-to-end.

### Requirements

**R** (≥ 4.0) with packages:

```bash
Rscript -e 'install.packages(c("MCMCglmm", "coda", "readxl", "dplyr", "tidyr", "ggplot2", "expm", "Matrix"))'
```

**Python** (≥ 3.9) with packages:

```bash
pip install numpy pandas matplotlib xlrd openpyxl
```

**Git LFS** (for large RData files in the C. elegans pipeline):

```bash
brew install git-lfs    # macOS
git lfs install
```

### Senecio pipeline

The Senecio pipeline estimates population-specific G-matrices for eight *S. lautus* populations, constructs the channel operator **S** = **G**^{1/2} **H** **G**^{1/2}, computes pairwise channel angles, and analyses the multigenerational transplant data that tests the framework's predictions about F1 vigour, F2 breakdown, and backcross recovery.

```bash
cd manuscript/analysis/senecio_pipeline/

# Install R packages (one-time)
Rscript requirements.R

# Quick test — short MCMC chains, ~5 minutes total
bash run_pipeline.sh --quick-test

# Full run — production MCMC chains, several hours
bash run_pipeline.sh
```

The pipeline runs in four stages:

| Stage | Script | What it does | Runtime |
|-------|--------|--------------|---------|
| 1 | `scripts/01_fit_G/population_G_D_snowball_analysis.R` | Fits population-specific **G** via MCMCglmm (8 pops × 2 trait blocks). Computes pairwise D and G geometry. | Hours (full) / ~5 min (quick) |
| 2 | `scripts/02_channel_operator/senecio_channel_operator_analysis.R` | Fits Lande–Arnold selection gradients per pop × env. Constructs **H**, computes **S**, extracts channel angles. Cross-validates with soil-based **H** proxy. | ~2 min |
| 3 | `scripts/03_figures/make_figures.py` | Summary figure: G anisotropy, pairwise g_max angles, S channel angles. | Seconds |
| 4 | `scripts/04_transplant/transplant_analysis.py` | Multigenerational fitness trajectory (F1 vigour → F2 breakdown → BC recovery). Extrinsic RI per generation. Seed-set trajectory. | Seconds |

Outputs land in `data/processed/` (Stages 1–2) and `results/` (Stages 3–4). To run individual stages:

```bash
bash run_pipeline.sh --stage 1        # G-matrix fitting only
bash run_pipeline.sh --stage 4        # transplant analysis only (no MCMC needed)
```

Key output files and the paper claims they support:

| Paper claim | Output file | Key column |
|---|---|---|
| G-matrix anisotropy per population | `vegetative_G_summary.csv`, `leaf_PCs_G_summary.csv` | `lambda1_over_lambda2` |
| Pairwise g_max angles (Schluter limit) | `pairwise_G_geometry_leafPC.csv` | `gmax_angle_deg` |
| S channel angles (full operator) | `S_channel_angles_LA.csv` | `S_channel_angle_deg` |
| F1 heterosis (+39% vs midparent) | `transplant_multigenerational_summary.csv` | `L10_mean` |
| F2 breakdown (−29% vs F1) | same | same |
| Extrinsic RI sign reversal (F1→F2) | `transplant_extrinsic_RI.csv` | `RI_extrinsic` |
| Seed-set F2 crash (0.61 → 0.12) | `seedset_trajectory.csv` | `mean_PropFilled` |
| Backcross recovery (83% of parental) | `transplant_fitness_by_generation.csv` | `L10` for BCHD in HEAD |

See `docs/expected_outputs.md` for validated reference numbers.

### C. elegans pipeline

The C. elegans pipeline tests the spectral-gap prediction using Mallard et al. (2023) data: 50 generations of adaptation to high salt in one ancestral and three evolved replicate populations.

```bash
cd manuscript/analysis/celegans_pipeline/

# Full run — deterministic, ~2 minutes
bash run_pipeline.sh
```

The pipeline parses posterior G-matrices from RData files, constructs a rank-1 surrogate **S** from the genetic selection gradient **β**, propagates uncertainty through the full MCMC posterior (10,000 draws × 4 populations × 2 environments), and produces a publication-ready violin figure and formatted results table.

Key output numbers:

| Population | λ₁/λ₂ | gmax → S rotation (°) | Spectral gap erosion |
|---|---|---|---|
| A6140 (ancestor) | 7.22 | 6.5° | — |
| GA150 | 1.76 | 27.8° | 84% |
| GA250 | 1.79 | 32.8° | 86% |
| GA450 | 3.20 | 15.2° | 49% |

See `docs/expected_outputs.md` in the C. elegans pipeline for full posterior intervals.

### Figures

The manuscript figures are composited in Keynote from individual PDF panels. To regenerate the panels:

```bash
# Conceptual figure (Fig. 1) — cones, grids, timeline
python3 scripts/figure_geometric_snowball_concept.py

# Simulation sweeps (Fig. 2) — geometric vs mutational snowball dynamics
python3 scripts/geometric_snowball_sims_checked.py

# Senecio summary (Fig. 5 panels) — run the Senecio pipeline Stage 3
cd manuscript/analysis/senecio_pipeline/
bash run_pipeline.sh --stage 3

# C. elegans figure — run the C. elegans pipeline
cd manuscript/analysis/celegans_pipeline/
bash run_pipeline.sh
```

Figure panels are written to `results/figures/` within each pipeline, and to the repository root `manuscript/figures/forKeynote/` for final compositing.

### Mathematical derivations

The Mathematica supplement (`manuscript/analysis/GeometricSnowballSupplement.nb`) contains all derivations:
the channel operator construction, the geometric overlap deficit, the Davis–Kahan perturbation bound, and the connection to Schluter's line of least resistance as a limiting case. The companion file `GeometricSnowballTools.wl` provides reusable functions.

The extended mathematical appendix (`extension/math_bootcamp_appendix.tex`, 61 pages) provides a pedagogical treatment of the linear algebra and matrix perturbation theory underlying the framework, written for quantitative biologists.

---

## Companion papers

- **Paper 1:** Additive channels in curved fitness landscapes — *GENETICS*, under review (2026).  
  Ortiz-Barrientos D & Cooper M. [bioRxiv 2026.03.21.713332](https://doi.org/10.1101/2026.03.21.713332)

- **Paper 2:** This manuscript — in preparation.

- **Paper 3:** The Hidden Geometry of Breeding Constraints — *Nature Plants*, accepted.

---

## Interactive tutorials

Pedagogical HTML apps accompanying the theory are in the repository and can be opened directly in any browser — no installation required. Each app teaches one concept from the framework through interactive visualisations with sliders, annotated equations (rendered with KaTeX), and step-by-step explanations.

Topics covered include:

- Channel nonlinearity: why hybrid cost scales as cos²θ
- The channel operator: how **G** and **H** jointly determine the adaptive direction
- Spectral gaps and rotational sensitivity

These apps were created with the assistance of Claude (Anthropic).

---

## Data sources and citations

No new empirical data were generated for this study. All data come from published sources:

| Dataset | Source | Location in repo |
|---|---|---|
| *S. lautus* common-garden morphology (Exp 1) | Walter et al. 2018 *Am Nat* | `data/amNat/` |
| *S. lautus* breeding-design morphology (Exp 2) | Walter et al. 2018 *Am Nat* | `data/amNat/` |
| *S. lautus* field transplant fitness | Walter et al. 2016 *Evolution*, 2020 *Evol Letters* | `data/evol/`, `data/evolLetters/` |
| *S. lautus* multigenerational transplant | Richards et al. | `data/evol2/` |
| *S. lautus* soil chemistry | Roda et al. 2013 *Evolution* | `data/` |
| *S. lautus* reproductive isolation crosses | James et al. 2025 *Mol Ecol* | `data/phyloandGenomics/` |
| *S. lautus* extrinsic RI (Cabarita) | Melo et al. 2014 *New Phyt* | `data/newPhyt/` |
| *C. elegans* G-matrices and selection gradients | Mallard et al. 2023 *eLife* | `data/elifeCelegans/` |
| Threespine stickleback supplementary | Schluter lab | `data/sticklebacks/` |

---

## Software and computation

| Tool | Version | Purpose |
|---|---|---|
| R | ≥ 4.0 | G-matrix estimation (MCMCglmm), Lande–Arnold regressions |
| Python | ≥ 3.9 | Figure generation, transplant analysis, pipeline orchestration |
| Mathematica | ≥ 13 | Symbolic derivations, supplement notebook |
| LaTeX | TeX Live 2024 | Extended appendix, tutorial documents |
| Git LFS | ≥ 3.0 | Large RData file storage (C. elegans pipeline) |

---

## License

Code: MIT — see `LICENSE`.

Manuscript text: © Ortiz-Barrientos, James & Burrage, all rights reserved until published.
