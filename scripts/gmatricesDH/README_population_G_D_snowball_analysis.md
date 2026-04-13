# Population-specific G and pairwise D analysis

This README goes with `population_G_D_snowball_analysis.R` and the two data files:

- `Data_Exp1_Dmatrix.csv`
- `Data_Exp2_Gvariance.csv`

The script is set up to estimate:

- **pairwise D** from experiment 1
- **population-specific G** from experiment 2

It does this in two trait blocks:

- **Vegetative traits:** `Height`, `MSL_W`, `SB`, `MSD`
- **Leaf traits:** reduced to **two axes** (`LeafPC1`, `LeafPC2`) using PCA on six leaf traits

## What the script is trying to do

The goal is to produce the quantities needed for the snowball / channel-style analyses without trying to fit a noisy full 10-trait G per population.

The workflow is:

1. Read the two CSV files.
2. Harmonise trait names between experiments.
3. Use experiment 1 to calculate population means and pairwise **D**.
4. Build a single leaf PCA basis from experiment 1.
5. Project experiment 2 individuals onto that same leaf basis.
6. Fit **population-specific G** matrices with `MCMCglmm` for:
   - a **4 × 4 vegetative G**
   - a **2 × 2 leaf-axis G**
7. Extract eigenvalues, eigenvectors, anisotropy summaries, and alignment summaries.
8. Save a compact results object plus human-readable CSV files.

## Input files expected in the working directory

The script assumes these files are in the folder where you run R:

- `population_G_D_snowball_analysis.R`
- `Data_Exp1_Dmatrix.csv`
- `Data_Exp2_Gvariance.csv`

## Expected columns in each file

### `Data_Exp1_Dmatrix.csv`

Required columns:

- `Ecotype`
- `Population`
- `VegHeight`
- `MSL_W`
- `SB`
- `MSD`
- `Area`
- `P2A2`
- `Circularity`
- `Nindents.Peri`
- `IndentWidth`
- `IndentDepth`

The script internally renames:

- `VegHeight` -> `Height`

### `Data_Exp2_Gvariance.csv`

Required columns:

- `Type`
- `Block`
- `sire`
- `dam`
- `Height`
- `MSL_W`
- `SB`
- `MSD`
- `Area`
- `P2A2`
- `Circularity`
- `Nindents.Peri`
- `IndentWidthMean`
- `IndentDepthMean`

The script internally renames:

- `IndentWidthMean` -> `IndentWidth`
- `IndentDepthMean` -> `IndentDepth`

## Trait blocks used in the analysis

### Vegetative block

These are fitted directly:

- `Height`
- `MSL_W`
- `SB`
- `MSD`

### Leaf block

The script starts from these six measured leaf traits:

- `Area`
- `P2A2`
- `Circularity`
- `Nindents.Peri`
- `IndentWidth`
- `IndentDepth`

It then builds a **global PCA basis from experiment 1** and keeps the first two axes:

- `LeafPC1`
- `LeafPC2`

Those two axes are what get used for the leaf-side G models.

## How D is estimated

The script does **not** build one multi-population D matrix in the Walter et al. ecotype MANOVA sense.

Instead, for each pair of breeding-design populations, it calculates:

\[
D_{ij} = (\mu_i - \mu_j)(\mu_i - \mu_j)^T
\]

So each pairwise D matrix is **rank 1**.

That means the biologically meaningful axis is the leading eigenvector only, which is just the standardised divergence vector between the two populations.

The script calculates pairwise D separately for:

- the 4 vegetative traits
- the 2 leaf PCA axes

## How G is estimated

G is estimated separately for each population in `Data_Exp2_Gvariance.csv`, where population is defined by `Block`.

For each population, the script fits a multivariate sire–dam breeding-design model using `MCMCglmm`:

- random terms: `us(trait):sire + us(trait):dam + us(trait):sire:dam`
- residual: `us(trait):units`

The additive genetic covariance matrix is extracted as:

\[
G = 4 \times \text{sire covariance}
\]

This follows the usual paternal half-sib logic used in the supplementary Walter code.

### Standardisation for G

By default:

```r
scale_within_population_for_G <- TRUE
```

So each trait is scaled **within population** before fitting G.

That means the resulting matrices are best interpreted for:

- **shape**
- **anisotropy**
- **relative orientation**
- **spectral summaries**

rather than absolute additive variance on the original raw scale.

## MCMC settings

At the top of the script you can change:

```r
quick_test <- FALSE
```

- `TRUE`: shorter chains for a quick dry run
- `FALSE`: longer chains matching the saved-sample depth used in the original setup of the script

Current settings are:

### Quick test

```r
MCMC_NITT   <- 260000
MCMC_BURNIN <- 60000
MCMC_THIN   <- 200
```

### Full run

```r
MCMC_NITT   <- 2100000
MCMC_BURNIN <- 100000
MCMC_THIN   <- 2000
```

## R packages required

Install these first:

```r
install.packages(c("MCMCglmm", "coda"))
```

The script loads:

- `MCMCglmm`
- `coda`

## How to run it

From R, in the folder containing the script and the two data files:

```r
source("population_G_D_snowball_analysis.R")
```

## Output directory

All results are written into:

```r
population_G_D_results/
```

This folder is created automatically if it does not already exist.

## Main output files

### Leaf PCA outputs

- `population_G_D_results/leaf_PCA_loadings.csv`
- `population_G_D_results/leaf_PCA_variance_explained.csv`

### Population means from experiment 1

- `population_G_D_results/Exp1_population_means_vegetative.csv`
- `population_G_D_results/Exp1_population_means_leafPCs.csv`

### Pairwise D outputs

- `population_G_D_results/pairwise_D_vegetative_summary.csv`
- `population_G_D_results/pairwise_D_leafPC_summary.csv`
- `population_G_D_results/pairwise_D_vegetative_matrices.rds`
- `population_G_D_results/pairwise_D_leafPC_matrices.rds`

### Population-specific G summaries

- `population_G_D_results/vegetative_G_summary.csv`
- `population_G_D_results/leaf_PCs_G_summary.csv`
- `population_G_D_results/vegetative_G_eigenvectors.csv`
- `population_G_D_results/leaf_PCs_G_eigenvectors.csv`
- `population_G_D_results/vegetative_G_diagnostics.csv`
- `population_G_D_results/leaf_PCs_G_diagnostics.csv`

### Pairwise geometry among G matrices

- `population_G_D_results/pairwise_G_geometry_vegetative.csv`
- `population_G_D_results/pairwise_G_geometry_leafPC.csv`

### Alignment of D with G

- `population_G_D_results/pairwise_D_alignment_to_G_vegetative.csv`
- `population_G_D_results/pairwise_D_alignment_to_G_leafPC.csv`

### Per-population posterior mean G matrices

Patterns like:

- `population_G_D_results/vegetative_Gmean_D1.csv`
- `population_G_D_results/vegetative_Gmean_D3.csv`
- `population_G_D_results/leaf_PCs_Gmean_D1.csv`
- `population_G_D_results/leaf_PCs_Gmean_D3.csv`

and so on for each breeding-design population.

### Compact R object

- `population_G_D_results/population_G_D_results_compact.rds`

## Key summary columns to look at

### In the G summary tables

Especially useful:

- `lambda1`
- `lambda2`
- `prop1`
- `prop2`
- `leading_plane_prop`
- `spectral_gap`
- `lambda1_over_lambda2`
- `isotropy12`

Practical reading:

- **high `lambda1_over_lambda2`** = stronger anisotropy
- **high `isotropy12`** = first two eigenvalues are more similar
- **high `leading_plane_prop`** = most variance concentrated in the first two axes
- **large `spectral_gap`** = stronger separation of `gmax` from the second axis

### In the pairwise G geometry tables

Useful columns:

- `gmax_angle_deg`
- `g2_angle_deg`
- `leading_plane_angle_deg`

These tell you how similar the G orientations are between populations.

### In the D-to-G alignment tables

Useful columns:

- `angle_to_gmax_deg`
- `angle_to_g2_deg`
- `leading_plane_capture`

These tell you how well the pairwise divergence vector aligns with each population's G.

## Important interpretation notes

1. **Pairwise D is rank 1 by construction.**
   So the first eigenvector is the only meaningful divergence axis for a given population pair.

2. **The leaf analysis is deliberately reduced-dimensional.**
   You are not fitting a 6-trait leaf G here; you are fitting a 2-axis G after PCA reduction.

3. **Population-specific G is ambitious for these data.**
   These models should be treated as useful exploratory estimates, but the MCMC diagnostics matter.

4. **When `scale_within_population_for_G = TRUE`, interpret G mainly as geometry.**
   That is, focus on eigenstructure, anisotropy, and alignment rather than raw-scale variance magnitude.

## Diagnostics to inspect before trusting the results

Open these files first:

- `population_G_D_results/vegetative_G_diagnostics.csv`
- `population_G_D_results/leaf_PCs_G_diagnostics.csv`

Things worth checking:

- lag-1 autocorrelation is not excessive
- effective sample sizes are reasonably large
- no populations show obviously pathological estimates
- posterior mean G matrices are biologically sensible

## Troubleshooting

### Error: missing package

Install:

```r
install.packages(c("MCMCglmm", "coda"))
```

### Error: missing columns

Check the input CSV headers match the column names listed above exactly.

### Error: cannot find files

Make sure you are running `source()` from the directory that contains:

- `Data_Exp1_Dmatrix.csv`
- `Data_Exp2_Gvariance.csv`
- `population_G_D_snowball_analysis.R`

### Run is taking too long

Set:

```r
quick_test <- TRUE
```

then confirm everything runs before switching back to the full chains.

## Suggested first pass

A sensible first pass is:

1. set `quick_test <- TRUE`
2. run the script once
3. inspect the diagnostics CSVs and the per-population G means
4. rerun with `quick_test <- FALSE` if the results look stable

