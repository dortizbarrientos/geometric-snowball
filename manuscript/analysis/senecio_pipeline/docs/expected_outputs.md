# Expected outputs — validated reference numbers

This file records the reference numbers from a validated full run of the
Senecio channel-operator pipeline. After running `bash run_pipeline.sh`
(full MCMC, not `--quick-test`), compare your outputs against these values.
Agreement to ~3 significant figures for eigenvalue ratios and within a few
degrees for principal angles confirms numerical parity.

> **Note:** MCMCglmm uses random starting values. Small run-to-run
> variation is expected. If you see larger discrepancies, check the
> diagnostics CSVs (`vegetative_G_diagnostics.csv`,
> `leaf_PCs_G_diagnostics.csv`) for low effective sample sizes or high
> lag-1 autocorrelation.

---

## How to populate this file

After your first validated run:

```bash
# From the pipeline root:
head -1 data/processed/population_G_D_results/leaf_PCs_G_summary.csv
cat data/processed/population_G_D_results/leaf_PCs_G_summary.csv
```

Copy the key columns into the tables below. The columns to record are
marked with `___` placeholders.

---

## Stage 1: population-specific G

### Leaf PCs G summary

Source: `data/processed/population_G_D_results/leaf_PCs_G_summary.csv`

| Population | Ecotype | λ₁ | λ₂ | λ₁/λ₂ | Leading plane prop. | Spectral gap |
|------------|---------|----|----|--------|---------------------|-------------|
| D1 | Dune | ___ | ___ | ___ | ___ | ___ |
| D3 | Dune | ___ | ___ | ___ | ___ | ___ |
| H1 | Head | ___ | ___ | ___ | ___ | ___ |
| H2 | Head | ___ | ___ | ___ | ___ | ___ |
| T1 | Table | ___ | ___ | ___ | ___ | ___ |
| T9 | Table | ___ | ___ | ___ | ___ | ___ |
| W2 | Wood | ___ | ___ | ___ | ___ | ___ |
| W4 | Wood | ___ | ___ | ___ | ___ | ___ |

### Vegetative G summary

Source: `data/processed/population_G_D_results/vegetative_G_summary.csv`

| Population | Ecotype | λ₁ | λ₂ | λ₁/λ₂ | Leading plane prop. | Spectral gap |
|------------|---------|----|----|--------|---------------------|-------------|
| D1 | Dune | ___ | ___ | ___ | ___ | ___ |
| D3 | Dune | ___ | ___ | ___ | ___ | ___ |
| H1 | Head | ___ | ___ | ___ | ___ | ___ |
| H2 | Head | ___ | ___ | ___ | ___ | ___ |
| T1 | Table | ___ | ___ | ___ | ___ | ___ |
| T9 | Table | ___ | ___ | ___ | ___ | ___ |
| W2 | Wood | ___ | ___ | ___ | ___ | ___ |
| W4 | Wood | ___ | ___ | ___ | ___ | ___ |

### Pairwise gmax angles (leaf PCs, selected Dune × Headland pairs)

Source: `data/processed/population_G_D_results/pairwise_G_geometry_leafPC.csv`

| Pop 1 | Pop 2 | gmax angle (°) | g₂ angle (°) | Leading plane angle (°) |
|-------|-------|-----------------|---------------|------------------------|
| D1 | H1 | ___ | ___ | ___ |
| D1 | H2 | ___ | ___ | ___ |
| D3 | H1 | ___ | ___ | ___ |
| D3 | H2 | ___ | ___ | ___ |

### MCMC diagnostics check

Source: `data/processed/population_G_D_results/leaf_PCs_G_diagnostics.csv`

Confirm for each population:
- [ ] Effective sample size > 200 for diagonal VCV terms
- [ ] Lag-1 autocorrelation < 0.3 for most parameters
- [ ] No warning messages during model fitting

---

## Stage 2: channel operator

### Lande–Arnold β gradients (leaf block, selected pop × env)

Source: `data/processed/channel_operator_results/lande_arnold_beta_all.csv`

| Pop | Env | Trait | β | R² |
|-----|-----|-------|---|----|
| D1 | DUNE | LeafPC1 | ___ | ___ |
| D1 | DUNE | LeafPC2 | ___ | ___ |
| H1 | HEAD | LeafPC1 | ___ | ___ |
| H1 | HEAD | LeafPC2 | ___ | ___ |

### S channel angles (Route A — Lande–Arnold, leaf block)

Source: `data/processed/channel_operator_results/S_channel_angles_LA.csv`

| Pop 1 | Pop 2 | Env | S channel angle (°) | S gap pop1 | S gap pop2 |
|-------|-------|-----|---------------------|------------|------------|
| D1 | H1 | DUNE | ___ | ___ | ___ |
| D1 | H1 | HEAD | ___ | ___ | ___ |
| D1 | H2 | DUNE | ___ | ___ | ___ |
| D1 | H2 | HEAD | ___ | ___ | ___ |
| D3 | H1 | DUNE | ___ | ___ | ___ |
| D3 | H1 | HEAD | ___ | ___ | ___ |

### S channel angles (Route B — soil proxy, leaf block)

Source: `data/processed/channel_operator_results/S_channel_angles_soil.csv`

| Pop 1 | Pop 2 | Env | S channel angle (°) |
|-------|-------|-----|---------------------|
| D1 | H1 | DUNE | ___ |
| D1 | H1 | HEAD | ___ |
| D1 | H2 | DUNE | ___ |
| D1 | H2 | HEAD | ___ |

### Cross-validation: Route A vs Route B agreement

Record the Spearman rank correlation between Route A and Route B S channel
angles across all Dune × Headland pairs. This should be positive if channel
angles are a real geometric feature rather than a regression artefact.

```
Spearman ρ (Route A vs Route B): ___
p-value: ___
```

### Comparison: gmax angle vs S channel angle

For each Dune × Headland pair in its home environment, compare the
Schluter-limit proxy (gmax angle from Stage 1) against the full channel
angle (S angle from Stage 2). A consistent **delta** indicates the
contribution of fitness curvature beyond G-alignment alone.

| Pop 1 | Pop 2 | gmax angle (°) | S channel angle (°) | Δ (S − gmax) |
|-------|-------|-----------------|---------------------|---------------|
| D1 | H1 | ___ | ___ | ___ |
| D1 | H2 | ___ | ___ | ___ |
| D3 | H1 | ___ | ___ | ___ |
| D3 | H2 | ___ | ___ | ___ |

---

## Populating this file automatically

After a successful run, you can extract numbers with:

```bash
# Stage 1 G summary
Rscript -e '
d <- read.csv("data/processed/population_G_D_results/leaf_PCs_G_summary.csv")
print(d[, c("population","ecotype","lambda1","lambda2","lambda1_over_lambda2","leading_plane_prop","spectral_gap")])
'

# Stage 2 Dune × Headland S angles
Rscript -e '
d <- read.csv("data/processed/channel_operator_results/S_channel_angles_LA.csv")
dh <- d[(grepl("^D", d$pop1) & grepl("^H", d$pop2)) | (grepl("^H", d$pop1) & grepl("^D", d$pop2)),]
print(dh[dh$block == "leaf", c("pop1","pop2","env","S_channel_angle_deg","S_gap_pop1","S_gap_pop2")])
'
```

Copy the output into the tables above and commit.
