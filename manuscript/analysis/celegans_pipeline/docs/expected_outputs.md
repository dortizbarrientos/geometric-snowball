# Expected outputs (reproducibility reference)

This document records the exact numerical outputs produced by the
pipeline when run against the supplied raw data. If you run
`bash run_pipeline.sh` on a clean clone and get different numbers,
something has gone wrong — either the raw data have been modified,
the software environment differs substantially, or there is a bug.

The pipeline is deterministic: there is no random seed, no bootstrap,
no Monte Carlo step. These numbers should reproduce exactly (to the
reported precision) on any modern machine.

---

## High-salt environment (NaCl) — primary analysis

Posterior median [95% credible interval] from 10,000 MCMC draws per
population.

| Population | λ₁/λ₂ | Spectral gap δ | Rotation (°) |
|---|---|---|---|
| A6140 | 7.22 [3.98, 12.22] | 0.819 [0.441, 1.306] | 6.5 [3.7, 12.1] |
| GA150 | 1.76 [1.12, 3.57]  | 0.134 [0.026, 0.364] | 27.8 [9.1, 84.6] |
| GA250 | 1.79 [1.12, 3.45]  | 0.118 [0.022, 0.301] | 32.8 [9.7, 85.5] |
| GA450 | 3.20 [1.44, 7.42]  | 0.418 [0.119, 0.936] | 15.2 [5.8, 47.7] |

Spectral-gap erosion relative to ancestor:

- GA150: 84%
- GA250: 86%
- GA450: 49%

## Low-salt environment (NGM) — secondary (context only)

| Population | λ₁/λ₂ | Rotation (°) |
|---|---|---|
| A6140 | 3.05 [1.97, 4.81] | 9.0 [5.4, 15.2] |
| GA150 | 2.22 [1.17, 5.06] | 24.5 [7.3, 80.9] |
| GA250 | 1.92 [1.14, 4.07] | 15.2 [4.4, 78.1] |
| GA450 | 2.86 [1.30, 6.60] | 20.9 [6.4, 73.3] |

The low-salt environment shows less dramatic spectral collapse than
high-salt, which is biologically expected: these worms were selected
in high salt, so the high-salt G matrix is where the signature of
directional selection is clearest.

---

## Sanity checks

Beyond matching the tables above, the pipeline should produce:

- `data/processed/beta_ancestor_highsalt_standardized.csv` with seven
  rows (SF, SB, FS, FB, BS, BF, Size) and β_SF = −0.229, β_Size = 0.439
  (posterior modes matching Mallard et al. Table 1).
- `data/processed/G_matrices_posterior_mode_long.csv` with exactly
  441 rows (9 matrices × 49 cells) and A6140 high-salt diagonal
  variances of 0.080, 0.117, 0.052, 0.156, 0.019, 0.128, 0.063.
- `data/processed/G_chains_high_salt_long.csv` of size ~70 MB with
  exactly 4 × 10,000 × 49 = 1,960,000 rows.
- `data/processed/celegans_posterior_perdraw.csv` with exactly 80,000
  rows (4 pops × 2 envs × 10,000 draws).
- `results/figures/celegans_posterior_violins.pdf` — two-panel violin
  figure (see `celegans_posterior_violins.png` for a preview).
- `results/tables/table_celegans_summary.docx` — six-column summary
  table matching the high-salt table above.

## If the numbers disagree

1. Verify that `data/raw/G_matrices_high_salt.RData` and
   `G_matrices_low_salt.RData` match the originals (file size 44 MB
   for each, SHA-256 sums in the Dryad repository).
2. Verify R version ≥ 4.0. The behaviour of `load()` and `matrix()`
   on this data has been stable across R 4.x.
3. Verify Python package versions against `requirements.txt`. In
   particular, numpy's `eigh` and the matrix-multiplication operator
   should produce identical results across recent versions.
4. Contact d.ortizbarrientos@uq.edu.au with the discrepancy and the
   console log from `bash run_pipeline.sh`.
