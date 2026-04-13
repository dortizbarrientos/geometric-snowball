Population-specific G and pairwise D analysis
============================================

Main files
----------
- population_G_D_snowball_analysis.R
- Data_Exp1_Dmatrix.csv
- Data_Exp2_Gvariance.csv

What the script does
--------------------
1. Uses Exp1 to calculate pairwise D among the breeding-design populations.
2. Uses Exp1 to build a 2-axis PCA basis for leaf traits.
3. Projects Exp2 offspring onto those leaf axes.
4. Fits population-specific G matrices with MCMCglmm:
   - 4 x 4 vegetative G for Height, MSL_W, SB, MSD
   - 2 x 2 leaf G for LeafPC1, LeafPC2
5. Saves summaries of eigenvalues, eigenvectors, anisotropy, G geometry,
   and D-to-G alignment.

How to run
----------
In R, from the folder containing the files above:

  source("population_G_D_snowball_analysis.R")

Packages required
-----------------
- MCMCglmm
- coda

Install with:

  install.packages(c("MCMCglmm", "coda"))

Important settings in the script
--------------------------------
- quick_test
- scale_within_population_for_G
- save_model_objects
- save_posterior_G

Output folder
-------------
All results are written to:

  population_G_D_results/

Key outputs
-----------
- leaf_PCA_loadings.csv
- leaf_PCA_variance_explained.csv
- Exp1_population_means_vegetative.csv
- Exp1_population_means_leafPCs.csv
- pairwise_D_vegetative_summary.csv
- pairwise_D_leafPC_summary.csv
- vegetative_G_summary.csv
- leaf_PCs_G_summary.csv
- vegetative_G_eigenvectors.csv
- leaf_PCs_G_eigenvectors.csv
- vegetative_G_diagnostics.csv
- leaf_PCs_G_diagnostics.csv
- pairwise_G_geometry_vegetative.csv
- pairwise_G_geometry_leafPC.csv
- pairwise_D_alignment_to_G_vegetative.csv
- pairwise_D_alignment_to_G_leafPC.csv
- population_G_D_results_compact.rds

Interpretation notes
--------------------
- Pairwise D is rank 1, so its leading eigenvector is the meaningful divergence axis.
- G is most useful here as a geometric object: focus on lambda1, lambda2,
  lambda1_over_lambda2, isotropy12, leading_plane_prop, and spectral_gap.
- With scale_within_population_for_G = TRUE, treat G mostly as shape/orientation,
  not raw-scale additive variance.

For a fuller version, use:
- README_population_G_D_snowball_analysis.md
