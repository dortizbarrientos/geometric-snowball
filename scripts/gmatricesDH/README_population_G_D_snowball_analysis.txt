Population-specific G and pairwise D analysis
===========================================

Files expected in the working directory
--------------------------------------
1. Data_Exp1_Dmatrix.csv
2. Data_Exp2_Gvariance.csv
3. population_G_D_snowball_analysis.R

What the script does
--------------------
1. Reads Exp1 and Exp2.
2. Uses Exp1 to compute pairwise D among the eight breeding-design populations.
3. Uses Exp1 to build a common 2-axis leaf PCA basis.
4. Projects Exp2 offspring onto that leaf basis.
5. Fits population-specific G matrices with MCMCglmm:
   - 4 x 4 vegetative G for Height, MSL_W, SB, MSD
   - 2 x 2 leaf G for LeafPC1, LeafPC2
6. Summarises eigenvalues, eigenvectors, anisotropy, and leading-plane variance.
7. Computes pairwise geometry among G matrices and alignment between D and G.

How to run
----------
Open R in the folder containing the files above and run:

  source("population_G_D_snowball_analysis.R")

Package requirements
--------------------
- MCMCglmm
- coda

Install with:

  install.packages(c("MCMCglmm", "coda"))

Important settings near the top of the script
--------------------------------------------
- quick_test: shorter MCMC chains for a dry run
- scale_within_population_for_G: TRUE focuses more on matrix shape / anisotropy
- save_model_objects: save fitted MCMCglmm model objects
- save_posterior_G: save posterior arrays of G matrices

Main outputs
------------
leaf_PCA_loadings.csv
leaf_PCA_variance_explained.csv
Exp1_population_means_vegetative.csv
Exp1_population_means_leafPCs.csv
pairwise_D_vegetative_summary.csv
pairwise_D_leafPC_summary.csv
pairwise_D_vegetative_matrices.rds
pairwise_D_leafPC_matrices.rds
vegetative_G_summary.csv
leaf_PCs_G_summary.csv
vegetative_G_eigenvectors.csv
leaf_PCs_G_eigenvectors.csv
pairwise_G_geometry_vegetative.csv
pairwise_G_geometry_leafPC.csv
pairwise_D_alignment_to_G_vegetative.csv
pairwise_D_alignment_to_G_leafPC.csv
population_G_D_results_compact.rds

Interpretation notes
--------------------
- Pairwise D is rank 1 by construction, so the biologically meaningful axis is
  the leading eigenvector (the standardised divergence vector).
- For G, useful shape summaries are:
  * lambda1_over_lambda2
  * isotropy12 (= lambda2 / lambda1)
  * leading_plane_prop
  * spectral_gap
- The D-to-G alignment tables report:
  * angle_to_gmax_deg
  * angle_to_g2_deg
  * leading_plane_capture
