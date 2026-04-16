################################################################################
# senecio_channel_operator_analysis.R
#
# "The Geometric Snowball of Speciation"
# Ortiz-Barrientos, James, Burrage
#
# PURPOSE
# -------
# This script estimates the full channel operator S = G^{1/2} H G^{1/2} for
# eight Senecio lautus populations (D1, D3, H1, H2, T1, T9, W2, W4) using:
#
#   (A) LANDE-ARNOLD SELECTION GRADIENTS (Option 3)
#       Family-mean morphology (glasshouse) linked to family-mean fitness
#       (field transplant) via shared full-sibling family IDs.
#       This estimates the fitness curvature matrix H in each environment
#       from the Lande-Arnold quadratic regression.
#
#   (B) SOIL-BASED ENVIRONMENTAL PROXY (Option 2)
#       The three biologically meaningful soil axes (Salinity, Nutrients,
#       Acidity) are used to construct a population-level environmental
#       distance matrix. The trait-environment relationship gives a
#       rank-1 approximation to H via the selection differential.
#       This serves as cross-validation for Option A.
#
# BIOLOGICAL BACKGROUND
# ---------------------
# In the geometric snowball framework, each population evolves a preferred
# direction of additive genetic response — its "genetic channel" — determined
# jointly by the G-matrix (available heritable variation) and the H matrix
# (local fitness curvature). When two populations adapt to different
# environments, their channels rotate apart. The angle between channels
# predicts hybrid fitness loss (I_geom ≈ k * theta^2).
#
# Estimating S directly, rather than relying on the Schluter-limit proxy,
# moves the empirical test from "does gmax alignment predict RI?" to
# "does full channel misalignment predict RI?" — a much stronger test.
#
# WORKFLOW
# --------
# 1.  Load G-matrices (from population_G_D_snowball_analysis.R outputs)
# 2.  Link glasshouse morphology to field fitness via family IDs
# 3.  Project leaf traits onto the global PCA basis
# 4.  Compute family-mean morphology and fitness per environment
# 5.  Run Lande-Arnold regression per population × environment:
#       (a) Directional gradients beta (linear selection)
#       (b) Quadratic gradients gamma (stabilising/disruptive selection)
# 6.  Construct H from gamma (H = -gamma, sign convention)
# 7.  Compute S = G^{1/2} H G^{1/2} for each population × environment
# 8.  Extract channel angles and compare to Schluter-limit (gmax) angles
# 9.  Build soil-based H proxy and cross-validate
# 10. Test whether S-based channel angles predict RI better than gmax angles
# 11. Save all outputs
#
# INPUT FILES NEEDED (in working directory or specify paths below)
# ---------------------------------------------------------------
# From the population_G_D_snowball_analysis.R outputs:
#   population_G_D_results/vegetative_Gmean_*.csv
#   population_G_D_results/leaf_PCs_Gmean_*.csv
#   population_G_D_results/leaf_PCA_loadings.csv
#
# From the project data directory:
#   Data_Exp2_Gvariance.csv     (glasshouse morphology)
#   Walter_Data_Exp2.csv        (field transplant fitness)
#   All_S__lautus_soil_data.xlsx (soil variables)
#
# Optional (for RI test):
#   CrossesDH-ReciprocalsAveraged.txt  (RI data)
#   allTraitsTheta.txt                  (PCVA angles)
#
# OUTPUT DIRECTORY
# ----------------
# All outputs written to: channel_operator_results/
#
# PACKAGES REQUIRED
# -----------------
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "expm", "Matrix"))
#
# USAGE
# -----
# Rscript senecio_channel_operator_analysis.R
# or in R: source("senecio_channel_operator_analysis.R")
#
# AUTHORS
# -------
# D. Ortiz-Barrientos, M.E. James, K. Burrage
# Code annotations: Claude (Anthropic), April 2026
################################################################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(expm)   # for matrix square root

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 0: PATHS AND SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
#
# HOW TO SET PATHS
# ----------------
# Option A (easiest): set your working directory to the folder that contains
#   Data_Exp2_Gvariance.csv, Walter_Data_Exp2.csv, and
#   All_S__lautus_soil_data.xlsx, and set G_RESULTS_DIR to wherever
#   population_G_D_snowball_analysis.R wrote its outputs.
#
#   In RStudio: Session → Set Working Directory → Choose Directory
#   In R:       setwd("/path/to/your/data/folder")
#
# Option B: set the full paths explicitly in the two variables below.
#
# EXAMPLE (macOS):
#   DATA_DIR      <- "/Users/uqdortiz/Desktop/snowball_data"
#   G_RESULTS_DIR <- "/Users/uqdortiz/Desktop/snowball_data/population_G_D_results"
#
# The script will tell you exactly which files it cannot find before stopping.

# ── Path overrides for pipeline layout ──────────────────────────────────────
# Original script assumed it ran from scripts/gmatricesDH/ with raw data two
# levels up. In this pipeline layout, raw data lives in data/raw/, Stage 1
# outputs live in data/processed/population_G_D_results/, and Stage 2 writes to
# data/processed/channel_operator_results/. These six lines are the only edits
# to the original scripts/gmatricesDH/senecio_channel_operator_analysis.R.
# The wrapper run_pipeline.sh always invokes Rscript from this stage directory,
# so "../.." resolves to the pipeline root.
PIPELINE_ROOT <- normalizePath(file.path("..", ".."), mustWork = FALSE)
G_RESULTS_DIR <- file.path(PIPELINE_ROOT, "data", "processed", "population_G_D_results")
MORPH_FILE    <- file.path(PIPELINE_ROOT, "data", "raw", "Data_Exp2_Gvariance.csv")
FITNESS_FILE  <- file.path(PIPELINE_ROOT, "data", "raw", "Walter_Data_Exp2.csv")
SOIL_FILE     <- file.path(PIPELINE_ROOT, "data", "raw", "All_S__lautus_soil_data.xlsx")
# ROOT is no longer used; kept for safety in case downstream code references it.
ROOT          <- PIPELINE_ROOT
# ─────────────────────────────────────────────────────────────────────────────

# Analysis settings
FITNESS_MEASURE <- "L10"   # "L10" = establishment (recommended), "Bud", "GE"
MIN_FAMILIES    <- 10      # minimum linked families per pop × env to run regression
RUN_QUADRATIC   <- TRUE    # also fit quadratic terms for gamma (H estimate)

# Output directory (pipeline layout: writes to data/processed/, not CWD)
OUT_DIR <- file.path(PIPELINE_ROOT, "data", "processed", "channel_operator_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== CHANNEL OPERATOR ANALYSIS ===\n")
cat("Working directory:", getwd(), "\n")
cat("Data directory:   ", normalizePath(".", mustWork = FALSE), "\n")
cat("G-matrix results: ", normalizePath(G_RESULTS_DIR, mustWork = FALSE), "\n\n")

# ── FILE EXISTENCE CHECK — stops early with clear diagnostics ────────────────
required_files <- list(
  DATA = c(MORPH_FILE, FITNESS_FILE, SOIL_FILE),
  G_MATRIX = c(
    file.path(G_RESULTS_DIR, "leaf_PCA_loadings.csv"),
    file.path(G_RESULTS_DIR, "leaf_PCs_Gmean_D1.csv"),
    file.path(G_RESULTS_DIR, "vegetative_Gmean_D1.csv"),
    file.path(G_RESULTS_DIR, "pairwise_G_geometry_leafPC.csv")
  )
)

missing <- Filter(function(f) !file.exists(f), unlist(required_files))

if (length(missing) > 0) {
  cat("\n*** MISSING FILES — cannot proceed ***\n")
  cat("The following required files were not found:\n")
  for (f in missing) cat("  MISSING:", f, "\n")
  cat("\nFIX: update MORPH_FILE, FITNESS_FILE, SOIL_FILE at the top of this script\n")
  cat("     to point to the correct directories, then re-run.\n\n")
  cat("Current MORPH_FILE  =", MORPH_FILE, "\n"); cat("Current FITNESS_FILE=", FITNESS_FILE, "\n"); cat("Current SOIL_FILE   =", SOIL_FILE, "\n")
  cat("Current G_RESULTS_DIR =", normalizePath(G_RESULTS_DIR, mustWork = FALSE), "\n")
  stop("Path configuration error — see messages above.", call. = FALSE)
}

cat("All required files found ✓\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: LOAD G-MATRICES
# ─────────────────────────────────────────────────────────────────────────────
# We load the posterior-mean G matrices estimated by population_G_D_snowball_analysis.R.
# Two trait blocks:
#   Vegetative: Height, MSL_W, SB, MSD  (4 × 4 matrices)
#   Leaf PCs:   LeafPC1, LeafPC2        (2 × 2 matrices)

cat("Loading G-matrices...\n")

populations <- c("D1", "D3", "H1", "H2", "T1", "T9", "W2", "W4")
ecotype_map <- c(D1="Dune", D3="Dune", H1="Head", H2="Head",
                 T1="Table", T9="Table", W2="Wood", W4="Wood")

veg_traits  <- c("Height", "MSL_W", "SB", "MSD")
leaf_traits <- c("LeafPC1", "LeafPC2")

load_G_matrix <- function(pop, block) {
  # Load the posterior-mean G matrix CSV for a given population and trait block
  # Returns a named matrix
  prefix <- ifelse(block == "veg", "vegetative", "leaf_PCs")
  fname  <- file.path(G_RESULTS_DIR, paste0(prefix, "_Gmean_", pop, ".csv"))
  if (!file.exists(fname)) {
    warning("G-matrix file not found: ", fname)
    return(NULL)
  }
  mat <- as.matrix(read.csv(fname, row.names = 1))
  # Ensure symmetry (minor numerical drift from MCMC)
  mat <- 0.5 * (mat + t(mat))
  traits <- if (block == "veg") veg_traits else leaf_traits
  rownames(mat) <- colnames(mat) <- traits
  return(mat)
}

G_veg  <- lapply(setNames(populations, populations),
                 function(p) load_G_matrix(p, "veg"))
G_leaf <- lapply(setNames(populations, populations),
                 function(p) load_G_matrix(p, "leaf"))

# Quick sanity: all matrices positive definite?
check_pd <- function(G, label) {
  ev <- tryCatch(eigen(G, only.values = TRUE)$values, error = function(e) NA)
  if (any(is.na(ev)) || any(ev <= 0)) {
    warning(label, ": G-matrix not positive definite. min eigenvalue = ",
            min(ev, na.rm = TRUE))
    return(FALSE)
  }
  return(TRUE)
}

cat("G-matrix PD check:\n")
for (p in populations) {
  ok_v <- if (!is.null(G_veg[[p]]))  check_pd(G_veg[[p]],  paste(p,"veg"))  else FALSE
  ok_l <- if (!is.null(G_leaf[[p]])) check_pd(G_leaf[[p]], paste(p,"leaf")) else FALSE
  cat(sprintf("  %s: veg=%s  leaf=%s\n", p,
              ifelse(ok_v, "OK", "WARN"), ifelse(ok_l, "OK", "WARN")))
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: LOAD AND LINK MORPHOLOGY-FITNESS DATA
# ─────────────────────────────────────────────────────────────────────────────
# Data_Exp2_Gvariance.csv: glasshouse morphology (same full-sib families as
#   used to estimate G matrices)
# Walter_Data_Exp2.csv:    field transplant fitness (same families, different
#   experiment — linked via Family ID)
#
# IMPORTANT: H1 families have a trailing "R" in the transplant file that must
# be stripped before matching.

cat("\nLoading morphology and fitness data...\n")

morph_raw <- read.csv(MORPH_FILE, stringsAsFactors = FALSE)
fit_raw   <- read.csv(FITNESS_FILE, stringsAsFactors = FALSE)

# Rename columns to match (IndentWidthMean → IndentWidth etc.)
morph_raw <- morph_raw %>%
  rename(IndentWidth = IndentWidthMean,
         IndentDepth = IndentDepthMean)

# Fix H1 family ID mismatch: strip trailing "R" from transplant data
# (H1-A01R in transplant ↔ H1-A01 in G-matrix glasshouse)
fit_raw$Family[fit_raw$TypeS == "H1"] <-
  sub("R$", "", fit_raw$Family[fit_raw$TypeS == "H1"])

# Restrict transplant data to parental populations matching G-matrix pops
fit_par <- fit_raw %>%
  filter(TypeS %in% populations)

# Find overlapping families
overlap_families <- intersect(
  unique(morph_raw$Family),
  unique(fit_par$Family)
)
cat(sprintf("  Overlapping families: %d across %d populations\n",
            length(overlap_families),
            length(unique(morph_raw$Family[morph_raw$Family %in% overlap_families]))))

# Filter both datasets to overlapping families only
morph_link <- morph_raw %>%
  filter(Family %in% overlap_families) %>%
  rename(Population = Block)

fit_link <- fit_par %>%
  filter(Family %in% overlap_families) %>%
  rename(Population = TypeS)

cat("  Families per population:\n")
print(table(morph_link$Population))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: COMPUTE LEAF PCA SCORES
# ─────────────────────────────────────────────────────────────────────────────
# Project the 6 leaf morphology traits in the G-matrix glasshouse data onto
# the global PCA basis built in population_G_D_snowball_analysis.R.
# This ensures our leaf axes are the same as those used in the G-matrix
# estimation — essential for S = G^{1/2} H G^{1/2} to make sense.

cat("\nProjecting leaf traits onto PCA basis...\n")

pca_loadings <- read.csv(file.path(G_RESULTS_DIR, "leaf_PCA_loadings.csv"),
                         stringsAsFactors = FALSE)

leaf_raw_traits <- c("Area", "P2A2", "Circularity",
                     "Nindents.Peri", "IndentWidth", "IndentDepth")

# Standardise leaf traits within the glasshouse data (same scaling as G-matrix fit)
# Note: we use the GLOBAL mean and SD across all populations (not within-pop)
# to preserve the between-population signal in PCA space.
leaf_data <- morph_link[, leaf_raw_traits]
leaf_scaled <- scale(leaf_data)  # global standardisation

# Extract PC loadings matrix (rows = traits, cols = PC1, PC2)
L <- as.matrix(pca_loadings[, c("LeafPC1_loading", "LeafPC2_loading")])
rownames(L) <- pca_loadings$trait

# Project: scores = scaled_data %*% L
leaf_scores <- leaf_scaled %*% L
colnames(leaf_scores) <- c("LeafPC1", "LeafPC2")

# Append to morphology data
morph_link$LeafPC1 <- leaf_scores[, "LeafPC1"]
morph_link$LeafPC2 <- leaf_scores[, "LeafPC2"]

cat("  LeafPC1 and LeafPC2 computed for all",
    nrow(morph_link), "individuals\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: COMPUTE FAMILY MEANS
# ─────────────────────────────────────────────────────────────────────────────
# The Lande-Arnold regression requires:
#   z_i = mean trait value of family i  (from glasshouse)
#   W_i = mean fitness of family i in each environment (from field)
#
# Using family means removes within-family environmental variance and focuses
# on among-family (additive genetic) variance, which is what G describes.

cat("\nComputing family means...\n")

all_morph_traits <- c(veg_traits, "LeafPC1", "LeafPC2")

fam_morph <- morph_link %>%
  group_by(Population, Family) %>%
  summarise(across(all_of(all_morph_traits), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

fam_fit <- fit_link %>%
  group_by(Population, Family, ENV) %>%
  summarise(fitness = mean(.data[[FITNESS_MEASURE]], na.rm = TRUE),
            n_seeds = n(),
            .groups = "drop") %>%
  filter(n_seeds >= 3)   # require at least 3 seeds per family-environment

cat(sprintf("  Morphology families: %d\n", nrow(fam_morph)))
cat(sprintf("  Fitness family-environments: %d\n", nrow(fam_fit)))

# Join morphology and fitness
fam_data <- fam_fit %>%
  inner_join(fam_morph, by = c("Population", "Family"))

cat("  Linked family-environment records:", nrow(fam_data), "\n")
cat("  Coverage:\n")
fam_data %>%
  group_by(Population, ENV) %>%
  summarise(n = n(), mean_fit = round(mean(fitness), 3),
            .groups = "drop") %>%
  print()

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: LANDE-ARNOLD SELECTION GRADIENT REGRESSION
# ─────────────────────────────────────────────────────────────────────────────
# For each population × environment combination, fit:
#
#   W_rel = a + z * beta + 0.5 * z^T * gamma * z + epsilon
#
# where:
#   W_rel  = relative fitness (W / mean(W)) — ensures mean = 1
#   z      = standardised trait vector (mean 0, SD 1 within pop × env)
#   beta   = directional selection gradients (p × 1 vector)
#   gamma  = matrix of nonlinear selection gradients (p × p symmetric)
#
# From Lande & Arnold (1983):
#   beta  = P^{-1} S    where S = cov(z, W_rel) — selection differential
#   gamma = P^{-1} C P^{-1}  where C = cov(z z^T, W_rel) — approximately
#
# In practice we use multiple regression which automatically handles
# the P^{-1} inversion:
#   lm(W_rel ~ z1 + z2 + ... + zp + z1^2 + z2^2 + z1:z2 + ...)
#
# IMPORTANT INTERPRETATION NOTE (for the biologist)
# -------------------------------------------------
# beta tells us the direction of steepest fitness ascent in trait space.
# In the geometric snowball framework, beta ≈ -H * d where d is the
# displacement of the population from its optimum. If the population is
# close to its optimum (d is small), beta is small and noisy.
#
# gamma gives the curvature of the fitness surface: gamma_ij < 0 means
# stabilising selection (fitness decreases away from optimum in direction ij);
# gamma_ij > 0 means disruptive selection.
#
# H = -gamma  (sign convention: H is positive definite at a fitness peak)
#
# With family-mean data (n ≈ 12-30 per pop × env):
#   - beta estimates are useful but have large standard errors
#   - gamma estimates are very noisy — treat as exploratory
#   - We report both but base H on the combination of approaches

cat("\n=== SECTION 5: LANDE-ARNOLD SELECTION GRADIENTS ===\n")

run_lande_arnold <- function(dat, traits, fitness_col = "fitness",
                             run_quad = TRUE, min_n = 10) {
  # Runs the Lande-Arnold regression for one population × environment block.
  #
  # Arguments:
  #   dat       : data.frame with columns for traits and fitness_col
  #   traits    : character vector of trait names to include
  #   fitness_col: name of the fitness column
  #   run_quad  : if TRUE, also fit quadratic terms for gamma
  #   min_n     : minimum sample size; returns NULL if fewer families
  #
  # Returns a list with:
  #   beta      : directional selection gradients (named vector)
  #   gamma     : nonlinear selection gradient matrix (if run_quad = TRUE)
  #   r2_linear : R^2 of the linear model
  #   r2_quad   : R^2 of the quadratic model (if run_quad = TRUE)
  #   n         : number of families
  #   converged : logical
  
  dat <- dat[complete.cases(dat[, c(traits, fitness_col)]), ]
  n   <- nrow(dat)
  
  if (n < min_n) {
    return(list(beta = NULL, gamma = NULL, n = n, converged = FALSE,
                message = paste("Too few families:", n)))
  }
  
  # Compute relative fitness (mean = 1 by construction)
  W_rel <- dat[[fitness_col]] / mean(dat[[fitness_col]])
  
  # Standardise traits within this pop × env block
  # This puts beta on the same scale as Lande-Arnold (1983):
  # beta_i is the change in W per SD change in trait i, holding others constant
  z_mat <- scale(dat[, traits])
  trait_means <- attr(z_mat, "scaled:center")
  trait_sds   <- attr(z_mat, "scaled:scale")
  
  # Guard against zero-variance traits
  zero_var <- trait_sds == 0
  if (any(zero_var)) {
    warning("Zero-variance traits removed: ", paste(traits[zero_var], collapse = ", "))
    traits  <- traits[!zero_var]
    z_mat   <- z_mat[, !zero_var, drop = FALSE]
  }
  
  # --- Linear model: W_rel ~ z (estimates beta) ---
  z_df     <- as.data.frame(z_mat)
  lin_form <- as.formula(paste("W_rel ~", paste(traits, collapse = " + ")))
  lin_fit  <- tryCatch(lm(lin_form, data = cbind(W_rel = W_rel, z_df)),
                       error = function(e) NULL)
  
  if (is.null(lin_fit)) {
    return(list(beta = NULL, gamma = NULL, n = n, converged = FALSE,
                message = "Linear model failed to converge"))
  }
  
  coefs     <- coef(lin_fit)
  beta_raw  <- coefs[traits]
  names(beta_raw) <- traits
  r2_linear <- summary(lin_fit)$r.squared
  
  # --- Quadratic model: W_rel ~ z + z^2 + z_i * z_j (estimates gamma) ---
  gamma_mat <- NULL
  r2_quad   <- NA
  
  if (run_quad && length(traits) >= 2) {
    # Build quadratic terms: z_i^2 (diagonal gamma) and z_i * z_j (off-diagonal)
    quad_terms <- list()
    quad_names <- c()
    
    for (i in seq_along(traits)) {
      # Diagonal: I(z_i^2)
      vname <- paste0(traits[i], "_sq")
      quad_terms[[vname]] <- z_mat[, i]^2
      quad_names <- c(quad_names, vname)
    }
    for (i in seq_along(traits)) {
      for (j in seq_along(traits)) {
        if (j > i) {
          # Off-diagonal: z_i * z_j
          vname <- paste0(traits[i], "_x_", traits[j])
          quad_terms[[vname]] <- z_mat[, i] * z_mat[, j]
          quad_names <- c(quad_names, vname)
        }
      }
    }
    quad_df  <- as.data.frame(quad_terms)
    quad_all <- cbind(W_rel = W_rel, z_df, quad_df)
    
    quad_form <- as.formula(paste("W_rel ~",
                                  paste(c(traits, quad_names), collapse = " + ")))
    quad_fit  <- tryCatch(lm(quad_form, data = quad_all),
                          error = function(e) NULL)
    
    if (!is.null(quad_fit)) {
      q_coefs   <- coef(quad_fit)
      r2_quad   <- summary(quad_fit)$r.squared
      p         <- length(traits)
      gamma_mat <- matrix(0, p, p, dimnames = list(traits, traits))
      
      # Fill diagonal: gamma_ii = 2 * coef(z_i^2)
      # (factor of 2 from the quadratic derivative convention)
      for (i in seq_along(traits)) {
        vname <- paste0(traits[i], "_sq")
        if (vname %in% names(q_coefs)) {
          gamma_mat[i, i] <- 2 * q_coefs[vname]
        }
      }
      # Fill off-diagonal: gamma_ij = gamma_ji = coef(z_i:z_j)
      for (i in seq_along(traits)) {
        for (j in seq_along(traits)) {
          if (j > i) {
            vname <- paste0(traits[i], "_x_", traits[j])
            if (vname %in% names(q_coefs)) {
              gamma_mat[i, j] <- q_coefs[vname]
              gamma_mat[j, i] <- q_coefs[vname]
            }
          }
        }
      }
    }
  }
  
  list(beta       = beta_raw,
       gamma      = gamma_mat,
       n          = n,
       converged  = TRUE,
       r2_linear  = r2_linear,
       r2_quad    = r2_quad,
       trait_means = trait_means,
       trait_sds   = trait_sds,
       message    = "OK")
}

# ── Run regressions for all pop × environment × trait-block combinations ─────

environments <- c("DUNE", "HEAD", "TABLE", "WOOD")

results_la <- list()

for (pop in populations) {
  for (env in environments) {
    for (block in c("veg", "leaf")) {
      traits <- if (block == "veg") veg_traits else leaf_traits
      dat_sub <- fam_data %>%
        filter(Population == pop, ENV == env) %>%
        as.data.frame()
      
      if (nrow(dat_sub) < MIN_FAMILIES) next
      
      res <- run_lande_arnold(dat_sub, traits,
                              fitness_col = "fitness",
                              run_quad = RUN_QUADRATIC,
                              min_n = MIN_FAMILIES)
      key <- paste(pop, env, block, sep = "_")
      results_la[[key]] <- c(res, list(pop = pop, env = env, block = block))
    }
  }
}

cat(sprintf("  Regressions run: %d\n", length(results_la)))
cat(sprintf("  Converged: %d\n",
            sum(sapply(results_la, function(x) isTRUE(x$converged)))))

# ── Summarise beta (directional gradients) ────────────────────────────────────
cat("\n--- Directional selection gradients (beta) summary ---\n")

beta_summary <- lapply(names(results_la), function(key) {
  r <- results_la[[key]]
  if (!isTRUE(r$converged) || is.null(r$beta)) return(NULL)
  data.frame(
    key   = key,
    pop   = r$pop,
    env   = r$env,
    block = r$block,
    n     = r$n,
    r2    = round(r$r2_linear, 3),
    trait = names(r$beta),
    beta  = round(r$beta, 4)
  )
}) %>% bind_rows()

cat("\n  Leaf PCs — selection gradients per pop × env:\n")
beta_leaf <- beta_summary %>% filter(block == "leaf")
print(
  beta_leaf %>%
    pivot_wider(names_from = c(env, trait), values_from = beta) %>%
    select(pop, n, r2, everything())
)

write.csv(beta_summary,
          file.path(OUT_DIR, "lande_arnold_beta_summary.csv"),
          row.names = FALSE)
cat("  Saved: lande_arnold_beta_summary.csv\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: CONSTRUCT H FROM GAMMA
# ─────────────────────────────────────────────────────────────────────────────
# The fitness curvature matrix H is related to gamma by:
#   H = -gamma
#
# This follows from the quadratic fitness approximation:
#   log W(z) = log W_0 - (1/2)(z - theta)^T H (z - theta)
#   => d^2 log W / dz_i dz_j = -H_ij
#   => gamma_ij ≈ -H_ij
#
# H must be positive definite for a fitness MAXIMUM to exist at the optimum.
# Negative definite gamma (all eigenvalues negative) implies H > 0 (a peak).
# If gamma is not negative definite, we project to the nearest PD matrix,
# which corresponds to "trimming" disruptive selection in directions where
# we expect stabilising selection near the optimum.
#
# IMPORTANT CAVEAT: With n = 12-30 families, gamma estimates have wide
# confidence intervals. We use H primarily for its EIGENVECTOR (direction
# of steepest curvature), not its magnitude. The beta-based rank-1
# approximation H_approx = beta * beta^T / ||beta||^2 is used as a
# fallback when gamma is not reliably negative definite.

cat("\n=== SECTION 6: CONSTRUCTING H MATRICES ===\n")

nearest_pd <- function(M, epsilon = 1e-4) {
  # Project a symmetric matrix to the nearest positive definite matrix.
  # Uses eigenvalue thresholding: replace negative eigenvalues with epsilon.
  ev  <- eigen(M, symmetric = TRUE)
  lam <- pmax(ev$values, epsilon)
  ev$vectors %*% diag(lam) %*% t(ev$vectors)
}

construct_H <- function(gamma_mat, beta_vec, method = "gamma") {
  # Construct the fitness curvature matrix H from gamma or beta.
  #
  # method = "gamma" : H = -gamma, projected to PD if needed
  # method = "beta"  : H ≈ beta %*% t(beta) / ||beta||^2 (rank-1)
  # method = "hybrid": use gamma eigenvectors, beta for magnitude scaling
  
  if (method == "beta" || is.null(gamma_mat)) {
    # Rank-1 approximation: dominant curvature direction = selection direction
    b  <- beta_vec / sqrt(sum(beta_vec^2))
    H  <- outer(b, b)
    attr(H, "method") <- "beta_rank1"
    return(H)
  }
  
  # H = -gamma (sign convention)
  H_raw  <- -gamma_mat
  H_sym  <- 0.5 * (H_raw + t(H_raw))
  
  ev     <- eigen(H_sym, symmetric = TRUE)
  is_pd  <- all(ev$values > 0)
  
  if (is_pd) {
    attr(H_sym, "method") <- "gamma_pd"
    return(H_sym)
  }
  
  # Not PD: project to nearest PD matrix
  H_pd <- nearest_pd(H_sym)
  attr(H_pd, "method") <- "gamma_projected"
  return(H_pd)
}

# Construct H for each pop × env × block
H_matrices <- list()

for (key in names(results_la)) {
  r <- results_la[[key]]
  if (!isTRUE(r$converged)) next
  
  gamma  <- r$gamma
  beta   <- r$beta
  
  # Decide method based on gamma quality.
  # Guard: gamma can contain NA/Inf when lm() is rank-deficient
  # (common with small n and many quadratic terms, e.g. veg block with n=12).
  gamma_ok <- !is.null(gamma) &&
    !any(is.na(gamma)) &&
    !any(is.infinite(gamma))
  
  if (gamma_ok) {
    ev_gamma <- eigen(gamma, symmetric = TRUE, only.values = TRUE)$values
    neg_def  <- all(ev_gamma < 0)   # gamma < 0 ↔ H > 0 (fitness peak)
    if (neg_def) {
      H <- construct_H(gamma, beta, method = "gamma")
    } else {
      # gamma not negative definite — project to nearest PD
      H <- construct_H(gamma, beta, method = "hybrid")
    }
  } else {
    # gamma unusable — fall back to rank-1 beta approximation
    H <- construct_H(NULL, beta, method = "beta")
  }
  
  H_matrices[[key]] <- list(H = H, pop = r$pop, env = r$env, block = r$block)
}

cat(sprintf("  H matrices constructed: %d\n", length(H_matrices)))

methods_used <- sapply(H_matrices, function(x) attr(x$H, "method"))
cat("  Methods used:\n")
print(table(methods_used))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: COMPUTE S = G^{1/2} H G^{1/2}
# ─────────────────────────────────────────────────────────────────────────────
# For each population × environment combination, compute the symmetric channel
# operator S and extract its leading eigenvector (the genetic channel direction
# in genetically whitened trait space).
#
# The channel angle between two populations in the SAME environment can then
# be compared to the PCVA-based angle to test whether incorporating H
# changes the prediction of RI.
#
# We also compute the angle between a population's home-environment channel
# and its away-environment channels — this is the rotation angle induced by
# adapting to different fitness surfaces.

cat("\n=== SECTION 7: CHANNEL OPERATOR S = G^{1/2} H G^{1/2} ===\n")

compute_channel_operator <- function(G, H) {
  # Compute S = G^{1/2} H G^{1/2}
  #
  # Steps:
  # 1. Compute G^{1/2} via spectral decomposition
  # 2. S = G_half %*% H %*% G_half
  # 3. Symmetrise S (guard against numerical drift)
  #
  # Both G and H must be positive definite.
  
  G_sym <- 0.5 * (G + t(G))
  H_sym <- 0.5 * (H + t(H))
  
  # G^{1/2} via eigendecomposition
  ev_G <- eigen(G_sym, symmetric = TRUE)
  lam  <- pmax(ev_G$values, 1e-12)  # numerical safety
  G_half <- ev_G$vectors %*% diag(sqrt(lam)) %*% t(ev_G$vectors)
  
  # S = G^{1/2} H G^{1/2}
  S     <- G_half %*% H_sym %*% G_half
  S_sym <- 0.5 * (S + t(S))
  
  # Eigendecomposition of S gives the channels
  ev_S <- eigen(S_sym, symmetric = TRUE)
  
  # Map leading eigenvector back to phenotype space: u = G^{1/2} v
  # (where v is the leading eigenvector of S in whitened space)
  v1_whitened  <- ev_S$vectors[, 1]
  u1_phenotype <- G_half %*% v1_whitened
  u1_phenotype <- u1_phenotype / sqrt(sum(u1_phenotype^2))
  
  list(
    S             = S_sym,
    eigenvalues   = ev_S$values,
    eigenvectors  = ev_S$vectors,           # whitened space
    gmax_whitened = v1_whitened,
    gmax_phenotype = c(u1_phenotype),        # phenotype space (for display)
    G_half        = G_half,
    spectral_gap  = ev_S$values[1] - ev_S$values[2],
    lambda1_over_lambda2 = ev_S$values[1] / max(ev_S$values[2], 1e-12)
  )
}

principal_angle <- function(u, v) {
  # Acute angle between two unit vectors (degrees)
  cos_theta <- abs(sum(u * v)) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
  cos_theta <- pmin(cos_theta, 1.0)
  acos(cos_theta) * 180 / pi
}

# Compute S for each combination
S_results <- list()

for (key in names(H_matrices)) {
  h_entry <- H_matrices[[key]]
  pop     <- h_entry$pop
  env     <- h_entry$env
  block   <- h_entry$block
  H       <- h_entry$H
  
  G <- if (block == "veg") G_veg[[pop]] else G_leaf[[pop]]
  if (is.null(G)) next
  
  # Rescale G and H to comparable scale
  # (G was estimated on standardised traits; H from regression on same scale)
  S_entry <- tryCatch(
    compute_channel_operator(G, H),
    error = function(e) {
      warning("S computation failed for ", key, ": ", e$message)
      NULL
    }
  )
  if (!is.null(S_entry)) {
    S_results[[key]] <- c(S_entry, list(pop = pop, env = env, block = block,
                                        H_method = attr(H, "method")))
  }
}

cat(sprintf("  S matrices computed: %d\n", length(S_results)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: CHANNEL ANGLE COMPARISONS
# ─────────────────────────────────────────────────────────────────────────────
# For each pair of populations in the SAME environment, compute the angle
# between their S-based channel directions (in whitened space).
# Compare to the Schluter-limit (gmax) angles from the G-matrix analysis.

cat("\n=== SECTION 8: CHANNEL ANGLE COMPARISONS ===\n")

pairwise_channel_angles <- function(S_list, env_filter, block_filter) {
  # Extract gmax directions for all pops in a given environment and block
  entries <- Filter(function(x) x$env == env_filter && x$block == block_filter,
                    S_list)
  if (length(entries) < 2) return(NULL)
  
  pops_here <- sapply(entries, function(x) x$pop)
  results   <- list()
  
  for (i in seq_along(pops_here)) {
    for (j in seq_along(pops_here)) {
      if (j <= i) next
      p1 <- pops_here[i]
      p2 <- pops_here[j]
      
      # Angle between leading eigenvectors in WHITENED space
      v1 <- entries[[i]]$gmax_whitened
      v2 <- entries[[j]]$gmax_whitened
      theta_S <- principal_angle(v1, v2)
      
      results[[paste(p1, p2, sep = "_")]] <- data.frame(
        env   = env_filter,
        block = block_filter,
        pop1  = p1, pop2 = p2,
        eco1  = ecotype_map[p1],
        eco2  = ecotype_map[p2],
        S_channel_angle_deg = round(theta_S, 2),
        S_gap_pop1  = round(entries[[i]]$spectral_gap, 4),
        S_gap_pop2  = round(entries[[j]]$spectral_gap, 4),
        S_lam_ratio_pop1 = round(entries[[i]]$lambda1_over_lambda2, 3),
        S_lam_ratio_pop2 = round(entries[[j]]$lambda1_over_lambda2, 3)
      )
    }
  }
  bind_rows(results)
}

# Compute angles for each environment × block
cat("\nS-based channel angles (Dune × Headland, leaf PCs):\n")
all_angle_results <- list()

for (env in environments) {
  for (blk in c("veg", "leaf")) {
    ang <- pairwise_channel_angles(S_results, env, blk)
    if (!is.null(ang)) {
      all_angle_results[[paste(env, blk, sep = "_")]] <- ang
    }
  }
}

S_angles <- bind_rows(all_angle_results)

# Print Dune × Headland comparisons
dh_S <- S_angles %>%
  filter((eco1 == "Dune" & eco2 == "Head") |
           (eco1 == "Head" & eco2 == "Dune")) %>%
  arrange(block, env, pop1, pop2)

print(dh_S[, c("block","env","pop1","pop2",
               "S_channel_angle_deg","S_gap_pop1","S_gap_pop2")])

write.csv(S_angles,
          file.path(OUT_DIR, "S_channel_angles_all.csv"),
          row.names = FALSE)
cat("\n  Saved: S_channel_angles_all.csv\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: COMPARE S-BASED VS SCHLUTER-LIMIT CHANNEL ANGLES
# ─────────────────────────────────────────────────────────────────────────────
# The Schluter-limit analysis (H = hI, gmax proxy) used angles between
# PCVA divergence vectors as a proxy for channel misalignment.
# The S-based analysis uses the full S = G^{1/2} H G^{1/2} eigenvectors.
#
# We test: do S-based channel angles differ from gmax angles?
# And: do S-based angles predict RI better than gmax angles?

cat("\n=== SECTION 9: S vs SCHLUTER-LIMIT COMPARISON ===\n")

# Load Schluter-limit gmax angles from G-matrix analysis
gmax_geom_leaf <- tryCatch(
  read.csv(file.path(G_RESULTS_DIR, "pairwise_G_geometry_leafPC.csv")),
  error = function(e) NULL
)

if (!is.null(gmax_geom_leaf)) {
  gmax_leaf_dh <- gmax_geom_leaf %>%
    filter((grepl("^D", pop1) & grepl("^H", pop2)) |
             (grepl("^H", pop1) & grepl("^D", pop2)))
  cat("\n  Schluter-limit (gmax) Dune×Headland angles (leaf PCs):\n")
  print(gmax_leaf_dh[, c("pop1","pop2","gmax_angle_deg","spectral_gap_pop1","spectral_gap_pop2")])
  
  cat("\n  S-based angles by home environment (leaf):\n")
  dh_S_leaf <- dh_S %>% filter(block == "leaf")
  print(dh_S_leaf)
  
  # For populations in their HOME environment, compare gmax vs S angles
  home_env <- c(D1 = "DUNE", D3 = "DUNE", H1 = "HEAD", H2 = "HEAD",
                T1 = "TABLE", T9 = "TABLE", W2 = "WOOD", W4 = "WOOD")
  
  cat("\n  Comparison: gmax angle vs S-based angle (home environment, leaf):\n")
  comp <- dh_S_leaf %>%
    filter(mapply(function(p1, p2, env)
      home_env[p1] == env | home_env[p2] == env,
      pop1, pop2, env)) %>%
    left_join(gmax_leaf_dh %>% select(pop1, pop2, gmax_angle_deg),
              by = c("pop1", "pop2")) %>%
    mutate(delta_angle = round(S_channel_angle_deg - gmax_angle_deg, 2))
  print(comp[, c("pop1","pop2","env","gmax_angle_deg",
                 "S_channel_angle_deg","delta_angle")])
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: SOIL-BASED H PROXY (OPTION 2)
# ─────────────────────────────────────────────────────────────────────────────
# Use the three soil composite axes (Salinity, Nutrients, Acidity) to build
# a population-level environmental characterisation.
#
# Strategy: the correlation between soil variables and trait means
# (across populations) gives the direction in trait space that aligns
# with the environmental gradient. This is the direction that
# natural selection is expected to act along — i.e., the leading
# direction of H.
#
# This is NOT a full H matrix (no magnitude, no off-diagonal curvature).
# It gives a rank-1 approximation H_soil = d_env %*% t(d_env) where
# d_env is the trait direction most aligned with the soil gradient.
#
# Biological justification: if soil salinity drives adaptive divergence,
# then H_salinity concentrates curvature in the trait directions most
# relevant to salinity adaptation. The angle between H_salinity eigenvectors
# across populations captures how different their fitness landscapes are —
# exactly what we need for S.

cat("\n=== SECTION 10: SOIL-BASED H PROXY ===\n")

# Load soil data
soil_raw <- read_excel(SOIL_FILE, sheet = "All soil data")

# Population ID map
pop_soil_map <- c(D1 = "D01", D3 = "D03", H1 = "H01", H2 = "H02",
                  T1 = "T1",  T9 = "T9",  W2 = "I5",  W4 = "I24")

salinity_vars  <- c("Electrical Conductivity @ 25_C", "Chloride",
                    "Exchangeable Sodium", "Exchangeable Magnesium",
                    "Exchangeable Calcium")
nutrient_vars  <- c("Total Nitrogen as N", "Total Organic Carbon",
                    "Bicarbonate Ext. P (Colwell)", "Cation Exchange Capacity")
acidity_vars   <- c("pH Value", "Exchangeable Aluminium", "Aluminium")

all_soil_vars <- c(salinity_vars, nutrient_vars, acidity_vars)

# Extract and clean soil data for our populations
soil_sub <- soil_raw %>%
  filter(Population %in% pop_soil_map) %>%
  mutate(G_pop = names(pop_soil_map)[match(Population, pop_soil_map)]) %>%
  select(G_pop, Ecotype, Population, all_of(all_soil_vars)) %>%
  mutate(across(all_of(all_soil_vars), as.numeric))

# Log-transform skewed variables and standardise
skewed_vars <- c(salinity_vars, "Aluminium", "Bicarbonate Ext. P (Colwell)",
                 "Total Nitrogen as N", "Total Organic Carbon")
soil_mat <- soil_sub[, all_soil_vars]
for (v in skewed_vars) {
  if (v %in% colnames(soil_mat)) {
    soil_mat[[v]] <- log1p(pmax(soil_mat[[v]], 0, na.rm = TRUE))
  }
}
soil_scaled <- as.data.frame(scale(soil_mat))
rownames(soil_scaled) <- soil_sub$G_pop

# PCA of each biological group
soil_pc1 <- function(vars, dat) {
  sub <- dat[, vars, drop = FALSE]
  sub <- sub[, apply(sub, 2, function(x) !all(is.na(x))), drop = FALSE]
  if (ncol(sub) == 0) return(rep(NA, nrow(dat)))
  sv  <- svd(sub, nu = 1, nv = 0)
  sv$u[, 1] * sv$d[1]
}

soil_composite <- data.frame(
  G_pop           = soil_sub$G_pop,
  Ecotype         = soil_sub$Ecotype,
  Salinity_PC1    = soil_pc1(salinity_vars, soil_scaled),
  Nutrients_PC1   = soil_pc1(nutrient_vars,  soil_scaled),
  Acidity_PC1     = soil_pc1(acidity_vars,   soil_scaled)
)
rownames(soil_composite) <- soil_composite$G_pop

cat("\nSoil composite axes per population:\n")
print(soil_composite[, c("Ecotype","Salinity_PC1","Nutrients_PC1","Acidity_PC1")])

# Soil-to-trait mapping: regress trait means on soil axes
# Population-mean morphology (standardised)
pop_morph_means <- fam_morph %>%
  group_by(Population) %>%
  summarise(across(all_of(all_morph_traits), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

# Standardise trait means across populations
pop_morph_std <- pop_morph_means %>%
  mutate(across(all_of(all_morph_traits),
                ~ scale(.)[,1]))

# Join with soil composites
soil_trait <- inner_join(
  soil_composite,
  pop_morph_std %>% rename(G_pop = Population),
  by = "G_pop"
)

cat("\nSoil → trait regressions (how well does each soil axis predict each trait?):\n")
soil_axes <- c("Salinity_PC1","Nutrients_PC1","Acidity_PC1")

soil_beta_mat <- matrix(NA, nrow = length(all_morph_traits),
                        ncol = length(soil_axes),
                        dimnames = list(all_morph_traits, soil_axes))

for (ax in soil_axes) {
  for (tr in all_morph_traits) {
    x <- soil_trait[[ax]]
    y <- soil_trait[[tr]]
    if (all(is.na(x)) || all(is.na(y))) next
    fit <- lm(y ~ x)
    soil_beta_mat[tr, ax] <- coef(fit)["x"]
  }
}

soil_beta_df <- as.data.frame(round(soil_beta_mat, 3))
cat("\n  Soil-trait regression coefficients (slope of trait ~ soil axis):\n")
print(soil_beta_df)

# The dominant direction of H in each environment can be approximated as the
# trait direction most aligned with the soil gradient in THAT environment.
# For simplicity, use salinity gradient (main Dune-Headland axis) as H's
# dominant eigenvector and construct rank-1 H_soil per population.

cat("\n  Building rank-1 H_soil per population (based on local soil gradient)...\n")

H_soil <- list()
env_soil_map <- c(DUNE = "Salinity_PC1", HEAD = "Salinity_PC1",
                  TABLE = "Acidity_PC1", WOOD = "Nutrients_PC1")

for (pop in populations) {
  for (env in environments) {
    ax <- env_soil_map[env]
    d_soil <- soil_beta_mat[, ax]  # trait direction aligned with env axis
    d_soil_complete <- d_soil[!is.na(d_soil)]
    if (length(d_soil_complete) == 0) next
    
    # Rank-1 H_soil for leaf PCs (project soil gradient onto leaf space)
    leaf_idx <- c("LeafPC1","LeafPC2")
    d_leaf   <- d_soil[leaf_idx]
    if (any(is.na(d_leaf))) next
    d_leaf_norm <- d_leaf / sqrt(sum(d_leaf^2))
    H_leaf_soil <- outer(d_leaf_norm, d_leaf_norm)
    
    key <- paste(pop, env, "leaf_soil", sep = "_")
    H_soil[[key]] <- list(H = structure(H_leaf_soil, method = "soil_rank1"),
                          pop = pop, env = env, block = "leaf")
    
    # Also for veg traits
    veg_idx    <- veg_traits
    d_veg      <- d_soil[veg_idx]
    if (any(is.na(d_veg))) next
    d_veg_norm <- d_veg / sqrt(sum(d_veg^2))
    H_veg_soil <- outer(d_veg_norm, d_veg_norm)
    
    key_v <- paste(pop, env, "veg_soil", sep = "_")
    H_soil[[key_v]] <- list(H = structure(H_veg_soil, method = "soil_rank1"),
                            pop = pop, env = env, block = "veg")
  }
}

# Compute S using soil-based H
S_soil_results <- list()
for (key in names(H_soil)) {
  h_e  <- H_soil[[key]]
  G    <- if (h_e$block %in% c("leaf","leaf_soil")) G_leaf[[h_e$pop]] else G_veg[[h_e$pop]]
  if (is.null(G)) next
  S_e <- tryCatch(compute_channel_operator(G, h_e$H),
                  error = function(e) NULL)
  if (!is.null(S_e)) {
    S_soil_results[[key]] <- c(S_e,
                               list(pop = h_e$pop, env = h_e$env, block = h_e$block))
  }
}

cat(sprintf("  S_soil matrices computed: %d\n", length(S_soil_results)))

# Pairwise channel angles from soil-based S
S_soil_angles_list <- list()
for (env in environments) {
  for (blk in c("leaf")) {
    block_filt <- paste0(blk, "_soil")
    ang <- pairwise_channel_angles(S_soil_results, env, block_filt)
    if (!is.null(ang)) S_soil_angles_list[[paste(env, blk, sep = "_")]] <- ang
  }
}
S_soil_angles <- bind_rows(S_soil_angles_list)

if (nrow(S_soil_angles) > 0) {
  cat("\n  Soil-based S channel angles (Dune × Headland, leaf):\n")
  dh_soil <- S_soil_angles %>%
    filter((eco1 == "Dune" & eco2 == "Head") |
             (eco1 == "Head" & eco2 == "Dune"))
  print(dh_soil[, c("env","pop1","pop2","S_channel_angle_deg")])
  write.csv(S_soil_angles,
            file.path(OUT_DIR, "S_soil_channel_angles.csv"),
            row.names = FALSE)
  cat("  Saved: S_soil_channel_angles.csv\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: SAVE ALL OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== SECTION 11: SAVING OUTPUTS ===\n")

# 1. Selection gradients
write.csv(beta_summary,
          file.path(OUT_DIR, "lande_arnold_beta_all.csv"),
          row.names = FALSE)

# 2. S channel angles (Lande-Arnold H)
write.csv(S_angles,
          file.path(OUT_DIR, "S_channel_angles_LA.csv"),
          row.names = FALSE)

# 3. S channel angles (soil H)
if (exists("S_soil_angles") && nrow(S_soil_angles) > 0) {
  write.csv(S_soil_angles,
            file.path(OUT_DIR, "S_channel_angles_soil.csv"),
            row.names = FALSE)
}

# 4. Soil composite table
write.csv(soil_composite,
          file.path(OUT_DIR, "soil_composite_indices.csv"),
          row.names = FALSE)

# 5. Soil → trait regression betas
write.csv(as.data.frame(soil_beta_mat, optional = TRUE),
          file.path(OUT_DIR, "soil_trait_regression_betas.csv"))

# 6. Compact RDS with all S matrices
saveRDS(list(S_LA    = S_results,
             S_soil  = S_soil_results,
             H_LA    = H_matrices,
             H_soil  = H_soil,
             G_veg   = G_veg,
             G_leaf  = G_leaf),
        file.path(OUT_DIR, "channel_operator_results.rds"))

cat("\nOutput files written to:", OUT_DIR, "\n")
cat(list.files(OUT_DIR, "\\.csv$"), sep = "\n")

cat("\n=== DONE ===\n")
cat("Next steps:\n")
cat("  1. Inspect lande_arnold_beta_all.csv — check sign and magnitude of beta\n")
cat("  2. Compare S_channel_angles_LA.csv vs S_channel_angles_soil.csv\n")
cat("  3. Cross-validate against pairwise_G_geometry_leafPC.csv (Schluter limit)\n")
cat("  4. Test whether S-based angles predict RI better than gmax angles\n")
cat("     (merge with CrossesDH-ReciprocalsAveraged.txt and run lm(RI ~ S_angle))\n")