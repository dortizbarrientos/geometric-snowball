# Population-specific G and pairwise D analyses for Senecio lautus
# -------------------------------------------------------------------
# This script is designed for the two uploaded data files:
#   1) Data_Exp1_Dmatrix.csv    -> pairwise D from common-garden means
#   2) Data_Exp2_Gvariance.csv  -> population-specific G from the breeding design
#
# Main choices in this script:
#   * Vegetative block: Height, MSL_W, SB, MSD
#   * Leaf block: reduce six leaf traits to two global axes (LeafPC1, LeafPC2)
#     using PCA on Exp1, then estimate a 2 x 2 G for each population.
#   * G is fitted separately for each breeding-design population.
#   * D is calculated pairwise among the populations present in Exp2 using
#       D_ij = (mu_i - mu_j) %*% t(mu_i - mu_j)
#     so the leading eigenvector is the standardized divergence vector.
#
# The code follows the Walter et al. workflow of estimating D from experiment 1
# and G from experiment 2, while adapting the G analysis to population-specific
# 4-trait vegetative matrices and 2-axis leaf matrices.
#
# IMPORTANT:
#   Population-level G estimates are more weakly replicated than the original
#   pooled-ecotype analysis, so inspect MCMC diagnostics before interpreting
#   biological patterns.
# -------------------------------------------------------------------

# ==========================
# User settings
# ==========================
# Path overrides for pipeline layout. Original script read CSVs from CWD and
# wrote outputs to "population_G_D_results/". In this pipeline layout raw data
# lives in data/raw/ and processed outputs live in data/processed/. These four
# lines are the only edits; all downstream analysis code is identical to the
# original scripts/gmatricesDH/population_G_D_snowball_analysis.R.
# The wrapper run_pipeline.sh always invokes Rscript from this stage directory,
# so "../.." resolves to the pipeline root.
PIPELINE_ROOT <- normalizePath(file.path("..", ".."), mustWork = FALSE)
exp1_file <- file.path(PIPELINE_ROOT, "data", "raw", "Data_Exp1_Dmatrix.csv")
exp2_file <- file.path(PIPELINE_ROOT, "data", "raw", "Data_Exp2_Gvariance.csv")
outdir    <- file.path(PIPELINE_ROOT, "data", "processed", "population_G_D_results")

# TRUE = shorter chains for a quick dry run.
# FALSE = chains matching the original paper's saved-sample depth.
quick_test <- FALSE

# If TRUE, standardise each trait within each population before fitting G.
# This is useful when the main goal is matrix shape / anisotropy rather than
# absolute genetic variance on the raw scale.
scale_within_population_for_G <- TRUE

# Save fitted MCMCglmm objects and posterior G arrays.
save_model_objects <- FALSE
save_posterior_G   <- FALSE

# ==========================
# Packages
# ==========================
required_pkgs <- c("MCMCglmm", "coda")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Please install the required packages first: ",
    paste(missing_pkgs, collapse = ", "),
    "\nExample:\ninstall.packages(c(\"MCMCglmm\", \"coda\"))"
  )
}

library(MCMCglmm)
library(coda)

# ==========================
# MCMC settings
# ==========================
if (quick_test) {
  MCMC_NITT   <- 260000
  MCMC_BURNIN <- 60000
  MCMC_THIN   <- 200
} else {
  MCMC_NITT   <- 2100000
  MCMC_BURNIN <- 100000
  MCMC_THIN   <- 2000
}

# ==========================
# Trait definitions
# ==========================
veg_traits  <- c("Height", "MSL_W", "SB", "MSD")
leaf_traits <- c("Area", "P2A2", "Circularity", "Nindents.Peri",
                 "IndentWidth", "IndentDepth")
leaf_axes   <- c("LeafPC1", "LeafPC2")

# ==========================
# Helper functions
# ==========================
msg <- function(...) cat(..., "\n")

safe_scale <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - m) / s)
}

symmetrise <- function(M) {
  (M + t(M)) / 2
}

make_psd <- function(M, tol = 1e-8) {
  M <- symmetrise(M)
  ee <- eigen(M, symmetric = TRUE)
  vals <- ee$values
  vals[vals < tol] <- tol
  ee$vectors %*% diag(vals, nrow = length(vals)) %*% t(ee$vectors)
}

orient_eigenvector <- function(v) {
  idx <- which.max(abs(v))
  if (length(idx) == 1 && v[idx] < 0) -v else v
}

angle_deg <- function(u, v) {
  u <- as.numeric(u)
  v <- as.numeric(v)
  uu <- sqrt(sum(u^2))
  vv <- sqrt(sum(v^2))
  if (uu == 0 || vv == 0) return(NA_real_)
  cc <- abs(sum(u * v) / (uu * vv))
  cc <- min(1, max(-1, cc))
  acos(cc) * 180 / pi
}

principal_angles_deg <- function(U, V) {
  sv <- svd(t(U) %*% V)$d
  sv[sv > 1] <- 1
  sv[sv < 0] <- 0
  acos(sv) * 180 / pi
}

leading_plane_capture <- function(v, U) {
  # Proportion of a unit vector v captured by the subspace spanned by U
  vv <- as.numeric(v)
  vv <- vv / sqrt(sum(vv^2))
  pr <- t(U) %*% vv
  as.numeric(sum(pr^2))
}

scale_within_group <- function(df, cols, group_col) {
  out <- df
  for (cc in cols) {
    out[[cc]] <- ave(df[[cc]], df[[group_col]], FUN = safe_scale)
  }
  out
}

build_leaf_pca_basis <- function(exp1_df, leaf_cols) {
  X <- scale(as.matrix(exp1_df[, leaf_cols]), center = TRUE, scale = TRUE)
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  list(
    center = attr(X, "scaled:center"),
    scale = attr(X, "scaled:scale"),
    rotation = pca$rotation[, 1:2, drop = FALSE],
    sdev = pca$sdev,
    var_explained = (pca$sdev^2) / sum(pca$sdev^2)
  )
}

project_leaf_scores <- function(df, leaf_cols, basis) {
  X <- as.matrix(df[, leaf_cols])
  X <- sweep(X, 2, basis$center, FUN = "-")
  X <- sweep(X, 2, basis$scale,  FUN = "/")
  scores <- X %*% basis$rotation
  scores <- as.data.frame(scores)
  names(scores) <- leaf_axes
  scores
}

write_named_matrix <- function(M, file) {
  out <- data.frame(Trait = rownames(M), M, row.names = NULL, check.names = FALSE)
  write.csv(out, file, row.names = FALSE)
}

make_pairwise_D <- function(mean_df, block_name = "block") {
  stopifnot("Population" %in% names(mean_df), "Ecotype" %in% names(mean_df))
  trait_cols <- setdiff(names(mean_df), c("Population", "Ecotype"))
  pop_names <- mean_df$Population
  mean_mat  <- as.matrix(mean_df[, trait_cols, drop = FALSE])
  rownames(mean_mat) <- pop_names

  pair_index <- utils::combn(pop_names, 2, simplify = FALSE)
  summary_rows <- vector("list", length(pair_index))
  D_list <- vector("list", length(pair_index))
  d_list <- vector("list", length(pair_index))

  for (ii in seq_along(pair_index)) {
    pp <- pair_index[[ii]]
    p1 <- pp[1]
    p2 <- pp[2]
    d  <- mean_mat[p1, ] - mean_mat[p2, ]
    dnorm <- sqrt(sum(d^2))
    u <- if (dnorm > 0) d / dnorm else rep(0, length(d))
    Dmat <- tcrossprod(d)
    key <- paste(p1, p2, sep = "__")

    row <- data.frame(
      block = block_name,
      pop1 = p1,
      pop2 = p2,
      ecotype1 = mean_df$Ecotype[match(p1, mean_df$Population)],
      ecotype2 = mean_df$Ecotype[match(p2, mean_df$Population)],
      delta_norm = dnorm,
      lambda1 = dnorm^2,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    delta_df <- stats::setNames(as.list(d), paste0("delta_", trait_cols))
    axis_df  <- stats::setNames(as.list(u), paste0("axis1_", trait_cols))
    summary_rows[[ii]] <- cbind(row, as.data.frame(delta_df), as.data.frame(axis_df),
                                check.names = FALSE)
    D_list[[key]] <- Dmat
    d_list[[key]] <- u
    names(D_list[[key]]) <- NULL
    names(d_list[[key]]) <- trait_cols
    rownames(D_list[[key]]) <- trait_cols
    colnames(D_list[[key]]) <- trait_cols
  }

  list(
    summary = do.call(rbind, summary_rows),
    matrices = D_list,
    unit_vectors = d_list,
    traits = trait_cols
  )
}

make_prior <- function(n_traits) {
  P <- diag(n_traits)
  list(
    R = list(V = P * 0.25, nu = n_traits - 0.998),
    G = list(
      G1 = list(V = P, nu = n_traits,
                alpha.mu = rep(0, n_traits), alpha.V = P),
      G2 = list(V = P, nu = n_traits,
                alpha.mu = rep(0, n_traits), alpha.V = P),
      G3 = list(V = P, nu = n_traits,
                alpha.mu = rep(0, n_traits), alpha.V = P)
    )
  )
}

fit_population_G <- function(df, population, trait_cols, trait_labels, block_name) {
  dat <- df[df$Block == population, c("Block", "Type", "sire", "dam", trait_cols), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  dat$sire <- factor(dat$sire)
  dat$dam  <- factor(dat$dam)

  safe_names <- paste0("T", seq_along(trait_cols))
  names(dat)[match(trait_cols, names(dat))] <- safe_names

  fixed_formula <- as.formula(
    paste0("cbind(", paste(safe_names, collapse = ", "), ") ~ trait - 1")
  )

  model <- MCMCglmm(
    fixed  = fixed_formula,
    random = ~ us(trait):sire + us(trait):dam + us(trait):sire:dam,
    rcov   = ~ us(trait):units,
    family = rep("gaussian", length(safe_names)),
    data   = dat,
    prior  = make_prior(length(safe_names)),
    nitt   = MCMC_NITT,
    burnin = MCMC_BURNIN,
    thin   = MCMC_THIN,
    pr     = TRUE,
    verbose = TRUE
  )

  list(
    model = model,
    data = dat,
    safe_names = safe_names,
    trait_labels = trait_labels,
    population = population,
    block_name = block_name
  )
}

extract_G_posterior <- function(model, n_traits) {
  # Because us(trait):sire is the first random term, the first n_traits^2 columns
  # correspond to the sire-level covariance matrix, matching the supplementary
  # code of Walter et al.
  vcv <- model$VCV
  n_keep <- n_traits^2
  if (ncol(vcv) < n_keep) {
    stop("VCV matrix does not contain enough columns to reconstruct the sire covariance matrix.")
  }

  G_array <- array(NA_real_, dim = c(n_traits, n_traits, nrow(vcv)))
  for (ii in seq_len(nrow(vcv))) {
    SireMat <- matrix(as.numeric(vcv[ii, 1:n_keep]), ncol = n_traits)
    Gmat <- 4 * symmetrise(SireMat)
    G_array[, , ii] <- make_psd(Gmat)
  }
  G_mean <- apply(G_array, c(1, 2), mean)
  list(G_array = G_array, G_mean = make_psd(G_mean))
}

summarise_G <- function(G_post, fit_obj) {
  G_array <- G_post$G_array
  G_mean  <- G_post$G_mean
  n_traits <- dim(G_array)[1]
  trait_labels <- fit_obj$trait_labels
  pop <- fit_obj$population
  block_name <- fit_obj$block_name

  rownames(G_mean) <- trait_labels
  colnames(G_mean) <- trait_labels

  ee_mean <- eigen(G_mean, symmetric = TRUE)
  ord <- order(ee_mean$values, decreasing = TRUE)
  vals <- ee_mean$values[ord]
  vecs <- ee_mean$vectors[, ord, drop = FALSE]
  vecs <- apply(vecs, 2, orient_eigenvector)
  if (is.vector(vecs)) vecs <- matrix(vecs, ncol = 1)
  colnames(vecs) <- paste0("g", seq_len(ncol(vecs)))
  rownames(vecs) <- trait_labels

  eig_samples <- t(vapply(seq_len(dim(G_array)[3]), function(i) {
    sort(eigen(G_array[, , i], symmetric = TRUE, only.values = TRUE)$values,
         decreasing = TRUE)
  }, numeric(n_traits)))

  n_keep <- min(2, n_traits)
  total_var <- sum(vals)
  lambda1 <- vals[1]
  lambda2 <- if (n_traits >= 2) vals[2] else NA_real_

  summary_df <- data.frame(
    block = block_name,
    population = pop,
    ecotype = unique(fit_obj$data$Type)[1],
    n_offspring = nrow(fit_obj$data),
    n_sire = nlevels(fit_obj$data$sire),
    n_dam = nlevels(fit_obj$data$dam),
    n_traits = n_traits,
    total_genetic_variance = total_var,
    lambda1 = lambda1,
    lambda1_lwr = stats::quantile(eig_samples[, 1], 0.025, na.rm = TRUE),
    lambda1_upr = stats::quantile(eig_samples[, 1], 0.975, na.rm = TRUE),
    lambda2 = lambda2,
    lambda2_lwr = if (n_traits >= 2) stats::quantile(eig_samples[, 2], 0.025, na.rm = TRUE) else NA_real_,
    lambda2_upr = if (n_traits >= 2) stats::quantile(eig_samples[, 2], 0.975, na.rm = TRUE) else NA_real_,
    prop1 = lambda1 / total_var,
    prop2 = if (n_traits >= 2) lambda2 / total_var else NA_real_,
    leading_plane_prop = sum(vals[1:n_keep]) / total_var,
    spectral_gap = if (n_traits >= 2) (lambda1 - lambda2) else NA_real_,
    lambda1_over_lambda2 = if (n_traits >= 2) (lambda1 / lambda2) else NA_real_,
    isotropy12 = if (n_traits >= 2) (lambda2 / lambda1) else NA_real_,
    stringsAsFactors = FALSE
  )

  loadings_rows <- vector("list", n_keep)
  for (kk in seq_len(n_keep)) {
    loadings_rows[[kk]] <- data.frame(
      block = block_name,
      population = pop,
      ecotype = unique(fit_obj$data$Type)[1],
      axis = paste0("g", kk),
      trait = trait_labels,
      loading = vecs[, kk],
      stringsAsFactors = FALSE
    )
  }

  ac_diag <- as.data.frame(coda::autocorr.diag(as.mcmc(fit_obj$model$VCV), lags = 1))
  ess <- coda::effectiveSize(as.mcmc(fit_obj$model$VCV))
  diagnostics <- data.frame(
    block = block_name,
    population = pop,
    ecotype = unique(fit_obj$data$Type)[1],
    parameter = colnames(fit_obj$model$VCV),
    lag1_autocorr = as.numeric(ac_diag[1, ]),
    effective_size = as.numeric(ess),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_df,
    loadings = do.call(rbind, loadings_rows),
    diagnostics = diagnostics,
    G_mean = G_mean,
    eig_values = vals,
    eig_vectors = vecs,
    G_array = G_array
  )
}

run_G_block_models <- function(df, populations, trait_cols, trait_labels,
                               block_name, outdir) {
  summary_rows <- list()
  loading_rows <- list()
  diagnostic_rows <- list()
  G_means <- list()
  eig_vectors <- list()
  eig_values <- list()

  for (pop in populations) {
    msg("Fitting", block_name, "G for population", pop)
    fit_obj <- fit_population_G(df, pop, trait_cols, trait_labels, block_name)
    if (save_model_objects) {
      saveRDS(fit_obj$model,
              file = file.path(outdir, paste0(block_name, "_model_", pop, ".rds")))
    }

    G_post <- extract_G_posterior(fit_obj$model, length(trait_cols))
    out <- summarise_G(G_post, fit_obj)

    summary_rows[[pop]] <- out$summary
    loading_rows[[pop]] <- out$loadings
    diagnostic_rows[[pop]] <- out$diagnostics
    G_means[[pop]] <- out$G_mean
    eig_vectors[[pop]] <- out$eig_vectors
    eig_values[[pop]] <- out$eig_values

    write_named_matrix(out$G_mean,
                       file = file.path(outdir, paste0(block_name, "_Gmean_", pop, ".csv")))

    if (save_posterior_G) {
      saveRDS(out$G_array,
              file = file.path(outdir, paste0(block_name, "_Gposterior_", pop, ".rds")))
    }

    rm(fit_obj, G_post, out)
    gc(verbose = FALSE)
  }

  summary_df <- do.call(rbind, summary_rows)
  loadings_df <- do.call(rbind, loading_rows)
  diagnostics_df <- do.call(rbind, diagnostic_rows)

  write.csv(summary_df,
            file = file.path(outdir, paste0(block_name, "_G_summary.csv")),
            row.names = FALSE)
  write.csv(loadings_df,
            file = file.path(outdir, paste0(block_name, "_G_eigenvectors.csv")),
            row.names = FALSE)
  write.csv(diagnostics_df,
            file = file.path(outdir, paste0(block_name, "_G_diagnostics.csv")),
            row.names = FALSE)
  saveRDS(G_means,
          file = file.path(outdir, paste0(block_name, "_Gmean_list.rds")))

  list(
    summary = summary_df,
    loadings = loadings_df,
    diagnostics = diagnostics_df,
    G_means = G_means,
    eig_vectors = eig_vectors,
    eig_values = eig_values,
    traits = trait_labels
  )
}

pairwise_G_geometry <- function(G_results, block_name) {
  pops <- names(G_results$eig_vectors)
  pairs <- utils::combn(pops, 2, simplify = FALSE)
  rows <- vector("list", length(pairs))

  for (ii in seq_along(pairs)) {
    pp <- pairs[[ii]]
    p1 <- pp[1]
    p2 <- pp[2]
    V1 <- G_results$eig_vectors[[p1]]
    V2 <- G_results$eig_vectors[[p2]]
    lam1 <- G_results$eig_values[[p1]]
    lam2 <- G_results$eig_values[[p2]]

    gmax_angle <- angle_deg(V1[, 1], V2[, 1])
    g2_angle <- if (ncol(V1) >= 2 && ncol(V2) >= 2) angle_deg(V1[, 2], V2[, 2]) else NA_real_
    leading_plane_angle <- if (ncol(V1) >= 2 && ncol(V2) >= 2) {
      max(principal_angles_deg(V1[, 1:2, drop = FALSE], V2[, 1:2, drop = FALSE]))
    } else {
      NA_real_
    }

    rows[[ii]] <- data.frame(
      block = block_name,
      pop1 = p1,
      pop2 = p2,
      gmax_angle_deg = gmax_angle,
      g2_angle_deg = g2_angle,
      leading_plane_angle_deg = leading_plane_angle,
      spectral_gap_pop1 = if (length(lam1) >= 2) lam1[1] - lam1[2] else NA_real_,
      spectral_gap_pop2 = if (length(lam2) >= 2) lam2[1] - lam2[2] else NA_real_,
      lambda1_over_lambda2_pop1 = if (length(lam1) >= 2) lam1[1] / lam1[2] else NA_real_,
      lambda1_over_lambda2_pop2 = if (length(lam2) >= 2) lam2[1] / lam2[2] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

pairwise_D_alignment_to_G <- function(D_result, G_results, block_name) {
  pair_keys <- names(D_result$unit_vectors)
  rows <- list()
  rr <- 1

  for (key in pair_keys) {
    bits <- strsplit(key, "__", fixed = TRUE)[[1]]
    p1 <- bits[1]
    p2 <- bits[2]
    d  <- D_result$unit_vectors[[key]]

    for (pop in c(p1, p2)) {
      V <- G_results$eig_vectors[[pop]]
      ang1 <- angle_deg(d, V[, 1])
      ang2 <- if (ncol(V) >= 2) angle_deg(d, V[, 2]) else NA_real_
      plane_cap <- if (ncol(V) >= 2) {
        leading_plane_capture(d, V[, 1:2, drop = FALSE])
      } else {
        1
      }

      rows[[rr]] <- data.frame(
        block = block_name,
        pair = key,
        focal_population = pop,
        other_population = if (pop == p1) p2 else p1,
        angle_to_gmax_deg = ang1,
        angle_to_g2_deg = ang2,
        leading_plane_capture = plane_cap,
        stringsAsFactors = FALSE
      )
      rr <- rr + 1
    }
  }
  do.call(rbind, rows)
}

# ==========================
# Load and harmonise data
# ==========================
if (!file.exists(exp1_file)) stop("Cannot find ", exp1_file)
if (!file.exists(exp2_file)) stop("Cannot find ", exp2_file)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

msg("Loading data ...")
exp1 <- read.csv(exp1_file, stringsAsFactors = FALSE)
exp2 <- read.csv(exp2_file, stringsAsFactors = FALSE)

# Harmonise names across experiments
exp1$Height      <- exp1$VegHeight
exp2$IndentWidth <- exp2$IndentWidthMean
exp2$IndentDepth <- exp2$IndentDepthMean

needed_exp1 <- c("Ecotype", "Population", veg_traits, leaf_traits)
needed_exp2 <- c("Type", "Block", "sire", "dam", veg_traits, leaf_traits)

missing_exp1 <- setdiff(needed_exp1, names(exp1))
missing_exp2 <- setdiff(needed_exp2, names(exp2))
if (length(missing_exp1) > 0) stop("Missing Exp1 columns: ", paste(missing_exp1, collapse = ", "))
if (length(missing_exp2) > 0) stop("Missing Exp2 columns: ", paste(missing_exp2, collapse = ", "))

breeding_pops <- sort(unique(exp2$Block))
msg("Breeding-design populations:", paste(breeding_pops, collapse = ", "))

exp1_sub <- exp1[exp1$Population %in% breeding_pops,
                 c("Ecotype", "Population", veg_traits, leaf_traits),
                 drop = FALSE]

# ==========================
# D analyses from Exp1
# ==========================
msg("Building global Exp1 standardisation and leaf PCA basis ...")

# Vegetative block for D: globally standardised across the Exp1 populations
exp1_veg_scaled <- as.data.frame(scale(exp1_sub[, veg_traits], center = TRUE, scale = TRUE))
exp1_veg_for_D <- cbind(exp1_sub[, c("Ecotype", "Population")], exp1_veg_scaled)
veg_means <- aggregate(exp1_veg_for_D[, veg_traits],
                       by = exp1_veg_for_D[c("Population", "Ecotype")],
                       FUN = mean)

# Leaf PCA basis: defined once on Exp1 across the same populations present in Exp2
leaf_basis <- build_leaf_pca_basis(exp1_sub, leaf_traits)
exp1_leaf_scores <- project_leaf_scores(exp1_sub, leaf_traits, leaf_basis)
exp1_leaf_for_D <- cbind(exp1_sub[, c("Ecotype", "Population")], exp1_leaf_scores)
leaf_means <- aggregate(exp1_leaf_for_D[, leaf_axes],
                        by = exp1_leaf_for_D[c("Population", "Ecotype")],
                        FUN = mean)

# Save leaf PCA basis
leaf_loading_df <- data.frame(
  trait = rownames(leaf_basis$rotation),
  LeafPC1_loading = leaf_basis$rotation[, 1],
  LeafPC2_loading = leaf_basis$rotation[, 2],
  stringsAsFactors = FALSE
)
leaf_var_df <- data.frame(
  PC = paste0("PC", seq_along(leaf_basis$var_explained)),
  variance_explained = leaf_basis$var_explained,
  cumulative_variance = cumsum(leaf_basis$var_explained),
  stringsAsFactors = FALSE
)
write.csv(leaf_loading_df, file = file.path(outdir, "leaf_PCA_loadings.csv"), row.names = FALSE)
write.csv(leaf_var_df, file = file.path(outdir, "leaf_PCA_variance_explained.csv"), row.names = FALSE)
write.csv(veg_means, file = file.path(outdir, "Exp1_population_means_vegetative.csv"), row.names = FALSE)
write.csv(leaf_means, file = file.path(outdir, "Exp1_population_means_leafPCs.csv"), row.names = FALSE)

# Pairwise D
msg("Computing pairwise D matrices ...")
D_veg  <- make_pairwise_D(veg_means,  block_name = "vegetative")
D_leaf <- make_pairwise_D(leaf_means, block_name = "leaf_PCs")
write.csv(D_veg$summary,  file = file.path(outdir, "pairwise_D_vegetative_summary.csv"), row.names = FALSE)
write.csv(D_leaf$summary, file = file.path(outdir, "pairwise_D_leafPC_summary.csv"), row.names = FALSE)
saveRDS(D_veg$matrices,   file = file.path(outdir, "pairwise_D_vegetative_matrices.rds"))
saveRDS(D_leaf$matrices,  file = file.path(outdir, "pairwise_D_leafPC_matrices.rds"))

# ==========================
# Prepare Exp2 for G analyses
# ==========================
msg("Projecting Exp2 leaf scores onto the Exp1 leaf basis ...")
exp2_leaf_scores <- project_leaf_scores(exp2, leaf_traits, leaf_basis)
exp2_work <- cbind(exp2, exp2_leaf_scores)

if (scale_within_population_for_G) {
  msg("Scaling Exp2 traits within population before fitting G ...")
  exp2_work <- scale_within_group(exp2_work, veg_traits,  "Block")
  exp2_work <- scale_within_group(exp2_work, leaf_axes,   "Block")
  veg_traits_for_G  <- veg_traits
  leaf_traits_for_G <- leaf_axes
} else {
  msg("Scaling Exp2 traits globally before fitting G ...")
  exp2_work[, veg_traits] <- lapply(exp2_work[, veg_traits, drop = FALSE], safe_scale)
  exp2_work[, leaf_axes]  <- lapply(exp2_work[, leaf_axes,  drop = FALSE], safe_scale)
  veg_traits_for_G  <- veg_traits
  leaf_traits_for_G <- leaf_axes
}

# ==========================
# Fit G models
# ==========================
msg("Starting vegetative G models ...")
veg_G_results <- run_G_block_models(
  df = exp2_work,
  populations = breeding_pops,
  trait_cols = veg_traits_for_G,
  trait_labels = veg_traits,
  block_name = "vegetative",
  outdir = outdir
)

msg("Starting leaf-axis G models ...")
leaf_G_results <- run_G_block_models(
  df = exp2_work,
  populations = breeding_pops,
  trait_cols = leaf_traits_for_G,
  trait_labels = leaf_axes,
  block_name = "leaf_PCs",
  outdir = outdir
)

# ==========================
# Pairwise geometry among G matrices
# ==========================
msg("Computing pairwise G geometry ...")
pairwise_G_veg  <- pairwise_G_geometry(veg_G_results,  block_name = "vegetative")
pairwise_G_leaf <- pairwise_G_geometry(leaf_G_results, block_name = "leaf_PCs")
write.csv(pairwise_G_veg,  file = file.path(outdir, "pairwise_G_geometry_vegetative.csv"), row.names = FALSE)
write.csv(pairwise_G_leaf, file = file.path(outdir, "pairwise_G_geometry_leafPC.csv"), row.names = FALSE)

# ==========================
# Alignment of pairwise D with population G
# ==========================
msg("Computing D-to-G alignment summaries ...")
align_veg  <- pairwise_D_alignment_to_G(D_veg,  veg_G_results,  block_name = "vegetative")
align_leaf <- pairwise_D_alignment_to_G(D_leaf, leaf_G_results, block_name = "leaf_PCs")
write.csv(align_veg,  file = file.path(outdir, "pairwise_D_alignment_to_G_vegetative.csv"), row.names = FALSE)
write.csv(align_leaf, file = file.path(outdir, "pairwise_D_alignment_to_G_leafPC.csv"), row.names = FALSE)

# ==========================
# Save a compact results object
# ==========================
results_object <- list(
  settings = list(
    quick_test = quick_test,
    scale_within_population_for_G = scale_within_population_for_G,
    MCMC_NITT = MCMC_NITT,
    MCMC_BURNIN = MCMC_BURNIN,
    MCMC_THIN = MCMC_THIN,
    populations = breeding_pops,
    veg_traits = veg_traits,
    leaf_traits = leaf_traits,
    leaf_axes = leaf_axes
  ),
  leaf_basis = leaf_basis,
  D_vegetative = D_veg,
  D_leaf = D_leaf,
  G_vegetative = list(
    summary = veg_G_results$summary,
    loadings = veg_G_results$loadings,
    pairwise_geometry = pairwise_G_veg,
    alignment = align_veg,
    G_means = veg_G_results$G_means
  ),
  G_leaf = list(
    summary = leaf_G_results$summary,
    loadings = leaf_G_results$loadings,
    pairwise_geometry = pairwise_G_leaf,
    alignment = align_leaf,
    G_means = leaf_G_results$G_means
  )
)
saveRDS(results_object, file = file.path(outdir, "population_G_D_results_compact.rds"))

msg("Done.")
msg("Outputs written to:", normalizePath(outdir))
