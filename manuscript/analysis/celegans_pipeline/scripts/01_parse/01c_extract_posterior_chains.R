#!/usr/bin/env Rscript
# ============================================================================
# 01c_extract_posterior_chains.R
# ============================================================================
#
# PURPOSE
# -------
# Read the MCMCglmm posterior chains from the Mallard et al. RData files
# and export them as language-neutral CSVs. After this step, the full
# analysis (scripts 02_analysis/*) can run entirely in Python without any
# R dependency.
#
# INPUT FILES (data/raw/)
# -----------------------
#   G_matrices_high_salt.RData    ~ 44 MB
#   G_matrices_low_salt.RData     ~ 44 MB
#
# Each file contains a single top-level object (VCV_mat_NaCl or
# VCV_mat_NGM) that is a list of 4 populations. Each population entry has:
#   Population     character          "A6140", "GA150", "GA250", "GA450"
#   N_measurement  integer            sample size
#   G1_mat         7x7 matrix         posterior-mode additive genetic G
#   G2_mat         7x7 matrix         posterior-mode date-block G
#   R_mat          7x7 matrix         posterior-mode residual covariance
#   VCV_Mat        mcmc 10000 x 147   the full posterior chain
#
# The 147 columns of VCV_Mat decompose as:
#   columns   1- 49  "pop_label"  -> additive genetic G  (7x7 stacked column-major)
#   columns  50- 98  "date_str"   -> date-block random effect
#   columns  99-147  "units"      -> residual covariance
#
# Only the first 49 columns (the additive genetic G) are needed for the
# channel-operator analysis. We export those 10,000 posterior draws per
# population to a long-format CSV.
#
# OUTPUTS (data/processed/)
# -------------------------
#   G_chains_high_salt_long.csv
#     Columns: population, draw, row_idx, col_idx, row_trait, col_trait, value
#     Rows: 4 pops * 10,000 draws * 49 cells = 1,960,000 rows (~120 MB)
#
#   G_chains_low_salt_long.csv    (same structure, low-salt environment)
#
# If the downstream analysis only ever needs one environment, comment out
# the other load call at the bottom of this script.

suppressPackageStartupMessages(library(stats))

args <- commandArgs(trailingOnly = TRUE)

HERE <- if (length(args) >= 1) args[1] else {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(sub("--file=", "", file_arg))
  } else getwd()
}

ROOT <- normalizePath(file.path(HERE, "..", ".."))
RAW  <- file.path(ROOT, "data", "raw")
OUT  <- file.path(ROOT, "data", "processed")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("SF", "SB", "FS", "FB", "BS", "BF", "Size")
P      <- length(TRAITS)        # 7
G_CHAIN_NCOLS <- P * P          # 49

# ----------------------------------------------------------------------
# GENERIC EXTRACTOR
# ----------------------------------------------------------------------
extract_chains <- function(rdata_path, list_name, out_csv) {
  cat("\nExtracting", rdata_path, "->", out_csv, "\n")

  env <- new.env()
  load(rdata_path, envir = env)
  if (!exists(list_name, envir = env)) {
    stop("Object '", list_name, "' not found in ", rdata_path)
  }
  vcv_list <- get(list_name, envir = env)

  # Prepare per-cell trait labels for long output (column-major order of a 7x7)
  row_idx <- rep(seq_len(P), times = P) - 1      # 0-indexed
  col_idx <- rep(seq_len(P), each  = P) - 1
  row_trait <- TRAITS[row_idx + 1]
  col_trait <- TRAITS[col_idx + 1]

  # Write CSV incrementally to avoid building a huge in-memory frame
  con <- file(out_csv, "w")
  writeLines("population,draw,row_idx,col_idx,row_trait,col_trait,value", con)

  for (i in seq_along(vcv_list)) {
    entry <- vcv_list[[i]]
    pop <- entry$Population
    vcv_mat <- entry$VCV_Mat
    n_draws <- nrow(vcv_mat)

    cat(sprintf("  %s: %d draws x %d cells\n", pop, n_draws, G_CHAIN_NCOLS))

    G_chain <- vcv_mat[, 1:G_CHAIN_NCOLS]  # additive genetic G only

    for (d in seq_len(n_draws)) {
      row_strs <- sprintf(
        "%s,%d,%d,%d,%s,%s,%.10g",
        pop, d,
        row_idx, col_idx,
        row_trait, col_trait,
        G_chain[d, ]
      )
      writeLines(row_strs, con)
    }
  }

  close(con)
  cat("  wrote", out_csv, "\n")
}

# ----------------------------------------------------------------------
# EXTRACT BOTH ENVIRONMENTS
# ----------------------------------------------------------------------
high_salt_path <- file.path(RAW, "G_matrices_high_salt.RData")
low_salt_path  <- file.path(RAW, "G_matrices_low_salt.RData")

if (file.exists(high_salt_path)) {
  extract_chains(
    rdata_path = high_salt_path,
    list_name  = "VCV_mat_NaCl",
    out_csv    = file.path(OUT, "G_chains_high_salt_long.csv")
  )
} else {
  cat("Skipping high salt: file not found\n")
}

if (file.exists(low_salt_path)) {
  # NOTE: object name in the low-salt file is "VCV_mat_NGM"
  extract_chains(
    rdata_path = low_salt_path,
    list_name  = "VCV_mat_NGM",
    out_csv    = file.path(OUT, "G_chains_low_salt_long.csv")
  )
} else {
  cat("Skipping low salt: file not found\n")
}

cat("\nAll posterior chains extracted to CSV.\n")
cat("Downstream scripts can now run without any R dependency.\n")
