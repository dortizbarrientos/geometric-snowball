# requirements.R
#
# Install all R packages required by the Senecio channel-operator pipeline.
#
# Usage: Rscript requirements.R
#
# Packages by stage:
#   Stage 1 (population_G_D_snowball_analysis.R): MCMCglmm, coda
#   Stage 2 (senecio_channel_operator_analysis.R): readxl, dplyr, tidyr,
#                                                  ggplot2, expm, Matrix
#
# Note: MCMCglmm has Matrix and coda as hard dependencies; we install them
# first explicitly to avoid subtle compile-time errors on fresh systems.

required_packages <- c(
  # Hard dependencies first
  "Matrix",
  "coda",
  # Stage 1
  "MCMCglmm",
  # Stage 2
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "expm"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)

if (length(to_install) == 0) {
  cat("All required R packages already installed.\n")
} else {
  cat("Installing missing packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

# Verification
missing <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing) > 0) {
  stop("Failed to install: ", paste(missing, collapse = ", "),
       "\nOn macOS, ensure Xcode command-line tools are installed first:\n",
       "  xcode-select --install")
} else {
  cat("\nAll packages installed successfully. Pipeline is ready.\n")
}
