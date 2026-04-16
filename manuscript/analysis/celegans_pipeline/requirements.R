# R dependencies for the C. elegans channel-operator pipeline
# =============================================================
# 
# Only BASE R is required.
# 
# The single R script in this pipeline (scripts/01_parse/01c_extract_posterior_chains.R)
# reads the MCMCglmm posterior chains from RData files and exports them to CSV.
# It uses only:
#   - load()          -- base R
#   - matrix()        -- base R
#   - writeLines()    -- base R
#   - sprintf()       -- base R
#
# No CRAN packages (including MCMCglmm itself) need to be installed,
# because we are not refitting the model -- we are reading already-fitted
# posterior chains stored as matrices.
#
# Tested on:  R 4.3.3 (2024-02-29) "Angel Food Cake"
# Expected compatible with: R >= 4.0
#
# To install R:
#   macOS:       brew install r        (or https://cran.r-project.org)
#   Ubuntu:      apt install r-base-core
#   Windows:     download from https://cran.r-project.org
