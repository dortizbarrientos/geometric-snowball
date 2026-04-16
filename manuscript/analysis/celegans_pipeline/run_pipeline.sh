#!/usr/bin/env bash
# ============================================================================
# run_pipeline.sh
# ============================================================================
#
# Execute the full C. elegans posterior analysis pipeline end-to-end.
# Run this from the repository root:
#
#   bash run_pipeline.sh
#
# The script stops on the first error. Each stage writes into data/processed
# or results/ and is idempotent: re-running is safe.
#
# STAGES
# ------
#   01a  Parse selection tables            (Python)
#   01b  Parse posterior-mode G matrices    (Python)
#   01c  Extract MCMC chains to CSV         (R; the only R stage)
#   02a  Main posterior analysis            (Python)
#   03a  Build figure                       (Python)
#   03b  Build Word + LaTeX summary table  (Python)

set -euo pipefail

# Ensure we are running from the repo root
cd "$(dirname "$0")"
ROOT="$(pwd)"
echo "Repo root: $ROOT"
echo

# ------------------------------------------------------------
# Stage 01: parsing
# ------------------------------------------------------------
echo "=== Stage 01a: parse selection tables ==="
python3 scripts/01_parse/01a_parse_selection_tables.py

echo
echo "=== Stage 01b: parse posterior-mode G matrices ==="
python3 scripts/01_parse/01b_parse_posterior_mode_matrices.py

echo
echo "=== Stage 01c: extract MCMC chains from RData (this takes ~30 s per environment) ==="
if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found. Install R (e.g. 'apt install r-base-core' or"
    echo "       see https://cran.r-project.org), then re-run the pipeline."
    exit 1
fi
Rscript scripts/01_parse/01c_extract_posterior_chains.R scripts/01_parse

# ------------------------------------------------------------
# Stage 02: analysis
# ------------------------------------------------------------
echo
echo "=== Stage 02a: main posterior analysis ==="
python3 scripts/02_analysis/02a_posterior_analysis.py

# ------------------------------------------------------------
# Stage 03: figures and tables
# ------------------------------------------------------------
echo
echo "=== Stage 03a: build figure ==="
python3 scripts/03_figures/03a_figure_posterior_violins.py

echo
echo "=== Stage 03b: build summary table ==="
python3 scripts/03_figures/03b_table_celegans_summary.py

# ------------------------------------------------------------
echo
echo "=== Pipeline complete ==="
echo "Outputs:"
echo "  results/figures/celegans_posterior_violins.{pdf,png}"
echo "  results/tables/celegans_posterior_summary.csv"
echo "  results/tables/table_celegans_summary.{docx,tex}"
echo "  data/processed/*.csv  (intermediate)"
