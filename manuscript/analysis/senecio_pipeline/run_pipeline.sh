#!/usr/bin/env bash
#
# Senecio channel-operator pipeline — end-to-end runner
# =====================================================
# Wraps the two validated R scripts from scripts/gmatricesDH/ in a clean
# reproducible pipeline plus a thin Python figure layer. Three stages:
#
#   Stage 1 (01_fit_G):
#     Reads Data_Exp1_Dmatrix.csv, Data_Exp2_Gvariance.csv
#     Fits population-specific G matrices for 8 populations (D1,D3,H1,H2,T1,
#     T9,W2,W4) in two trait blocks (4x4 vegetative, 2x2 leaf PCs).
#     Computes pairwise D and pairwise G geometry.
#     Outputs to data/processed/population_G_D_results/
#
#   Stage 2 (02_channel_operator):
#     Reads Stage 1 outputs plus Walter_Data_Exp2.csv and soil data.
#     Fits Lande-Arnold selection gradients (beta, gamma) per pop x env.
#     Constructs H from quadratic gradients.
#     Computes S = G^{1/2} H G^{1/2} and extracts channel angles.
#     Cross-validates with soil-based rank-1 H proxy.
#     Outputs to data/processed/channel_operator_results/
#
#   Stage 3 (03_figures):
#     Reads Stage 1 and 2 CSVs. Produces publication-quality summary figure
#     and formatted table. Outputs to results/figures/ and results/tables/.
#
# USAGE
#   bash run_pipeline.sh                 # full run, hours
#   bash run_pipeline.sh --quick-test    # short MCMC chains, ~5 min
#   bash run_pipeline.sh --stage 1       # run only Stage 1
#   bash run_pipeline.sh --stage 2       # run only Stage 2
#   bash run_pipeline.sh --stage 3       # figures only (Stages 1+2 must exist)
#   bash run_pipeline.sh --quick-test --stage 1
#
# REQUIREMENTS
#   R >= 4.0 with packages MCMCglmm, coda, readxl, dplyr, tidyr, ggplot2, expm
#   Python >= 3.9 with numpy, pandas, matplotlib (Stage 3 only)
#   Install via: Rscript -e 'source("requirements.R")'
#                pip install -r requirements.txt
#
# ------------------------------------------------------------------------------

set -euo pipefail

# Pipeline root = directory containing this script
PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_ROOT"

# Default arguments
QUICK_TEST=0
STAGES="1 2 3"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --quick-test)
      QUICK_TEST=1
      shift
      ;;
    --stage)
      STAGES="$2"
      shift 2
      ;;
    -h|--help)
      sed -n '2,40p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

echo "=============================================================="
echo " Senecio channel-operator pipeline"
echo "=============================================================="
echo " Pipeline root : $PIPELINE_ROOT"
echo " Quick test    : $([ $QUICK_TEST -eq 1 ] && echo 'YES (short MCMC)' || echo 'NO (full MCMC)')"
echo " Stages to run : $STAGES"
echo "=============================================================="

# Sanity check: required raw data present
REQUIRED_FILES=(
  "data/raw/Data_Exp1_Dmatrix.csv"
  "data/raw/Data_Exp2_Gvariance.csv"
  "data/raw/Walter_Data_Exp2.csv"
  "data/raw/All_S__lautus_soil_data.xlsx"
)
MISSING=0
for f in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "MISSING: $f" >&2
    MISSING=$((MISSING + 1))
  fi
done
if [[ $MISSING -gt 0 ]]; then
  echo "Required raw data missing ($MISSING files). See data/raw/README.md." >&2
  exit 2
fi
echo "All raw data files present ✓"

# Ensure processed-data directory exists
mkdir -p data/processed

# Helper: run a stage with a temporary quick_test patch if requested
run_stage_1() {
  echo ""
  echo "=============================================================="
  echo " Stage 1: population-specific G and pairwise D"
  echo "=============================================================="
  cd "$PIPELINE_ROOT/scripts/01_fit_G"
  if [[ $QUICK_TEST -eq 1 ]]; then
    # In-place regex swap of quick_test <- FALSE to TRUE for this run,
    # then restore on exit.
    cp population_G_D_snowball_analysis.R population_G_D_snowball_analysis.R.bak
    sed -i.tmp 's/^quick_test <- FALSE$/quick_test <- TRUE/' population_G_D_snowball_analysis.R
    rm -f population_G_D_snowball_analysis.R.tmp
    trap 'mv -f population_G_D_snowball_analysis.R.bak population_G_D_snowball_analysis.R 2>/dev/null || true' EXIT
  fi
  Rscript population_G_D_snowball_analysis.R
  if [[ $QUICK_TEST -eq 1 ]]; then
    mv -f population_G_D_snowball_analysis.R.bak population_G_D_snowball_analysis.R
    trap - EXIT
  fi
  cd "$PIPELINE_ROOT"
  echo "Stage 1 complete ✓"
}

run_stage_2() {
  echo ""
  echo "=============================================================="
  echo " Stage 2: channel operator S = G^{1/2} H G^{1/2}"
  echo "=============================================================="
  cd "$PIPELINE_ROOT/scripts/02_channel_operator"
  Rscript senecio_channel_operator_analysis.R
  cd "$PIPELINE_ROOT"
  echo "Stage 2 complete ✓"
}

run_stage_3() {
  echo ""
  echo "=============================================================="
  echo " Stage 3: summary figures and tables"
  echo "=============================================================="
  cd "$PIPELINE_ROOT"
  python3 scripts/03_figures/make_figures.py
  echo "Stage 3 complete ✓"
}

# Dispatch
for STAGE in $STAGES; do
  case "$STAGE" in
    1) run_stage_1 ;;
    2) run_stage_2 ;;
    3) run_stage_3 ;;
    *) echo "Unknown stage: $STAGE (expected 1, 2, or 3)" >&2; exit 1 ;;
  esac
done

echo ""
echo "=============================================================="
echo " Pipeline complete"
echo "=============================================================="
echo " Outputs:"
echo "   data/processed/population_G_D_results/      (Stage 1)"
echo "   data/processed/channel_operator_results/    (Stage 2)"
echo "   results/figures/                            (Stage 3)"
echo "   results/tables/                             (Stage 3)"
echo ""
echo " Key files to inspect:"
echo "   - population_G_D_results/vegetative_G_summary.csv"
echo "   - population_G_D_results/pairwise_G_geometry_leafPC.csv"
echo "   - channel_operator_results/lande_arnold_beta_all.csv"
echo "   - channel_operator_results/S_channel_angles_LA.csv"
echo ""
echo " See docs/expected_outputs.md for validated reference numbers."
echo "=============================================================="
