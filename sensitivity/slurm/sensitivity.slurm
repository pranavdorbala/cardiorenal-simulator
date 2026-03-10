#!/bin/bash
#SBATCH --job-name=sensitivity-analysis
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sensitivity_%j.out
#SBATCH --error=logs/sensitivity_%j.err
#
# Run sensitivity analysis (Morris + Sobol) for the Hallow model.
#
# Usage:
#   sbatch sensitivity/slurm/sensitivity.slurm
#
#   # Morris only (faster):
#   METHOD=morris sbatch sensitivity/slurm/sensitivity.slurm
#
#   # Override parameters:
#   MORRIS_R=30 SOBOL_N=2048 sbatch sensitivity/slurm/sensitivity.slurm
# ============================================================================

set -euo pipefail

# ---- Configuration ----
METHOD="${METHOD:-full}"
MORRIS_R="${MORRIS_R:-20}"
SOBOL_N="${SOBOL_N:-1024}"
TOP_K="${TOP_K:-50}"
T_HOURS="${T_HOURS:-168}"
SEED="${SEED:-42}"
LEVELS="${LEVELS:-4}"

PROJECT_DIR="${PROJECT_DIR:-$(dirname $(dirname $(realpath $0)))}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/results}"

echo "============================================"
echo "Sensitivity Analysis"
echo "============================================"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "CPUs:         ${SLURM_CPUS_PER_TASK}"
echo "Method:       ${METHOD}"
echo "Morris r:     ${MORRIS_R}"
echo "Sobol N:      ${SOBOL_N}"
echo "Top-k:        ${TOP_K}"
echo "Duration:     ${T_HOURS}h"
echo "Output:       ${OUTPUT_DIR}"
echo "============================================"

mkdir -p logs "${OUTPUT_DIR}"

# Ensure C library is compiled
cd "${PROJECT_DIR}/.."
if [ ! -f "hallow_rhs.so" ]; then
    echo "C shared library not found - running build.sh"
    bash build.sh
fi

echo ""
echo "Launching sensitivity analysis..."
echo ""

python3 -m sensitivity.run_analysis \
    --method "${METHOD}" \
    --r "${MORRIS_R}" \
    --sobol_n "${SOBOL_N}" \
    --top_k "${TOP_K}" \
    --t_hours "${T_HOURS}" \
    --levels "${LEVELS}" \
    --seed "${SEED}" \
    --n_workers "${SLURM_CPUS_PER_TASK}" \
    --output_dir "${OUTPUT_DIR}"

echo ""
echo "============================================"
echo "Sensitivity analysis complete"
echo "Results: ${OUTPUT_DIR}"
echo "============================================"
