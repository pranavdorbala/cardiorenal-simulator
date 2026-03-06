#!/bin/bash
#SBATCH --job-name=node-cr-datagen
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/datagen_%j.out
#SBATCH --error=logs/datagen_%j.err
#
# Generate training data for the Neural ODE cardiorenal surrogate.
#
# This runs the C-compiled mechanistic solver across a Latin Hypercube
# sample of the parameter space. CPU-only (no GPU needed).
#
# Usage:
#   sbatch generate_data.slurm
#
#   # Override sample count:
#   N_SAMPLES=10000 sbatch generate_data.slurm
# ============================================================================

set -euo pipefail

# ---- Configuration ----
N_SAMPLES="${N_SAMPLES:-5000}"
T_HOURS="${T_HOURS:-168}"
DT_OUTPUT="${DT_OUTPUT:-1.0}"
SEED="${SEED:-42}"
TIMEOUT="${TIMEOUT:-120}"

PROJECT_DIR="${PROJECT_DIR:-$(dirname $(dirname $(realpath $0)))}"
OUTPUT_PATH="${OUTPUT_PATH:-${PROJECT_DIR}/training_data.h5}"

echo "============================================"
echo "Training Data Generation"
echo "============================================"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "CPUs:       ${SLURM_CPUS_PER_TASK}"
echo "Samples:    ${N_SAMPLES}"
echo "Duration:   ${T_HOURS}h per trajectory"
echo "Output:     ${OUTPUT_PATH}"
echo "============================================"

mkdir -p logs

# Load modules (adjust for your cluster)
# module purge
# module load python/3.10 gcc/11.3

# Ensure C library is compiled
cd "${PROJECT_DIR}/.."
if [ ! -f "${PROJECT_DIR}/../hallow_rhs.so" ]; then
    echo "C shared library not found — running build.sh"
    bash build.sh
fi

echo ""
echo "Launching data generation..."
echo ""

python3 -m neural_surrogate.generate_training_data \
    --n_samples "${N_SAMPLES}" \
    --t_hours "${T_HOURS}" \
    --dt_output "${DT_OUTPUT}" \
    --seed "${SEED}" \
    --output "${OUTPUT_PATH}"

echo ""
echo "============================================"
echo "Data generation complete"
echo "Output: ${OUTPUT_PATH}"
echo "============================================"
