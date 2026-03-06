#!/bin/bash
#SBATCH --job-name=node-cr-train
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/train_%j.out
#SBATCH --error=logs/train_%j.err
#SBATCH --signal=B:SIGUSR1@120
#SBATCH --requeue
#
# Neural ODE Cardiorenal Surrogate — Training on HPC
#
# Usage:
#   # Single-node, 4 GPUs:
#   sbatch train.slurm
#
#   # Override defaults:
#   sbatch --gpus-per-node=2 --time=48:00:00 train.slurm
#
#   # Multi-node (2 nodes × 4 GPUs = 8 GPUs):
#   sbatch --nodes=2 train.slurm
#
#   # CPU-only (no GPUs):
#   sbatch --partition=cpu --gpus-per-node=0 train.slurm
#
# The script automatically handles:
#   - Preemption: saves checkpoint on SIGUSR1, job is requeued
#   - Resume: detects existing checkpoints and resumes training
#   - Multi-node: sets up MASTER_ADDR/PORT for DDP
# ============================================================================

set -euo pipefail

# ---- Configuration (modify these or override via environment) ----
DATA_PATH="${DATA_PATH:-/path/to/training_data.h5}"
EPOCHS="${EPOCHS:-300}"
BATCH_SIZE="${BATCH_SIZE:-64}"
LR="${LR:-1e-3}"
LATENT_DIM="${LATENT_DIM:-24}"
CURRICULUM="${CURRICULUM:-1}"       # 1=enable, 0=disable
WINDOW_MIN="${WINDOW_MIN:-8}"
WINDOW_MAX="${WINDOW_MAX:-48}"
SAVE_EVERY="${SAVE_EVERY:-10}"

# ---- Paths ----
PROJECT_DIR="${PROJECT_DIR:-$(dirname $(dirname $(realpath $0)))}"
RUN_DIR="${RUN_DIR:-${PROJECT_DIR}/runs/slurm_${SLURM_JOB_ID}}"

# ---- Environment setup ----
echo "============================================"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Job Name:     ${SLURM_JOB_NAME}"
echo "Nodes:        ${SLURM_JOB_NUM_NODES}"
echo "GPUs/node:    ${SLURM_GPUS_ON_NODE:-0}"
echo "CPUs/task:    ${SLURM_CPUS_PER_TASK}"
echo "Partition:    ${SLURM_JOB_PARTITION}"
echo "Node list:    ${SLURM_JOB_NODELIST}"
echo "Project dir:  ${PROJECT_DIR}"
echo "Run dir:      ${RUN_DIR}"
echo "Data:         ${DATA_PATH}"
echo "============================================"

# Create output directories
mkdir -p "${RUN_DIR}/checkpoints" "${RUN_DIR}/logs" "logs"

# Load modules (adjust for your HPC cluster)
# Uncomment and modify the lines relevant to your system:
#
# --- Generic module system ---
# module purge
# module load python/3.10
# module load cuda/12.1
# module load gcc/11.3
#
# --- Conda ---
# source activate cardiorenal
#
# --- Singularity/Apptainer ---
# CONTAINER="/path/to/pytorch_latest.sif"

# Verify Python and PyTorch
echo ""
echo "Environment:"
python3 --version
python3 -c "import torch; print(f'PyTorch {torch.__version__}, CUDA available: {torch.cuda.is_available()}')"
python3 -c "import torch; print(f'CUDA devices: {torch.cuda.device_count()}') if torch.cuda.is_available() else print('CPU mode')"
echo ""

# ---- DDP Configuration ----
# Master address: first node in the allocation
MASTER_ADDR=$(scontrol show hostnames "${SLURM_JOB_NODELIST}" | head -n1)
MASTER_PORT=$((29500 + SLURM_JOB_ID % 1000))
export MASTER_ADDR MASTER_PORT

NGPUS_PER_NODE="${SLURM_GPUS_ON_NODE:-1}"
NNODES="${SLURM_JOB_NUM_NODES}"
WORLD_SIZE=$((NNODES * NGPUS_PER_NODE))

echo "DDP: MASTER=${MASTER_ADDR}:${MASTER_PORT}, WORLD_SIZE=${WORLD_SIZE}"

# ---- Build curriculum flag ----
CURRICULUM_FLAG=""
if [ "${CURRICULUM}" = "1" ]; then
    CURRICULUM_FLAG="--curriculum --window_min ${WINDOW_MIN} --window_max ${WINDOW_MAX}"
fi

# ---- Check for existing checkpoint (auto-resume) ----
RESUME_FLAG=""
if [ -d "${RUN_DIR}/checkpoints" ]; then
    CKPT_COUNT=$(find "${RUN_DIR}/checkpoints" -name "*.pt" 2>/dev/null | wc -l)
    if [ "${CKPT_COUNT}" -gt 0 ]; then
        echo "Found ${CKPT_COUNT} existing checkpoint(s) — resuming"
        RESUME_FLAG="--resume auto"
    fi
fi

# ---- Launch training ----
echo ""
echo "Launching training..."
echo "============================================"

cd "${PROJECT_DIR}/.."

if [ "${WORLD_SIZE}" -gt 1 ]; then
    # Multi-GPU: use torchrun for elastic launch
    srun torchrun \
        --nnodes="${NNODES}" \
        --nproc_per_node="${NGPUS_PER_NODE}" \
        --rdzv_id="${SLURM_JOB_ID}" \
        --rdzv_backend=c10d \
        --rdzv_endpoint="${MASTER_ADDR}:${MASTER_PORT}" \
        -m neural_surrogate.train_hpc \
        --data "${DATA_PATH}" \
        --epochs "${EPOCHS}" \
        --batch_size "${BATCH_SIZE}" \
        --lr "${LR}" \
        --latent_dim "${LATENT_DIM}" \
        --run_dir "${RUN_DIR}" \
        --save_every "${SAVE_EVERY}" \
        --num_workers "${SLURM_CPUS_PER_TASK}" \
        ${CURRICULUM_FLAG} \
        ${RESUME_FLAG}
else
    # Single GPU or CPU
    python3 -m neural_surrogate.train_hpc \
        --data "${DATA_PATH}" \
        --epochs "${EPOCHS}" \
        --batch_size "${BATCH_SIZE}" \
        --lr "${LR}" \
        --latent_dim "${LATENT_DIM}" \
        --run_dir "${RUN_DIR}" \
        --save_every "${SAVE_EVERY}" \
        --num_workers "${SLURM_CPUS_PER_TASK}" \
        ${CURRICULUM_FLAG} \
        ${RESUME_FLAG}
fi

EXIT_CODE=$?

echo ""
echo "============================================"
echo "Training exited with code ${EXIT_CODE}"
echo "Run directory: ${RUN_DIR}"
echo "============================================"

exit ${EXIT_CODE}
