#!/bin/bash
#SBATCH --job-name=hBN_compile_pt2
#SBATCH --time=00:15:00             # AOTInductor compilation takes a few minutes
#SBATCH --mem=12G                   # AOTInductor compilation can be memory intensive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8           # Request 8 cores for parallel C++ code generation
#SBATCH --output=slurm-%x-%j.out

set -e

# 1. Load exact Triton software environment modules
module purge
module load mamba
module load triton/2025.1-gcc
module load pair_nequip_allegro/19-9-2025
source activate coeff-env

# 2. Configure PyTorch compiler threads for compilation speed
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 3. Dynamically locate the best.ckpt file in the latest training run
# PyTorch Lightning places checkpoints in: lightning_logs/version_X/checkpoints/best.ckpt
LATEST_VERSION_DIR=$(ls -td lightning_logs/version_*/ 2>/dev/null | head -1)

if [ -n "$LATEST_VERSION_DIR" ]; then
    # Search recursively in the latest directory for any .ckpt files (typically best.ckpt)
    CHECKPOINT=$(find "$LATEST_VERSION_DIR" -name "best.ckpt" | head -1)
fi

# Fallback: check if best.ckpt is directly in the working directory
if [ -z "$CHECKPOINT" ] || [ ! -f "$CHECKPOINT" ]; then
    if [ -f "best.ckpt" ]; then
        CHECKPOINT="best.ckpt"
    else
        echo "Error: Could not find 'best.ckpt' in 'lightning_logs/' or the current folder."
        echo "Please make sure your training has started and saved a checkpoint."
        exit 1
    fi
fi

echo "========================================================================="
echo " Deploying Allegro Model with AOTInductor Compilation"
echo " Source Checkpoint: $CHECKPOINT"
echo " Target Device:    CPU"
echo "========================================================================="

OUTPUT_MODEL="hBN_model_cpu.nequip.pt2"

# 4. Execute the official nequip-compile command using AOTInductor target pair_allegro
nequip-compile \
    "$CHECKPOINT" \
    "$OUTPUT_MODEL" \
    --device cpu \
    --mode aotinductor \
    --target pair_allegro

echo ""
echo "========================================================================="
echo " SUCCESS! Highly parallelized AOTInductor model created."
echo " Deployed Binary File: $OUTPUT_MODEL"
echo ""
echo "You can now load this optimized model inside your LAMMPS input script:"
echo "  pair_style allegro"
echo "  pair_coeff * * $OUTPUT_MODEL B N"
echo "========================================================================="
