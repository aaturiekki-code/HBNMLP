#!/bin/bash
#SBATCH --job-name=hBN_allegro_cpu
#SBATCH --time=24:00:00             
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16          # Request 16 CPU cores for parallel math
#SBATCH --output=slurm-%x-%j.out

set -e

# 1. Load Triton software environment modules
module purge
module load mamba
module load triton/2025.1-gcc
module load pair_nequip_allegro/19-9-2025
source activate coeff-env

# 2. Configure PyTorch to utilize the requested CPU cores efficiently
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=== CPU Environment Information ==="
echo "Running on $SLURM_CPUS_PER_TASK CPU cores."
python3 -c "import torch; print('PyTorch Threads:', torch.get_num_threads())"
echo "==================================="

# 3. Verify the compiled dataset path (aligned with YAML config)
TARGET_XYZ="compiled_dataset/dataset.xyz"
if [ ! -f "$TARGET_XYZ" ]; then
    echo "Error: Could not find '$TARGET_XYZ'"
    echo "Please ensure you have run create_features_and_labels.py first!"
    exit 1
fi

# 4. Run training using the explicit --config flag
echo "Starting CPU-only Allegro ML training..."
nequip-train --config-name allegro_config_cpu.yaml

echo "Allegro training completed successfully!"
