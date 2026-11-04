#!/bin/bash
#SBATCH --job-name=hBN_allegro_gpu
#SBATCH --time=04:00:00             # GPU training is fast, 4 hours is plenty
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1                # Request 1 physical GPU dynamically on Triton
#SBATCH --cpus-per-task=4           # CPU cores for fast parallel data loading
#SBATCH --output=slurm-%x-%j.out

set -e

# Load Triton software environment modules
module purge
module load mamba                   # Triton-supported Mamba package manager

# Activate your custom conda environment containing PyTorch + CUDA + Allegro
source activate coeff-env

# Verify CUDA visibility inside Slurm task allocation
echo "=== CUDA Device Information ==="
nvidia-smi
python3 -c "import torch; print('CUDA Available:', torch.cuda.is_available())"
echo "==============================="

# Convert raw numpy labels to XYZ if not done already
if [ ! -f "dataset.xyz" ]; then
    echo "Converting raw arrays to Extended XYZ format..."
    python3 convert_to_xyz.py
fi

# Launch training run using the GPU config
echo "Starting GPU-accelerated Allegro ML training..."
python3 -m nequip.scripts.train --config-name allegro_config_gpu.yaml

echo "Allegro training completed successfully!"
