#!/bin/bash
#SBATCH --job-name=hBN_minimization
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16          # Multi-threading parallel cores
#SBATCH --output=slurm-%x-%j.out

set -e

# Load Triton modules and activate environment
module purge
module load mamba
module load triton/2025.1-gcc
module load pair_nequip_allegro/19-9-2025
source activate coeff-env

# Set OpenMP threads to match CPUs requested
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=== [1/2] Running LAMMPS minimization ==="
srun lmp -sf omp -pk omp $SLURM_CPUS_PER_TASK -in fixedcell.script > fixedcell.log

#echo "=== [2/2] Generating Box-Enabled XYZ File ==="
python3 ../pyscripts/data_to_poscar.py

echo "All steps completed successfully!"
