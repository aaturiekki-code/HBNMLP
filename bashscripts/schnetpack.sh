#!/bin/bash
#SBATCH --job-name=shnetpack
#SBATCH --time=00:01:00
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-%x-%j.out

source activate schnetpack

srun spktrain experiment=qm9_atomwise