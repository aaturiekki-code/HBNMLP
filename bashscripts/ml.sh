#!/bin/bash
#SBATCH --job-name=model
#SBATCH --time=02:00:00
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

mamba run -n coeff-env python3 model.py