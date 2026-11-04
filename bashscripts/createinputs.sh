#!/bin/bash
#SBATCH --job-name=inputs
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --constraint=milan
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

mamba run -n coeff-env python3 pyscripts/batch_create.py --parent-dir "QEinputs/twist0n3" --distance-range 3.2 3.7 4 --scale-range 1.42 1.45 4 --slide-x-range 0 1 4 --slide-y-range 0 1 4 --mode bilayer --bilayer-stack AA --vacuum 30 --dev-xy 0.05 --dev-z 0.05 --n1 3 --n2 3 --slurm-time 24:00:00 --slurm-mem 50G
