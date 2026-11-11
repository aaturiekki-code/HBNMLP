#!/bin/bash
#SBATCH --job-name=restart
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

mamba run -n coeff-env python3 pyscripts/restart_unfinished.py \
                                --dir QEinputs/Twist_m0_2_sweep \
                                --workflow scf-only --dev-xy 0 --dev-z 0 \
                                --regenerate --force \
                                --run
