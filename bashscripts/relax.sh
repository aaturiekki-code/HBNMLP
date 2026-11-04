#!/bin/bash
#SBATCH --job-name=relax
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

NAME="relax"

module restore qemod

cd QEinputs
srun pw.x -np 16 -inp aa.relax.in > ../QEoutputs/aa.relax.out
