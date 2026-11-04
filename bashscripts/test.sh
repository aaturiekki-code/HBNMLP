#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=00:00:10
#SBATCH --mem=1G
#SBATCH --constraint=milan
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

export WORKDIR="${LOCAL_SCRATCH:-/tmp}/qe_$SLURM_JOB_ID"
mkdir -p "$WORKDIR/out"