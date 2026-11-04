#!/bin/bash
#SBATCH --job-name=hBN_orchestrator
#SBATCH --time=00:15:00             # Runs quickly to generate directories
#SBATCH --mem=4G                    # Low footprint
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1           # Only coordinates other jobs
#SBATCH --output=slurm-%x-%j.out

set -e

# Default parameter settings
N=1
M0=""

# Help menu description
usage() {
  echo "Usage: sbatch bashscripts/$(basename "$0") [options]"
  echo "Options:"
  echo "  -n INT    Generate untwisted bilayer sweeps with NxN supercell expansion (default: 1)"
  echo "  -m INT    Generate twisted moiré bilayer sweep with index m0 (switches mode to 'twist')"
  echo "  -h        Show this help menu"
  exit 1
}

# Parse command line flags passed through sbatch
while getopts "n:m:h" opt; do
  case "$opt" in
    n) N=$OPTARG ;;
    m) M0=$OPTARG ;;
    h|*) usage ;;
  esac
done

# Load modules required for running Python/Mamba configurations
module purge
module load mamba
source activate coeff-env

# Set absolute parent project directory paths
PARENT_DIR="/scratch/work/riekkia3/Special_Assignment"
SCRIPT_DIR="$PARENT_DIR/bashscripts"
PYTHON_DIR="$PARENT_DIR/pyscripts"

# Path to the executable batch creator
BATCH_CREATE="$PYTHON_DIR/batch_create.py"

# Ensure the output sibling directory exists
mkdir -p "$PARENT_DIR/QEinputs"

# -------------------------------------------------------------------------
# Dynamic SLURM Resource Allocation Logic
# -------------------------------------------------------------------------
NTASKS=16

if [ -n "$M0" ]; then
  # Twisted moiré (e.g., m0=1 containing 28 atoms)
  TIME_LIMIT="24:00:00"
  MEMORY="50G"
  echo "Setting SLURM limits for twisted moiré sweep: --slurm-time=$TIME_LIMIT --slurm-mem=$MEMORY"
else
  # Untwisted bilayers based on N
  if [ "$N" -eq 1 ]; then
    TIME_LIMIT="00:30:00"
    MEMORY="10G"
  elif [ "$N" -eq 2 ]; then
    TIME_LIMIT="05:00:00"
    MEMORY="20G"
  elif [ "$N" -eq 3 ]; then
    TIME_LIMIT="24:00:00"
    MEMORY="50G"
  else
    # Fallback default for any larger grids
    TIME_LIMIT="24:00:00"
    MEMORY="64G"
  fi
  echo "Setting SLURM limits for ${N}x${N} untwisted sweep: --slurm-time=$TIME_LIMIT --slurm-mem=$MEMORY"
fi

# -------------------------------------------------------------------------
# Execution Mode: Twisted Moiré Bilayer (if -m was specified)
# -------------------------------------------------------------------------
if [ -n "$M0" ]; then
  echo "========================================================================="
  echo " Running Generation for Twisted Moiré Bilayer (m0 = $M0)"
  echo " Target Directory:  $PARENT_DIR/QEinputs/Twist_m0_${M0}_sweep"
  echo "========================================================================="

  echo -e "\n---> Generating Twisted Stacking Sweep (m0=$M0, 1x1 supercell)..."
  python3 "$BATCH_CREATE" \
    --parent-dir "$PARENT_DIR/QEinputs/Twist_m0_${M0}_sweep" \
    --distance-range 3.0 3.8 9 \
    --scale-range 1.41 1.47 5 \
    --mode "twist" \
    --m0 "$M0" \
    --vacuum 30 \
    --n1 1 --n2 1 \
    --dev-xy 0.04 --dev-z 0.04 --seed 555 \
    --slurm-time "$TIME_LIMIT" \
    --slurm-mem "$MEMORY" \
    --slurm-cores "$NTASKS"

  echo "Queuing Twisted Sweep jobs in Slurm..."
  cd "$PARENT_DIR/QEinputs/Twist_m0_${M0}_sweep"
  sbatch submit_all.sh
  cd "$PARENT_DIR"

# -------------------------------------------------------------------------
# Execution Mode: Untwisted Bilayer sweeps (if -n was specified or default)
# -------------------------------------------------------------------------
else
  echo "========================================================================="
  echo " Running Generation for Untwisted Bilayers (Supercell size: ${N}x${N})"
  echo " Target Directory:  $PARENT_DIR/QEinputs"
  echo "========================================================================="

  # 1. AA' Stacking Sweep
  echo -e "\n---> [1/4] Generating Scaled AA' Stacking Sweep (${N}x${N})..."
  python3 "$BATCH_CREATE" \
    --parent-dir "$PARENT_DIR/QEinputs/AA_prime_sweep_${N}x${N}" \
    --distance-range 3.0 3.8 9 \
    --scale-range 1.41 1.47 5 \
    --mode "bilayer" \
    --bilayer-stack "AA" \
    --mirror 1 \
    --vacuum 30 \
    --n1 "$N" --n2 "$N" \
    --dev-xy 0.04 --dev-z 0.04 --seed 111 \
    --slurm-time "$TIME_LIMIT" \
    --slurm-mem "$MEMORY" \
    --slurm-cores "$NTASKS"

  echo "Queuing AA' ${N}x${N} Sweep jobs in Slurm..."
  cd "$PARENT_DIR/QEinputs/AA_prime_sweep_${N}x${N}"
  sbatch submit_all.sh
  cd "$PARENT_DIR"

  # 2. AB Stacking Sweep
  echo -e "\n---> [2/4] Generating Scaled AB Stacking Sweep (${N}x${N})..."
  python3 "$BATCH_CREATE" \
    --parent-dir "$PARENT_DIR/QEinputs/AB_sweep_${N}x${N}" \
    --distance-range 3.0 3.8 9 \
    --scale-range 1.41 1.47 5 \
    --mode "bilayer" \
    --bilayer-stack "AB" \
    --mirror 0 \
    --vacuum 30 \
    --n1 "$N" --n2 "$N" \
    --dev-xy 0.04 --dev-z 0.04 --seed 222 \
    --slurm-time "$TIME_LIMIT" \
    --slurm-mem "$MEMORY" \
    --slurm-cores "$NTASKS"

  echo "Queuing AB ${N}x${N} Sweep jobs in Slurm..."
  cd "$PARENT_DIR/QEinputs/AB_sweep_${N}x${N}"
  sbatch submit_all.sh
  cd "$PARENT_DIR"

  # 3. AA Stacking Sweep
  echo -e "\n---> [3/4] Generating Scaled AA Stacking Sweep (${N}x${N})..."
  python3 "$BATCH_CREATE" \
    --parent-dir "$PARENT_DIR/QEinputs/AA_sweep_${N}x${N}" \
    --distance-range 3.0 3.8 9 \
    --scale-range 1.41 1.47 5 \
    --mode "bilayer" \
    --bilayer-stack "AA" \
    --mirror 0 \
    --vacuum 30 \
    --n1 "$N" --n2 "$N" \
    --dev-xy 0.04 --dev-z 0.04 --seed 333 \
    --slurm-time "$TIME_LIMIT" \
    --slurm-mem "$MEMORY" \
    --slurm-cores "$NTASKS"

  echo "Queuing AA ${N}x${N} Sweep jobs in Slurm..."
  cd "$PARENT_DIR/QEinputs/AA_sweep_${N}x${N}"
  sbatch submit_all.sh
  cd "$PARENT_DIR"

  # 4. Stacking-Boundary Sliding Sweep
  echo -e "\n---> [4/4] Generating Sliding Path Sweep (${N}x${N})..."
  python3 "$BATCH_CREATE" \
    --parent-dir "$PARENT_DIR/QEinputs/Sliding_sweep_${N}x${N}" \
    --slide-x-range 0.0 0.5 6 \
    --slide-y-range 0.0 0.5 6 \
    --distance-range 3.0 3.8 3 \
    --mode "bilayer" \
    --bilayer-stack "AA" \
    --mirror 1 \
    --vacuum 30 \
    --n1 "$N" --n2 "$N" \
    --dev-xy 0.04 --dev-z 0.04 --seed 444 \
    --slurm-time "$TIME_LIMIT" \
    --slurm-mem "$MEMORY" \
    --slurm-cores "$NTASKS"

  echo "Queuing Sliding ${N}x${N} Sweep jobs in Slurm..."
  cd "$PARENT_DIR/QEinputs/Sliding_sweep_${N}x${N}"
  sbatch submit_all.sh
  cd "$PARENT_DIR"
fi

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
echo -e "\n========================================================================="
echo " Generation and Submission Phase Complete!"
echo " All corresponding jobs have been queued in Slurm."
echo "========================================================================="
