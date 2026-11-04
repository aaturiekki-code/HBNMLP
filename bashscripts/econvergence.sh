#!/usr/bin/env bash
#SBATCH --job-name=qe_ecut_conv
#SBATCH --output=qe_ecut_conv_%j.out
#SBATCH --error=qe_ecut_conv_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G


set -euo pipefail

# ecut_convergence.sbatch
# Generate QE inputs varying ecutwfc, run pw.x with srun, and test energy convergence in Ry/atom.

# -------- User-configurable defaults (override via CLI args or env) --------
ECUT_START=${1:-30}                 # Ry
ECUT_END=${2:-100}                  # Ry
ECUT_STEP=${3:-5}                   # Ry
THRESHOLD_RY_PER_ATOM=${4:-1e-5}    # Ry/atom convergence threshold

ECRATIO=${ECRATIO:-4}               # ecutrho = ECRATIO * ecutwfc
PW_CMD=${PW_CMD:-pw.x}              # pw.x in PATH, or set to full path
INPUT_DIR=${INPUT_DIR:-"../QEinputs"}
TEMPLATE=${TEMPLATE:-"aa.in"}
RUN_ROOT=${RUN_ROOT:-"../runs_ecutwfc"}
WORKDIR_NAME=${WORKDIR_NAME:-"ecut_sweep"}
VERBOSE=${VERBOSE:-0}

module restore qemod

# Parallel setup
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
NTASKS=${SLURM_NTASKS:-1}

# Scratch base for QE outdir (use local node storage if available)
SCRATCH_BASE="${SLURM_TMPDIR:-$PWD}"

# Constants
RY2EV=13.605693009

echo "=== QE ecutwfc convergence job on Triton ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "NTASKS: ${NTASKS}, OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "Using PW_CMD: ${PW_CMD}"

# Checks
if [[ ! -f "${INPUT_DIR}/${TEMPLATE}" ]]; then
  echo "Error: template input not found at ${INPUT_DIR}/${TEMPLATE}" >&2
  exit 1
fi

if ! command -v "${PW_CMD}" >/dev/null 2>&1; then
  echo "Error: pw.x command not found. Load a QE module or set PW_CMD." >&2
  exit 1
fi

mkdir -p "${RUN_ROOT}"
WORKDIR="${RUN_ROOT}/${WORKDIR_NAME}_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-local}"
mkdir -p "${WORKDIR}"
echo "Working directory: ${WORKDIR}"

# Extract nat (number of atoms) from template for per-atom conversion
NAT=$(awk '
  BEGIN { in_sys=0; nat="" }
  /^\s*&SYSTEM/ { in_sys=1 }
  in_sys && $0 ~ /\// { in_sys=0 }
  in_sys && $0 ~ /nat/ {
    match($0, /nat[[:space:]]*=[[:space:]]*([0-9]+)/, a)
    if (a[1] != "") { print a[1]; exit }
  }
' "${INPUT_DIR}/${TEMPLATE}")

if [[ -z "${NAT}" ]]; then
  echo "Warning: could not parse nat from template. Defaulting to nat=1 for convergence calc."
  NAT=1
fi

echo "Parsed nat = ${NAT}"

# Copy template for provenance
cp "${INPUT_DIR}/${TEMPLATE}" "${WORKDIR}/template.in"

# Helper: generate a run-specific input by replacing ecutwfc/ecutrho/outdir/prefix
gen_input() {
  local ecutwfc="$1" ecutrho="$2" runid="$3"
  local infile="${WORKDIR}/aa_ec${ecutwfc}.in"
  local outdir_path="/tmp/"

  mkdir -p "${outdir_path}"

  sed -E \
    -e "s/^([[:space:]]*ecutwfc[[:space:]]*=[[:space:]]*).*/\1${ecutwfc},/" \
    -e "s/^([[:space:]]*ecutrho[[:space:]]*=[[:space:]]*).*/\1${ecutrho},/" \
    -e "s|^([[:space:]]*outdir[[:space:]]*=[[:space:]]*).*|\1'${outdir_path//\//\\/}',|" \
    -e "s/^([[:space:]]*prefix[[:space:]]*=[[:space:]]*).*/\1'aa_ec${ecutwfc}',/" \
    "${WORKDIR}/template.in" > "${infile}"

  echo "${infile}"
}

echo "# ecutwfc_sweep: start=${ECUT_START} end=${ECUT_END} step=${ECUT_STEP} ratio=${ECRATIO} threshold=${THRESHOLD_RY_PER_ATOM} Ry/atom"
echo "ecutwfc_ry, ecutrho_ry, total_energy_ry, total_energy_eV, delta_ry_per_atom, delta_meV_per_atom" | tee "${WORKDIR}/results.csv" > /dev/null

prev_energy_ry=""
prev_energy_ev=""
converged_at=""
converged_delta_ry=""

for (( ecut=${ECUT_START}; ecut<=ECUT_END; ecut+=ECUT_STEP )); do
  ecutrho=$(( ecut * ECRATIO ))
  runid="ec${ecut}"

  infile=$(gen_input "${ecut}" "${ecutrho}" "${runid}")
  outfile="${WORKDIR}/aa_ec${ecut}.out"

  if [[ "${VERBOSE}" -eq 1 ]]; then
    echo "Running ecutwfc=${ecut} Ry, ecutrho=${ecutrho} Ry ..."
  fi

  # Run QE with srun (MPI aware)
  srun --ntasks="${NTASKS}" --cpus-per-task="${OMP_NUM_THREADS}" "${PW_CMD}" -in "${infile}" > "${outfile}" 2>&1

  # Parse total energy from output: line like "!    total energy =   -XXX.XXXXXXXX Ry"
  energy_ry=$(grep -m1 "total energy" "${outfile}" | awk '{print $(NF-1)}' || true)

  if [[ -z "${energy_ry}" ]]; then
    echo "Error: could not parse total energy from ${outfile}"
    echo "See ${outfile} for details."
    exit 1
  fi

  energy_ev=$(awk -v e="${energy_ry}" -v f="${RY2EV}" 'BEGIN { printf("%.10f", e*f) }')

  delta_ry_per_atom="NA"
  delta_mev_per_atom="NA"
  if [[ -n "${prev_energy_ry}" ]]; then
    delta_ry_per_atom=$(awk -v eprev="${prev_energy_ry}" -v ecur="${energy_ry}" -v nat="${NAT}" '
      BEGIN {
        d = (ecur - eprev) / nat;
        if (d < 0) d = -d;
        printf("%.8e", d);
      }')
    delta_mev_per_atom=$(awk -v dry="${delta_ry_per_atom}" -v f="${RY2EV}" '
      BEGIN { printf("%.6f", dry * f * 1000.0) }')
  fi

  echo "${ecut}, ${ecutrho}, ${energy_ry}, ${energy_ev}, ${delta_ry_per_atom}, ${delta_mev_per_atom}" | tee -a "${WORKDIR}/results.csv" > /dev/null

  # Convergence check using Ry/atom (avoid set -e abort by testing within if)
  if [[ "${delta_ry_per_atom}" != "NA" ]]; then
    if awk -v d="${delta_ry_per_atom}" -v t="${THRESHOLD_RY_PER_ATOM}" 'BEGIN { exit !(d <= t) }'; then
      converged_at="${ecut}"
      converged_delta_ry="${delta_ry_per_atom}"
      # Stop early on first convergence
      break
    fi
  fi

  prev_energy_ry="${energy_ry}"
  prev_energy_ev="${energy_ev}"
done

echo
if [[ -n "${converged_at}" ]]; then
  echo "Converged at ecutwfc = ${converged_at} Ry with delta ≈ ${converged_delta_ry} Ry/atom."
else
  echo "Did not reach the specified threshold within the range."
fi

echo "Results saved to: ${WORKDIR}/results.csv"
echo "All inputs/outputs saved to: ${WORKDIR}"

# Optionally archive outputs
# tar -C "${RUN_ROOT}" -czf "${WORKDIR}.tar.gz" "$(basename "${WORKDIR}")"
