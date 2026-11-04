#!/usr/bin/env bash
#SBATCH --job-name=qe_ecrho_conv
#SBATCH --output=qe_ecrho_conv_%j.out
#SBATCH --error=qe_ecrho_conv_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

# ecutrho_convergence.sbatch
# Sweep ecutrho (charge density cutoff) while keeping ecutwfc fixed,
# run pw.x with srun, and test energy convergence in Ry/atom.

# -------- User-configurable defaults (override via CLI args or env) --------
# Usage examples:
#   sbatch scripts/ecutrho_convergence.sbatch 120 320 20 1e-5
#   ECUTWFC=60 sbatch scripts/ecutrho_convergence.sbatch 180 360 30 5e-6
RHO_START=${1:-160}                 # ecutrho start (Ry)
RHO_END=${2:-320}                   # ecutrho end (Ry)
RHO_STEP=${3:-20}                   # ecutrho step (Ry)
THRESHOLD_RY_PER_ATOM=${4:-1e-5}    # Ry/atom convergence threshold

PW_CMD=${PW_CMD:-pw.x}              # pw.x in PATH, or set to full path
INPUT_DIR=${INPUT_DIR:-"../QEinputs"}
TEMPLATE=${TEMPLATE:-"aa.in"}
RUN_ROOT=${RUN_ROOT:-"../runs_ecrho"}
WORKDIR_NAME=${WORKDIR_NAME:-"ecrho_sweep"}
VERBOSE=${VERBOSE:-0}

# If set, override ecutwfc in the generated inputs
# Example: ECUTWFC=60 sbatch ...
ECUTWFC_OVERRIDE=${ECUTWFC:-}

module restore qemod

# Parallel setup
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
NTASKS=${SLURM_NTASKS:-1}

# Scratch base for QE outdir (use local node storage if available)
SCRATCH_BASE="${SLURM_TMPDIR:-$PWD}"

# Constants
RY2EV=13.605693009

echo "=== QE ecutrho convergence job on Triton ==="
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

# Extract current ecutwfc from template (to report and sanity check)
ECUTWFC_TEMPLATE=$(awk '
  BEGIN { in_sys=0; val="" }
  /^\s*&SYSTEM/ { in_sys=1 }
  in_sys && $0 ~ /\// { in_sys=0 }
  in_sys {
    match($0, /ecutwfc[[:space:]]*=[[:space:]]*([0-9.]+)/, a)
    if (a[1] != "") { val=a[1] }
  }
  END { if (val!="") print val }
' "${INPUT_DIR}/${TEMPLATE}" || true)

if [[ -n "${ECUTWFC_OVERRIDE}" ]]; then
  echo "Overriding ecutwfc to ${ECUTWFC_OVERRIDE} Ry for all runs."
  ECUTWFC_USED="${ECUTWFC_OVERRIDE}"
else
  ECUTWFC_USED="${ECUTWFC_TEMPLATE:-""}"
  if [[ -z "${ECUTWFC_USED}" ]]; then
    echo "Warning: could not parse ecutwfc from template; continuing without override."
  else
    echo "Using template ecutwfc = ${ECUTWFC_USED} Ry."
  fi
fi

# Copy template for provenance
cp "${INPUT_DIR}/${TEMPLATE}" "${WORKDIR}/template.in"

# Helper: generate a run-specific input by replacing:
# - ecutrho (swept)
# - optionally ecutwfc (if ECUTWFC_OVERRIDE is set)
# - prefix and outdir
gen_input() {
  local ecrho="$1" runid="$2"
  local infile="${WORKDIR}/aa_ecrho${ecrho}.in"
  local outdir_path="${SCRATCH_BASE}/tmp_${runid}"

  mkdir -p "${outdir_path}"

  # Build sed script dynamically to only replace needed fields
  local sedscript=()
  sedscript+=(-e "s/^([[:space:]]*ecutrho[[:space:]]*=[[:space:]]*).*/\\1${ecrho},/")
  if [[ -n "${ECUTWFC_OVERRIDE}" ]]; then
    sedscript+=(-e "s/^([[:space:]]*ecutwfc[[:space:]]*=[[:space:]]*).*/\\1${ECUTWFC_OVERRIDE},/")
  fi
  sedscript+=(-e "s|^([[:space:]]*outdir[[:space:]]*=[[:space:]]*).*|\\1'${outdir_path//\//\\/}',|")
  sedscript+=(-e "s/^([[:space:]]*prefix[[:space:]]*=[[:space:]]*).*/\\1'aa_ecrho${ecrho}',/")

  sed -E "${sedscript[@]}" "${WORKDIR}/template.in" > "${infile}"

  echo "${infile}"
}

# Build ecutrho sequence (allow integer or float step via seq)
build_rho_list() {
  local start="$1" end="$2" step="$3"
  # seq is commonly available; format to avoid FP artifacts
  seq -f "%.6f" "${start}" "${step}" "${end}" | sed 's/[.]000000$//'  # trim .000000 for integers
}

RHO_LIST=$(build_rho_list "${RHO_START}" "${RHO_END}" "${RHO_STEP}")
echo "# ecutrho sweep: start=${RHO_START} end=${RHO_END} step=${RHO_STEP} threshold=${THRESHOLD_RY_PER_ATOM} Ry/atom"
echo "ecutrho_ry, ecutwfc_ry, total_energy_ry, total_energy_eV, delta_ry_per_atom, delta_meV_per_atom" | tee "${WORKDIR}/results.csv" > /dev/null

prev_energy_ry=""
prev_energy_ev=""
converged_at=""
converged_delta_ry=""

for ecrho in ${RHO_LIST}; do
  # Optional sanity: ensure ecutrho >= ecutwfc if we know ecutwfc
  if [[ -n "${ECUTWFC_USED}" ]]; then
    awk -v ecr="${ecrho}" -v ecw="${ECUTWFC_USED}" 'BEGIN { if (ecr < ecw) exit 1 }' || {
      echo "Warning: ecutrho (${ecrho} Ry) < ecutwfc (${ECUTWFC_USED} Ry); QE may reset/complain."
    }
  fi

  runid="ecrho${ecrho}"
  infile=$(gen_input "${ecrho}" "${runid}")
  outfile="${WORKDIR}/aa_ecrho${ecrho}.out"

  if [[ "${VERBOSE}" -eq 1 ]]; then
    echo "Running ecutrho=${ecrho} Ry ..."
  fi

  # Run QE with srun (MPI aware)
  srun --ntasks="${NTASKS}" --cpus-per-task="${OMP_NUM_THREADS}" "${PW_CMD}" -in "${infile}" > "${outfile}" 2>&1

  # Parse total energy from output: line like "!    total energy =   -XXX.XXXXXXXX Ry"
  energy_ry=$(grep -m1 '^!.*total energy' "${outfile}" | awk '{print $(NF-1)}' || true)

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

  echo "${ecrho}, ${ECUTWFC_USED:-NA}, ${energy_ry}, ${energy_ev}, ${delta_ry_per_atom}, ${delta_mev_per_atom}" | tee -a "${WORKDIR}/results.csv" > /dev/null

  # Convergence check using Ry/atom (avoid set -e abort by testing within if)
  if [[ "${delta_ry_per_atom}" != "NA" ]]; then
    if awk -v d="${delta_ry_per_atom}" -v t="${THRESHOLD_RY_PER_ATOM}" 'BEGIN { exit !(d <= t) }'; then
      converged_at="${ecrho}"
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
  echo "Converged at ecutrho = ${converged_at} Ry with delta ≈ ${converged_delta_ry} Ry/atom."
else
  echo "Did not reach the specified threshold within the range."
fi

echo "Results saved to: ${WORKDIR}/results.csv"
echo "All inputs/outputs saved to: ${WORKDIR}"

# Optionally archive outputs
# tar -C "${RUN_ROOT}" -czf "${WORKDIR}.tar.gz" "$(basename "${WORKDIR}")"
