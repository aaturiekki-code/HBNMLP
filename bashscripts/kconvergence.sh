#!/usr/bin/env bash
#SBATCH --job-name=qe_kpt_conv
#SBATCH --output=qe_kpt_conv_%j.out
#SBATCH --error=qe_kpt_conv_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G


set -euo pipefail

# kpt_convergence.sbatch
# Generate QE inputs varying k-point mesh (k x k x 1), run pw.x with srun,
# and test energy convergence in Ry/atom.

# -------- User-configurable defaults (override via CLI args or env) --------
# Pass a whitespace- or comma-separated list of k-values as the first arg, e.g.:
#   sbatch scripts/kpt_convergence.sbatch "6 8 10 12" 1e-5
# or  sbatch scripts/kpt_convergence.sbatch 6,8,10,12 1e-5
K_LIST_ARG=${1:-"6 8 10 12"}
THRESHOLD_RY_PER_ATOM=${2:-1e-5}    # Ry/atom convergence threshold

PW_CMD=${PW_CMD:-pw.x}              # pw.x in PATH, or set to full path
INPUT_DIR=${INPUT_DIR:-"../QEinputs"}
TEMPLATE=${TEMPLATE:-"aa.in"}
RUN_ROOT=${RUN_ROOT:-"../runs_kpt"}
WORKDIR_NAME=${WORKDIR_NAME:-"kpt_sweep"}
VERBOSE=${VERBOSE:-0}

module restore qemod

# Parallel setup
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
NTASKS=${SLURM_NTASKS:-1}

# Scratch base for QE outdir (use local node storage if available)
SCRATCH_BASE="${SLURM_TMPDIR:-$PWD}"

# Constants
RY2EV=13.605693009

echo "=== QE k-point convergence job on Triton ==="
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

# Normalize K list: replace commas with spaces, squeeze spaces
K_LIST=$(echo "${K_LIST_ARG}" | tr ',' ' ' | xargs)
if [[ -z "${K_LIST}" ]]; then
  echo "Error: empty k-point list." >&2
  exit 1
fi
echo "k-point meshes to test (k x k x 1): ${K_LIST}"

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

# Helper: generate a run-specific input by replacing:
# - prefix
# - outdir (to a unique per-run scratch)
# - K_POINTS automatic mesh line to "k k 1 0 0 0"
gen_input() {
  local k="$1" runid="$2"
  local infile="${WORKDIR}/aa_k${k}.in"
  local outdir_path="${SCRATCH_BASE}/tmp_${runid}"

  mkdir -p "${outdir_path}"

  awk -v k="${k}" -v outdir="${outdir_path}" '
    function repl_qe_str(val) {
      # surround with quotes
      return "\x27" val "\x27";
    }
    BEGIN { replace_next_k_line=0 }
    # Replace prefix and outdir in &CONTROL (match anywhere to be robust)
    /^[[:space:]]*prefix[[:space:]]*=/ {
      sub(/=[[:space:]]*.*/, "= " repl_qe_str("aa_k" k) ",")
    }
    /^[[:space:]]*outdir[[:space:]]*=/ {
      sub(/=[[:space:]]*.*/, "= " repl_qe_str(outdir) ",")
    }
    # Handle K_POINTS automatic block: replace the next line with new mesh
    /^K_POINTS[[:space:]]+automatic/ {
      print
      replace_next_k_line=1
      next
    }
    replace_next_k_line==1 {
      # Print our replacement mesh and skip original line
      printf("%d %d %d 0 0 0\n", k, k, 1)
      replace_next_k_line=0
      next
    }
    { print }
  ' "${WORKDIR}/template.in" > "${infile}"

  echo "${infile}"
}

echo "# k-point sweep: meshes=${K_LIST} threshold=${THRESHOLD_RY_PER_ATOM} Ry/atom"
echo "k, kx, ky, kz, total_energy_ry, total_energy_eV, delta_ry_per_atom, delta_meV_per_atom" | tee "${WORKDIR}/results.csv" > /dev/null

prev_energy_ry=""
prev_energy_ev=""
converged_at=""
converged_delta_ry=""

for k in ${K_LIST}; do
  if ! [[ "${k}" =~ ^[0-9]+$ ]]; then
    echo "Skipping invalid k value: ${k}" >&2
    continue
  fi

  runid="k${k}"
  infile=$(gen_input "${k}" "${runid}")
  outfile="${WORKDIR}/aa_k${k}.out"

  if [[ "${VERBOSE}" -eq 1 ]]; then
    echo "Running k-mesh = ${k} x ${k} x 1 ..."
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

  echo "${k}, ${k}, ${k}, 1, ${energy_ry}, ${energy_ev}, ${delta_ry_per_atom}, ${delta_mev_per_atom}" | tee -a "${WORKDIR}/results.csv" > /dev/null

  # Convergence check using Ry/atom (avoid set -e abort by testing within if)
  if [[ "${delta_ry_per_atom}" != "NA" ]]; then
    if awk -v d="${delta_ry_per_atom}" -v t="${THRESHOLD_RY_PER_ATOM}" 'BEGIN { exit !(d <= t) }'; then
      converged_at="${k}"
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
  echo "Converged at k-mesh = ${converged_at} x ${converged_at} x 1 with delta ≈ ${converged_delta_ry} Ry/atom."
else
  echo "Did not reach the specified threshold within the provided k-mesh list."
fi

echo "Results saved to: ${WORKDIR}/results.csv"
echo "All inputs/outputs saved to: ${WORKDIR}"

# Optionally archive outputs
# tar -C "${RUN_ROOT}" -czf "${WORKDIR}.tar.gz" "$(basename "${WORKDIR}")"
