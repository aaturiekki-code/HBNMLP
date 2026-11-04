#!/usr/bin/env bash
#
# SBATCH directives (adjust as needed for Triton)
#SBATCH --job-name=qe_diel_sweep
#SBATCH --output=qe_diel_sweep_%A_%a.out
#SBATCH --error=qe_diel_sweep_%A_%a.err
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --export=ALL
# Submit with:
#   TOTAL=$((NA*NS*NL))    # you compute this
#   sbatch --array=0-$((TOTAL-1)) scripts/dielectric_sweep_array.sbatch "A_SCALES" "SCALES" "SLIDES"

set -euo pipefail
module restore qemod

# Inputs: three lists (positional or env), separators: space/comma/semicolon
A_SCALE_LIST_RAW=${1:-${A_SCALE_LIST:-"1.00"}}
SCALE_LIST_RAW=${2:-${SCALE_LIST:-"1.00"}}
SLIDE_LIST_RAW=${3:-${SLIDE_LIST:-"0.00:0.00"}}

# Tools/paths
PW_CMD=${PW_CMD:-pw.x}
PH_CMD=${PH_CMD:-ph.x}
INPUT_DIR=${INPUT_DIR:-"../QEinputs"}
TEMPLATE=${TEMPLATE:-"aa.in"}
PH_TEMPLATE=${PH_TEMPLATE:-"../QEinputs/aa.ph.in"}

RUN_ROOT=${RUN_ROOT:-"../runs_interlayer"}
WORKDIR_NAME=${WORKDIR_NAME:-"dielectric_sweep"}
VERBOSE=${VERBOSE:-0}

# Normalization constants (Å)
C_VAC=${C_VAC:-30.0}
CORE_PAD_A=${CORE_PAD_A:-1.17}

# Per-task feature/label generator (DScribe 2.x)
TASK_PY=${TASK_PY:-"../pyscripts/create_features_and_labels_task.py"}

# Normalize list separators: space/comma/semicolon → spaces; trim
norm_list() { echo "$1" | tr ',;' ' ' | xargs; }
A_SCALE_LIST=$(norm_list "${A_SCALE_LIST_RAW}")
SCALE_LIST=$(norm_list "${SCALE_LIST_RAW}")
SLIDE_LIST=$(norm_list "${SLIDE_LIST_RAW}")

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
NTASKS=${SLURM_NTASKS:-1}
RY2EV=13.605693009

# Checks
[[ -f "${TASK_PY}" ]] || { echo "Error: per-task Python script ${TASK_PY} not found" >&2; exit 1; }
[[ -f "${INPUT_DIR}/${TEMPLATE}" ]] || { echo "Error: template ${INPUT_DIR}/${TEMPLATE} not found" >&2; exit 1; }
command -v "${PW_CMD}" >/dev/null 2>&1 || { echo "Error: pw.x not found" >&2; exit 1; }
command -v "${PH_CMD}" >/dev/null 2>&1 || { echo "Error: ph.x not found" >&2; exit 1; }

# Shared base directory for entire array job
BASEDIR="${RUN_ROOT}/${WORKDIR_NAME}_${SLURM_ARRAY_JOB_ID:-local}"
mkdir -p "${BASEDIR}"

# Build parameter arrays
mapfile -t A_SCALES < <(for a in ${A_SCALE_LIST}; do echo "$a"; done)
mapfile -t ZSCALES  < <(for s in ${SCALE_LIST};   do echo "$s"; done)
mapfile -t SLIDES   < <(for p in ${SLIDE_LIST};   do echo "$p"; done)

NA=${#A_SCALES[@]}
NS=${#ZSCALES[@]}
NL=${#SLIDES[@]}
(( NA>0 && NS>0 && NL>0 )) || { echo "Error: one of the lists is empty (NA=${NA}, NS=${NS}, NL=${NL})" >&2; exit 1; }

# TOTAL compute tasks
TOTAL=$(( NA * NS * NL ))
TASK=${SLURM_ARRAY_TASK_ID:-0}
if (( TASK >= TOTAL )); then
  echo "Info: array task ${TASK} >= TOTAL ${TOTAL}; nothing to do."
  exit 0
fi

# Map TASK → (ia,is,il)
prod_sl=$(( NS * NL ))
ia=$(( TASK / prod_sl ))
rem=$(( TASK % prod_sl ))
is=$(( rem / NL ))
il=$(( rem % NL ))

A_SCALE="${A_SCALES[$ia]}"
SCALE="${ZSCALES[$is]}"
pair="${SLIDES[$il]}"
dx="${pair%%:*}"
dy="${pair##*:}"

# Per-task working directory
WORKDIR="${BASEDIR}/task_${TASK}"
mkdir -p "${WORKDIR}"

# Copy templates for this task
cp "${INPUT_DIR}/${TEMPLATE}" "${WORKDIR}/template.in"
if [[ -f "${PH_TEMPLATE}" ]]; then
  cp "${PH_TEMPLATE}" "${WORKDIR}/template.ph.in"
fi

# Parse nat and CELL_PARAMETERS; compute |c|
read -r NAT AX AY AZ BX BY BZ CX CY CZ C_NORM <<<"$(
  awk '
    BEGIN{in_sys=0; in_cell=0; ncell=0}
    /^\s*&SYSTEM/ {in_sys=1}
    in_sys && /^[[:space:]]*\// {in_sys=0}
    in_sys && /nat[[:space:]]*=/ {match($0,/nat[[:space:]]*=[[:space:]]*([0-9]+)/,a); if(a[1]!="") nat=a[1]}
    /^CELL_PARAMETERS/ {in_cell=1; ncell=0; next}
    in_cell && ncell<3 {
      ncell++;
      if(ncell==1){ax=$1; ay=$2; az=$3}
      else if(ncell==2){bx=$1; by=$2; bz=$3}
      else{cx=$1; cy=$2; cz=$3; in_cell=0}
      next
    }
    END {
      c=(cx==""?0:sqrt(cx*cx+cy*cy+cz*cz));
      printf("%s %s %s %s %s %s %s %s %s %s %.12f\n", nat, ax,ay,az, bx,by,bz, cx,cy,cz, c)
    }' "${WORKDIR}/template.in"
)"
[[ -n "${NAT}" && "${C_NORM}" != "0" ]] || { echo "Error: failed to parse nat or CELL_PARAMETERS" >&2; exit 1; }

# Compute initial interlayer separation (fractional), split around z=0.5
read -r SEP0_FRAC <<<"$(
  awk -v nat="${NAT}" '
    BEGIN{in_atoms=0; n=0; sum_l=0; sum_u=0; cnt_l=0; cnt_u=0}
    /^ATOMIC_POSITIONS[[:space:]]*\(crystal\)/{in_atoms=1; n=0; next}
    in_atoms && n<nat {
      n++; z=$4+0.0;
      if(z<0.5){sum_l+=z; cnt_l++} else {sum_u+=z; cnt_u++}
      if(n==nat) in_atoms=0
    }
    END{
      if(cnt_l==0 || cnt_u==0){print 0}else{print (sum_u/cnt_u)-(sum_l/cnt_l)}
    }' "${WORKDIR}/template.in"
)"

# Derive separation and thickness (Å)
awk -v v="${SCALE}" 'BEGIN{exit !(v+0==v)}' || SCALE="1.0"
sepA=$(awk -v s="${SCALE}" -v s0="${SEP0_FRAC}" -v c="${C_NORM}" 'BEGIN{printf("%.6f", s*s0*c)}')
thickness=$(awk -v sep="${sepA}" -v pad="${CORE_PAD_A}" 'BEGIN{printf("%.6f", sep + 2.0*pad)}')

OUTDIR_TASK="${WORKDIR}/tmp_aS${A_SCALE}_zS${SCALE}_dx${dx}_dy${dy}"
mkdir -p "${OUTDIR_TASK}"

# Generators
gen_pw_input() {
  local infile="${WORKDIR}/aa_aS${A_SCALE}_zS${SCALE}_dx${dx}_dy${dy}.in"
  local outdir="${OUTDIR_TASK}"
  awk -v aS="${A_SCALE}" -v SCALE="${SCALE}" -v dx="${dx}" -v dy="${dy}" -v outdir="${outdir}" -v nat="${NAT}" \
      -v ax="${AX}" -v ay="${AY}" -v az="${AZ}" -v bx="${BX}" -v by="${BY}" -v bz="${BZ}" -v cx="${CX}" -v cy="${CY}" -v cz="${CZ}" '
    BEGIN{in_ctrl=0; in_cell=0; in_atoms=0; n=0; cell_rows=0}
    /^\s*&CONTROL/ {in_ctrl=1}
    in_ctrl && /^[[:space:]]*\// {in_ctrl=0}
    in_ctrl && /^[[:space:]]*prefix[[:space:]]*=/ {sub(/=[[:space:]]*.*/,"= \"" "aa_aS" aS "_zS" SCALE "_dx" dx "_dy" dy "\" ,"); print; next}
    in_ctrl && /^[[:space:]]*outdir[[:space:]]*=/ {sub(/=[[:space:]]*.*/,"= \"" outdir "\" ,"); print; next}
    in_ctrl && /^[[:space:]]*pseudo_dir[[:space:]]*=/ {sub(/=[[:space:]]*.*/,"= '\''../pseudos/'\'',"); print; next}

    /^CELL_PARAMETERS/ {
      print; in_cell=1; cell_rows=0;
      printf("% .12f % .12f % .12f\n", ax*aS, ay*aS, az*aS);
      printf("% .12f % .12f % .12f\n", bx*aS, by*aS, bz*aS);
      printf("% .12f % .12f % .12f\n", cx, cy, cz);
      cell_rows=3; next
    }
    in_cell && cell_rows>0 { cell_rows--; if(cell_rows==0) in_cell=0; next }

    /^ATOMIC_POSITIONS[[:space:]]*\(crystal\)/ {print; in_atoms=1; n=0; next}
    in_atoms {
      if(n>=nat){in_atoms=0}
      if(!in_atoms){print; next}
      n++
      label=$1; x=$2; y=$3; z=$4
      if(z==""){print; next}
      znew = 0.5 + SCALE * ((z+0.0) - 0.5)
      while(znew<0.0) znew+=1.0
      while(znew>=1.0) znew-=1.0
      xnew = x + 0.0; ynew = y + 0.0
      if((z+0.0) > 0.5) { xnew = x + dx; ynew = y + dy }
      while(xnew<0.0) xnew+=1.0
      while(xnew>=1.0) xnew-=1.0
      while(ynew<0.0) ynew+=1.0
      while(ynew>=1.0) ynew-=1.0
      printf("%-3s %20.12f %20.12f %20.12f", label, xnew, ynew, znew)
      for(i=5;i<=NF;i++) printf(" %s", $i)
      printf("\n")
      if(n==nat){in_atoms=0}
      next
    }
    {print}
  ' "${WORKDIR}/template.in" > "${infile}"
  echo "${infile}"
}

gen_ph_input() {
  local infile="${WORKDIR}/aa_aS${A_SCALE}_zS${SCALE}_dx${dx}_dy${dy}.ph.in"
  local outdir="${OUTDIR_TASK}"
  if [[ -f "${WORKDIR}/template.ph.in" ]]; then
    awk -v aS="${A_SCALE}" -v SCALE="${SCALE}" -v dx="${dx}" -v dy="${dy}" -v outdir="${outdir}" '
      BEGIN{in_ph=0}
      /^\s*&INPUTPH/ {in_ph=1}
      in_ph && /^[[:space:]]*\// {in_ph=0}
      in_ph && /^[[:space:]]*prefix[[:space:]]*=/ {sub(/=[[:space:]]*.*/,"= \"" "aa_aS" aS "_zS" SCALE "_dx" dx "_dy" dy "\" ,"); print; next}
      in_ph && /^[[:space:]]*outdir[[:space:]]*=/ {sub(/=[[:space:]]*.*/,"= \"" outdir "\" ,"); print; next}
      {print}
    ' "${WORKDIR}/template.ph.in" > "${infile}"
  else
    cat > "${infile}" <<EOF
&INPUTPH
  prefix = 'aa_aS${A_SCALE}_zS${SCALE}_dx${dx}_dy${dy}',
  outdir = '${outdir}',
  tr2_ph = 1.0d-12,
  epsil = .true.,
  ldisp = .true.,
  nq1 = 1,
  nq2 = 1,
  nq3 = 1,
/
0.0 0.0 0.0
EOF
  fi
  echo "${infile}"
}

pw_in=$(gen_pw_input)
ph_in=$(gen_ph_input)
pw_out="${WORKDIR}/aa_aS${A_SCALE}_zS${SCALE}_dx${dx}_dy${dy}.pw.out"
ph_out="${WORKDIR}/aa_aS${A_SCALE}_zS${SCALE}_dx${dx}_dy${dy}.ph.out"

# Run pw.x
if [[ "${VERBOSE}" -eq 1 ]]; then
  echo "pw.x: aS=${A_SCALE}, zS=${SCALE}, dx=${dx}, dy=${dy}, sep≈${sepA} Å, t=${thickness} Å"
fi
srun --ntasks="${NTASKS}" --cpus-per-task="${OMP_NUM_THREADS}" "${PW_CMD}" -in "${pw_in}" > "${pw_out}" 2>&1

# Run ph.x (allow ph to fail without stopping the task)
if [[ "${VERBOSE}" -eq 1 ]]; then echo "ph.x: dielectric tensor at Gamma"; fi
srun --ntasks="${NTASKS}" --cpus-per-task="${OMP_NUM_THREADS}" "${PH_CMD}" -in "${ph_in}" > "${ph_out}" 2>&1 || true

# Per-task SOAP features + labels (NPZ + CSV) in task directory
mamba run -n coeff-env python3 "${TASK_PY}" \
  --out "${WORKDIR}/features_labels_task_${TASK}.npz" \
  --csv "${WORKDIR}/features_labels_task_${TASK}.csv" \
  --l_max "${SOAP_LMAX:-8}" \
  --n_max "${SOAP_NMAX:-8}" \
  --r_cut "${SOAP_RCUT:-6.0}" \
  --sigma "${SOAP_SIGMA:-0.5}" \
  ${SOAP_PERIODIC:+--periodic} \
  --aggregate "${SOAP_AGGREGATE:-mean}" \
  --stress-unit "${SOAP_STRESS_UNIT:-kbar}"

touch "${WORKDIR}/.done"
