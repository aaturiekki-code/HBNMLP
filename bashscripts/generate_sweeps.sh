#!/bin/bash
#SBATCH --job-name=hBN_orchestrator
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

set -e

N=1
M0=""
SWEEP=all
MIRROR=both
TAG=""
OFFSET=0
DEV_XY=0.04
DEV_Z=0.04
WORKFLOW=auto

usage() {
    cat <<EOF
Usage: sbatch bashscripts/$(basename "$0") [options]

-n INT    Untwisted NxN expansion
-m INT    Twisted m0 index
-t NAME   all, aa-prime, ab, aa, slide-prime, slide-aa
-p VALUE  Twisted mirror: 0, 1, both
-v TAG    Additional realization tag
-s INT    Seed offset
-x FLOAT  In-plane deviation in Angstrom
-z FLOAT  Out-of-plane deviation in Angstrom
-w MODE   auto, scf-only, electronic-only, full-response
-h        Help

Defaults:
  N=1, sweep=all, mirror=both, deviations=0.04 Angstrom
  auto: untwisted/m0=1 electronic-only; m0=2 SCF-only

Existing folders are skipped. No calculation jobs are submitted.
EOF
}

while getopts "n:m:t:p:v:s:x:z:w:h" option; do
    case "$option" in
        n) N="$OPTARG" ;;
        m) M0="$OPTARG" ;;
        t) SWEEP="$OPTARG" ;;
        p) MIRROR="$OPTARG" ;;
        v) TAG="$OPTARG" ;;
        s) OFFSET="$OPTARG" ;;
        x) DEV_XY="$OPTARG" ;;
        z) DEV_Z="$OPTARG" ;;
        w) WORKFLOW="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

shift $((OPTIND - 1))
if [ "$#" -ne 0 ]; then
    echo "Unexpected arguments: $*" >&2
    exit 1
fi

[[ "$N" =~ ^[1-9][0-9]*$ ]] || {
    echo "Invalid N" >&2; exit 1;
}
[[ -z "$M0" || "$M0" =~ ^[1-9][0-9]*$ ]] || {
    echo "Invalid m0" >&2; exit 1;
}
[[ "$OFFSET" =~ ^(0|[1-9][0-9]*)$ ]] || {
    echo "Invalid seed offset" >&2; exit 1;
}
[[ -z "$TAG" || "$TAG" =~ ^[A-Za-z0-9_-]+$ ]] || {
    echo "Invalid variant tag" >&2; exit 1;
}

case "$SWEEP" in
    all|aa-prime|ab|aa|slide-prime|slide-aa) ;;
    *) echo "Invalid sweep" >&2; exit 1 ;;
esac

case "$MIRROR" in
    0|1|both) ;;
    *) echo "Invalid mirror" >&2; exit 1 ;;
esac

case "$WORKFLOW" in
    auto|scf-only|electronic-only|full-response) ;;
    *) echo "Invalid workflow" >&2; exit 1 ;;
esac

if [[ -n "$M0" && "$SWEEP" != all ]]; then
    echo "-t applies only to untwisted sweeps" >&2
    exit 1
fi

if [[ -z "$M0" && "$MIRROR" != both ]]; then
    echo "Use -p with -m" >&2
    exit 1
fi

if [[ -z "$TAG" && "$OFFSET" != 0 ]]; then
    echo "Use -v when changing seeds" >&2
    exit 1
fi

if [[ "$WORKFLOW" == auto ]]; then
    if [[ "$M0" == 2 ]]; then
        WORKFLOW=scf-only
    elif [[ -z "$M0" || "$M0" == 1 ]]; then
        WORKFLOW=electronic-only
    else
        echo "Choose -w explicitly for other m0 values" >&2
        exit 1
    fi
fi

module purge
module load mamba
source activate coeff-env

python3 - "$DEV_XY" "$DEV_Z" <<'PY'
import math
import sys

for value in sys.argv[1:]:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise SystemExit("Deviations must be finite and nonnegative")
PY

ROOT="/scratch/work/riekkia3/Special_Assignment"
BATCH="$ROOT/pyscripts/batch_create.py"
mkdir -p "$ROOT/QEinputs"

if [[ -n "$M0" ]]; then
    TIME="24:00:00"
    MEM="50G"
else
    case "$N" in
        1) TIME="00:30:00"; MEM="10G" ;;
        2) TIME="05:00:00"; MEM="20G" ;;
        3) TIME="24:00:00"; MEM="50G" ;;
        *) TIME="24:00:00"; MEM="64G" ;;
    esac
fi

COMMON=(
    --vacuum 30
    --dev-xy "$DEV_XY" --dev-z "$DEV_Z"
    --slurm-time "$TIME" --slurm-mem "$MEM" --slurm-cores 16
    "--$WORKFLOW"
)

[[ -z "$TAG" ]] || COMMON+=(--variant-tag "$TAG")

stacking() {
    local name="$1" stack="$2" mirror="$3" seed="$4"

    python3 "$BATCH" \
        --parent-dir "$ROOT/QEinputs/${name}_${N}x${N}" \
        --distance-range 3.0 3.8 9 \
        --scale-range 1.41 1.47 5 \
        --mode bilayer --bilayer-stack "$stack" --mirror "$mirror" \
        --n1 "$N" --n2 "$N" --seed "$((seed + OFFSET))" \
        "${COMMON[@]}"
}

sliding() {
    local name="$1" mirror="$2" seed="$3"

    python3 "$BATCH" \
        --parent-dir "$ROOT/QEinputs/${name}_${N}x${N}" \
        --slide-x-range 0.0 0.5 6 \
        --slide-y-range 0.0 0.5 6 \
        --distance-range 3.0 3.8 3 \
        --mode bilayer --bilayer-stack AA --mirror "$mirror" \
        --n1 "$N" --n2 "$N" --seed "$((seed + OFFSET))" \
        "${COMMON[@]}"
}

twist() {
    local mirror="$1" folder seed

    if [[ "$mirror" -eq 0 ]]; then
        folder="Twist_m0_${M0}_sweep"
        seed=555
    else
        folder="Twist_m0_${M0}_mirror1_sweep"
        seed=666
    fi

    python3 "$BATCH" \
        --parent-dir "$ROOT/QEinputs/$folder" \
        --distance-range 3.0 3.8 9 \
        --scale-range 1.41 1.47 5 \
        --mode twist --m0 "$M0" --mirror "$mirror" \
        --n1 1 --n2 1 --seed "$((seed + OFFSET))" \
        "${COMMON[@]}"
}

echo "Workflow for new configurations: $WORKFLOW"

if [[ -n "$M0" ]]; then
    case "$MIRROR" in
        0) twist 0 ;;
        1) twist 1 ;;
        both) twist 0; twist 1 ;;
    esac
else
    case "$SWEEP" in
        all)
            stacking AA_prime_sweep AA 1 111
            stacking AB_sweep AB 0 222
            stacking AA_sweep AA 0 333
            sliding Sliding_sweep 1 444
            sliding Sliding_AA_sweep 0 777
            ;;
        aa-prime) stacking AA_prime_sweep AA 1 111 ;;
        ab) stacking AB_sweep AB 0 222 ;;
        aa) stacking AA_sweep AA 0 333 ;;
        slide-prime) sliding Sliding_sweep 1 444 ;;
        slide-aa) sliding Sliding_AA_sweep 0 777 ;;
    esac
fi

echo "Generation finished. No calculation jobs submitted."
