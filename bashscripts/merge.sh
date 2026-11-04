#!/usr/bin/env bash
#
# SBATCH directives (adjust as needed)
#SBATCH --job-name=qe_merge_all
#SBATCH --output=qe_merge_all_%j.out
#SBATCH --error=qe_merge_all_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --export=ALL

set -euo pipefail
module restore qemod

# Config (relative to bashscripts/)
PARENT_DIR=${PARENT_DIR:-"../runs_interlayer"}
MERGE_ALL_PY=${MERGE_ALL_PY:-"../pyscripts/merge_all_to_one.py"}

# Optional overrides
OUT_NPZ=${OUT_NPZ:-"${PARENT_DIR}/combined_features_labels_all.npz"}
OUT_CSV=${OUT_CSV:-"${PARENT_DIR}/combined_features_labels_all.csv"}
INCLUDE_PATTERNS=${INCLUDE:-}     # e.g., 'dielectric_sweep_*'
EXCLUDE_PATTERNS=${EXCLUDE:-}     # e.g., '*tmp*'
OVERWRITE_FLAG=${OVERWRITE:-0}
DRY_RUN_FLAG=${DRY_RUN:-0}

[[ -d "${PARENT_DIR}" ]] || { echo "Error: parent dir not found: ${PARENT_DIR}" >&2; exit 1; }
[[ -f "${MERGE_ALL_PY}" ]] || { echo "Error: merge script not found: ${MERGE_ALL_PY}" >&2; exit 1; }

echo "Merging ALL per-task NPZ/CSV into ONE NPZ/CSV"
echo "Parent:    ${PARENT_DIR}"
echo "Out NPZ:   ${OUT_NPZ}"
echo "Out CSV:   ${OUT_CSV}"
[[ -n "${INCLUDE_PATTERNS}" ]] && echo "Include:   ${INCLUDE_PATTERNS}"
[[ -n "${EXCLUDE_PATTERNS}" ]] && echo "Exclude:   ${EXCLUDE_PATTERNS}"
[[ "${OVERWRITE_FLAG}" = "1" ]] && echo "Overwrite: enabled" || echo "Overwrite: disabled"
[[ "${DRY_RUN_FLAG}" = "1" ]] && echo "Dry run:   enabled" || echo "Dry run:   disabled"

args=( --parent-dir "${PARENT_DIR}" --out-npz "${OUT_NPZ}" --out-csv "${OUT_CSV}" )
# shellcheck disable=SC2206
if [[ -n "${INCLUDE_PATTERNS}" ]]; then args+=( --include ${INCLUDE_PATTERNS} ); fi
# shellcheck disable=SC2206
if [[ -n "${EXCLUDE_PATTERNS}" ]]; then args+=( --exclude ${EXCLUDE_PATTERNS} ); fi
if [[ "${OVERWRITE_FLAG}" = "1" ]]; then args+=( --overwrite ); fi
if [[ "${DRY_RUN_FLAG}" = "1" ]]; then args+=( --dry-run ); fi

mamba run -n coeff-env python3 "${MERGE_ALL_PY}" "${args[@]}"

echo "Global merge completed."
