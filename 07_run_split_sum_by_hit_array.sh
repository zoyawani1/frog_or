#!/bin/bash
#SBATCH --job-name=step_1_filt_split_by_hit
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-189
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_1_chunks/logs/split_by_hit_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_1_chunks/logs/split_by_hit_%A_%a.err

set -euo pipefail

module purge
module load R

SUM_DIR="/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20"
CHUNK_BASE_DIR="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_1_chunks"
LOG_DIR="${CHUNK_BASE_DIR}/logs"
SCRIPT_PATH="/project/stuckert/users/wani/f2/scripts/06_split_sum_by_hit.R"

mkdir -p "${CHUNK_BASE_DIR}"
mkdir -p "${LOG_DIR}"

SUM_LIST="${CHUNK_BASE_DIR}/sum_files.list"

if [[ ! -f "${SUM_LIST}" ]]; then
    find "${SUM_DIR}" -maxdepth 1 -type f -name "*.sum" | sort > "${SUM_LIST}"
fi

SUM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SUM_LIST}")

if [[ -z "${SUM_FILE:-}" ]]; then
    echo "ERROR: No sum file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BASENAME=$(basename "${SUM_FILE}")
GENOME_ID="${BASENAME%_OR_query_1e-20.sum}"

echo "=== 07_run_split_sum_by_hit_array.sh ==="
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "SUM_FILE=${SUM_FILE}"
echo "BASENAME=${BASENAME}"
echo "GENOME_ID=${GENOME_ID}"
echo "CHUNK_BASE_DIR=${CHUNK_BASE_DIR}"
echo "LOG_DIR=${LOG_DIR}"
echo "SCRIPT_PATH=${SCRIPT_PATH}"

Rscript "${SCRIPT_PATH}" "${SUM_FILE}" "${GENOME_ID}" "${CHUNK_BASE_DIR}"

echo "DONE"
