#!/bin/bash
#SBATCH --job-name=step2_prune_chunks
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_2_prunning/logs/prune_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_2_prunning/logs/prune_%A_%a.err

set -euo pipefail

module purge
module load R

STEP2_BASE="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_2_prunning"
MASTER_LIST="${STEP2_BASE}/chunk_files_master"
SCRIPT_PATH="/project/stuckert/users/wani/f2/scripts/08_prune_chunk_with_granges.R"
GENOME_BASE="/project/stuckert/users/wani/f2/anuran_genomes/ncbi_dataset/data"

LOG_DIR="${STEP2_BASE}/logs"
SUM_OUT_BASE="${STEP2_BASE}/cleaned_sum"
BED_OUT_BASE="${STEP2_BASE}/bed_files"

mkdir -p "${LOG_DIR}" "${SUM_OUT_BASE}" "${BED_OUT_BASE}"

CHUNK_SUM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MASTER_LIST}")

if [[ -z "${CHUNK_SUM_FILE:-}" ]]; then
    echo "ERROR: no chunk file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

GENOME_ID=$(basename "$(dirname "${CHUNK_SUM_FILE}")")
GENOME_DIR="${GENOME_BASE}/${GENOME_ID}"

if [[ ! -d "${GENOME_DIR}" ]]; then
    echo "ERROR: genome directory not found: ${GENOME_DIR}"
    exit 1
fi

FAI_PATH=$(find "${GENOME_DIR}" -maxdepth 1 -type f -name "*.fai" | head -n 1)

if [[ -z "${FAI_PATH:-}" ]]; then
    echo "ERROR: no .fai file found for genome ${GENOME_ID}"
    exit 1
fi

GENOME_SUM_OUT_DIR="${SUM_OUT_BASE}/${GENOME_ID}"
GENOME_BED_OUT_DIR="${BED_OUT_BASE}/${GENOME_ID}"

mkdir -p "${GENOME_SUM_OUT_DIR}" "${GENOME_BED_OUT_DIR}"

echo "=== 09_run_prune_chunk_array.sh ==="
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "CHUNK_SUM_FILE=${CHUNK_SUM_FILE}"
echo "GENOME_ID=${GENOME_ID}"
echo "GENOME_DIR=${GENOME_DIR}"
echo "FAI_PATH=${FAI_PATH}"
echo "GENOME_SUM_OUT_DIR=${GENOME_SUM_OUT_DIR}"
echo "GENOME_BED_OUT_DIR=${GENOME_BED_OUT_DIR}"
echo "SCRIPT_PATH=${SCRIPT_PATH}"

Rscript "${SCRIPT_PATH}" \
  "${CHUNK_SUM_FILE}" \
  "${FAI_PATH}" \
  "${GENOME_ID}" \
  "${GENOME_SUM_OUT_DIR}" \
  "${GENOME_BED_OUT_DIR}"

echo "DONE"
