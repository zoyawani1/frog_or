#!/bin/bash
#SBATCH --job-name=getfasta_ORs
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_3_getfasta/logs/getfasta_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_3_getfasta/logs/getfasta_%A_%a.err

set -euo pipefail

module purge
module load BEDTools

BED_DIR="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_2_prunning/final_bed_files"
GENOME_BASE="/project/stuckert/users/wani/f2/anuran_genomes/ncbi_dataset/data"
OUT_BASE="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_3_getfasta"
LOG_DIR="${OUT_BASE}/logs"

mkdir -p "$OUT_BASE" "$LOG_DIR"

BED_LIST="${OUT_BASE}/combined_bed_files.list"

BED_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BED_LIST")

if [[ -z "${BED_FILE:-}" ]]; then
    echo "ERROR: No BED file found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BED_BASENAME=$(basename "$BED_FILE")
GENOME_ID="${BED_BASENAME%_combined.bed}"

GENOME_DIR="${GENOME_BASE}/${GENOME_ID}"
GENOME_FASTA=$(find "$GENOME_DIR" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | head -n 1)

if [[ -z "${GENOME_FASTA:-}" ]]; then
    echo "ERROR: No genome FASTA found for ${GENOME_ID} in ${GENOME_DIR}"
    exit 1
fi

OUT_DIR="${OUT_BASE}/${GENOME_ID}"
mkdir -p "$OUT_DIR"

OUT_FASTA="${OUT_DIR}/${GENOME_ID}_olfacUniqueBlastHits.stranded.fasta"

echo "BED_FILE=$BED_FILE"
echo "GENOME_ID=$GENOME_ID"
echo "GENOME_FASTA=$GENOME_FASTA"
echo "OUT_FASTA=$OUT_FASTA"

bedtools getfasta -s -name \
    -bed "$BED_FILE" \
    -fi "$GENOME_FASTA" \
    -fo "$OUT_FASTA"

echo "DONE"
