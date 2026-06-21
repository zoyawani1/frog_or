#!/bin/bash
#SBATCH --job-name=einsi_refiltered
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_7_einsi_refiltered/logs/einsi2_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_7_einsi_refiltered/logs/einsi2_%A_%a.err

set -euo pipefail

source /project/stuckert/software/anaconda3/etc/profile.d/conda.sh
conda activate trees

which einsi
which mafft

FILTERED_LIST="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_7_einsi_refiltered/filtered_fasta_files.list"
OUT_BASE="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_7_einsi_refiltered"

FILTERED_FASTA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILTERED_LIST")

if [[ -z "${FILTERED_FASTA:-}" ]]; then
    echo "ERROR: No filtered FASTA found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$FILTERED_FASTA" ]]; then
    echo "ERROR: Filtered FASTA does not exist: $FILTERED_FASTA"
    exit 1
fi

GENOME_ID=$(basename "$(dirname "$FILTERED_FASTA")")
OUT_DIR="${OUT_BASE}/${GENOME_ID}"
mkdir -p "$OUT_DIR"

REALIGNED="${OUT_DIR}/${GENOME_ID}_step_4_result_mafftAlignment_removed.fasta"

echo "GENOME_ID=$GENOME_ID"
echo "FILTERED_FASTA=$FILTERED_FASTA"
echo "REALIGNED=$REALIGNED"

einsi --preservecase --thread 8 --inputorder "$FILTERED_FASTA" > "$REALIGNED"

echo "DONE"
