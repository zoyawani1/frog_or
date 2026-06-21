#!/bin/bash
#SBATCH --job-name=einsi_ORs
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_5_einsi/logs/einsi_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_5_einsi/logs/einsi_%A_%a.err

set -euo pipefail

source /project/stuckert/software/anaconda3/etc/profile.d/conda.sh
conda activate trees

which einsi
which mafft

ORF_LIST="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_5_einsi/orf_fasta_files.list"
HUMAN_REF="/project/stuckert/users/wani/f2/scripts/Human_OR2J3.fasta"

ORF_FASTA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ORF_LIST")

if [[ -z "${ORF_FASTA:-}" ]]; then
    echo "ERROR: No ORF FASTA found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$ORF_FASTA" ]]; then
    echo "ERROR: ORF FASTA does not exist: $ORF_FASTA"
    exit 1
fi

if [[ ! -f "$HUMAN_REF" ]]; then
    echo "ERROR: Human reference does not exist: $HUMAN_REF"
    exit 1
fi

OUT_DIR=$(dirname "$ORF_FASTA")
ORF_BASENAME=$(basename "$ORF_FASTA")
GENOME_ID=$(basename "$OUT_DIR")

ADDED_HUMAN="${OUT_DIR}/${GENOME_ID}_ORF_longThan_810_bp_addedHuman_OR2J3.fasta"
ALIGNMENT="${OUT_DIR}/${GENOME_ID}_step_4_result_mafftAlignment.fasta"

echo "GENOME_ID=$GENOME_ID"
echo "ORF_FASTA=$ORF_FASTA"
echo "ADDED_HUMAN=$ADDED_HUMAN"
echo "ALIGNMENT=$ALIGNMENT"

cat "$HUMAN_REF" "$ORF_FASTA" > "$ADDED_HUMAN"

einsi --preservecase --thread 8 --inputorder "$ADDED_HUMAN" > "$ALIGNMENT"

echo "DONE"
