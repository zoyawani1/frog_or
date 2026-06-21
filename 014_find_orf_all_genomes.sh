#!/bin/bash
#SBATCH --job-name=find_ORF_all_genomes
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_4_orfs/logs/orf_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_4_orfs/logs/orf_%A_%a.err

set -euo pipefail

FASTA_LIST="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_4_orfs/stranded_fasta_files.list"
PERL_SCRIPT="/project/stuckert/users/wani/f2/scripts/013_find_orf.pl"
EXPECT_LENGTH=810

FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FASTA_LIST")

if [[ -z "${FASTA_FILE:-}" ]]; then
    echo "ERROR: No FASTA file found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$FASTA_FILE" ]]; then
    echo "ERROR: FASTA file does not exist: $FASTA_FILE"
    exit 1
fi

OUT_DIR=$(dirname "$FASTA_FILE")

echo "FASTA_FILE=$FASTA_FILE"
echo "OUT_DIR=$OUT_DIR"
echo "EXPECT_LENGTH=$EXPECT_LENGTH"

cd "$OUT_DIR"

perl "$PERL_SCRIPT" "$FASTA_FILE" "$EXPECT_LENGTH"

echo "DONE"
