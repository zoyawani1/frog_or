#!/usr/bin/env bash
#SBATCH --job-name=anuran_blastdb
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-1
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/blastdb/logs/blastdb_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/blastdb/logs/blastdb_%A_%a.err

set -euo pipefail

# Load modules
module load SAMtools/1.18-GCC-12.3.0
module load BLAST+/2.14.1-gompi-2023a

# -----------------------------
# Absolute paths (IMPORTANT)
# -----------------------------

DATA_DIR="/project/stuckert/users/wani/f2/anuran_genomes/ncbi_dataset/data"
LIST="$DATA_DIR/genomes.fna.list"

BLASTDB_DIR="/project/stuckert/users/wani/f2/anuran_genomes/blastdb"
LOG_DIR="$BLASTDB_DIR/logs"

mkdir -p "$BLASTDB_DIR" "$LOG_DIR"

# -----------------------------
# Select genome for this task
# -----------------------------

fna="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")"

if [[ -z "${fna:-}" ]]; then
  echo "ERROR: No genome found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

# Assembly name comes from folder
asm_dir="$(dirname "$fna")"
asm_name="$(basename "$asm_dir")"

echo "===================================="
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Assembly: $asm_name"
echo "Genome: $fna"
echo "Start time: $(date)"
echo "===================================="

# -----------------------------
# Index + BLAST DB
# -----------------------------

samtools faidx "$fna"

makeblastdb \
  -in "$fna" \
  -out "$BLASTDB_DIR/$asm_name" \
  -dbtype nucl \
  -hash_index

echo "Finished at: $(date)"
