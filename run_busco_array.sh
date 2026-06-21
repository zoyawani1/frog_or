#!/bin/bash
#SBATCH --job-name=busco_frogs
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-225
#SBATCH --output=/project/stuckert/users/wani/f2/busco/busco_genome_results/logs/busco_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/busco/busco_genome_results/logs/busco_%A_%a.err

source ~/.bashrc
conda activate /project/stuckert/users/wani/f2/busco/busco_env

# 👉 UPDATED: list is in scripts folder
LIST="/project/stuckert/users/wani/f2/scripts/genome_fasta_list.txt"

BUSCO_DIR="/project/stuckert/users/wani/f2/busco"
OUTDIR="${BUSCO_DIR}/busco_genome_results"
LOGDIR="${OUTDIR}/logs"
LINEAGE="tetrapoda_odb12"

# Make sure directories exist
mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [ -z "$GENOME" ]; then
  echo "No genome found for array task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

ACC=$(basename "$(dirname "$GENOME")")

echo "Running BUSCO"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Genome: $GENOME"
echo "Accession: $ACC"

busco \
  -i "$GENOME" \
  -l "$LINEAGE" \
  -m genome \
  -o "$ACC" \
  --out_path "$OUTDIR" \
  -c "$SLURM_CPUS_PER_TASK"

echo "Finished BUSCO for $ACC"
