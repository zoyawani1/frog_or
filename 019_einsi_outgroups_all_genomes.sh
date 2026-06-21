#!/bin/bash
#SBATCH --job-name=einsi_outgroups
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_9_einsi_outgroups/logs/einsi_out_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_9_einsi_outgroups/logs/einsi_out_%A_%a.err

set -euo pipefail

source /project/stuckert/software/anaconda3/etc/profile.d/conda.sh
conda activate trees

SCRIPT_DIR="/project/stuckert/users/wani/f2/scripts"
PICKEDM_LIST="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_9_einsi_outgroups/pickedM_files.list"
OUT_BASE="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_9_einsi_outgroups"

OUTGROUPS="${SCRIPT_DIR}/niimura.genbank.outgroupSeqs.fasta"
ZEBRA="${SCRIPT_DIR}/zebra_query.fasta"

PICKEDM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PICKEDM_LIST")

if [[ ! -f "$PICKEDM" ]]; then
  echo "ERROR: pickedM file missing: $PICKEDM"
  exit 1
fi

if [[ ! -f "$OUTGROUPS" ]]; then
  echo "ERROR: missing outgroups: $OUTGROUPS"
  exit 1
fi

if [[ ! -f "$ZEBRA" ]]; then
  echo "ERROR: missing zebra query: $ZEBRA"
  exit 1
fi

GENOME_ID=$(basename "$(dirname "$PICKEDM")")
OUT_DIR="${OUT_BASE}/${GENOME_ID}"
mkdir -p "$OUT_DIR"

COMBINED="${OUT_DIR}/step_5_result.pickedM.wOutgroups.wRepresentatives.fasta"
ALIGNMENT="${OUT_DIR}/step_6_result_mafftAlignment.wOutgroups.fasta"

echo "GENOME_ID=$GENOME_ID"
echo "PICKEDM=$PICKEDM"
echo "COMBINED=$COMBINED"
echo "ALIGNMENT=$ALIGNMENT"

cat "$OUTGROUPS" "$ZEBRA" "$PICKEDM" > "$COMBINED"

einsi --preservecase --thread 8 --inputorder "$COMBINED" > "$ALIGNMENT"

echo "DONE"
