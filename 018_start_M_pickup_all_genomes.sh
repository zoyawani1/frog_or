#!/bin/bash
#SBATCH --job-name=start_M_pickup
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_8_startM/logs/startM_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_8_startM/logs/startM_%A_%a.err

set -euo pipefail

PYTHON="/home/zwani2/miniconda3/envs/orfilter/bin/python"
SCRIPT_DIR="/project/stuckert/users/wani/f2/scripts"
DOMAIN_FILE="${SCRIPT_DIR}/HumanORTMD.txt"
REALIGNED_LIST="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_8_startM/realigned_filtered_files.list"
OUT_BASE="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_8_startM"

REALIGNED_FASTA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$REALIGNED_LIST")

if [[ -z "${REALIGNED_FASTA:-}" ]]; then
    echo "ERROR: No realigned FASTA for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$REALIGNED_FASTA" ]]; then
    echo "ERROR: FASTA does not exist: $REALIGNED_FASTA"
    exit 1
fi

GENOME_ID=$(basename "$(dirname "$REALIGNED_FASTA")")
OUT_DIR="${OUT_BASE}/${GENOME_ID}"
mkdir -p "$OUT_DIR"

# Copy realigned fasta into step 8 folder with the exact filename Perl expects
cp "$REALIGNED_FASTA" "${OUT_DIR}/${GENOME_ID}_step_4_result_mafftAlignment_removed.fasta"

echo "GENOME_ID=$GENOME_ID"
echo "REALIGNED_FASTA=$REALIGNED_FASTA"
echo "OUT_DIR=$OUT_DIR"

output=$("$PYTHON" "${SCRIPT_DIR}/4_filter_pseudogene.py" -get_tm1_pos \
  --query_file_path "${OUT_DIR}/${GENOME_ID}_step_4_result_mafftAlignment_removed.fasta" \
  --domain_file_path "$DOMAIN_FILE")

echo "$output"

tm1_pos=$(echo "$output" | grep -o '[0-9]*' | tail -1)

echo "TM1_POS=$tm1_pos"

perl "${SCRIPT_DIR}/5_Start_M_pick_up.120.pl" \
  "$OUT_DIR" \
  "$tm1_pos"

sed -i 's/>//g' "${OUT_DIR}/step_5_result.pickedM.CoordinatesOfBestStart.txt"
sed -i 's/Human_OR2J3[[:space:]]*//g' "${OUT_DIR}/step_5_result.pickedM.CoordinatesOfBestStart.txt"

echo "DONE"
