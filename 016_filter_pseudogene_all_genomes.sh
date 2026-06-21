#!/bin/bash
#SBATCH --job-name=filter_pseudo
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_6_pseudogene/logs/pseudo_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_6_pseudogene/logs/pseudo_%A_%a.err

set -euo pipefail

PYTHON="/home/zwani2/miniconda3/envs/orfilter/bin/python"

SCRIPT_DIR="/project/stuckert/users/wani/f2/scripts"
PSEUDO_SCRIPT="${SCRIPT_DIR}/4_filter_pseudogene.py"
ALIGN_LIST="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_3_getfasta/mafft_files.list"
DOMAIN_FILE="${SCRIPT_DIR}/HumanORTMD.txt"

OUT_BASE="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_6_pseudogene"
LOG_DIR="${OUT_BASE}/logs"

mkdir -p "$LOG_DIR"

ALIGNMENT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ALIGN_LIST")

if [[ -z "${ALIGNMENT:-}" ]]; then
    echo "ERROR: No alignment found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$ALIGNMENT" ]]; then
    echo "ERROR: alignment file does not exist: $ALIGNMENT"
    exit 1
fi

if [[ ! -f "$PYTHON" ]]; then
    echo "ERROR: Python executable does not exist: $PYTHON"
    exit 1
fi

if [[ ! -f "$PSEUDO_SCRIPT" ]]; then
    echo "ERROR: pseudogene script does not exist: $PSEUDO_SCRIPT"
    exit 1
fi

if [[ ! -f "$DOMAIN_FILE" ]]; then
    echo "ERROR: domain file does not exist: $DOMAIN_FILE"
    exit 1
fi

GENOME_ID=$(basename "$(dirname "$ALIGNMENT")")
OUT_DIR="${OUT_BASE}/${GENOME_ID}"
mkdir -p "$OUT_DIR"

echo "PYTHON=$PYTHON"
"$PYTHON" -c "import pandas, numpy; print('pandas/numpy OK')"

echo "GENOME_ID=$GENOME_ID"
echo "ALIGNMENT=$ALIGNMENT"
echo "DOMAIN_FILE=$DOMAIN_FILE"
echo "OUT_DIR=$OUT_DIR"

cd "$SCRIPT_DIR"

"$PYTHON" "$PSEUDO_SCRIPT" \
  --save_fasta True \
  --save_data True \
  --query_file_path "$ALIGNMENT" \
  --domain_file_path "$DOMAIN_FILE" \
  --output_path "$OUT_DIR"

echo "DONE"
