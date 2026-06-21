#!/bin/bash
#SBATCH --job-name=clustalw_tree
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_10_clustalw/logs/clustalw_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_10_clustalw/logs/clustalw_%A_%a.err

set -eo pipefail

source /home/zwani2/miniconda3/etc/profile.d/conda.sh
conda activate orf_env

which clustalw

ALIGN_LIST="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_10_clustalw/outgroup_alignment_files.list"
OUT_BASE="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_10_clustalw"

ALIGNMENT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ALIGN_LIST")

if [[ -z "${ALIGNMENT:-}" ]]; then
    echo "ERROR: No alignment found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$ALIGNMENT" ]]; then
    echo "ERROR: Alignment file does not exist: $ALIGNMENT"
    exit 1
fi

GENOME_ID=$(basename "$(dirname "$ALIGNMENT")")
OUT_DIR="${OUT_BASE}/${GENOME_ID}"
mkdir -p "$OUT_DIR"

OUT_TREE="${OUT_DIR}/step_6_result_mafftAlignment.wOutgroups.phb"

echo "GENOME_ID=$GENOME_ID"
echo "ALIGNMENT=$ALIGNMENT"
echo "OUT_TREE=$OUT_TREE"
date

#changing bootstrap from 1000 to 500 due to errors from high sequence counts
clustalw \
  -CLUSTERING=NJ \
  -bootstrap=100 \
  -KIMURA \
  -TOSSGAPS \
  -BOOTLABELS=node \
  -quiet \
  -infile="$ALIGNMENT" \
  -outfile="$OUT_TREE"

echo "DONE"
date
