#!/bin/bash
#SBATCH --job-name=test_pseudo
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_6_pseudogene/logs/test_%j.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_6_pseudogene/logs/test_%j.err

set -euo pipefail

OUTDIR=/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_6_pseudogene/GCA_002284835.2

mkdir -p "$OUTDIR"

echo "Starting pseudogene filtering"
date

/home/zwani2/miniconda3/envs/orfilter/bin/python \
/project/stuckert/users/wani/f2/scripts/4_filter_pseudogene.py \
  --save_fasta True \
  --save_data True \
  --query_file_path /project/stuckert/users/wani/f2/anuran_genomes/filtering_step_3_getfasta/GCA_002284835.2/GCA_002284835.2_step_4_result_mafftAlignment.fasta \
  --domain_file_path /project/stuckert/users/wani/f2/scripts/HumanORTMD.txt \
  --output_path "$OUTDIR"

echo "DONE"
date
