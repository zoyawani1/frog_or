#!/bin/bash
#SBATCH --job-name=tblastn_coqui
#SBATCH --output=tblastn_coqui.out
#SBATCH --error=tblastn_coqui.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

module purge
# activate your env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate orf_env

tblastn \
  -evalue 1e-20 \
  -query /project/stuckert/users/wani/f2/scripts/OR_query.fasta \
  -db /project/stuckert/users/wani/f2/genomes/Eleutherodactylus_coqui/ncbi_dataset/data/GCF_035609145.1/Eleutherodactylus_coqui \
  -out /project/stuckert/users/wani/f2/blast_results/test_coqui.out \
  -outfmt 5 \
  -max_target_seqs 20000 \
  -num_threads 4
