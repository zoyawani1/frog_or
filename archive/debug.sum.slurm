#!/bin/bash
#SBATCH --job-name=debug_one
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filter_sum_logs/debug_two_%j.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filter_sum_logs/debug_two_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=1

module load R/4.4.1-foss-2022b

Rscript /project/stuckert/users/wani/f2/scripts/2_Filter_ORs_debug.R /project/stuckert/users/wani/f2/anuran_genomes/ncbi_dataset/data/GCA_004786255.1 GCA_004786255.1_Pads_1.0_genomic.fna.fai







#
#echo "Using Rscript: $(which Rscript)"
#Rscript --version

#SCRIPT="/project/stuckert/users/wani/f2/scripts/2_Filter_ORs.R"
#GENOME_DIR="/project/stuckert/users/wani/f2/anuran_genomes/ncbi_dataset/data/GCA_002284835.2"
#FAI_NAME="GCA_002284835.2_RCv2.1_genomic.fna.fai"

#echo "SCRIPT path: $SCRIPT"
#echo "Genome dir: $GENOME_DIR"
#echo "FAI name: $FAI_NAME"

#Rscript "$SCRIPT" "$GENOME_DIR" "$FAI_NAME"
