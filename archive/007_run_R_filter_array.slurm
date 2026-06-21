#!/bin/bash
#SBATCH --job-name=filter_ORs_array
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/filter_ORs_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/filter_ORs_%A_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=400G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-189%10

module load R/4.4.1-foss-2022b
#which Rscript
#Rscript -e "library(GenomicRanges); print('GenomicRanges loaded successfully')"





mkdir -p /project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs

GENOME_LIST="/project/stuckert/users/wani/f2/scripts/genome_dirs.txt"
SCRIPT="/project/stuckert/users/wani/f2/scripts/06_Filter_ORs_debug.R"
MISSING_FAI="/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/missing_fai_files_array.txt"
FAILED_GENOMES="/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/failed_filter_genomes_array.txt"

genome_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENOME_LIST")

if [ -z "$genome_dir" ]; then
    echo "ERROR: No genome directory found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

genome_id=$(basename "$genome_dir")

fai_file=$(find "$genome_dir" -maxdepth 1 -name "*.fna.fai" | head -n 1)

if [ -z "$fai_file" ]; then
    echo "WARNING: No .fna.fai file found for $genome_id"
    echo "$genome_id" >> "$MISSING_FAI"
    exit 0
fi

fai_name=$(basename "$fai_file")

echo "========================================"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Genome: $genome_id"
echo "Genome dir: $genome_dir"
echo "FAI file: $fai_name"
echo "Using Rscript: $(which Rscript)"
echo "========================================"

if Rscript "$SCRIPT" "$genome_dir" "$fai_name"; then
    echo "Finished genome: $genome_id"
else
    echo "ERROR: Rscript failed for $genome_id"
    echo "$genome_id" >> "$FAILED_GENOMES"
    exit 1
fi
