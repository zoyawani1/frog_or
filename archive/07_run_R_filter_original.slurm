#!/bin/bash
#SBATCH --job-name=filter_ORs_with_checks
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/filter_sum_logs/1_filter_ORs_%j.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/filter_sum_logs/1_filter_ORs_%j.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=1


DATA_DIR="/project/stuckert/users/wani/f2/anuran_genomes/ncbi_dataset/data"
SCRIPT="/project/stuckert/users/wani/f2/scripts/06_Filter_ORs_debug.R"
MISSING_FAI="/project/stuckert/users/wani/f2/anuran_genomes/filter_sum_logs/missing_fai_files.txt"
FAILED_GENOMES="/project/stuckert/users/wani/f2/anuran_genomes/filter_sum_logs/failed_filter_genomes.txt"

for genome_dir in "$DATA_DIR"/GCA_* "$DATA_DIR"/GCF_*; do
    [ -d "$genome_dir" ] || continue

    genome_id=$(basename "$genome_dir")

    fai_file=$(find "$genome_dir" -maxdepth 1 -name "*.fna.fai" | head -n 1)

    if [ -z "$fai_file" ]; then
        echo "WARNING: No .fna.fai file found for $genome_id"
        echo "$genome_id" >> "$MISSING_FAI"
        continue
    fi

    fai_name=$(basename "$fai_file")

    echo "========================================"
    echo "Running genome: $genome_id"
    echo "Genome dir: $genome_dir"
    echo "FAI file: $fai_name"

    if Rscript "$SCRIPT" "$genome_dir" "$fai_name"; then
        echo "Finished genome: $genome_id"
    else
        echo "ERROR: Rscript failed for $genome_id"
        echo "$genome_id" >> "$FAILED_GENOMES"
    fi

    echo "========================================"
done
