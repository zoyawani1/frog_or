#!/bin/bash
#SBATCH --job-name=filter_ORs_array
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/filter_ORs_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/filter_ORs_%A_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-189%10

module load R/4.4.1-foss-2022b

mkdir -p /project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs
mkdir -p /project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20/chunks

GENOME_LIST="/project/stuckert/users/wani/f2/scripts/genome_dirs.txt"
FILTER_SCRIPT="/project/stuckert/users/wani/f2/scripts/06_Filter_ORs_debug.R"
CHUNK_SCRIPT="/project/stuckert/users/wani/f2/scripts/make_hit_chunks.R"

MISSING_FAI="/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/missing_fai_files_array.txt"
MISSING_SUM="/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/missing_sum_files_array.txt"
FAILED_GENOMES="/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/failed_filter_genomes_array.txt"
FAILED_CHUNKS="/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/failed_filter_chunks_array.txt"

BED_DIR="/project/stuckert/users/wani/f2/anuran_genomes/bed_files"
BLAST_DIR="/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20"

# Change this if you want more or fewer chunks
N_CHUNKS=10

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
sum_file="${BLAST_DIR}/${genome_id}_OR_query_1e-20.sum"

if [ ! -f "$sum_file" ]; then
    echo "WARNING: No .sum file found for $genome_id"
    echo "$genome_id" >> "$MISSING_SUM"
    exit 0
fi

chunk_dir="${BLAST_DIR}/chunks/${genome_id}"
mkdir -p "$chunk_dir"

echo "========================================"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Genome: $genome_id"
echo "Genome dir: $genome_dir"
echo "FAI file: $fai_file"
echo "SUM file: $sum_file"
echo "Chunk dir: $chunk_dir"
echo "N_CHUNKS: $N_CHUNKS"
echo "Using Rscript: $(which Rscript)"
echo "========================================"

###############################################################################
# 1. Split the full .sum file by Hit into balanced chunk files
###############################################################################
echo
echo "=== STEP 1: make chunk files ==="

if ! Rscript "$CHUNK_SCRIPT" "$fai_file" "$sum_file" "$chunk_dir" "$N_CHUNKS"; then
    echo "ERROR: make_hit_chunks.R failed for $genome_id"
    echo "$genome_id" >> "$FAILED_GENOMES"
    exit 1
fi

###############################################################################
# 2. Run filter script on each chunk
###############################################################################
echo
echo "=== STEP 2: run filter script on chunks ==="

shopt -s nullglob
chunk_files=("${chunk_dir}"/chunk_*.sum)
shopt -u nullglob

if [ ${#chunk_files[@]} -eq 0 ]; then
    echo "ERROR: No chunk files were produced for $genome_id"
    echo "$genome_id" >> "$FAILED_GENOMES"
    exit 1
fi

for chunk in "${chunk_files[@]}"; do
    label=$(basename "$chunk" .sum)

    echo
    echo "Running chunk: $label"
    echo "Chunk file: $chunk"

    if ! Rscript "$FILTER_SCRIPT" "$genome_dir" "$fai_name" "$chunk" "$label"; then
        echo "ERROR: Filter script failed for $genome_id / $label"
        echo "${genome_id}    ${label}" >> "$FAILED_CHUNKS"
        echo "$genome_id" >> "$FAILED_GENOMES"
        exit 1
    fi
done

###############################################################################
# 3. Merge chunk outputs back together
###############################################################################
echo
echo "=== STEP 3: merge chunk outputs ==="

final_bed="${BED_DIR}/${genome_id}_complete.cleaned.bed"
final_sum="${BED_DIR}/${genome_id}_complete.cleaned.sum"
final_info="${BED_DIR}/${genome_id}_complete.pullout.info.txt"

rm -f "$final_bed" "$final_sum" "$final_info"

# Merge BED
shopt -s nullglob
bed_files=("${BED_DIR}/${genome_id}"_chunk_*_cleaned.bed)
sum_files=("${BED_DIR}/${genome_id}"_chunk_*_cleaned.sum)
info_files=("${BED_DIR}/${genome_id}"_chunk_*_pullout.info.txt)
shopt -u nullglob

if [ ${#bed_files[@]} -eq 0 ]; then
    echo "ERROR: No chunk BED files found for $genome_id"
    echo "$genome_id" >> "$FAILED_GENOMES"
    exit 1
fi

cat "${bed_files[@]}" | sort -k1,1 -k2,2n | uniq > "$final_bed"

# Merge SUM, keeping one header
first=1
for f in "${sum_files[@]}"; do
    if [ $first -eq 1 ]; then
        cat "$f" > "$final_sum"
        first=0
    else
        tail -n +2 "$f" >> "$final_sum"
    fi
done

# Merge INFO, keeping one header
first=1
for f in "${info_files[@]}"; do
    if [ $first -eq 1 ]; then
        cat "$f" > "$final_info"
        first=0
    else
        tail -n +2 "$f" >> "$final_info"
    fi
done

echo "Merged BED : $final_bed"
echo "Merged SUM : $final_sum"
echo "Merged INFO: $final_info"

echo
echo "Finished genome: $genome_id"
