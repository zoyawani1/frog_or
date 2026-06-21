#!/bin/bash
#SBATCH --job-name=filter_sum_1e-20
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20/log1e-20/filter_sum_1e-20_%j.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e^20/log1e-20/filter_sum_1e-20_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1

INPUT_DIR="/project/stuckert/users/wani/f2/anuran_genomes/blast_results"
OUTPUT_DIR="/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/log1e-20"

echo "Starting evalue filtering job"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

for f in "$INPUT_DIR"/*_OR_query.sum; do
    [ -f "$f" ] || continue

    base=$(basename "$f" .sum)
    out="$OUTPUT_DIR/${base}_1e-20.sum"

    echo "Processing: $f"
    echo "Output: $out"

    awk '
    BEGIN { FS=OFS="/" }
    NR==1 {
        print
        for (i=1; i<=NF; i++) {
            if ($i == "E-value") ecol=i
        }
        if (!ecol) {
            print "ERROR: no E-value column found in " FILENAME > "/dev/stderr"
            exit 1
        }
        next
    }
    ($ecol+0) <= 1e-20 { print }
    ' "$f" > "$out"

    if [ $? -ne 0 ]; then
        echo "ERROR processing $f" >&2
    else
        echo "Finished: $out"
    fi
done

echo "Done"
