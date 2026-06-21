#!/bin/bash
set -euo pipefail

BASE_DIR="/project/stuckert/users/wani/f2/anuran_genomes/filtering_step_2_prunning"
BED_DIR="${BASE_DIR}/bed_files"
OUT_DIR="${BASE_DIR}/final_bed_files"

mkdir -p "$OUT_DIR"

echo "Combining cleaned BED files per genome..."

for genome_dir in "$BED_DIR"/GC[AF]_*; do
    [ -d "$genome_dir" ] || continue

    genome=$(basename "$genome_dir")
    out_file="${OUT_DIR}/${genome}_combined.bed"

    echo "Processing $genome"

    find "$genome_dir" -type f -name "*_cleaned.bed" -print0 | \
      xargs -0 cat | \
      sort -k1,1 -k2,2n > "$out_file"

    echo "Wrote $out_file"
done

echo "DONE"
