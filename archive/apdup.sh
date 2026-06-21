#!/bin/bash

set -euo pipefail

DUP_DIR="./duplicates"
LOG_DIR="/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20/log1e-20"
LOG_FILE="${LOG_DIR}/duplicate_gca_archived.log"

mkdir -p "$DUP_DIR"
mkdir -p "$LOG_DIR"

echo "Run started: $(date)" > "$LOG_FILE"
echo "Searching for GCA/GCF duplicate sum files..." >> "$LOG_FILE"
echo >> "$LOG_FILE"

shopt -s nullglob

for gcf_file in GCF_*_OR_query_1e-20.sum; do
    accession_part="${gcf_file#GCF_}"
    matching_gca="GCA_${accession_part}"

    if [[ -f "$matching_gca" ]]; then
        mv "$matching_gca" "$DUP_DIR/"
        echo "Moved: $matching_gca -> $DUP_DIR/" | tee -a "$LOG_FILE"
    fi
done

echo >> "$LOG_FILE"
echo "Run finished: $(date)" >> "$LOG_FILE"
