#!/usr/bin/env bash

set -euo pipefail

# -------------------------------
# Configuration
# -------------------------------
TAXON_ID=8342
OUTDIR="/project/stuckert/users/wani/f2/anuran_genomes"
LOGDIR="${OUTDIR}/logs"
ASSEMBLY_LEVELS="complete,chromosome,scaffold,contig"

# -------------------------------
# Create output directory
# -------------------------------
mkdir -p "${OUTDIR}"
mkdir -p "${LOGDIR}"
# -------------------------------
# Download genomes
# -------------------------------
echo "Downloading anuran genomes into ${OUTDIR}/ ..."

datasets download genome taxon "${TAXON_ID}" \
  --assembly-level "${ASSEMBLY_LEVELS}" \
  --include genome \
  --filename "${OUTDIR}/anura_genomes.zip"

# -------------------------------
# Extract
# -------------------------------
echo "Extracting genomes..."
unzip -o "${OUTDIR}/anura_genomes.zip" -d "${OUTDIR}"
rm "${OUTDIR}/anura_genomes.zip"

echo "Done!"
echo "Genomes are in: ${OUTDIR}/ncbi_dataset/data/"
