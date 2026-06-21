#!/bin/bash
#SBATCH --job-name=tblastn_ORs
#SBATCH --time=12-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --array=1-225
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/blast_results/logs/tblastn_%A_%a.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/blast_results/logs/tblastn_%A_%a.err

module load BLAST+/2.14.1-gompi-2023a

# Paths
DBDIR="/project/stuckert/users/wani/f2/anuran_genomes/blastdb"
LIST="${DBDIR}/db_roots.list"

SCRIPT_DIR="/project/stuckert/users/wani/f2/scripts"
PERL_SCRIPT="${SCRIPT_DIR}/1_ab_tblastN_ORs.pl"
QUERY="${SCRIPT_DIR}/OR_query.fasta"

OUTDIR="/project/stuckert/users/wani/f2/anuran_genomes/blast_results"
LOGDIR="${OUTDIR}/logs"

# Ensure directories exist
mkdir -p "${OUTDIR}"
mkdir -p "${LOGDIR}"

# Select database for this array task
DBROOT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}")
DB="${DBDIR}/${DBROOT}"

# Run BLAST via Perl wrapper
perl "${PERL_SCRIPT}" "${DB}" "${QUERY}" "${OUTDIR}" "${DBROOT}"

