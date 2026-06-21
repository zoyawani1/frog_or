#!/usr/bin/env bash
#SBATCH --job-name=anuran_dl
#SBATCH --output=/project/stuckert/users/wani/f2/anuran_genomes/logs/anuran_dl_%j.out
#SBATCH --error=/project/stuckert/users/wani/f2/anuran_genomes/logs/anuran_dl_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=6-00:00:00

set -euo pipefail

BASE="/project/stuckert/users/wani/f2/anuran_genomes"

mkdir -p "$BASE"
mkdir -p "$BASE/logs"

cd /project/stuckert/users/wani/f2/scripts
bash download_anuran_genomes.sh
