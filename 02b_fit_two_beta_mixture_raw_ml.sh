#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH -c 4
#SBATCH --job-name=two_beta_raw_ml
#SBATCH -o logs/02b_fit_two_beta_mixture_raw_ml_%j.txt
#SBATCH -e logs/02b_fit_two_beta_mixture_raw_ml_%j.txt
#SBATCH --time=24:00:00

set -euo pipefail

BATCH="${1:-m6A_2OmeA_I_genome}"
SCRIPT_DIR=/dcs04/hicks/data/sparthib/drna_retinalrs/code/04_modkit_stats/01_ml_probability

mkdir -p logs
module load conda_R/4.4

 echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "Batch: ${BATCH}"

Rscript "${SCRIPT_DIR}/02b_fit_two_beta_mixture_raw_ml.R" \
    "${SCRIPT_DIR}/config_${BATCH}.R"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
