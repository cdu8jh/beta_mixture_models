#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH --job-name=two_beta_plot
#SBATCH -o logs/03b_plot_two_beta_mixture_raw_ml_%j.txt
#SBATCH -e logs/03b_plot_two_beta_mixture_raw_ml_%j.txt
#SBATCH --time=01:00:00

set -euo pipefail
BATCH="${1:-m6A_2OmeA_I_genome}"
SCRIPT_DIR=/dcs04/hicks/data/sparthib/drna_retinalrs/code/04_modkit_stats/01_ml_probability
OUT_DIR=/dcs04/hicks/data/sparthib/drna_retinalrs/processed_data/04_pileup/01_ml_probability/${BATCH}/02b_two_beta_mixture_raw_ml

mkdir -p logs
module load conda_R/4.4
Rscript "${SCRIPT_DIR}/03b_plot_two_beta_mixture_raw_ml.R" "${OUT_DIR}"
