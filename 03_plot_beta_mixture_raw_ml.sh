#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH --job-name=beta_mix_plot
#SBATCH -o logs/03_plot_beta_mixture_raw_ml_%j.txt
#SBATCH -e logs/03_plot_beta_mixture_raw_ml_%j.txt
#SBATCH --time=01:00:00

BATCH="${1:-m6A_2OmeA_I_genome}"

SCRIPT_DIR=/dcs04/hicks/data/sparthib/drna_retinalrs/code/04_modkit_stats/01_ml_probability

OUT_DIR=/dcs04/hicks/data/sparthib/drna_retinalrs/processed_data/04_pileup/01_ml_probability/${BATCH}/02_beta_mixture_raw_ml

mkdir -p logs

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "Batch: ${BATCH}"
echo "Output dir: ${OUT_DIR}"

module load conda_R/4.4

Rscript \
    "${SCRIPT_DIR}/03_plot_beta_mixture_raw_ml.R" \
    "${OUT_DIR}"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"