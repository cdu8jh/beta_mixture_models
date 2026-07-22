#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH -c 4
#SBATCH --job-name=beta_mix_raw_ml
#SBATCH -o logs/02_fit_beta_mixture_raw_ml_%j.txt
#SBATCH -e logs/02_fit_beta_mixture_raw_ml_%j.txt
#SBATCH --time=24:00:00

BATCH="${1:-m6A_2OmeA_I_genome}"

SCRIPT_DIR=/dcs04/hicks/data/sparthib/drna_retinalrs/code/04_modkit_stats/01_ml_probability

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Batch: ${BATCH}"
echo "****"

mkdir -p logs

module load conda_R/4.4

Rscript "${SCRIPT_DIR}/02_fit_beta_mixture_raw_ml.R" \
    "${SCRIPT_DIR}/config_${BATCH}.R"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"