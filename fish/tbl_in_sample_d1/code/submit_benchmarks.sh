#!/bin/bash
# tbl:in-sample_d1 / tbl:out-sample -- stage 1: the six benchmark models, 50 random splits.
#   cd RUM_replication/fish/tbl_in_sample_d1 && sbatch code/submit_benchmarks.sh
# One array task per split. Each also writes that split's EM input files, so this
# must finish before code/submit_em.sh is submitted.
#SBATCH --job-name=inout_bm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=day
#SBATCH --time=12:00:00
#SBATCH --array=1-50
#SBATCH --requeue
#SBATCH --output=output/logs/bm_%A_%a.out
#SBATCH --error=output/logs/bm_%A_%a.err

module load R/4.4.1-foss-2022b
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

cd /home/hc654/RUM_replication/fish/tbl_in_sample_d1
export FISH_DATA=../../data/fish/fish_individual.csv
export OUT_DIR=output
mkdir -p output/raw output/logs output/splits output/em_input

Rscript --vanilla code/in_out_fit.R $SLURM_ARRAY_TASK_ID
