#!/bin/bash
# Table 3 (tbl:subs) -- maximal substitution, DVD, greedy algorithm.
#   cd RUM_replication/dvd/tbl_subs && sbatch code/submit_substitute.sh
#
# One array task per off-diagonal cell of the 4x4 table: case 1..12 selects the
# (chosen, removed) pair, so 12 tasks fill the whole table.
# Arguments are  case  dim  fe :
#   dim = 1  linear mixed logit (what the table reports)
#   fe  = 0  no fixed effects  (with free fixed effects the quantity is always 1)
# Aggregate with:  Rscript code/collect.R
#
#SBATCH --job-name=DVD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --array=1-12
#SBATCH --output=output/logs/subs_%A_%a.out
#SBATCH --error=output/logs/subs_%A_%a.err

module load miniconda
conda deactivate
conda activate approx
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

cd /home/hc654/RUM_replication/dvd/tbl_subs

# where the patched script reads characteristics and writes results
export DVD_DATA=../../data/dvd/dvd_J4_K2.csv
export OUT_DIR=output/raw
mkdir -p output/raw output/logs

python -u code/greedy_substitute.py $SLURM_ARRAY_TASK_ID 1 0
