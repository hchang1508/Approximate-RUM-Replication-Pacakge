#!/bin/bash
# tbl:fish_subs -- maximal substitution, fishing data, greedy.
#   cd RUM_replication/fish/tbl_fish_subs && sbatch code/submit_substitute.sh
# One array task per off-diagonal cell (case 1..12). Args: case dim fe.
#SBATCH --job-name=F_subs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --array=1-12
#SBATCH --requeue
#SBATCH --output=output/logs/subs_%A_%a.out
#SBATCH --error=output/logs/subs_%A_%a.err

module load miniconda
conda deactivate
conda activate approx
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

cd /home/hc654/RUM_replication/fish/tbl_fish_subs
export FISH_DATA=../../data/fish/fish_J4_K2.csv
export OUT_DIR=output/raw
mkdir -p output/raw output/logs

python -u code/greedy_substitute.py $SLURM_ARRAY_TASK_ID 1 0
