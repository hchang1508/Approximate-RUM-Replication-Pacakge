#!/bin/bash
# tbl:in-sample_d1 / tbl:out-sample -- stage 2: the "our method" row, 50 random splits.
#   cd RUM_replication/fish/tbl_in_sample_d1 && sbatch code/submit_em.sh
# Consumes output/em_input/, which code/submit_benchmarks.sh produces. Submit only
# after that array has finished.
#SBATCH --job-name=inout_em
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=day
#SBATCH --time=04:00:00
#SBATCH --array=1-50
#SBATCH --requeue
#SBATCH --output=output/logs/em_%A_%a.out
#SBATCH --error=output/logs/em_%A_%a.err

module load miniconda
conda deactivate
conda activate approx
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

cd /home/hc654/RUM_replication/fish/tbl_in_sample_d1
export EM_INPUT_DIR=output/em_input
export OUT_DIR=output/raw
mkdir -p output/raw output/logs

python -u code/approx_rum_EM_fish_inout.py $SLURM_ARRAY_TASK_ID
