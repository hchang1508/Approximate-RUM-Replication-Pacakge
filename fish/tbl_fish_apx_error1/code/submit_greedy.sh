#!/bin/bash
# Greedy algorithm for tbl:fish_apx_error1 (fish, no fixed effects).
#
#   cd RUM_replication/fish/tbl_fish_apx_error1 && sbatch code/submit_greedy.sh
#
# One array task per ranking; each runs order=1 (linear) and order=2 (quadratic).
# Writes output/raw/Fish_gd_<order>_<index>.txt. Aggregate with code/collect.R.
#
#SBATCH --job-name=F1_greedy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --array=0-23
#SBATCH --requeue
#SBATCH --output=output/logs/greedy_%A_%a.out
#SBATCH --error=output/logs/greedy_%A_%a.err

set -euo pipefail

cd /home/hc654/RUM_replication/fish/tbl_fish_apx_error1
mkdir -p output/raw output/logs

PYTHON="${PYTHON:-/gpfs/gibbs/project/narita/hc654/conda_envs/approx/bin/python}"

export OUT_DIR=output/raw
export FISH_DATA="${FISH_DATA:-../../data/fish/fish_J4_K2.csv}"
export SEED_BASE="${SEED_BASE:-20240808}"

IDX="${SLURM_ARRAY_TASK_ID:-0}"
echo "host        : $(hostname)"
echo "ranking idx : $IDX"
echo "data        : $FISH_DATA"
echo "seed base   : $SEED_BASE"
echo "python      : $PYTHON"

for ORDER in 1 2; do
  echo "=== greedy, ranking $IDX, order $ORDER ==="
  "$PYTHON" -u code/approx_rum_GREEDY_vertices_fish.py "$IDX" "$ORDER" "$FISH_DATA"
done

echo "done"
