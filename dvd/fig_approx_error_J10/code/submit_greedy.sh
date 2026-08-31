#!/bin/bash
# Figure A.3 (fig:approx_error_J10) -- approximation errors for |J|=10, K=3.
#
#   cd RUM_replication/dvd/fig_approx_error_J10 && sbatch code/submit_greedy.sh
#
# One array task per randomly drawn ranking, 100 in total, matching the figure's caption
# ("100 randomly generated rankings"). Task i draws its ranking with np.random.seed(i),
# so the 100 rankings are distinct and the whole figure is reproducible.
#
# Only order=1 (linear, K=3) is needed: the figure reports linear mixed-logit errors.
# Raw results land in output/raw/Large_DVD_gd_1_<index>.txt, each holding
# [pi(0..9), approximation_error].
# Aggregate and plot with:  Rscript code/collect_and_plot.R
#
# This is the heaviest of the DVD exhibits: |cal_D| = C(10,2)+C(10,3)+1 = 166 menus, so
# each task is far slower than an |J|=4 task. The original used partition=week with a
# 7-day limit; that is kept.
#
#SBATCH --job-name=FA3_greedy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --array=0-99
#SBATCH --output=output/logs/greedy_%A_%a.out
#SBATCH --error=output/logs/greedy_%A_%a.err

set -euo pipefail

cd /home/hc654/RUM_replication/dvd/fig_approx_error_J10
mkdir -p output/raw output/logs

PYTHON="${PYTHON:-/gpfs/gibbs/project/narita/hc654/conda_envs/approx/bin/python}"

export OUT_DIR=output/raw
export DVD_DATA="${DVD_DATA:-../../data/dvd/dvd_J10_K3.csv}"

IDX="${SLURM_ARRAY_TASK_ID:-0}"
echo "host        : $(hostname)"
echo "draw idx    : $IDX   (ranking drawn with np.random.seed($IDX))"
echo "data        : $DVD_DATA"
echo "python      : $PYTHON"

echo "=== greedy, |J|=10, draw $IDX, order 1 ==="
"$PYTHON" -u code/approx_rum_GREEDY_vertices_DVD_large.py "$IDX" 1 "$DVD_DATA"

echo "done"
