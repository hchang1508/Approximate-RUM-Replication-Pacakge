#!/bin/bash
# EM at parity for the two apx_error1 vertices exhibits. Results and logs go to scratch;
# this script lives in the repo.
#
#   sbatch --export=ALL,DATASET=fish code/submit_em_vertices.sh
#   sbatch --export=ALL,DATASET=dvd  code/submit_em_vertices.sh
#
# OUT_DIR and --output both resolve on scratch, the log directory is created before
# submission, and python is the conda interpreter by absolute path.
#
# An identical copy lives in dvd/tbl_apx_error1/code/. DATASET selects the exhibit, so
# either copy can run either dataset; keep the two files in sync.
#
# One array task per ranking; each runs order=1 (linear) and order=2 (quadratic).
# Writes $OUT_DIR/<PREFIX>_em_<order>_<index>.txt. Aggregate with the exhibit's collect.R.
#
#SBATCH --job-name=emvert
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --array=0-23
#SBATCH --requeue
#SBATCH --output=/vast/palmer/scratch/narita/hc654/RUM_replication/logs_em_vertices/em_%A_%a.out
#SBATCH --error=/vast/palmer/scratch/narita/hc654/RUM_replication/logs_em_vertices/em_%A_%a.err

set -euo pipefail

DATASET="${DATASET:?set DATASET=fish or DATASET=dvd via --export}"
HOME_REPL=/home/hc654/RUM_replication
SCRATCH_REPL=/vast/palmer/scratch/narita/hc654/RUM_replication

case "$DATASET" in
  fish)
    EXHIBIT=fish/tbl_fish_apx_error1
    SCRIPT=code/approx_rum_EM_vertices_fish.py
    export FISH_DATA="${FISH_DATA:-../../data/fish/fish_J4_K2.csv}"
    DATA="$FISH_DATA"
    ;;
  dvd)
    EXHIBIT=dvd/tbl_apx_error1
    SCRIPT=code/approx_rum_EM_vertices_DVD.py
    export DVD_DATA="${DVD_DATA:-../../data/dvd/dvd_J4_K2.csv}"
    DATA="$DVD_DATA"
    ;;
  *) echo "DATASET must be fish or dvd, got '$DATASET'" >&2; exit 2 ;;
esac

cd "$HOME_REPL/$EXHIBIT"

export OUT_DIR="${OUT_DIR:-$SCRATCH_REPL/$EXHIBIT/output/raw_em_tight}"
mkdir -p "$OUT_DIR"

export SEED_BASE="${SEED_BASE:-20240808}"
export EM_STARTS="${EM_STARTS:-10}"
export EM_MIXTURE="${EM_MIXTURE:-18}"
export EM_MAX_ITER="${EM_MAX_ITER:-1000}"
export EM_TOL="${EM_TOL:-1e-8}"
export EM_STANDARDIZE="${EM_STANDARDIZE:-1}"

# Adds a floor to every mixture weight, keeping zero-weight components alive
# (approx_rum_EM_vertices_fish.py:61). At 0 the fish run returns NaN for rankings 3
# (1>3>4>2) and 11 (2>4>3>1) -- a degenerate component at M=18.
#
# WARNING: the committed fish table was produced by job 60161730 with a NONZERO value that
# was passed on the sbatch command line and is recorded nowhere -- not in this script, not
# in any log, and the batch script is no longer retrievable from SLURM. Running as-is
# reproduces the NaNs. Pin the production value here before relying on this script.
export EM_LAMBDA_FLOOR="${EM_LAMBDA_FLOOR:-0}"

PYTHON="${PYTHON:-/gpfs/gibbs/project/narita/hc654/conda_envs/approx/bin/python}"

IDX="${SLURM_ARRAY_TASK_ID:-0}"
echo "host        : $(hostname)"
echo "dataset     : $DATASET"
echo "exhibit     : $EXHIBIT"
echo "ranking idx : $IDX"
echo "data        : $DATA"
echo "out dir     : $OUT_DIR"
echo "seed base   : $SEED_BASE"
echo "EM starts   : $EM_STARTS"
echo "mixtures    : $EM_MIXTURE"
echo "max iter    : $EM_MAX_ITER"
echo "tolerance   : $EM_TOL"
echo "standardise : $EM_STANDARDIZE"
echo "lambda floor: $EM_LAMBDA_FLOOR"
echo "python      : $PYTHON"

for ORDER in ${ORDERS:-1 2}; do
  echo "=== EM, ranking $IDX, order $ORDER ==="
  "$PYTHON" -u "$SCRIPT" "$IDX" "$ORDER" "$DATA"
done

echo "done"
