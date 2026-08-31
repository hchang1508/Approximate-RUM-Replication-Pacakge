#!/bin/bash
# Greedy algorithm for tbl:apx_error2 (DVD, fixed-effects grid).
#
#   cd RUM_replication/dvd/tbl_apx_error2 && bash code/submit_greedy.sh <case> <degree>
#   ... code/submit_greedy.sh all 1          # all six cases at degree 1
#
# One array task per entry of to_do_list_GREEDY/<case>_d<degree>_case_to_do.csv, so the
# array addresses remaining work rather than grid points directly. Rebuild the lists with
# code/make_todo_list_GREEDY.R before submitting, and not while an array is still queued.
#
# Writes output/raw/DVD_hh_gd_<case>_<degree>_<gridindex>.txt. Aggregate with code/collect.R.

set -euo pipefail

CASE="${1:?usage: submit_greedy.sh <case|all> <degree>}"
DEGREE="${2:?usage: submit_greedy.sh <case|all> <degree>}"
case "$DEGREE" in 1|2) ;; *) echo "degree must be 1 or 2" >&2; exit 1 ;; esac
case "$CASE" in all) CASES="1 2 3 4 5 6" ;; [1-6]) CASES="$CASE" ;;
  *) echo "case must be 1-6 or all" >&2; exit 1 ;; esac

cd /home/hc654/RUM_replication/dvd/tbl_apx_error2

TODO_DIR="${TODO_DIR:-to_do_list_GREEDY}"
PYTHON="${PYTHON:-/gpfs/gibbs/project/narita/hc654/conda_envs/approx/bin/python}"
PARTITION="${PARTITION:-scavenge}"
DRY_RUN="${DRY_RUN:-0}"

export DVD_DATA="${DVD_DATA:-../../data/dvd/dvd_J4_K2.csv}"
export TODO_DIR
export OUT_DIR="${OUT_DIR:-output/raw}"
mkdir -p "$OUT_DIR" output/logs

for c in $CASES; do
  todo=$TODO_DIR/${c}_d${DEGREE}_case_to_do.csv
  [ -f "$todo" ] || { echo "missing $todo -- run code/make_todo_list_GREEDY.R first" >&2; exit 1; }
  n=$(wc -l < "$todo")
  [ "$n" -gt 0 ] || { echo "case $c d$DEGREE: nothing to do"; continue; }

  CMD=(sbatch
    --job-name=T2gd_c${c}d${DEGREE}
    --partition=${PARTITION}
    --time=1-00:00:00
    --ntasks=1 --cpus-per-task=1 --nodes=1
    --array=0-$((n - 1))
    --requeue
    --output=output/logs/gd_c${c}d${DEGREE}_%A_%a.out
    --error=output/logs/gd_c${c}d${DEGREE}_%A_%a.err
    --export=ALL,DVD_DATA=${DVD_DATA},TODO_DIR=${TODO_DIR},OUT_DIR=${OUT_DIR}
    --wrap="$PYTHON -u code/approx_rum_GREEDY_hhmix_DVD.py \$SLURM_ARRAY_TASK_ID ${c} ${DEGREE}"
  )
  if [ "$DRY_RUN" = "1" ]; then printf '%q ' "${CMD[@]}"; echo
  else "${CMD[@]}"; fi
  echo "case $c d$DEGREE: $n tasks -> $OUT_DIR"
done
