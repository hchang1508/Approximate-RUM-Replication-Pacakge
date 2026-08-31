#!/bin/bash
# Half-half EM at a tightened tolerance. One replicate = one random start over the full
# 9261-point grid; run several replicates to raise the start count.
#
#   bash submit_em_hhmix.sh <dvd|fish> <degree> <replicate>
#
# <replicate> is either a single number or an inclusive range:
#   ... fish 2 3      replicate 3 only
#   ... fish 2 2-10   replicates 2 through 10
#
# Output goes to output/raw_em_tight/rep<N>/, one directory per replicate, so the
# collector can take the minimum over exactly the replicates it is pointed at.

set -euo pipefail

DS="${1:?usage: submit_em_hhmix.sh <dvd|fish> <degree> <replicate>}"
DEGREE="${2:?usage: submit_em_hhmix.sh <dvd|fish> <degree> <replicate>}"
REPSPEC="${3:?usage: submit_em_hhmix.sh <dvd|fish> <degree> <replicate>}"

case "$DS" in
  dvd)  EX=dvd/tbl_apx_error2;        PY=approx_rum_EM_hhmix_DVD.py;  TAG=T2emT ;;
  fish) EX=fish/tbl_fish_apx_error2;  PY=approx_rum_EM_hhmix_fish.py; TAG=F2emT ;;
  *)    echo "dataset must be dvd or fish" >&2; exit 1 ;;
esac
case "$DEGREE" in
  1) WALL=2:00:00 ;;   # ~21 s/fit, ~0.15 h per chunk
  2) WALL=6:00:00 ;;   # ~250 s/fit, ~1.74 h per chunk
  *) echo "degree must be 1 or 2" >&2; exit 1 ;;
esac

# <replicate> is either "N" or "A-B". Reject anything else before splitting: a bare "-",
# a dangling end ("2-"), or a second dash ("1-2-3") would otherwise split into nonsense.
case "$REPSPEC" in
  ''|*[!0-9-]*|-*|*-|*-*-*) echo "replicate must be a number or a range like 2-10" >&2; exit 1 ;;
esac
REP_FIRST="${REPSPEC%%-*}"
REP_LAST="${REPSPEC##*-}"
if [ "$REP_FIRST" -lt 1 ] || [ "$REP_LAST" -lt "$REP_FIRST" ]; then
  echo "replicate range must start at 1 or more and not run backwards" >&2; exit 1
fi

# Overridable so the file runs unedited from a checkout elsewhere.
REPO="${RUM_REPO:-/home/hc654/RUM_replication}"
SCRATCH="${RUM_SCRATCH:-/vast/palmer/scratch/narita/hc654/RUM_replication}"
PYTHON="${RUM_PYTHON:-/gpfs/gibbs/project/narita/hc654/conda_envs/approx/bin/python}"

CHUNK=25
NCHUNK=371                      # ceil(9261/25); the last chunk holds 11 points

for REP in $(seq "$REP_FIRST" "$REP_LAST"); do
  SEED_BASE=$(( 20240808 + 1000000 * REP ))
  OUT_DIR=$SCRATCH/$EX/output/raw_em_tight/rep$REP
  LOG_ROOT=$SCRATCH/$EX/output/logs_em_tight/rep$REP

  # A replicate that already holds results is being recomputed, not extended: same seed,
  # same filenames. Warn rather than skip, so the choice stays with the caller.
  if [ -d "$OUT_DIR" ] && [ -n "$(ls -A "$OUT_DIR" 2>/dev/null)" ]; then
    echo "WARNING: $OUT_DIR already holds $(ls "$OUT_DIR" | wc -l) files -- resubmitting overwrites them" >&2
  fi
  mkdir -p "$OUT_DIR"

  for CASE in 1 2 3 4 5 6; do
    LOGDIR=$LOG_ROOT/c${CASE}d${DEGREE}
    mkdir -p "$LOGDIR"          # SLURM opens --output before the job script runs

    sbatch \
      --job-name=${TAG}_r${REP}_c${CASE}d${DEGREE} \
      --partition=day \
      --time=${WALL} \
      --ntasks=1 --cpus-per-task=1 --nodes=1 --mem=5G \
      --array=1-${NCHUNK} \
      --requeue \
      --output=${LOGDIR}/em_%A_%a.out \
      --error=${LOGDIR}/em_%A_%a.err \
      --export=ALL,OUT_DIR=${OUT_DIR},EM_STARTS=1,EM_TOL=1e-8,EM_STANDARDIZE=1,EM_CHUNK=${CHUNK},SEED_BASE=${SEED_BASE} \
      --wrap="cd $REPO/$EX && $PYTHON -u code/$PY \$SLURM_ARRAY_TASK_ID ${CASE} ${DEGREE}"
  done

  echo "$DS d$DEGREE replicate $REP: SEED_BASE=$SEED_BASE  chunk=$CHUNK  wall=$WALL  arrays=6 x $NCHUNK tasks  out=$OUT_DIR"
done
