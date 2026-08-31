#!/usr/bin/env Rscript
# Rebuild the greedy to-do lists from whatever has already finished.
#
#   Rscript code/make_todo_list_GREEDY.R          # run from the exhibit directory
#
# Diffs the set of completed output files against the full grid and writes the remainder,
# so a resubmission picks up only what is missing. Paths resolve relative to this exhibit.
#
# Input : output/raw/DVD_hh_gd_<case>_<d>_<gridindex>.txt   (one file per finished cell)
# Output: to_do_list_GREEDY/<case>_d<d>_case_to_do.csv       (one grid index per line)
#
# The greedy script reads its to-do list and does `index = int(to_do[index])`, so a SLURM
# array index addresses a position in the REMAINING work, not a grid point directly.
# That is what makes the resubmit loop work: rerun this script, then submit an array
# sized to the new list. code/resubmit_greedy.sh does both in one step.

RAW_DIR  <- Sys.getenv("RAW_DIR", "output/raw")
TODO_DIR <- Sys.getenv("TODO_DIR", "to_do_list_GREEDY")
CASES    <- as.integer(strsplit(Sys.getenv("CASES", "1,2,3,4,5,6"), ",")[[1]])
DEGREES  <- as.integer(strsplit(Sys.getenv("DEGREES", "1,2"), ",")[[1]])
GRID_MAX <- as.integer(Sys.getenv("GRID_MAX", "9260"))     # 21^3 - 1

dir.create(TODO_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(RAW_DIR, pattern = "^DVD_hh_gd_[0-9]+_[0-9]+_[0-9]+\\.txt$")
cat(sprintf("raw dir  : %s   (%d finished files)\n", RAW_DIR, length(files)))
cat(sprintf("todo dir : %s\n", TODO_DIR))
cat(sprintf("grid     : 0..%d  (%d points)\n\n", GRID_MAX, GRID_MAX + 1))

# DVD_hh_gd_<case>_<d>_<gridindex>.txt
parts <- regmatches(files, regexec("^DVD_hh_gd_([0-9]+)_([0-9]+)_([0-9]+)\\.txt$", files))
done <- data.frame(
  case = as.integer(sapply(parts, `[`, 2)),
  d    = as.integer(sapply(parts, `[`, 3)),
  grid = as.integer(sapply(parts, `[`, 4))
)

total_left <- 0L
for (d in DEGREES) {
  for (case in CASES) {
    finished <- sort(unique(done$grid[done$case == case & done$d == d]))
    todo <- setdiff(0:GRID_MAX, finished)
    f <- file.path(TODO_DIR, sprintf("%d_d%d_case_to_do.csv", case, d))
    write.table(todo, f, sep = ",", col.names = FALSE, row.names = FALSE)
    total_left <- total_left + length(todo)
    cat(sprintf("  case %d d%d : %5d done, %5d to do  -> %s\n",
                case, d, length(finished), length(todo), f))
  }
}
cat(sprintf("\n%d cells remaining across %d case/degree combinations\n",
            total_left, length(CASES) * length(DEGREES)))
if (total_left == 0L) cat("nothing left -- aggregate with code/collect.R\n")
