#!/usr/bin/env Rscript
# Aggregate the per-grid-point greedy and EM results into tbl:apx_error2.
#
#   Rscript code/collect.R [raw_dir] [out_csv]
#
# raw_dir holds greedy/ and em_tight/rep<N>/, as in
# $SCRATCH/RUM_replication/results_final/dvd/tbl_apx_error2. The exhibit's own
# output/{raw,raw_em_tight} layout is accepted too.
#
# Each cell is the minimum over the 9261 fixed-effects grid points, and for EM over the
# replicate directories as well. Errors are rescaled by sqrt(11); printed values are
# truncated at the third decimal, not rounded. See README.md.

args    <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else "output/raw"
out_csv <- if (length(args) >= 2) args[2] else "output/tbl_apx_error2.csv"

N_ALT      <- 4
COUNT      <- choose(N_ALT, 2) + choose(N_ALT, 3) + 1
FACTOR     <- COUNT / sqrt(COUNT)
GRID       <- 21^3
CASES      <- 1:6
DEGREES    <- 1:2

PREF <- list(
  c("1>3>4>2", "2>4>3>1"),
  c("1>3>2>4", "4>2>3>1"),
  c("1>4>3>2", "2>3>4>1"),
  c("2>4>1>3", "3>1>4>2"),
  c("2>1>4>3", "3>4>1>2"),
  c("3>1>2>4", "4>2>1>3")
)

pick <- function(...) { for (p in c(...)) if (dir.exists(p)) return(p); NA_character_ }
greedy_dir <- pick(file.path(raw_dir, "greedy"), file.path(raw_dir, "raw"), raw_dir)
em_root    <- pick(file.path(raw_dir, "em_tight"), file.path(raw_dir, "raw_em_tight"))
em_reps    <- if (is.na(em_root)) character(0) else
              sort(list.dirs(em_root, recursive = FALSE))

if (is.na(greedy_dir)) stop("no greedy directory under ", raw_dir)
if (!length(em_reps))  stop("no em_tight replicate directories under ", raw_dir)

scan_files <- function(dir, glob, prog, fs = NULL) {
  cmd <- sprintf("find %s -maxdepth 1 -name %s -exec awk %s%s {} +",
                 shQuote(dir), shQuote(glob),
                 if (is.null(fs)) "" else paste0("-F", shQuote(fs), " "),
                 shQuote(prog))
  v <- suppressWarnings(as.numeric(system(cmd, intern = TRUE, ignore.stderr = TRUE)))
  if (!length(v)) return(numeric(0))
  v[!is.na(v)]
}

read_greedy <- function(case, degree) {
  scan_files(greedy_dir, sprintf("DVD_hh_gd_%d_%d_*.txt", case, degree), "FNR==5")
}

read_em <- function(case, degree) {
  unlist(lapply(em_reps, function(rp)
    scan_files(rp, sprintf("hh_DVD_%dd%d_EM_error_*.csv", case, degree),
               "FNR>1{print $2}", fs = ",")), use.names = FALSE)
}

cat(sprintf("greedy dir    %s\n", greedy_dir))
cat(sprintf("em replicates %d  (%s)\n", length(em_reps),
            paste(basename(em_reps), collapse = " ")))
cat(sprintf("scale         d = raw * count/sqrt(|D|) = raw * %.4f\n", FACTOR))
cat("EM config     tight: EM_STARTS=1 per replicate, EM_TOL=1e-8, EM_STANDARDIZE=1\n\n")

rows <- list()
cov  <- list()
for (case in CASES) {
  cell <- list(case = case, ranking = PREF[[case]][1], ranking_2 = PREF[[case]][2])
  for (degree in DEGREES) {
    g <- read_greedy(case, degree)
    e <- read_em(case, degree)
    key <- if (degree == 1) "linear" else "quadratic"
    cell[[paste0(key, "_greedy")]] <- if (length(g)) min(g) * FACTOR else NA_real_
    cell[[paste0(key, "_em")]]     <- if (length(e)) min(e) * FACTOR else NA_real_
    cov[[length(cov) + 1]] <- data.frame(
      case = case, degree = degree,
      greedy_pts = length(g), greedy_pct = 100 * length(g) / GRID,
      em_rows = length(e), em_pct = 100 * length(e) / (GRID * length(em_reps)))
  }
  rows[[length(rows) + 1]] <- as.data.frame(cell, stringsAsFactors = FALSE)
}

tab <- do.call(rbind, rows)
val_cols <- c("linear_greedy", "linear_em", "quadratic_greedy", "quadratic_em")

cov <- do.call(rbind, cov)
cat("=== coverage ===\n")
print(cov, row.names = FALSE, digits = 4)
if (any(cov$greedy_pts < GRID))
  cat(sprintf("\nnote: %d of 12 greedy cells are short of the %d-point grid\n",
              sum(cov$greedy_pts < GRID), GRID))

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(tab, out_csv, row.names = FALSE, quote = FALSE)

trunc3 <- function(x) floor(x * 1000) / 1000
printed <- tab
printed[val_cols] <- lapply(tab[val_cols], function(x) sprintf("%.3f", trunc3(x)))

cat("\n=== Table 2 (truncated at three decimals) ===\n")
print(printed[c("case", "ranking", "ranking_2", val_cols)], row.names = FALSE)

cat("\n=== full precision ===\n")
print(tab[c("case", "ranking", val_cols)], row.names = FALSE, digits = 4)

cat(sprintf("\nwrote %s\n", out_csv))
