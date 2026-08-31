#!/usr/bin/env Rscript
# Aggregate the per-task greedy and EM results into tbl:apx_error1.
#
#   Rscript code/collect.R [raw_dir] [out_csv]
#
# raw_dir is either a single directory holding DVD_gd_<order>_<index>.txt and
# DVD_em_<order>_<index>.txt, or a parent holding greedy/ and em_tight/. Each file
# holds five numbers: pi(0), pi(1), pi(2), pi(3), raw_error.
#
# Errors are rescaled by count/sqrt(|D|) = sqrt(11); see README.md, "Normalisation".

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else "output/raw"
out_csv <- if (length(args) >= 2) args[2] else "output/tbl_apx_error1.csv"

N_ALT  <- 4
N_RANK <- factorial(N_ALT)
COUNT  <- choose(N_ALT, 2) + choose(N_ALT, 3) + 1
CARD_D <- COUNT
FACTOR <- COUNT / sqrt(CARD_D)

split_layout <- dir.exists(file.path(raw_dir, "greedy")) &&
                dir.exists(file.path(raw_dir, "em_tight"))
alg_dir <- function(alg) {
  if (!split_layout) return(raw_dir)
  file.path(raw_dir, if (alg == "gd") "greedy" else "em_tight")
}

cat(sprintf("layout        %s\n",
            if (split_layout) "split (greedy/ + em_tight/)" else "flat (single raw dir)"))
cat(sprintf("code divisor  count = |cal_D| = %d  (pairs %d + triples %d + full set 1)\n",
            COUNT, choose(N_ALT, 2), choose(N_ALT, 3)))
cat(sprintf("paper norm    |D| = count = %d\n", CARD_D))
cat(sprintf("scale         d = raw * count/sqrt(|D|) = raw * %.4f\n", FACTOR))

read_one <- function(alg, order, index) {
  f <- file.path(alg_dir(alg), sprintf("DVD_%s_%d_%d.txt", alg, order, index))
  if (!file.exists(f)) return(NULL)
  v <- scan(f, quiet = TRUE)
  if (length(v) != N_ALT + 1) {
    warning(sprintf("%s: expected %d numbers, got %d", f, N_ALT + 1, length(v)))
    return(NULL)
  }
  list(ranking = v[1:N_ALT], raw = v[N_ALT + 1])
}

lab <- function(r) paste(r + 1, collapse = ">")

rows <- list()
missing <- character(0)
for (index in 0:(N_RANK - 1)) {
  cells <- list()
  ranking <- NULL
  for (order in 1:2) {
    for (alg in c("gd", "em")) {
      got <- read_one(alg, order, index)
      key <- paste0(ifelse(order == 1, "linear", "quadratic"), "_",
                    ifelse(alg == "gd", "greedy", "em"))
      if (is.null(got)) {
        missing <- c(missing, sprintf("%s order=%d index=%d", alg, order, index))
        cells[[key]] <- NA_real_
      } else {
        cells[[key]] <- got$raw
        if (is.null(ranking)) ranking <- got$ranking
      }
    }
  }
  if (is.null(ranking)) next
  rows[[length(rows) + 1]] <- c(list(ranking = lab(ranking)), cells)
}

if (!length(rows)) stop("no result files found in ", raw_dir,
                        " -- submit code/submit_greedy.sh and code/submit_em.sh first")

tab <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
val_cols <- c("linear_greedy", "linear_em", "quadratic_greedy", "quadratic_em")

tab_scaled <- tab
tab_scaled[val_cols] <- tab[val_cols] * FACTOR

tab_scaled$representable <- as.integer(
  apply(tab_scaled[c("linear_greedy", "linear_em")], 1,
        function(z) all(is.na(z) | z < 1e-3)))

ord <- order(tab_scaled$representable, -tab_scaled$linear_greedy)
tab_scaled <- tab_scaled[ord, c("ranking", "representable", val_cols)]

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(format(tab_scaled, digits = 4, trim = TRUE), out_csv,
          row.names = FALSE, quote = FALSE)

cat(sprintf("\ncollected %d of %d rankings", nrow(tab_scaled), N_RANK))
if (length(missing)) cat(sprintf("  (%d cells missing)", length(missing)))
cat(sprintf("\nscale applied  : count/sqrt(|D|) (x%.4f)\n", FACTOR))

cat("\n=== Table 1 ===\n")
print(tab_scaled, row.names = FALSE, digits = 4)

cat(sprintf("\nwrote %s\n", out_csv))
