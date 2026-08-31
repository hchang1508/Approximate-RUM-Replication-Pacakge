#!/usr/bin/env Rscript
# Aggregate the per-task results of submit_greedy.sh and draw Figure A.3.
#
#   Rscript code/collect_and_plot.R        # run from the exhibit directory
#
# Input : output/raw/Large_DVD_gd_1_<index>.txt, each holding 11 numbers:
#         pi(0..9) then raw_error
# Output: output/fig_approx_error_J10.csv   one row per draw
#         output/approximation_error.png    the histogram
#
# Rescaling. The Python script reports raw = ||rho - rho_hat|| / count, dividing by the
# NUMBER of choice sets rather than its square root. The paper's metric is
# d = ||rho - rho_hat|| / sqrt(|D|), so in general
#         d = raw * count / sqrt(|D|).
#
# Singleton menus are never evaluated, so D is exactly the family the code enumerates,
#         count = |D| = choose(10,2) + choose(10,3) + 1 = 166,
# and count/sqrt(|D|) reduces to sqrt(166). The original extractor's sqrt(cardD) is
# therefore correct for this figure.
#
# Same convention as the |J|=4 tables, which set |D| = count = 11 and rescale by sqrt(11);
# only the menu count differs. See ../../tbl_apx_error1/code/collect.R.

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else "output/raw"
out_csv <- if (length(args) >= 2) args[2] else "output/fig_approx_error_J10.csv"
out_png <- if (length(args) >= 3) args[3] else "output/approximation_error.png"

N_ALT <- 10
N_DRAW <- as.integer(Sys.getenv("N_DRAW", "100"))
COUNT  <- choose(N_ALT, 2) + choose(N_ALT, 3) + 1        # 166 = the code's divisor
CARD_D <- COUNT                                           # 166 = |D| for this exercise
SCALE  <- COUNT / sqrt(CARD_D)                            # = sqrt(166)

cat(sprintf("count = |D| = %d  (pairs %d + triples %d + full set 1)\n",
            COUNT, choose(N_ALT, 2), choose(N_ALT, 3)))
cat(sprintf("scale : d = raw * count/sqrt(|D|) = raw * %.4f\n", SCALE))

recs <- list()
missing <- integer(0)
for (index in 0:(N_DRAW - 1)) {
  f <- file.path(raw_dir, sprintf("Large_DVD_gd_1_%d.txt", index))
  if (!file.exists(f)) { missing <- c(missing, index); next }
  v <- scan(f, quiet = TRUE)
  if (length(v) != N_ALT + 1) {
    warning(sprintf("%s: expected %d numbers, got %d", f, N_ALT + 1, length(v)))
    next
  }
  recs[[length(recs) + 1]] <- data.frame(
    draw    = index,
    ranking = paste(v[1:N_ALT] + 1, collapse = ">"),
    error   = v[N_ALT + 1] * SCALE,
    stringsAsFactors = FALSE
  )
}

if (!length(recs)) stop("no result files found in ", raw_dir,
                        " -- submit code/submit_greedy.sh first")

tab <- do.call(rbind, recs)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(tab, out_csv, row.names = FALSE, quote = FALSE)

cat(sprintf("\ncollected %d of %d draws", nrow(tab), N_DRAW))
if (length(missing)) {
  cat(sprintf("  -- MISSING %d: %s", length(missing),
              paste(head(missing, 20), collapse = ",")))
  if (length(missing) > 20) cat(" ...")
}
cat("\n\n=== approximation errors ===\n")
print(summary(tab$error))
cat(sprintf("\nall errors non-negligible (> 0.1)? %s\n", all(tab$error > 0.1)))
cat(sprintf("range: %.3f to %.3f\n", min(tab$error), max(tab$error)))
cat("(the appendix text describes a range of 0.3 to 0.7)\n")

png(out_png, width = 1600, height = 1200, res = 220)
hist(tab$error,
     breaks = 20,
     main = sprintf("Approximation errors for %d randomly generated rankings", nrow(tab)),
     xlab = "Approximation error",
     ylab = "Count",
     col = "grey80", border = "white")
box()
invisible(dev.off())

cat(sprintf("\nwrote %s\n", out_csv))
cat(sprintf("wrote %s\n", out_png))
