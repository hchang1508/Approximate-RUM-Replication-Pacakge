#!/usr/bin/env Rscript
# Aggregate the 12 per-cell results of submit_substitute.sh into Table 3.
#
#   Rscript code/collect.R                 # run from the exhibit directory
#
# Input : output/raw/subs_d<dim>_fe<fe>_<chosen>_<removed>.txt
#         one file per off-diagonal cell; np.savetxt of greedy_result[1][-1:], i.e. the
#         final value of  1 - sqrt(subs_diff)  from the greedy iteration.
# Output: output/tbl_subs.csv          the 4x4 table, rows = dropped j, cols = chosen l
#
# No rescaling. Unlike the approximation-error tables, this quantity is already a
# difference of choice probabilities in [0,1] -- see eq (eq:dist.2) -- so the |D|
# normalisation that applies to the distance metric is not relevant here.

DIM <- as.integer(Sys.getenv("DIM", "1"))     # 1 = linear mixed logit, as the table reports
FE  <- as.integer(Sys.getenv("FE", "0"))      # 0 = no fixed effects
raw_dir <- Sys.getenv("RAW_DIR", "output/raw")
out_csv <- Sys.getenv("OUT_CSV", "output/tbl_subs.csv")

N_ALT <- 4
# case -> (chosen x, removed x_r), exactly as branched in greedy_substitute.py
CASES <- list(
  c(1, 0), c(2, 0), c(3, 0),
  c(0, 1), c(2, 1), c(3, 1),
  c(0, 2), c(1, 2), c(3, 2),
  c(0, 3), c(1, 3), c(2, 3)
)

cat(sprintf("dim = %d (%s), fe = %d\n", DIM,
            ifelse(DIM == 1, "linear", "quadratic"), FE))

tab <- matrix(NA_real_, N_ALT, N_ALT,
              dimnames = list(paste0("drop_", 1:N_ALT), paste0("choose_", 1:N_ALT)))
missing <- character(0)
for (k in seq_along(CASES)) {
  x <- CASES[[k]][1]; xr <- CASES[[k]][2]
  f <- file.path(raw_dir, sprintf("subs_d%d_fe%d_%d_%d.txt", DIM, FE, x, xr))
  if (!file.exists(f)) { missing <- c(missing, basename(f)); next }
  v <- scan(f, quiet = TRUE)
  # rows are the dropped alternative, columns the one substituted to; +1 for 1-based labels
  tab[xr + 1, x + 1] <- v[length(v)]
}

if (all(is.na(tab))) stop("no result files in ", raw_dir,
                          " -- submit code/submit_substitute.sh first")

cat(sprintf("collected %d of 12 cells", 12 - length(missing)))
if (length(missing)) cat(sprintf("  -- MISSING: %s", paste(missing, collapse = ", ")))
cat("\n\n=== Table 3: maximal substitution ===\n")
print(round(tab, 4), na.print = "-")

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
df <- data.frame(drop_j = 1:N_ALT, tab, row.names = NULL, check.names = FALSE)
names(df) <- c("drop_j", paste0("choose_", 1:N_ALT))
for (j in 1:N_ALT) df[j, j + 1] <- NA          # diagonal is not defined
write.csv(format(df, digits = 4, trim = TRUE, na.encode = FALSE), out_csv,
          row.names = FALSE, quote = FALSE, na = "")
cat(sprintf("\nwrote %s\n", out_csv))
