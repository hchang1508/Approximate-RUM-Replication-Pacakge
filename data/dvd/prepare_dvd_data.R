#!/usr/bin/env Rscript
# Prepare the DVD analysis datasets from the source table.
#
#   Rscript prepare_dvd_data.R
#
# Input :  farias2009_table2.csv   (all 15 products, raw units, as printed in the source)
# Output:  dvd_J4_K2.csv           |J|=4,  K=2  -- main text Section 5
#          dvd_J10_K3.csv          |J|=10, K=3  -- appendix sec:J10K2
#
# Scaling, applied uniformly:
#          price              / 100
#          avg_price_per_disc / 10
#          total_helpful_votes / 1000
# Per-column rescaling is innocuous for the theory -- a diagonal positive rescaling of
# the characteristics is absorbed into beta -- but it is reproduced here so that output
# is directly comparable across exhibits.
#
# Every downstream script reads its characteristics from the CSVs this script writes;
# none hard-codes them.

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(here) || !length(here)) here <- "."
setwd(here)

SCALE <- c(price = 100, avg_price_per_disc = 10, total_helpful_votes = 1000)

raw <- read.csv("farias2009_table2.csv")
stopifnot(nrow(raw) == 15, all(names(SCALE) %in% names(raw)))
# The paper selects "the first four / first ten rows", which are the most expensive.
stopifnot(!is.unsorted(rev(raw$price)))

scale_cols <- function(d) {
  for (nm in names(SCALE)) d[[nm]] <- d[[nm]] / SCALE[[nm]]
  d
}

emit <- function(d, chars, file) {
  out <- d[, c("product_id", chars)]
  names(out)[1] <- "alt_id"
  write.csv(out, file, row.names = FALSE, quote = FALSE)
  cat(sprintf("  %-26s %d x %d\n", file, nrow(out), length(chars)))
  out
}

K2 <- c("avg_price_per_disc", "total_helpful_votes")
K3 <- c("avg_price_per_disc", "total_helpful_votes", "price")

cat("writing datasets:\n")
j4  <- emit(scale_cols(raw[1:4, ]),  K2, "dvd_J4_K2.csv")
j10 <- emit(scale_cols(raw[1:10, ]), K3, "dvd_J10_K3.csv")

cat("\nsanity checks:\n")
# With K=2 and |J|=4 the set can never be affinely independent, which is exactly the
# failure the paper's application is built on. Assert it rather than assume it.
r4 <- qr(cbind(1, as.matrix(j4[, K2])))$rank      # affinely independent <=> rank == |J|
cat(sprintf("  |J|=4   rank[1 X] = %d, need 4 -> %s\n",
            r4, ifelse(r4 == 4, "independent", "NOT independent, K = 2 < 3 = |J|-1")))
stopifnot(r4 == 3)
r10 <- qr(cbind(1, as.matrix(j10[, K3])))$rank
cat(sprintf("  |J|=10  rank[1 X] = %d of 10 -> NOT independent (K = 3 < 9)\n", r10))
stopifnot(r10 == 4)
cat("\ndone\n")
