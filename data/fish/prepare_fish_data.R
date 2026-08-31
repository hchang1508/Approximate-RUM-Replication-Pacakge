#!/usr/bin/env Rscript
# Prepare the fishing-site choice datasets.
#
#   Rscript prepare_fish_data.R
#
# Source: the `Fishing` dataset from the R package `mlogit` (Croissant 2020), originally
# Thomson & Crooke (1991); see also Herriges & Kling (1999) and Cameron & Trivedi (2005)
# p.464. 1182 respondents choosing among 4 fishing modes, each described by a price and a
# catch rate (per hour fished, for major species by mode, summed over the respondent's
# targeted species).
#
# Output
#   fish_individual.csv   4728 rows = 1182 individuals x 4 alternatives, long format
#   fish_aggregated.csv   4 rows: characteristics averaged over individuals
#   fish_empirical_cp.csv 4 rows: empirical choice probabilities (market shares)
#   fish_J4_K2.csv        4 rows: the analysis input for the approximation-error and
#                         substitution exhibits -- price/100 and catch, in the estimation
#                         code's units (see "Units" below)
#
# Alternatives are numbered 1=beach, 2=boat, 3=charter, 4=pier. NOTE that mlogit's own
# column order is beach, pier, boat, charter, so the ordering is set explicitly here and
# every output file carries an alt_id alongside the name.
#
# Price is written POSITIVE and UNSCALED. Downstream code may want a different convention
# and should apply it explicitly.

suppressMessages(library(mlogit))

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(here) || !length(here)) here <- "."
setwd(here)

ALT <- c("beach", "boat", "charter", "pier")   # alt_id 1..4

data("Fishing", package = "mlogit")

## ---- long format, one row per (individual, alternative) --------------------------
long <- do.call(rbind, lapply(seq_along(ALT), function(i) {
  m <- ALT[i]
  data.frame(
    chid   = seq_len(nrow(Fishing)),
    alt_id = i,
    alt    = m,
    chosen = as.integer(as.character(Fishing$mode) == m),
    price  = Fishing[[paste0("price.", m)]],
    catch  = Fishing[[paste0("catch.", m)]],
    row.names = NULL
  )
}))
long <- long[order(long$chid, long$alt_id), ]
stopifnot(nrow(long) == 4 * nrow(Fishing),
          all(tapply(long$chosen, long$chid, sum) == 1))
write.csv(long, "fish_individual.csv", row.names = FALSE, quote = FALSE)

## ---- characteristics averaged over individuals -----------------------------------
agg <- data.frame(
  alt_id = seq_along(ALT),
  alt    = ALT,
  price  = sapply(ALT, function(m) mean(Fishing[[paste0("price.", m)]])),
  catch  = sapply(ALT, function(m) mean(Fishing[[paste0("catch.", m)]])),
  row.names = NULL
)
write.csv(agg, "fish_aggregated.csv", row.names = FALSE, quote = FALSE)

## ---- empirical choice probabilities ----------------------------------------------
cp <- data.frame(
  alt_id      = seq_along(ALT),
  alt         = ALT,
  choice_prob = sapply(ALT, function(m) mean(as.character(Fishing$mode) == m)),
  row.names = NULL
)
stopifnot(abs(sum(cp$choice_prob) - 1) < 1e-12)
write.csv(cp, "fish_empirical_cp.csv", row.names = FALSE, quote = FALSE)

## ---- analysis input: price/100, catch unscaled ------------------------------------
# The estimation scripts work in these units. Dividing price by 100 is a per-column
# rescaling, absorbed into beta, so it changes no theoretical quantity; it is applied so
# the two characteristics are on comparable numeric scales.
ana <- data.frame(
  alt_id = agg$alt_id,
  price_scaled = agg$price / 100,
  catch = agg$catch
)
write.csv(ana, "fish_J4_K2.csv", row.names = FALSE, quote = FALSE)

## ---- report ----------------------------------------------------------------------
cat("\n=== aggregated characteristics ===\n")
print(agg, row.names = FALSE, digits = 8)
cat("\n=== analysis input (fish_J4_K2.csv) ===\n")
print(ana, row.names = FALSE, digits = 8)
cat("\n=== empirical choice probabilities ===\n")
print(cp, row.names = FALSE, digits = 6)
cat("\n=== written ===\n")
for (f in c("fish_individual.csv", "fish_aggregated.csv", "fish_empirical_cp.csv",
            "fish_J4_K2.csv"))
  cat(sprintf("  %-24s %s rows\n", f, format(nrow(read.csv(f)), big.mark = ",")))
cat("\ndone\n")
