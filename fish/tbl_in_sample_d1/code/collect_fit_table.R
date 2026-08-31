#!/usr/bin/env Rscript
# Shared aggregator for tbl:in-sample_d1 and tbl:out-sample. Both tables come from the
# same 50 splits, so both are built here; each exhibit's code/collect.R is a two-line
# wrapper that sets `args` and sources this file.
#
#   args[1] raw_dir  directory holding benchmarks_<r>.csv and em_<r>.csv
#   args[2] out_csv  where to write the table
#   args[3] sample   "in" or "out"
#
# Input : raw_dir/benchmarks_<r>.csv  from code/in_out_fit.R
#         raw_dir/em_<r>.csv          from code/approx_rum_EM_fish_inout.py
#         both long: split, spec, sample, beach, boat, charter, pier, loss
# Output: out_csv, one row per model:
#         model, beach, beach_sd, ..., pred_error, pred_error_sd
#         means and standard deviations across the splits, matching the paper's
#         "In parentheses are standard deviations obtained by repeating the same
#         analyses 50 times."
#
# No rescaling is applied. Both sources already report the paper's metric, the plain
# Euclidean ||rho_hat - rho||_2: R via dist.prob on the six *_aggre specifications,
# Python via l2_distance with a single choice set (count = 1). The x2 correction in
# this exhibit's README concerns the older saved RDS, which used dist.1 = ||.||_2/2,
# and the /sqrt(n.alt) in run.alt.regressions.macro applies only to the *_micro
# specifications, which neither table reports.

if (!exists('args')) args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else 'output/raw'
out_csv <- if (length(args) >= 2) args[2] else 'output/fit_table.csv'
sample_ <- if (length(args) >= 3) args[3] else 'in'
stopifnot(sample_ %in% c('in', 'out'))

# paper row order and labels; the names on the left are what the raw files carry
ROWS <- c(our_method            = 'our_method',
          mnl_aggre             = 'mnl',
          nest_charter_aggre    = 'nest_charter',
          nest_boat_aggre       = 'nest_boat',
          random_coef_aggre     = 'mixed_logit_lognormal',
          mnl_FE_aggre          = 'mnl_FE',
          random_coef_FE_aggre  = 'mixed_logit_lognormal_FE')
ALTS <- c('beach', 'boat', 'charter', 'pier')

files <- c(list.files(raw_dir, pattern = '^benchmarks_[0-9]+\\.csv$', full.names = TRUE),
           list.files(raw_dir, pattern = '^em_[0-9]+\\.csv$',         full.names = TRUE))
if (length(files) == 0) {
  stop(sprintf('no results in %s -- run code/submit_benchmarks.sh and code/submit_em.sh first',
               raw_dir))
}

res <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
res <- res[res$sample == sample_, ]
res$spec[res$spec %in% names(ROWS)] <- ROWS[res$spec[res$spec %in% names(ROWS)]]

# report what is actually there before averaging anything
n_by_spec <- table(factor(res$spec, levels = unname(ROWS)))
cat(sprintf('%s-sample: %d split(s) per model\n', sample_, max(n_by_spec)))
for (m in unname(ROWS)) cat(sprintf('  %-28s %3d\n', m, n_by_spec[[m]]))
if (length(unique(n_by_spec[n_by_spec > 0])) > 1) {
  cat('WARNING: unequal split counts across models -- some array tasks are missing\n')
}
if (max(n_by_spec) < 50) {
  cat(sprintf('WARNING: the paper averages 50 splits; this table averages %d\n', max(n_by_spec)))
}

tab <- do.call(rbind, lapply(unname(ROWS), function(m) {
  r <- res[res$spec == m, ]
  if (nrow(r) == 0) return(NULL)
  vals <- c(unlist(lapply(ALTS, function(a) c(mean(r[[a]]), sd(r[[a]])))),
            mean(r$loss), sd(r$loss))
  data.frame(model = m, t(vals), n = nrow(r), stringsAsFactors = FALSE)
}))
colnames(tab) <- c('model',
                   as.vector(rbind(ALTS, paste0(ALTS, '_sd'))),
                   'pred_error', 'pred_error_sd', 'n')

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(tab, out_csv, row.names = FALSE)
cat(sprintf('\nwrote %s\n\n', out_csv))
print(format(tab[, -ncol(tab)], digits = 3), row.names = FALSE)
