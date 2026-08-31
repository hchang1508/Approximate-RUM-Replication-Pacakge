#!/usr/bin/env Rscript
# tbl:in-sample_d1 and tbl:out-sample -- the six benchmark models, one random split per task.
#
#   Rscript code/in_out_fit.R <split_index>      # run from the exhibit directory
#
# Characteristics are built from ../../data/fish/fish_individual.csv. Fish_d1.csv
# holds mode, alt, chid, X1.0 = -price and X0.1 = catch, with X1.0 then divided by
# 10; both conventions are reproduced below and asserted against the aggregated file.
#
# For split r this script
#   1. splits the 1,182 individuals into K = 2 halves with set.seed(r),
#   2. estimates every benchmark specification on half k = 1's training set and
#      predicts on both halves,
#   3. writes the aggregated characteristics and choice probabilities of that same
#      split to output/em_input/, which is what code/approx_rum_EM_fish_inout.py
#      (the "our method" row) consumes.
#
# Only the six *_aggre specifications are estimated.
#
# Output: output/raw/benchmarks_<r>.csv   long, one row per (spec, sample):
#           spec, sample in {in, out}, beach, boat, charter, pier, loss
#         output/splits/data.{train,test}.<r>.<k>.csv
#         output/em_input/data.{train,test}[.CP]<r>.1.csv
#
# The metric is dist.prob = plain Euclidean ||rho_hat - rho||_2.

suppressPackageStartupMessages({
  library(dplyr)
  library(fastDummies)
  library(mlogit)
  library(Rsolnp)
})

code_dir <- dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
if (is.na(code_dir) || !nzchar(code_dir)) code_dir <- 'code'
source(file.path(code_dir, 'Fish_functions.R'))
source(file.path(code_dir, 'Fish_functions_0528.R'))

args  <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])
if (is.na(index)) stop('usage: Rscript code/in_out_fit.R <split_index>')

fish_data <- Sys.getenv('FISH_DATA', '../../data/fish/fish_individual.csv')
out_dir   <- Sys.getenv('OUT_DIR', 'output')
split_dir <- file.path(out_dir, 'splits/')
em_dir    <- file.path(out_dir, 'em_input/')
raw_dir   <- file.path(out_dir, 'raw')
for (d in c(split_dir, em_dir, raw_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

set.seed(index)


######################################################
############## Preparation ###########################
######################################################

d         <- 1                     # degree of polynomials
pol.label <- c('X1.0', 'X0.1')     # characteristics to be used
K         <- 2                     # splits sample half-and-half

# Load characteristics and create the universal choice set.
# Column ORDER matters: do.empirical.analysis addresses characteristics positionally
# as dat[, c(4, 5)], so it must be mode, alt, chid, X1.0, X0.1.
fish <- read.csv(fish_data, stringsAsFactors = FALSE)
fish <- fish[order(fish$chid, fish$alt_id), ]
dat  <- data.frame(
  mode = fish$chosen,
  alt  = fish$alt,
  chid = fish$chid,
  X1.0 = -fish$price,              # price enters with a negative sign, as in Fish_d1.csv
  X0.1 = fish$catch,
  stringsAsFactors = FALSE
)
stopifnot(identical(colnames(dat)[4:5], pol.label))
stopifnot(identical(unique(dat$alt), c('beach', 'boat', 'charter', 'pier')))

dat[, 4] <- dat[, 4] / 10          # normalize

CS        <- unique(dat$alt)       # list of alternatives
nalt      <- length(CS)
indiv.idx <- unique(dat$chid)      # list of individuals
cat(sprintf('split %d: %d individuals x %d alternatives\n', index, length(indiv.idx), nalt))


######################################################
############## Benchmark models ######################
######################################################

specify.scale <- NULL
output <- do.empirical.analysis(dat = dat, K, pol.label, specify.scale, index, split_dir)
result <- output$results

spec <- c('mnl_aggre', 'nest_charter_aggre', 'nest_boat_aggre', 'random_coef_aggre',
          'mnl_FE_aggre', 'random_coef_FE_aggre')

### We only use the results of the first train-test split.
rows <- do.call(rbind, lapply(spec, function(method) {
  rbind(
    data.frame(split = index, spec = method, sample = 'in',
               t(setNames(as.numeric(result[[1]]$`in-sample_pred`[[method]]), CS)),
               loss = as.numeric(result[[1]]$`in-sample_loss`[[method]]),
               stringsAsFactors = FALSE),
    data.frame(split = index, spec = method, sample = 'out',
               t(setNames(as.numeric(result[[1]]$`out-of-sample_pred`[[method]]), CS)),
               loss = as.numeric(result[[1]]$`out-of-sample_loss`[[method]]),
               stringsAsFactors = FALSE)
  )
}))
colnames(rows) <- c('split', 'spec', 'sample', CS, 'loss')

write.csv(rows, file.path(raw_dir, sprintf('benchmarks_%d.csv', index)), row.names = FALSE)


######################################################
############## Inputs for the EM ("our method") ######
######################################################
# Body of make_files_for_EM.R, applied to the k = 1 split of this task.

data_temp_train <- output$samples[[1]][['train']]
data_temp_test  <- output$samples[[1]][['test']]

emp_CP_train <- compute.emp.CP(data_temp_train)
emp_CP_test  <- compute.emp.CP(data_temp_test)

data_temp_train <- gen.macro.data(data_temp_train, pol.label)
data_temp_test  <- gen.macro.data(data_temp_test,  pol.label)

write.csv(data_temp_train, paste0(em_dir, 'data.train.',    index, '.1', '.csv'))
write.csv(emp_CP_train,    paste0(em_dir, 'data.train.CP',  index, '.1', '.csv'))
write.csv(data_temp_test,  paste0(em_dir, 'data.test.',     index, '.1', '.csv'))
write.csv(emp_CP_test,     paste0(em_dir, 'data.test.CP',   index, '.1', '.csv'))

cat(sprintf('split %d done\n', index))
