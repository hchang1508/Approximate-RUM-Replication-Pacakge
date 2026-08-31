#!/usr/bin/env Rscript
# tbl:in-sample_d1 -- aggregate the 50 splits into Table "In-Sample Fit".
#
#   Rscript code/collect.R            # run from the exhibit directory
#
# Output: output/tbl_in_sample_d1.csv
args <- c('output/raw', 'output/tbl_in_sample_d1.csv', 'in')
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])),
                 'collect_fit_table.R'))
