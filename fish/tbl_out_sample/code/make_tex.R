#!/usr/bin/env Rscript
# tbl:out-sample -- render output/tbl_out_sample.csv as the Out-of-Sample Fit table.
#
#   Rscript code/make_tex.R            # run from the exhibit directory
#
# There is no separate computation for this table; the renderer lives with the run, in
# ../tbl_in_sample_d1, exactly as code/collect.R does.
#
# Input : output/tbl_out_sample.csv
# Output: table_latex/tbl_out_sample.tex
args <- c('output/tbl_out_sample.csv', 'table_latex/tbl_out_sample.tex', 'out')
source('../tbl_in_sample_d1/code/make_fit_tex.R')
