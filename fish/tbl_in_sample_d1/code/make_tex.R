#!/usr/bin/env Rscript
# tbl:in-sample_d1 -- render output/tbl_in_sample_d1.csv as the In-Sample Fit table.
#
#   Rscript code/make_tex.R            # run from the exhibit directory
#
# Input : output/tbl_in_sample_d1.csv
# Output: table_latex/tbl_in_sample_d1.tex
args <- c('output/tbl_in_sample_d1.csv', 'table_latex/tbl_in_sample_d1.tex', 'in')
source('code/make_fit_tex.R')
