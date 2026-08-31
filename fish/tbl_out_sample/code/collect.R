#!/usr/bin/env Rscript
# tbl:out-sample -- aggregate the 50 splits into Table "Out-of-Sample Fit".
#
#   Rscript code/collect.R            # run from the exhibit directory
#
# There is no separate computation for this table. Both fit tables come from the same
# run, which lives in ../tbl_in_sample_d1: its code/in_out_fit.R estimates every model
# on the training half and predicts on both halves, and its
# code/approx_rum_EM_fish_inout.py does the same for "our method". Submit the two
# arrays there, then run this.
#
# Input : ../tbl_in_sample_d1/output/raw
# Output: output/tbl_out_sample.csv
args <- c('../tbl_in_sample_d1/output/raw', 'output/tbl_out_sample.csv', 'out')
source('../tbl_in_sample_d1/code/collect_fit_table.R')
