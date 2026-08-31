# tbl:out-sample — Out-of-Sample Fit (Fish)

- **Paper**: online appendix `sec:insamp`
- Same design, models, and 50 splits as `../tbl_in_sample_d1` — read that README first.
- **Metric**: l2 norm ||rho_hat_test - rho_test||_2, predictions formed from
  train-sample estimates applied to test-sample characteristics.

## What the table shows
Our method (0.052) is on par with the two fixed-effects models (0.048) and clearly ahead
of the multinomial and nested-boat specifications (0.058), but nested logit / charter
(0.050) edges it. 

## Running it

There is no separate computation here. Both fit tables come out of the same run, in
`../tbl_in_sample_d1`: `code/in_out_fit.R` estimates each model on the training half
and predicts on both halves, and `code/approx_rum_EM_fish_inout.py` does the same for
"our method". This directory holds only the collector.

    cd ../tbl_in_sample_d1
    sbatch code/submit_benchmarks.sh     # then wait
    sbatch code/submit_em.sh
    cd ../tbl_out_sample && Rscript code/collect.R    # -> output/tbl_out_sample.csv

`code/collect.R` reads `../tbl_in_sample_d1/output/raw` and wraps the shared
aggregator `../tbl_in_sample_d1/code/collect_fit_table.R`.

## Where this table's numbers come from

For the six benchmark models: `run.alt.regressions.macro` in
`Fish_functions_0528.R`, which already predicted on the test half
(`pred.out` / `loss.out`).

For "our method": `approx_rum_EM_fish_inout.py`, which evaluates the mixture estimated
on the training half at the test half's average characteristics and scores it under the
same l2 metric, mirroring what the benchmarks do.

## Values

The full 50 splits were run on 2026-08-09. `output/tbl_out_sample.csv` holds the table;
`table_latex/tbl_out_sample.tex` is the formatted version.

| model | pred_error | sd |
|---|---|---|
| mnl_FE | 0.0476 | 0.0221 |
| mixed_logit_lognormal_FE | 0.0480 | 0.0169 |
| nest_charter | 0.0498 | 0.0198 |
| our_method | 0.0516 | 0.0216 |
| mixed_logit_lognormal | 0.0576 | 0.0184 |
| mnl | 0.0580 | 0.0191 |
| nest_boat | 0.0580 | 0.0191 |


### `table_latex/tbl_out_sample.tex`

The formatted table, holding nothing but the `table` environment — no preamble, no
comments — so it can be `\input` straight into the manuscript. It needs `multirow`,
`graphicx` for `\scalebox`, and `threeparttable` for the `tablenotes` block, and its note
refers to `\ref{tbl:out-sample}` and `\ref{eq:dist.1}`.

Point estimates and standard deviations are formatted to three decimals from
`output/tbl_out_sample.csv`; that file carries the unrounded figures, and the per-split
results behind them are in `../tbl_in_sample_d1/output/raw/`. Regenerate with
`Rscript code/collect.R` then `Rscript code/make_tex.R`; the file is generated, not
hand-edited.
