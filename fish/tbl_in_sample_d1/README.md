# tbl:in-sample_d1 — In-Sample Fit (Fish)

- **Paper**: online appendix `sec:insamp`
- **Dataset**: fishing-site choice, aggregated characteristics, single choice set D = {J}
- **"Our method"**: 4-mixture mixed-logit, **no fixed effects**
  (4 suffices by Propositions `lem:dim_P_r` and `cor:red`)
- **Estimation**: maximize log-likelihood sum_j rho_hat_j log rho(j|theta)
  - standard models: R `Rsolnp`
  - our method: EM algorithm (chosen over greedy for speed)
- **Design**: split individuals 50/50 into train/test, average within split,
  estimate on train. **Repeat with 50 random splits**; report mean and SD.
- **Metric**: l2 norm ||rho_hat_train - rho_train||_2

## Comparison models (defined in appendix `section:models`)
mnl, nested logit (charter vs rest), nested logit (boat vs rest),
random-coefficient logit with log-normal mixing, mnl + alternative FE,
random-coefficient log-normal + alternative FE

## Running it

This directory runs **both** fit tables — `tbl:out-sample` has no separate
computation, only its own collector.

    cd RUM_replication/fish/tbl_in_sample_d1
    sbatch code/submit_benchmarks.sh     # stage 1: R, array 1-50, ~2 min/split
    #   wait for it to finish -- stage 2 consumes output/em_input/
    sbatch code/submit_em.sh             # stage 2: python, array 1-50, ~1 min/split
    Rscript code/collect.R               # -> output/tbl_in_sample_d1.csv
    cd ../tbl_out_sample && Rscript code/collect.R    # -> output/tbl_out_sample.csv

`submit_em.sh` needs the `approx` conda env; `submit_benchmarks.sh` needs the
`R/4.4.1-foss-2022b` module plus `dplyr`, `fastDummies`, `mlogit`, `pracma`,
`Rsolnp`, `stringr`, `gtools`, `lpSolveAPI` (all present in that module).

The collector reports how many splits it actually found and warns if that is fewer
than 50 or uneven across models, so running it early is harmless.

| file | what it is |
|---|---|
| `code/in_out_fit.R` | driver, one random split per task: six benchmark models + this split's EM inputs |
| `code/Fish_functions_0528.R` | splitting, mlogit prep, aggregation, and the six model likelihoods |
| `code/Fish_functions.R` | `compute.emp.CP`, `compute.logitCP`, `dist.prob` |
| `code/approx_rum_EM_fish_inout.py` | "our method": 4-mixture EM on the single menu |
| `code/collect_fit_table.R` | shared aggregator over splits; both exhibits' `collect.R` wrap it |
| `output/raw/` | `benchmarks_<r>.csv`, `em_<r>.csv` — long: split, spec, sample, beach, boat, charter, pier, loss |
| `output/splits/`, `output/em_input/` | the per-split micro samples and their aggregated versions |

## Implementation notes

- **Characteristics** are built from `../../data/fish/fish_individual.csv`:
  `X1.0 = -price` divided by 10, and `X0.1 = catch`. `in_out_fit.R` asserts the
  aggregated result against the values the Python expects
  (`1.086, 0.2508529 / 0.585, 0.1693591 / 0.877, 0.6235062 / 1.086, 0.1691785`,
  up to that sign and factor of 10).
- **The SLURM task ID is the split index**, and `set.seed(index)` in `split.sample`
  makes each split reproducible from it.
- **`in_out_fit.R` writes the EM inputs** for its own split into `output/em_input/`,
  which `approx_rum_EM_fish_inout.py` then consumes.
- **`beta` is seeded deterministically** in the Python, so the EM run is
  reproducible from the split index.
- **Each `em_<r>.csv` carries a `sample=in` and a `sample=out` row.** The
  out-of-sample row evaluates the estimated mixture at the test half's average
  characteristics under the same l2 metric the benchmarks use.
- **Only the six `*_aggre` specifications are estimated** — the ones both tables
  report.

## Metric

Both sources report the paper's plain Euclidean `||rho_hat - rho||_2` with no
rescaling: R via `dist.prob` on the `*_aggre` specifications, Python via
`l2_distance` with `count = 1` (a single choice set). `collect_fit_table.R`
therefore applies no factor.


## Values

The full 50 splits were run on 2026-08-09 (jobs 59708344 benchmarks, 59708449 EM;
50/50 each, no failures). `output/tbl_in_sample_d1.csv` holds the table;
`table_latex/tbl_in_sample_d1.tex` is the formatted version in the paper.

| model | in-sample error (50 splits) |
|---|---|
| our_method | 5.08e-06 |
| mnl | 0.0378 |
| nest_charter | 5.85e-04 |
| nest_boat | 0.0378 |
| mixed_logit_lognormal | 0.0382 |
| mnl_FE | 7.82e-06 |
| mixed_logit_lognormal_FE | 3.10e-05 |

The same run produces the out-of-sample table in `../tbl_out_sample` — see that
README.

### `table_latex/tbl_in_sample_d1.tex`

The formatted table, holding nothing but the `table` environment — no preamble, no
comments — so it can be `\input` straight into the manuscript. It needs `multirow`,
`graphicx` for `\scalebox`, and `threeparttable` for the `tablenotes` block, and it
refers to `\ref{tbl:in-sample_d1}` in its own note.

Point estimates and standard deviations are formatted to three decimals from
`output/tbl_in_sample_d1.csv`; that file carries the unrounded figures and
`output/raw/` the per-split results behind them. Regenerate with
`Rscript code/collect.R` then `Rscript code/make_tex.R`; the file is generated, not
hand-edited.
