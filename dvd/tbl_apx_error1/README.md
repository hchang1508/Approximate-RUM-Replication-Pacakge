# tbl:apx_error1 — Approximation errors to rho^pi (DVD, no fixed effects)

Paper: main text, Section 5.1 (`sec:without_fixed`).

## Running it

    sbatch code/submit_greedy.sh                                    # 24 array tasks
    sbatch --export=ALL,DATASET=dvd code/submit_em_vertices.sh      # 24 array tasks
    Rscript code/collect.R                                     # -> output/tbl_apx_error1.csv
    Rscript code/make_tex.R                                    # -> table_latex/tbl_apx_error1.tex

`code/submit_em_vertices.sh` exports every EM setting explicitly, so this reproduces the
committed CSV with no further arguments. Results and SLURM logs go to scratch; `DATASET`
selects the exhibit, so the same script serves `dvd` and `fish`.

`collect.R` takes an optional raw directory and output path. It accepts either a single
directory holding both `DVD_gd_*` and `DVD_em_*`, or a parent holding `greedy/` and
`em_tight/`, and prints which layout it detected:

    Rscript code/collect.R \
      /vast/palmer/scratch/narita/hc654/RUM_replication/results_final/dvd/tbl_apx_error1

Both solvers take `argv = index, order, data`:

- greedy: `code/approx_rum_GREEDY_vertices_DVD.py`
- EM: `code/approx_rum_EM_vertices_DVD.py`

## Design

- **Dataset**: DVD, `../../data/dvd/dvd_J4_K2.csv`, first 4 rows of `farias2009_table2.csv`
- **Characteristics**: avg. price per disc, total helpful votes — so |J| = 4, K = 2
- **Fixed effects**: none (eta = 0)
- **Choice sets**: pairs, triples, and the full set, built in `cal_D`. Singletons are never
  constructed, so |D| = 11 rather than 2^J - 1 = 15. A one-item menu pins the choice
  probability to 1 under every model, so it contributes no squared error and nothing to
  `dim P_r + 1 = 1 + sum_D (|D|-1) = 18`.

| col | model | algorithm |
|---|---|---|
| (1) | linear mixed-logit (K=2) | greedy |
| (2) | linear mixed-logit (K=2) | EM |
| (3) | quadratic mixed-logit | greedy |
| (4) | quadratic mixed-logit | EM |

`order=1` is linear in characteristics, `order=2` quadratic.

## Settings

Greedy runs 1000 iterations, fixed at the call site in
`approx_rum_GREEDY_vertices_DVD.py`. EM settings are overridable by environment variable and
all exported explicitly by `code/submit_em_vertices.sh`:

| variable | value | meaning |
|---|---|---|
| `EM_STARTS` | 10 | random initial points; the reported error is the minimum over them |
| `EM_MIXTURE` | 18 | mixture components, = `dim P_r + 1` for \|J\|=4 |
| `EM_MAX_ITER` | 1000 | iteration cap per fit |
| `EM_TOL` | 1e-8 | absolute change in error between iterations at which EM stops |
| `EM_STANDARDIZE` | 1 | centre and scale features after any order-2 expansion |
| `SEED_BASE` | 20240808 | per-task seed is `SEED_BASE + 1000*order + index` |

Standardising is a reparametrisation — utilities shift by a constant within each choice
set, which cancels in the logit, and beta rescales — so the attainable error is unchanged
and only the iteration path differs.

These match `fish/tbl_fish_apx_error1`, so the two EM columns are directly comparable.

## Normalisation — sqrt(11)

`l2_distance` returns

    raw = ||rho - rho_hat||_2 / count,     count = |cal_D| = 11

dividing by the number of choice sets, not its square root. The paper's metric divides by
`sqrt(|D|)`, and the simulation never evaluates a singleton menu, so |D| = count = 11.
`collect.R` therefore rescales by

    d = raw * count/sqrt(|D|) = raw * 11/sqrt(11) = raw * sqrt(11) = raw * 3.3166

This is the only scaling applied, here and in every other |J|=4 exhibit.
