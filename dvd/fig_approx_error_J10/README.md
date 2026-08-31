# fig:approx_error_J10 — Approximation errors, |J|=10 and K=3 (DVD)

Paper: online appendix `sec:J10K2`.
Artifact: `MS_submission/approximation_error.png`.

## Running it

    sbatch code/submit_greedy.sh                  # 100 array tasks, partition=week

    module load R/4.4.1-foss-2022b                # Rscript is not on the bare PATH
    Rscript code/collect_and_plot.R               # -> output/fig_approx_error_J10.csv
                                                  #    output/approximation_error.png

Each task writes `output/raw/Large_DVD_gd_1_<index>.txt`, eleven numbers: `pi(0..9)` then
the raw error. Task `i` draws its ranking with `np.random.seed(i)`, so the figure is
reproducible bin-for-bin.

The raw per-task output lives in the consolidated results tree on scratch.
`collect_and_plot.R` takes an optional raw directory, output CSV, and output PNG; point
it there:

    Rscript code/collect_and_plot.R \
      /vast/palmer/scratch/narita/hc654/RUM_replication/results_final/dvd/fig_approx_error_J10/raw

Solver: `code/approx_rum_GREEDY_vertices_DVD_large.py`, taking `argv = index, order, data`.
Only `order=1` is needed — the figure reports linear mixed-logit errors.

## Design

- **Dataset**: DVD, first 10 rows of `../../data/dvd/farias2009_table2.csv`, prepared as
  `../../data/dvd/dvd_J10_K3.csv`
- **Alternatives**: |J| = 10
- **Characteristics**: avg. price per disc, total helpful votes, and price, so K = 3
- **Choice sets**: pairs, triples, and the full set only — quartets up are intractable at
  |J|=10 — so |D| = C(10,2) + C(10,3) + 1 = 166
- **Sample**: 100 rankings drawn uniformly from all 10!

At 166 menus per task this is the heaviest DVD exhibit, hence `partition=week` and the
7-day limit; an |J|=4 task is far cheaper.

## Normalisation — sqrt(166)

`l2_distance` returns

    raw = ||rho - rho_hat||_2 / count,     count = |cal_D| = 166

dividing by the number of choice sets, not its square root. The paper's metric divides by
`sqrt(|D|)`, and singleton menus are never evaluated, so |D| = count = 166.
`collect_and_plot.R` therefore rescales by

    d = raw * count/sqrt(|D|) = raw * 166/sqrt(166) = raw * sqrt(166) = raw * 12.8841

Same convention as the |J|=4 exhibits, which set |D| = count = 11 and rescale by sqrt(11);
only the menu count differs. See `../tbl_apx_error1/code/collect.R`.
