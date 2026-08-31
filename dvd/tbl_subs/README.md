# tbl:subs — Maximal substitution of the linear mixed-logit models (DVD)

    sbatch code/submit_substitute.sh     # 12 array tasks, one per table cell
    Rscript code/collect.R               # -> output/tbl_subs.csv
    Rscript code/make_tex.R              # -> table_latex/tbl_subs.tex

Both readers take the raw directory from `RAW_DIR`, defaulting to `output/raw`. That
default is right once the submit step has run here. The twelve files behind the committed
table live in the consolidated results tree on scratch; to build from there, set `RAW_DIR`
explicitly:

    R=/vast/palmer/scratch/narita/hc654/RUM_replication/results_final/dvd/tbl_subs/raw
    RAW_DIR=$R Rscript code/collect.R
    RAW_DIR=$R Rscript code/make_tex.R

Without it both fail with `missing output/raw/subs_d1_fe0_1_0.txt`. `OUT_CSV` and
`OUT_TEX` likewise override the destinations, and `DIM` / `FE` select the model — the
table reports `DIM=1`, `FE=0`.

`make_tex.R` reads full precision from `RAW_DIR` rather than the 4-decimal CSV, so the
printed value is truncated once rather than twice.

- **Paper**: main text, Section 5.2 (`sec:substitution`)
- **Dataset**: DVD rows 1-4, K=2 — `../../data/dvd/dvd_J4_K2.csv`
- **Quantity**: eq (`eq:dist.2`), sup over P_ml(0) of [ rho(J\{j}, l) - rho(J, l) ]
- **Algorithm**: greedy (the EM algorithm cannot be adapted to this objective)
- **Fixed effects**: none. With free fixed effects the quantity is always 1.
- Alternatives 1-4 are the four most expensive DVDs.

## How the 12 tasks map to the table

`case` (the array index, 1..12) selects an ordered pair via the branches in
`greedy_substitute.py`: `x` is the alternative substituted *to*, `x_r` the one removed.
In the table, rows are the dropped alternative and columns the chosen one, so cell
(drop j, choose l) comes from the task with `x_r = j-1`, `x = l-1`. The diagonal is
undefined and left blank. `collect.R` holds the same map and it has been checked against
the twelve `elif case==` branches.

Arguments are `case dim fe`; the submit script passes `dim=1` (linear mixed logit, which
is what the table reports) and `fe=0`.

## Values

`output/tbl_subs.csv` holds the table at four decimals; `RAW_DIR` holds the same values at
full precision; `table_latex/tbl_subs.tex` is the formatted version.

| | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| **drop 1** | – | 0.999 | 0.138 | 0.999 |
| **drop 2** | 0.999 | – | 0.997 | 0.104 |
| **drop 3** | 0.214 | 0.999 | – | 0.999 |
| **drop 4** | 0.999 | 0.288 | 0.999 | – |

The 1<->3 and 2<->4 pairs are the ones the linear model cannot generate. Section 5.2 and
the assortment example rest on that structure, not on the magnitudes.

Values are **truncated** at the third decimal (`floor(x * 1000) / 1000`,
`make_tex.R:32`) and read from `RAW_DIR`, not from the 4-decimal CSV, so they are cut once
rather than twice. Rerun `make_tex.R` when the numbers move; the file is generated, not
hand-edited.
