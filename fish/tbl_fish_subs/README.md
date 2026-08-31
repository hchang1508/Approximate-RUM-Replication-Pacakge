# tbl:fish_subs — Maximal substitution of the linear mixed-logit models (Fish)

    sbatch code/submit_substitute.sh     # 12 array tasks, one per off-diagonal cell
    Rscript code/collect.R               # -> output/tbl_fish_subs.csv
    Rscript code/make_tex.R              # -> table_latex/tbl_fish_subs.tex

Both readers take the raw directory from `RAW_DIR`, defaulting to `output/raw`. That
default is right once the submit step has run here. The twelve files behind the committed
table live in the consolidated results tree on scratch; to build from there, set `RAW_DIR`
explicitly:

    R=/vast/palmer/scratch/narita/hc654/RUM_replication/results_final/fish/tbl_fish_subs/raw
    RAW_DIR=$R Rscript code/collect.R
    RAW_DIR=$R Rscript code/make_tex.R

Without it both fail with `missing output/raw/subs_d1_fe0_1_0.txt`. `OUT_CSV` and
`OUT_TEX` likewise override the destinations, and `DIM` / `FE` select the model — the
table reports `DIM=1`, `FE=0`.

`make_tex.R` reads full precision from `RAW_DIR` rather than the 4-decimal CSV, so the
printed value is truncated once rather than twice.

- **Paper**: online appendix `sec:fish_substitution`
- **Dataset / characteristics**: same as `tbl_fish_apx_error1` — `../../data/fish/fish_J4_K2.csv`, |J|=4, K=2
- **Quantity**: eq (`eq:dist.2`), sup over P_ml(0) of [ rho(J\{j}, l) - rho(J, l) ]
- **Algorithm**: greedy (the EM algorithm cannot be adapted to this objective)
- **Fixed effects**: none. With free fixed effects the quantity is always 1.
- Alternatives: 1 = beach, 2 = boat, 3 = charter, 4 = pier

## How the 12 tasks map to the table

`case` (the array index, 1..12) selects an ordered pair via the branches in
`greedy_substitute.py`: `x` is the alternative substituted *to*, `x_r` the one removed.
In the table, rows are the dropped alternative and columns the chosen one, so cell
(drop j, choose l) comes from the task with `x_r = j-1`, `x = l-1`. The diagonal is
undefined. `collect.R` holds the same map.

Arguments are `case dim fe`; the submit script passes `dim=1` (linear mixed logit, which
is what the table reports) and `fe=0`.


## Values

`output/tbl_fish_subs.csv` holds the table at four decimals;
`table_latex/tbl_fish_subs.tex` is the formatted version, truncated to three:

| | 1 beach | 2 boat | 3 charter | 4 pier |
|---|---|---|---|---|
| **1 beach** | – | 0.119 | 0.997 | 0.997 |
| **2 boat** | 0.316 | – | 0.999 | 0.997 |
| **3 charter** | 0.997 | 0.999 | – | 0.285 |
| **4 pier** | 0.994 | 0.997 | 0.137 | – |

The row-wise minima put the constrained pairs at **beach <-> boat** and
**charter <-> pier**: substitution within those pairs is what the linear mixed-logit
model cannot generate. That structure, not the magnitudes, is what the appendix argues
from.

### `table_latex/tbl_fish_subs.tex`

The formatted table, holding the `table` environment and one leading `% Generated <date>`
comment line and nothing else — no preamble — so it can be `\input` straight into the
manuscript. It needs `diagbox` and `graphicx` for `\scalebox`, and its note refers to
`\ref{eq:dist.2}`.

Values are **truncated** at the third decimal (`floor(x * 1000) / 1000`, `make_tex.R:32`)
and read from `RAW_DIR`, not from the 4-decimal CSV, so they are cut once rather than
twice. Rerun `make_tex.R` when the numbers move; the file is generated, not hand-edited.
