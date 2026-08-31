# tbl:apx_error2 — Approximation errors to a half-half mixture (DVD, with fixed effects)

Paper: main text, Section 5.1 (also labelled `tab:representability_macro`).

For each of six rankings `pi`, this table reports how closely a mixed-logit model can
approximate the half-half mixture `(1/2) rho^pi + (1/2) rho^pi-bar`, where `pi-bar` reverses
`pi`. Four numbers per ranking: linear and quadratic models, each fitted by a greedy
algorithm and by EM.

## Requirements

- Python 3 with `numpy`, `scipy`, `pandas`, `joblib`
- R (base only)
- SLURM, for the two array submissions. Each cell is a 9261-point grid search, so this is
  not a laptop job — see "Scale" below.
- `collect.R` shells out to `find` and `awk`.

## Running it

Four steps, in order. Run everything from this directory.

    Rscript code/make_todo_list_GREEDY.R          # -> to_do_list_GREEDY/
    bash code/submit_greedy.sh all 1              # degree 1; repeat with 2
    bash code/submit_em.sh dvd 1 1-5              # degree 1, replicates 1-5; repeat with 2
    Rscript code/collect.R                        # -> output/tbl_apx_error2.csv
    Rscript code/make_tex.R                       # -> table_latex/tbl_apx_error2.tex

**Step 1 must precede every greedy submission.** The greedy script resolves its array index
through the to-do list (`index = int(to_do[index])`), so a SLURM array index addresses a
position in the *remaining* work, not a grid point directly. `make_todo_list_GREEDY.R` diffs
the finished output files against the full grid and writes the remainder, which is what makes
the job resumable: rerun it after a partial run and resubmit to pick up only what is missing.

Never rebuild the to-do lists while an array is queued or running. The list is read at task
runtime, so rebuilding it mid-flight shifts every index under the running tasks.

`collect.R` reads `output/raw` and `output/raw_em_tight` by default. Pass a different root
as the first argument, or set `GREEDY_DIR` and `EM_ROOT` to point the two halves at separate
trees:

    Rscript code/collect.R <root>            # expects <root>/greedy and <root>/em_tight/rep<N>
    GREEDY_DIR=... EM_ROOT=... Rscript code/collect.R

It takes about 90 seconds over a complete set.

## Scale

Per degree, the greedy step is 6 cases x 9261 grid points = 55,566 independent tasks. The EM
step is 6 x 371 chunks = 2,226 tasks per replicate, and the table uses five replicates. Both
degrees together are roughly 111,000 greedy tasks and 22,000 EM tasks.

Both steps are embarrassingly parallel and fully resumable, so a smaller allocation only
costs wall-clock. `collect.R` takes the minimum over the grid points it finds and prints
per-cell coverage, so it can be run at any point to see where the arrays stand.

## Design

- **Dataset**: DVD, `data/dvd/dvd_J4_K2.csv` — average price per disc and total helpful
  votes for the four most expensive products, so |J| = 4 and K = 2.
- **Target**: `(1/2) rho^pi + (1/2) rho^pi-bar` for each of the six cases below.
- **Fixed effects**: grid search over {-10,...,10}^3, with the fourth alternative's effect
  normalised to zero, giving 21^3 = 9261 points per cell. **Each reported cell is the
  minimum over that grid.** Reported errors are therefore grid-search errors; the true
  infimum over unrestricted fixed effects is weakly smaller.
- **Degrees**: `degree=1` is linear in the characteristics, `degree=2` quadratic.

| case | pi | pi-bar |
|---|---|---|
| 1 | 1>3>4>2 | 2>4>3>1 |
| 2 | 1>3>2>4 | 4>2>3>1 |
| 3 | 1>4>3>2 | 2>3>4>1 |
| 4 | 2>4>1>3 | 3>1>4>2 |
| 5 | 2>1>4>3 | 3>4>1>2 |
| 6 | 3>1>2>4 | 4>2>1>3 |

These are the six rankings that a linear mixed-logit cannot represent on this dataset.

## Settings

Greedy runs 1000 iterations, fixed at the call site in `approx_rum_GREEDY_hhmix_DVD.py`.

EM is controlled by environment variables, all exported by `code/submit_em.sh`:

| variable | value | meaning |
|---|---|---|
| `EM_STARTS` | 1 | random initial points per task — see below |
| `EM_MAX_ITER` | 1000 | iteration cap per fit |
| `EM_TOL` | 1e-8 | absolute change in error at which EM stops |
| `EM_STANDARDIZE` | 1 | centre and scale features after any degree-2 expansion |
| `EM_CHUNK` | 25 | grid points per array task; 371 chunks tile the grid |
| `SEED_BASE` | `20240808 + 1000000*<rep>` | per-task seed is `SEED_BASE + 100000*degree + 1000*case + index` |
| `EM_RNG_SEED` | derived | seeds the generator behind the initial points; set it to reproduce one draw |

M = 18 mixture components, `dim P_r + 1` for |J| = 4.

**Restarts come from replicates, not from `EM_STARTS`.** A task that carries all its starts
internally writes output only if it runs to completion, which is a poor fit for a 9261-point
grid. So each replicate is one random start over the whole grid, run as its own array under
its own `SEED_BASE` and written to its own directory, and the collector takes the minimum
over whichever replicates are present. Replicates 1..k of a 5-replicate run are exactly a
valid k-start result.

Standardising is a reparametrisation — utilities shift by a constant within each choice set,
which cancels in the logit, and beta rescales — so the attainable error is unchanged and only
the iteration path differs.

## Output formats

Greedy writes one file per grid point:

    output/raw/DVD_hh_gd_<case>_<degree>_<gridindex>.txt

Five numbers: three fixed effects, the normalised zero, then **the error on line 5**.

EM writes one CSV per chunk, under one directory per replicate:

    output/raw_em_tight/rep<N>/hh_DVD_<case>d<degree>_EM_error_<start>_<end>.csv

Each has a `,0` **header row that must be skipped** — read as data it parses as a zero and
every minimum comes out exactly 0. With `EM_STARTS=1` there is one data row per grid point.

`collect.R` writes `output/tbl_apx_error2.csv`: one row per case, with `linear_greedy`,
`linear_em`, `quadratic_greedy`, `quadratic_em`. A cell with no finished grid points is `NA`.

## Normalisation — sqrt(11)

The solvers report `raw = ||rho - rho_hat||_2 / count`, dividing by the number of choice sets
rather than its square root. The paper's metric divides by `sqrt(|D|)`, and singleton menus
are never evaluated, so `|D| = count = 11`. `collect.R` therefore rescales by

    d = raw * count/sqrt(|D|) = raw * 11/sqrt(11) = raw * sqrt(11) = raw * 3.3166

This is the only scaling applied. Printed values are truncated at the third decimal, not
rounded.

## Source

- greedy: `code/approx_rum_GREEDY_hhmix_DVD.py` (argv: to-do index, case, degree)
- EM: `code/approx_rum_EM_hhmix_DVD.py` (argv: chunk index, case, degree)

Both take an optional fourth argument giving the data file, defaulting to `DVD_DATA`.
