# Fishing-site choice data

    Rscript prepare_fish_data.R

Requires the R package `mlogit` (present in the cluster's `R/4.4.1` module).

## Files

| file | contents | used by |
|---|---|---|
| `fish_individual.csv` | 4,728 rows = 1,182 individuals × 4 alternatives, long format | `tbl:in-sample_d1`, `tbl:out-sample` (they split individuals) |
| `fish_aggregated.csv` | 4 rows, characteristics averaged over individuals | `tbl:fish_apx_error1`, `tbl:fish_apx_error2`, `tbl:fish_subs` |
| `fish_empirical_cp.csv` | 4 rows, empirical choice probabilities | both fit tables |

## Source

The `Fishing` dataset from the R package `mlogit` (Croissant 2020) — originally
Thomson & Crooke (1991); see also Herriges & Kling (1999) and Cameron & Trivedi (2005)
p.464. It is loaded from the package rather than a checked-in copy so the provenance is
unambiguous.

1,182 respondents each choose one of four fishing modes. The package stores the data wide,
with one column per (field, mode) pair: `price.beach`, `catch.beach`, `price.boat`, and so
on, plus the chosen `mode` and `income`.

Each alternative is described by two characteristics:

- **price** — the cost of that fishing mode
- **catch** — the catch rate, per hour fished, for major species by mode, summed over each
  respondent's targeted species

## Construction

**Alternative numbering.** Alternatives are numbered `1=beach, 2=boat, 3=charter, 4=pier`.
The package's own column order is beach, pier, boat, charter, so the numbering is set
explicitly in the script and every output file carries an `alt_id` alongside the name.

**`fish_individual.csv`** — the wide data reshaped to long, one row per
(individual, alternative): `chid`, `alt_id`, `alt`, `chosen`, `price`, `catch`. `chosen` is
1 for the mode the respondent picked and 0 otherwise, so it sums to 1 within each `chid`.
Sorted by `chid` then `alt_id`.

**`fish_aggregated.csv`** — `price` and `catch` averaged over all 1,182 individuals, one
row per alternative. This is the input for the exhibits that work from aggregated
characteristics.

| alt_id | alt | price | catch |
|---|---|---|---|
| 1 | beach | 103.422005 | 0.24101134 |
| 2 | boat | 55.256570 | 0.17121464 |
| 3 | charter | 84.379244 | 0.62936794 |
| 4 | pier | 103.422005 | 0.16222369 |

`price` is equal for beach and pier because both are the shore-fishing price.

**`fish_empirical_cp.csv`** — the share of respondents choosing each mode:

    beach 0.113367   boat 0.353638   charter 0.382403   pier 0.150592

## Units and conventions

`price` is written **positive** and **unscaled**, and `catch` in its original units. No
sign flip, normalisation or polynomial expansion is applied. Exhibit code should set
whatever convention it needs explicitly.

## Checks

The script asserts that the reshape is well formed: 4 rows per individual, and exactly one
chosen alternative per individual, with the empirical choice probabilities summing to 1.
