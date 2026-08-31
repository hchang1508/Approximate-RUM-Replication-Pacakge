# DVD data

    Rscript prepare_dvd_data.R

## Files

| file | contents | generated |
|---|---|---|
| `farias2009_table2.csv` | source table, all 15 products, raw units | no — transcribed by hand |
| `dvd_J4_K2.csv` | \|J\|=4, K=2 — main text Section 5 | yes |
| `dvd_J10_K3.csv` | \|J\|=10, K=3 — appendix `sec:J10K2` | yes |

## Source

`farias2009_table2.csv` is Table 2 of Farias, Jagabathula & Shah (2009),
<https://arxiv.org/pdf/0910.0063>, reproduced as `tab:amazon_dvd` in the paper; the
underlying data is from Rusmevichientong et al. (2010). It lists 15 DVD products with three fields: `price`, `avg_price_per_disc` and
`total_helpful_votes`.

## Construction

**Product selection.** Rows are ordered by descending `price`, so the paper's "first four
rows" and "first ten rows" are the four and ten most expensive products.

| dataset | products | characteristics | K |
|---|---|---|---|
| `dvd_J4_K2.csv` | 1–4 | `avg_price_per_disc`, `total_helpful_votes` | 2 |
| `dvd_J10_K3.csv` | 1–10 | `avg_price_per_disc`, `total_helpful_votes`, `price` | 3 |

**Scaling.** Each field is divided by a fixed constant:

    price / 100        avg_price_per_disc / 10        total_helpful_votes / 1000

This is a diagonal positive rescaling of the characteristics, which is absorbed into beta
and so leaves every theoretical quantity unchanged. It is applied to keep the
characteristics on comparable numeric scales.

**Output format.** One row per alternative. First column `alt_id` (the product ID),
remaining columns the K characteristics, in the order given in the table above.

## Property of the data

`prepare_dvd_data.R` checks affine independence directly, since it is the condition the
paper's main theorem turns on. The characteristics are affinely independent exactly when
`[1 X]` has rank \|J\|:

| dataset | rank `[1 X]` | needed | affinely independent? |
|---|---|---|---|
| `dvd_J4_K2.csv` | 3 | 4 | no — K = 2 < 3 = \|J\|−1 |
| `dvd_J10_K3.csv` | 4 | 10 | no — K = 3 < 9 = \|J\|−1 |

Both fail, which is what the application is built to illustrate. The script asserts these
ranks, so a change to the source table or the selection cannot pass unnoticed.
