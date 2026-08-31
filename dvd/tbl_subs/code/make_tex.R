#!/usr/bin/env Rscript
# Render the maximal-substitution table as table_latex/tbl_subs.tex.
#
#   Rscript code/make_tex.R
#
# Reads full precision from RAW_DIR, not the 4-decimal output/tbl_subs.csv, so the
# printed value is truncated once rather than twice. Values are truncated at the third
# decimal, not rounded.

DIM     <- as.integer(Sys.getenv("DIM", "1"))
FE      <- as.integer(Sys.getenv("FE", "0"))
raw_dir <- Sys.getenv("RAW_DIR", "output/raw")
out_tex <- Sys.getenv("OUT_TEX", "table_latex/tbl_subs.tex")

N_ALT <- 4
CASES <- list(
  c(1, 0), c(2, 0), c(3, 0),
  c(0, 1), c(2, 1), c(3, 1),
  c(0, 2), c(1, 2), c(3, 2),
  c(0, 3), c(1, 3), c(2, 3)
)

tab <- matrix(NA_real_, N_ALT, N_ALT)
for (k in seq_along(CASES)) {
  x <- CASES[[k]][1]; xr <- CASES[[k]][2]
  f <- file.path(raw_dir, sprintf("subs_d%d_fe%d_%d_%d.txt", DIM, FE, x, xr))
  if (!file.exists(f)) stop("missing ", f)
  v <- scan(f, quiet = TRUE)
  tab[xr + 1, x + 1] <- v[length(v)]
}

fmt <- function(x) if (is.na(x)) "-" else sprintf("%.3f", floor(x * 1000) / 1000)
row_tex <- function(i) sprintf("%d & %s & %s & %s & %s \\\\", i,
                               fmt(tab[i, 1]), fmt(tab[i, 2]), fmt(tab[i, 3]), fmt(tab[i, 4]))

tex <- c(
  sprintf("%% Generated %s", format(Sys.Date(), "%Y-%m-%d")),
  "\\begin{table}[ht]",
  "\\begin{center}",
  "  \\caption{Maximal substitution  of the linear mixed-logit models}\\label{tbl:subs}",
  "\\scalebox{0.9}{",
  "\\begin{tabular}{|c|c|c|c|c|}",
  " \\hline",
  "\\diagbox[width=10em, trim=l]{$j$ (drop)}{$l$ (choose)} & 1 & 2 & 3 & 4 \\\\",
  "  \\hline",
  row_tex(1), "\\hline",
  row_tex(2), "\\hline",
  row_tex(3), "\\hline",
  row_tex(4), "  \\hline",
  "\\end{tabular}",
  "}",
  "\\end{center}",
  "\\begin{scriptsize}",
  paste0("\\textit{Note}: The numbers in the table show the value of (\\ref{eq:dist.2}) ",
         "for each $j, l \\in \\{1,2,3,4\\}$ s.t. $j\\neq l$. All numbers are truncated to ",
         "three decimal places."),
  "\\end{scriptsize}",
  "\\end{table}")

dir.create(dirname(out_tex), showWarnings = FALSE, recursive = TRUE)
writeLines(tex, out_tex)
cat(sprintf("wrote %s\n", out_tex))
