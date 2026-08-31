#!/usr/bin/env Rscript
# Render output/tbl_fish_apx_error1.csv as table_latex/tbl_fish_apx_error1.tex.
#
#   Rscript code/make_tex.R [in_csv] [out_tex]
#
# Values are truncated at the third decimal, not rounded; see README.md.

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "output/tbl_fish_apx_error1.csv"
out_tex <- if (length(args) >= 2) args[2] else "table_latex/tbl_fish_apx_error1.tex"

tab <- read.csv(in_csv, stringsAsFactors = FALSE)
for (c in names(tab)[-1]) tab[[c]] <- as.numeric(tab[[c]])
val_cols <- c("linear_greedy", "linear_em", "quadratic_greedy", "quadratic_em")

ORDER <- c("1>2>3>4","1>2>4>3","1>3>2>4","1>4>2>3","2>1>3>4","2>1>4>3",
           "3>2>4>1","3>4>1>2","3>4>2>1","4>2>3>1","4>3>1>2","4>3>2>1")

fmt <- function(x) if (is.na(x)) "---" else sprintf("%.3f", floor(x * 1000) / 1000)
texify <- function(l) paste(sprintf("\\pi(%s)", strsplit(l, ">")[[1]]), collapse = "> ")

u <- tab[match(ORDER, tab$ranking), ]
stopifnot(!any(is.na(u$ranking)), all(u$representable == 0))
body <- vapply(seq_len(nrow(u)), function(i)
  sprintf("$%s$ &%s&%s& %s & %s \\\\", texify(u$ranking[i]),
          fmt(u$linear_greedy[i]), fmt(u$linear_em[i]),
          fmt(u$quadratic_greedy[i]), fmt(u$quadratic_em[i])), character(1))

r <- tab[tab$representable == 1, ]
summ <- sprintf("Linearly Representable Rankings & %s & %s & %s & %s \\\\",
                fmt(max(r$linear_greedy,    na.rm = TRUE)),
                fmt(max(r$linear_em,        na.rm = TRUE)),
                fmt(max(r$quadratic_greedy, na.rm = TRUE)),
                fmt(max(r$quadratic_em,     na.rm = TRUE)))

# Verbatim from the committed table; edit here to change the printed note.
note <- paste0(
  "\\textit{Note}: The numbers in the table show the approximation errors for each ",
  "$\\rho^{\\pi}$, where each ranking $\\pi$ is defined in the leftmost column. Alternative ",
  "numbers 1, 2, 3, 4 denote beach, boat, charter, and pier, respectively. For each ranking, ",
  "columns (1) and (2) show the approximation errors of the linear mixed-logit models ",
  "computed by the greedy algorithm and the EM algorithm, respectively. Columns (3) and (4) ",
  "show the approximation errors of the quadratic mixed-logit models calculated by each ",
  "algorithm. All numbers are truncated to three decimal places. For the greedy algorithm we ",
  "set the number of iterations to 1000. For the EM algorithm we set the number of random ",
  "initial points to 10 and the number of mixture components to 18.")

tex <- c(
  sprintf("%% Generated %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "\\begin{table}[ht]",
  "\\begin{center}",
  "  \\caption{Approximation errors to  $\\rho^{\\pi}$}\\label{tbl:fish_apx_error1}",
  "\\scalebox{0.9}{",
  "\\begin{tabular}{|c|c|c|c|c|}",
  "\\hline",
  "\\multirow{2}{*}{Ranking $\\pi$ } &\\multicolumn{2}{c|}{Linear mixed-logit}&\\multicolumn{2}{c|}{Quadratic mixed-logit} \\\\ \\cline{2-5}",
  " & Greedy  & EM    & Greedy  & EM \\\\",
  "   &  (1) &  (2) &(3) &(4) \\\\",
  "\\hline",
  "Linearly Unrepresentable Rankings &&&  &\\\\",
  body,
  "\\hline",
  summ,
  "\\hline",
  "\\end{tabular}",
  "}",
  "\\end{center}",
  "\\begin{scriptsize}",
  note,
  "\\end{scriptsize}",
  "\\end{table}")

dir.create(dirname(out_tex), showWarnings = FALSE, recursive = TRUE)
writeLines(tex, out_tex)
cat(sprintf("wrote %s (%d data rows)\n", out_tex, nrow(u)))
