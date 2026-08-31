#!/usr/bin/env Rscript
# Render output/tbl_apx_error2.csv as table_latex/tbl_apx_error2.tex.
#
#   Rscript code/make_tex.R [in_csv] [out_tex]
#
# Values are truncated at the third decimal, not rounded; see README.md.

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "output/tbl_apx_error2.csv"
out_tex <- if (length(args) >= 2) args[2] else "table_latex/tbl_apx_error2.tex"

tab <- read.csv(in_csv, stringsAsFactors = FALSE)
trunc3 <- function(x) sprintf("%.3f", floor(x * 1000) / 1000)

texify <- function(lbl) {
  paste(sprintf("\\pi(%s)", strsplit(lbl, ">")[[1]]), collapse = "> ")
}

body <- vapply(seq_len(nrow(tab)), function(i) {
  sprintf("$%s$ &%s& %s& {%s} & {%s} \\\\", texify(tab$ranking[i]),
          trunc3(tab$linear_greedy[i]),    trunc3(tab$linear_em[i]),
          trunc3(tab$quadratic_greedy[i]), trunc3(tab$quadratic_em[i]))
}, character(1))

# Verbatim from the committed table; edit here to change the printed note.
note <- "\\textit{Note:} The numbers in the table show the approximation errors to $\\frac{1}{2}\\rho^{\\pi}+ \\frac{1}{2}\\rho^{\\pi^-}$, where $\\pi$ is defined in the leftmost column. All numbers are truncated to three decimal places. For the greedy algorithm we set the number of iterations to 1000. For the EM algorithm we set the number of random initial points to 5. We optimize over fixed effects from -10 to 10 with a step size of 1. "

tex <- c(
  sprintf("%% Generated %s", format(Sys.Date(), "%Y-%m-%d")),
  "\\begin{table}[ht]",
  "\\caption{ Approximation errors to random utility models $\\frac{1}{2}\\rho^{\\pi}+\\frac{1}{2}\\rho^{\\pi^-}$}\\label{tbl:apx_error2}",
  "\\label{tab:representability_macro}",
  "\\begin{center}",
  "\\scalebox{0.9}{",
  "\\begin{tabular}{|c|c|c|c|c|}",
  "\\hline",
  "\\multirow{2}{*}{Ranking $\\pi$ } &\\multicolumn{2}{c|}{Linear mixed-logit}&\\multicolumn{2}{c|}{Quadratic mixed-logit} \\\\ \\cline{2-5}",
  " & Greedy  & EM    & Greedy  & EM \\\\ ",
  "   &  (1) &  (2) &(3) &(4) \\\\  ",
  "\\hline",
  "Linearly unrepresentable rankings &&&&\\\\",
  body,
  "\\hline",
  "Linearly representable rankings & 0.000 & 0.000 & {0.000}& {0.000}\\\\",
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
cat(sprintf("wrote %s (%d data rows)\n", out_tex, nrow(tab)))
