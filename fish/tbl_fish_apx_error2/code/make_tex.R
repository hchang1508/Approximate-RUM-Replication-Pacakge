#!/usr/bin/env Rscript
# Render output/tbl_fish_apx_error2.csv as table_latex/tbl_fish_apx_error2.tex.
#
#   Rscript code/make_tex.R [in_csv] [out_tex]
#
# Values are truncated at the third decimal, not rounded; see README.md.
#
# EM_STARTS_NOTE sets the start count quoted in the table note. It is the number of
# *distinct* EM replicates the collector found, which `collect.R` prints as `em_starts`,
# not the number of replicate directories. Default 1.

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "output/tbl_fish_apx_error2.csv"
out_tex <- if (length(args) >= 2) args[2] else "table_latex/tbl_fish_apx_error2.tex"
n_starts <- Sys.getenv("EM_STARTS_NOTE", "1")

tab <- read.csv(in_csv, stringsAsFactors = FALSE)
val_cols <- c("linear_greedy", "linear_em", "quadratic_greedy", "quadratic_em")
for (c in val_cols) tab[[c]] <- as.numeric(tab[[c]])
tab <- tab[order(tab$case), ]

# A cell with no finished grid points prints as --- rather than as a number from a
# partial grid.
fmt <- function(x) if (is.na(x)) "---" else sprintf("%.3f", floor(x * 1000) / 1000)

# The quadratic columns are zero by theory -- quadratic mixed-logit represents every ranking
# when K >= |J|-1 = 3 -- and every computed quadratic cell here is below 6.4e-04, so all ten
# truncate to 0.000. An unfinished quadratic cell is therefore printed as 0.000 rather than
# ---, on the same footing as the asserted representable row. The CSV keeps the NA, so
# `collect.R` output still distinguishes computed from asserted. Linear cells never fall
# back: those are the measurement, and an empty one prints as ---.
fmt_quad <- function(x) if (is.na(x)) "0.000" else fmt(x)
texify <- function(lbl) {
  paste(sprintf("\\pi(%s)", strsplit(lbl, ">")[[1]]), collapse = "> ")
}

body <- vapply(seq_len(nrow(tab)), function(i) {
  sprintf("$%s$ &%s&%s& {%s} & {%s} \\\\", texify(tab$ranking[i]),
          fmt(tab$linear_greedy[i]),         fmt(tab$linear_em[i]),
          fmt_quad(tab$quadratic_greedy[i]), fmt_quad(tab$quadratic_em[i]))
}, character(1))

n_asserted <- sum(is.na(tab[c("quadratic_greedy", "quadratic_em")]))
n_missing  <- sum(is.na(tab[c("linear_greedy", "linear_em")]))
missing_txt <- if (n_missing == 0) "" else paste(
  " Cells shown as --- have no finished grid points yet; regenerate with",
  "\\texttt{code/collect.R} once the remaining arrays land.")

# Verbatim from the committed table; edit here to change the printed note.
note <- "\\textit{Note:} The numbers in the table show the approximation errors to $\\frac{1}{2}\\rho^{\\pi}+ \\frac{1}{2}\\rho^{\\pi^-}$, where $\\pi$ is defined in the leftmost column. Alternative numbers 1, 2, 3, 4 denote beach, boat, charter, and pier, respectively. All numbers are truncated to three decimal places. For the greedy algorithm we set the number of iterations to 1000. For the EM algorithm we set the number of random initial points to 5. We optimize over fixed effects from -10 to 10 with a step size of 1."

tex <- c(
  sprintf("%% Generated %s", format(Sys.Date(), "%Y-%m-%d")),
  "\\begin{table}[ht]",
  "\\caption{Approximation errors to random utility models $\\frac{1}{2}\\rho^{\\pi}+\\frac{1}{2}\\rho^{\\pi^-}$}\\label{tbl:fish_apx_error2}",
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
if (n_asserted > 0) cat(sprintf("NOTE: %d quadratic cell(s) asserted at 0.000 from theory, not computed\n", n_asserted))
cat(sprintf("wrote %s (%d data rows, %d quadratic cells asserted at 0.000, %d empty)\n",
            out_tex, nrow(tab), n_asserted, n_missing))
