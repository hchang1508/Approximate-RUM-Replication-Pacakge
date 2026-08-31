#!/usr/bin/env Rscript
# Shared renderer for the two fit tables. Both come from the same run, so both are
# rendered here; see collect_fit_table.R, which shares its inputs the same way.
#
# Sourced with:  args <- c(in_csv, out_tex, mode)   mode = "in" | "out"
#
# Values are truncated at the third decimal, not rounded.

stopifnot(exists("args"), length(args) >= 3)
in_csv <- args[1]; out_tex <- args[2]; mode <- args[3]
stopifnot(mode %in% c("in", "out"))

tab <- read.csv(in_csv, stringsAsFactors = FALSE)

MODELS <- c("our_method", "mnl", "nest_charter", "nest_boat",
            "mixed_logit_lognormal", "mnl_FE", "mixed_logit_lognormal_FE")
LABELS <- list(
  our_method               = c("Our method", ""),
  mnl                      = c("Multinomial logit", ""),
  nest_charter             = c("Nested logit", "(charter and others)"),
  nest_boat                = c("Nested logit", "(boat and others)"),
  mixed_logit_lognormal    = c("Mixed-logit with", "log normal distribution"),
  mnl_FE                   = c("Multinomial logit", "with fixed effects"),
  mixed_logit_lognormal_FE = c("Mixed-logit with log normal", "distribution and fixed effects")
)
VAL <- c("beach", "boat", "charter", "pier", "pred_error")
SD  <- paste0(VAL, "_sd")

stopifnot(setequal(tab$model, MODELS))
fmt <- function(x) sprintf("%.3f", floor(x * 1000) / 1000)

body <- unlist(lapply(MODELS, function(m) {
  r <- tab[tab$model == m, ]
  c(sprintf("%s & %s \\\\", LABELS[[m]][1],
            paste(vapply(VAL, function(v) fmt(r[[v]]), character(1)), collapse = " & ")),
    sprintf("%s & %s \\\\", LABELS[[m]][2],
            paste(vapply(SD, function(v) sprintf("(%s)", fmt(r[[v]])), character(1)), collapse = " & ")),
    "\\hline")
}), use.names = FALSE)

n_rep <- unique(tab$n)
if (mode == "in") {
  caption <- "In-Sample Fit"; label <- "tbl:in-sample_d1"
  note <- paste0(
    "Table \\ref{", label, "} summarizes the in-sample fit of different models. The row ",
    "``our method'' presents choice probabilities predicted by the four-mixture ",
    "mixed-logit model and the prediction error. The remaining rows present in-sample ",
    "predicted choice probabilities and prediction errors obtained by standard models. ",
    "In parentheses are standard deviations obtained by repeating the same analyses ",
    n_rep, " times. All numbers are truncated to three decimal places.")
} else {
  caption <- "Out-of-Sample Fit"; label <- "tbl:out-sample"
  note <- paste0(
    "Table \\ref{", label, "} summarizes the out-of-sample fit of different models. The row ",
    "``our method'' presents choice probabilities predicted by the four-mixture ",
    "mixed-logit model and the prediction error (\\ref{eq:dist.1}). The remaining rows ",
    "present out-of-sample predicted choice probabilities and prediction errors obtained ",
    "by standard models. In parentheses are standard deviations obtained by repeating the ",
    "same analyses ", n_rep, " times. All numbers are truncated to three decimal places.")
}

tex <- c(
  sprintf("%% Generated %s", format(Sys.Date(), "%Y-%m-%d")),
  "\\begin{table}[ht]",
  sprintf("\\caption{%s}", caption),
  "\\begin{center}",
  "\\scalebox{0.85}{",
  "\\begin{tabular}{|c|c|c|c|c|c|}",
  "  \\hline",
  "  \\multirow{2}{*}{Model } &\\multicolumn{4}{c|}{Choice probabilities}&\\multicolumn{1}{c|}{Prediction error} \\\\ \\cline{2-6}",
  " & Beach  & Boat    & Charter  & Pier & \\\\",
  "   &  (1) &  (2) & (3) & (4) & (5)\\\\",
  "  \\hline",
  body,
  "\\end{tabular}",
  "}",
  "\\end{center}",
  sprintf("\\label{%s}", label),
  "\\begin{tablenotes}",
  "\\item \\begin{scriptsize} \\textit{Note}:",
  note,
  "\\end{scriptsize}",
  "\\end{tablenotes}",
  "\\end{table}")

dir.create(dirname(out_tex), showWarnings = FALSE, recursive = TRUE)
writeLines(tex, out_tex)
cat(sprintf("wrote %s (%d models)\n", out_tex, length(MODELS)))
