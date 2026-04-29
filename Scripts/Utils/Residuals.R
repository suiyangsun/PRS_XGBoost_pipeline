#!/usr/bin/env Rscript
############################################
# Residualization script for PRS
# Author: Yang Sui
# Updated: 2026-04-29
############################################
suppressMessages({
  library(docopt)
  library(data.table)
})
doc <- "
Extract residuals from linear regression.
Usage:
  Residuals.R -f <formula> -t <title>
Options:
  -f <formula>   Regression formula (e.g. PRS ~ PC1 + PC2) (required)
  -t <title>     Output column name for residuals (required)
"
opts <- docopt(doc)
reg.formula <- opts$f
title <- opts$t
# -----------------------------
# Check required arguments
# -----------------------------
if (is.null(opts$f)) {
  stop("ERROR: -f (regression formula) is required")
}
if (is.null(opts$t)) {
  stop("ERROR: -t (output column name) is required")
}
# -----------------------------
# Read data
# -----------------------------
df <- fread("stdin")
if (nrow(df) == 0) {
  stop("ERROR: input is empty")
}
# -----------------------------
# Check formula variables
# -----------------------------
vars <- all.vars(as.formula(reg.formula))
if (!all(vars %in% colnames(df))) {
  missing_vars <- vars[!vars %in% colnames(df)]
  stop(paste("ERROR: variables not found in input data:", paste(missing_vars, collapse=", ")))
}
# -----------------------------
# Remove NA
# -----------------------------
n_before <- nrow(df)
df <- na.omit(df)
n_after <- nrow(df)
cat(sprintf("[Residuals] Samples removed due to NA: %d\n", n_before - n_after), file=stderr())
# -----------------------------
# Fit model
# -----------------------------
cat("[Residuals] Fitting model:", reg.formula, "\n", file=stderr())
fit <- lm(as.formula(reg.formula), data=df)
# -----------------------------
# Extract residuals
# -----------------------------
df[[title]] <- residuals(fit)
cat(sprintf("[Residuals] Residuals added as column: %s\n", title), file=stderr())
# -----------------------------
# Output
# -----------------------------
fwrite(df, file="", sep="\t")
cat("[Residuals] Done!\n", file=stderr())