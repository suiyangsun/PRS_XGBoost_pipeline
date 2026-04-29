#!/usr/bin/env Rscript
"
Apply saved incremental XGBoost models to TEST and output scores.

Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
Last update time: 2026-04-29

Reads TEST from stdin. Loads an RDS created by the incremental training script
(--save-model). For each saved model (feature subset), predicts probabilities
and outputs columns xgb_score_<k>, where <k> is the number of features used.

Usage:
  WeightedScore_XGB_test_score.R --models <RDS> [-o OUTCOME]
  WeightedScore_XGB_test_score.R (-h | --help)

Options:
  --models <RDS>   Path to incremental models RDS (required).
  -o OUTCOME       Outcome column name in TEST (optional; if found, it is included).
" -> doc

suppressMessages(library(docopt))
suppressMessages(library(data.table))
suppressMessages(library(xgboost))

opts <- docopt(doc)

# ---- args
if (is.null(opts$models) || !nzchar(opts$models)) {
  stop("--models <RDS> is required.")
}

# ---- load models RDS
inc <- readRDS(opts$models)
if (is.null(inc$models)) stop("Models list not found in RDS.")
models_list <- inc$models
available_idxs <- which(vapply(models_list, function(x) !is.null(x) && !is.null(x$model), logical(1)))
if (length(available_idxs) == 0) stop("No non-NULL models found in RDS.")

feature_order <- inc$feature_order
saved_outcome <- inc$outcome_col

# ---- read TEST from stdin
df <- tryCatch(
  suppressWarnings(data.table::fread("stdin", data.table = FALSE)),
  error = function(e) {
    suppressWarnings(read.table(file("stdin"), header=TRUE, check.names=FALSE, sep="", quote="", comment.char=""))
  }
)
if (nrow(df) == 0) stop("No input rows read from stdin.")
colnames(df) <- trimws(colnames(df))

IID_col <- colnames(df)[1]

# ---- optional outcome inclusion
outcome_col <- NULL
if (!is.null(opts$o) && nzchar(opts$o)) {
  if (opts$o %in% colnames(df)) {
    outcome_col <- opts$o
  } else {
    idx <- which(tolower(colnames(df)) == tolower(trimws(opts$o)))
    if (length(idx) == 1) outcome_col <- colnames(df)[idx]
  }
} else if (!is.null(saved_outcome) && nzchar(saved_outcome)) {
  if (saved_outcome %in% colnames(df)) {
    outcome_col <- saved_outcome
  } else {
    idx <- which(tolower(colnames(df)) == tolower(trimws(saved_outcome)))
    if (length(idx) == 1) outcome_col <- colnames(df)[idx]
  }
}
if (!is.null(outcome_col) && !outcome_col %in% colnames(df)) outcome_col <- NULL

# ---- prepare output
out <- data.frame(IID = df[[IID_col]], stringsAsFactors = FALSE)
if (!is.null(outcome_col)) out$outcome <- df[[outcome_col]]

# ---- score each saved model (adds one score column per model)
missing_any <- FALSE
scored_count <- 0L

for (k in available_idxs) {
  entry <- models_list[[k]]
  booster  <- entry$model
  req_feats <- entry$features

  # ensure required features present in TEST
  missing_feats <- setdiff(req_feats, colnames(df))
  if (length(missing_feats) > 0) {
    missing_any <- TRUE
    cat(sprintf("WARN: Skipping model index %s (needs %d features); missing in TEST: %s\n",
                as.character(k), length(req_feats), paste(missing_feats, collapse=", ")), file=stderr())
    next
  }

  # build matrix in the exact order trained
  X <- as.matrix(df[, req_feats, drop=FALSE])

  # predict using all trees (model already trained to its own best rounds)
  preds <- predict(booster, xgb.DMatrix(X))

  # name column by feature count (e.g., xgb_score_2, xgb_score_3, ...)
  col_name <- paste0("xgb_score_", length(req_feats))
  out[[col_name]] <- as.numeric(preds)
  scored_count <- scored_count + 1L
}

if (scored_count == 0L) {
  stop("No models could be scored on TEST (features missing?). Check warnings in stderr.")
}
if (missing_any) {
  cat("NOTE: Some models were skipped due to missing TEST features.\n", file=stderr())
}

# ---- write scores to stdout
write.table(out, "", row.names=FALSE, col.names=TRUE, quote=FALSE)

