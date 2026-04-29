#!/usr/bin/env Rscript
"
Tune XGBoost hyperparameters using random search CV (imbalance-safe) and save best params.

Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
Last update time: 2026-04-29

Usage:
  WeightedScore_XGB_train_Tuned_IB.R --save-model <SAVE> -o OUTCOME [--seed SEED] [--n-iter N] [--nfold FOLD] [--nrounds NROUND] [--early-stop EARLY] [--nthread NTHREAD]

Options:
  --save-model <SAVE>   Path to write tuned parameter RDS (required).
  -o OUTCOME            Outcome column name.
  --seed SEED           Random seed [default: 123].
  --n-iter N            Number of random parameter sets [default: 100].
  --nfold FOLD          Number of folds [default: 5].
  --nrounds NROUND      Max boosting rounds [default: 1000].
  --early-stop EARLY    Early stopping rounds [default: 20].
  --nthread NTHREAD     Number of parallel threads for XGBoost [default: 4].


" -> doc

suppressMessages(library(docopt))
suppressMessages(library(data.table))
suppressMessages(library(xgboost))

opts <- docopt(doc)


# -----------------------------
# Robust argument parsing
# -----------------------------

get_opt <- function(optlist, keys, default) {
  for (k in keys) {
    if (!is.null(optlist[[k]]) && length(optlist[[k]]) == 1 && nzchar(as.character(optlist[[k]]))) {
      v <- suppressWarnings(as.numeric(optlist[[k]]))
      if (!is.na(v)) return(v)
      # if not numeric (e.g., for strings), return original
      return(optlist[[k]])
    }
  }
  default
}

# Use robust getter for hyphenated flags
N_ITER     <- get_opt(opts, c("n-iter","n_iter","n.iter"), 100)
EARLY_STOP <- get_opt(opts, c("early-stop","early_stop","early.stop"), 20)

# Non-hyphenated are fine as-is, but you can still harden:
N_FOLD     <- get_opt(opts, c("nfold"), 5)
N_ROUNDS   <- get_opt(opts, c("nrounds"), 1000)

# nthread input (default 4)
N_THREAD   <- get_opt(opts, c("nthread","n-thread","n.thread"), 4)
N_THREAD   <- as.integer(ifelse(is.na(as.numeric(N_THREAD)), 4, N_THREAD))
if (is.na(N_THREAD) || N_THREAD < 1) N_THREAD <- 1L

# Seed too:
SEED       <- get_opt(opts, c("seed"), 123)
set.seed(SEED)

# Outcome name (string)
outcome_col_name <- opts$o
if (is.null(outcome_col_name) || !nzchar(outcome_col_name)) {
  stop("You must provide -o OUTCOME (the outcome column name).")
}

# # docopt converts --save-model to save_model (dash -> underscore)
save_path <- opts$save_model
if (is.null(save_path) || !is.character(save_path) || length(save_path) != 1 || !nzchar(save_path)) {
  stop("❌ --save-model must be a non-empty file path string.")
}


# Debug: what docopt gave us
cat("DEBUG keys:\n", paste(names(opts), collapse=", "), "\n", file=stderr())
cat(sprintf("Parsed -> n_iter=%s, early_stop=%s, nfold=%s, nrounds=%s, seed=%s, nthread=%s\n",
            N_ITER, EARLY_STOP, N_FOLD, N_ROUNDS, SEED, N_THREAD), file=stderr())
cat("DEBUG: opts$save_model =", save_path, "\n", file=stderr())

# -----------------------------
# Read input (stdin) robustly
# -----------------------------
df <- na.omit(read.table(file("stdin"), header = TRUE, check.names = FALSE, sep = "", quote = "", comment.char = ""))

# Outcome name must be present
outcome_col_name <- opts$o
if (is.null(outcome_col_name) || !nzchar(outcome_col_name)) {
  stop("You must provide -o OUTCOME (the outcome column name).")
}


if (!outcome_col_name %in% colnames(df)) {
  stop(sprintf("Outcome column '%s' not found in input data.", outcome_col_name))
}

# Basic columns
IID_column <- colnames(df)[1]
PGS_columns <- setdiff(colnames(df), c(IID_column, outcome_col_name))
if (length(PGS_columns) == 0) stop("No feature columns found after excluding IID and outcome.")

# -----------------------------
# Map outcome to strict 0/1
# -----------------------------

y_raw <- df[[outcome_col_name]]
if (is.numeric(y_raw)) {
  vals <- sort(unique(y_raw))
  if (all(vals %in% c(0,1))) {
    y <- as.integer(y_raw)
  } else if (length(vals) == 2) {
    y <- as.integer(y_raw == max(vals))  # larger -> 1
  } else {
    stop("Numeric outcome must be binary (0/1) or exactly two distinct values.")
  }
} else {
  fac <- factor(y_raw)
  if (nlevels(fac) != 2) stop("Outcome must have exactly two classes for binary:logistic.")
  y <- as.integer(fac) - 1L
}

# -----------------------------
# Class stats & nfold guard
# -----------------------------
tab <- table(y)
if (length(tab) != 2) stop("Outcome must contain both classes (0 and 1).")
pos <- as.integer(tab["1"]); neg <- as.integer(tab["0"])
if (is.na(pos) || is.na(neg)) stop("Outcome mapping failed; need classes 0 and 1.")
max_k <- min(pos, neg)
if (N_FOLD > max_k) {
  warning(sprintf("nfold=%d > smallest class size=%d; reducing nfold to %d.", N_FOLD, max_k, max_k))
  N_FOLD <- max_k
}
scale_pos_weight <- if (pos > 0) neg / pos else 1

# -----------------------------
# Explicit stratified folds
# -----------------------------
build_stratified_folds <- function(y, k, seed=123) {
  pos_idx <- which(y == 1)
  neg_idx <- which(y == 0)
  set.seed(seed)
  pos_idx <- sample(pos_idx); neg_idx <- sample(neg_idx)
  rr <- function(idx, k) {
    buckets <- vector("list", k)
    for (i in seq_along(idx)) {
      buckets[[ ((i-1) %% k) + 1 ]] <- c(buckets[[ ((i-1) %% k) + 1 ]], idx[i])
    }
    buckets
  }
  p_b <- rr(pos_idx, k); n_b <- rr(neg_idx, k)
  lapply(seq_len(k), function(i) c(p_b[[i]], n_b[[i]]))
}
# folds_list <- build_stratified_folds(y, N_FOLD, seed = get_num(opts$seed,123))
folds_list <- build_stratified_folds(y, N_FOLD, seed = SEED)

# -----------------------------
# Data & search grid
# -----------------------------
x_all <- as.matrix(df[, PGS_columns, drop = FALSE])
dtrain <- xgb.DMatrix(data = x_all, label = y)

# ---------- Random search grid ----------
full_grid <- expand.grid(
  eta              = c(0.01, 0.03, 0.05, 0.1, 0.2),
  max_depth        = c(3, 5, 7, 9),
  subsample        = c(0.6, 0.8, 1.0),
  colsample_bytree = c(0.6, 0.8, 1.0),
  min_child_weight = c(1, 3, 5),
  gamma            = c(0, 1, 5),
  lambda           = c(1, 3, 5),
  alpha            = c(0, 1, 5)
)
N_ITER <- min(as.integer(N_ITER), nrow(full_grid))
search_idx <- sample.int(nrow(full_grid), N_ITER)
search_grid <- full_grid[search_idx, , drop = FALSE]

# -----------------------------
# Header
# -----------------------------
header_text <- paste0(
  sprintf("Performing random search over %d parameter combinations...\n", N_ITER),
  sprintf("Using %d-fold CV (stratified), up to %d rounds, early stopping = %d, nthread = %d\n",
          N_FOLD, N_ROUNDS, EARLY_STOP, N_THREAD),
  sprintf("Positives=%d, Negatives=%d, PosRate=%.4f, scale_pos_weight=%.3f\n\n",
          pos, neg, pos/(pos+neg), scale_pos_weight)
)
cat(header_text)


# -----------------------------
# CV runner (fallback without ES on error)
# -----------------------------
run_cv_safely <- function(params, dtrain, N_ROUNDS, folds_list, EARLY_STOP) {
  res <- try({
    xgb.cv(
      params = params,
      data = dtrain,
      nrounds = N_ROUNDS,
      folds = folds_list,
      early_stopping_rounds = EARLY_STOP,
      maximize = TRUE,
      verbose = 0
    )
  }, silent = TRUE)
  if (inherits(res, "try-error")) {
    # fallback without early stopping; choose best manually later
    res <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = N_ROUNDS,
      folds = folds_list,
      maximize = TRUE,
      verbose = 0
    )
    attr(res, "no_es") <- TRUE
  } else {
    attr(res, "no_es") <- FALSE
  }
  res
}

# -----------------------------
# Random search loop (optimize AUCPR)
# -----------------------------

best_score   <- -Inf   # PR-AUC target
best_params  <- list()
best_nrounds <- 0

for (i in seq_len(nrow(search_grid))) {
  p <- search_grid[i,]
  cat(sprintf(
    "[%d/%d] eta=%.3f depth=%d subsample=%.2f colsample=%.2f min_child=%d gamma=%.2f lambda=%.2f alpha=%.2f\n",
    i, nrow(search_grid), p$eta, p$max_depth, p$subsample, p$colsample_bytree,
    p$min_child_weight, p$gamma, p$lambda, p$alpha
  ))

  params <- list(
    objective        = "binary:logistic",
    eval_metric      = c("aucpr","auc"),  # tune on PR-AUC; also log ROC AUC
    eta              = p$eta,
    max_depth        = p$max_depth,
    subsample        = p$subsample,
    colsample_bytree = p$colsample_bytree,
    min_child_weight = p$min_child_weight,
    gamma            = p$gamma,
    lambda           = p$lambda,
    alpha            = p$alpha,
    scale_pos_weight = scale_pos_weight,
    nthread          = N_THREAD
  )

  cv <- run_cv_safely(params, dtrain, N_ROUNDS, folds_list, EARLY_STOP)

  # PR-AUC column name in evaluation_log
  if (!"test_aucpr_mean" %in% names(cv$evaluation_log)) {
    warning("No test_aucpr_mean in evaluation_log; skipping combo.")
    next
  }
  pr <- cv$evaluation_log$test_aucpr_mean
  # Skip if all NA/Inf
  if (all(!is.finite(pr))) {
    warning(sprintf("All PR-AUC values NA for combo %d — skipping.", i))
    next
  }
  finite_idx <- which(is.finite(pr))
  best_idx <- finite_idx[ which.max(pr[finite_idx]) ]
  score <- pr[best_idx]

  best_iter <- cv$best_iteration
  if (is.null(best_iter) || !is.finite(best_iter) || best_iter <= 0) best_iter <- best_idx

  roc_score <- NA_real_
  if ("test_auc_mean" %in% names(cv$evaluation_log)) {
    ra <- cv$evaluation_log$test_auc_mean
    if (is.finite(ra[best_iter])) roc_score <- ra[best_iter]
  }

  cat(sprintf("   -> AUCPR = %.5f at round %d", score, best_iter))
  if (is.finite(roc_score)) cat(sprintf("  |  ROC AUC = %.5f", roc_score))
  if (isTRUE(attr(cv, "no_es"))) cat("  |  (no early stopping fallback)")
  cat("\n")

  if (!is.na(score) && score > best_score) {
    best_score   <- score
    best_params  <- params
    best_nrounds <- best_iter
  }
}

# ---------- SAFE SAVE BLOCK ----------
# Normalize path & ensure directory exists
save_path <- normalizePath(save_path, mustWork = FALSE)
dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)

cat("SAVE CHECK:\n", file=stderr())
cat("  getwd()        = ", getwd(), "\n", file=stderr())
cat("  save_path      = ", save_path, "\n", file=stderr())
cat("  typeof(path)   = ", typeof(save_path), " length=", length(save_path), "\n", file=stderr())
cat("  write access?  = ", file.access(dirname(save_path), 2), " (0 means OK)\n", file=stderr())

stopifnot(is.character(save_path), length(save_path) == 1, nzchar(save_path))

# Probe write to catch path issues early
probe_path <- sub("\\.rds$", ".probe.rds", save_path)
ok <- tryCatch({ saveRDS(list(ok=TRUE), file=probe_path); TRUE },
               error=function(e) { cat("PROBE ERROR:", conditionMessage(e), "\n", file=stderr()); FALSE })
if (!ok) stop("Probe write failed, cannot save to: ", probe_path)

# Save tuned results (RDS)
tuned_obj <- list(
  best_params       = best_params,
  best_nrounds      = best_nrounds,
  best_aucpr        = best_score,
  features_used     = PGS_columns,
  nfold             = N_FOLD,
  nrounds           = N_ROUNDS,
  early_stop        = EARLY_STOP,
  n_iter            = N_ITER,
  nthread           = N_THREAD,   
  seed              = SEED,
  outcome_col       = outcome_col_name,
  pos               = pos,
  neg               = neg,
  scale_pos_weight  = scale_pos_weight
)

tryCatch({
  saveRDS(tuned_obj, file = save_path)
  cat("Saved model to: ", save_path, "\n", file=stderr())
}, error=function(e) {
  cat("SAVE ERROR: ", conditionMessage(e), "\n", file=stderr())
  stop("saveRDS failed for path: ", save_path)
})

# ---------- Write .summary.txt ----------
header_text <- paste0(
  sprintf("Performing random search over %d parameter combinations...\n", N_ITER),
  sprintf("Using %d-fold CV (stratified), up to %d rounds, early stopping = %d\n",
          N_FOLD, N_ROUNDS, EARLY_STOP),
  sprintf("Positives=%d, Negatives=%d, PosRate=%.4f, scale_pos_weight=%.3f\n\n",
          pos, neg, pos/(pos+neg), scale_pos_weight)
)

summary_text <- paste0(
  "\n====================\n",
  sprintf("Best AUCPR: %.5f\n", best_score),
  sprintf("Best nrounds: %d\n", best_nrounds),
  "Best parameters:\n",
  paste(capture.output(print(best_params)), collapse = "\n"),
  "\nRun settings:\n",
  sprintf("Outcome: %s\nSeed: %d\nN_iter: %d\nN_fold: %d\nN_rounds(max): %d\nEarly_stop: %d\n",
          outcome_col_name, SEED, N_ITER, N_FOLD, N_ROUNDS, EARLY_STOP),  
  sprintf("Pos: %d  Neg: %d  PosRate: %.6f  scale_pos_weight: %.3f\n",
          pos, neg, pos/(pos+neg), scale_pos_weight),
  "Features used: ",
  paste(PGS_columns, collapse = ", "),
  "\n====================\n"
)
summary_file <- sub("\\.rds$", ".summary.txt", save_path)
writeLines(paste0(header_text, summary_text), con = summary_file)
cat(summary_text)


