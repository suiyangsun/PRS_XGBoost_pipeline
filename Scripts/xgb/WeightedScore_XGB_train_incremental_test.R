#!/usr/bin/env Rscript
"
Train incremental XGBoost models on TRAIN using tuned params.
For p PGS features, this fits p-1 models: i = 2..p. Each subset (first i features)
gets its own best nrounds via stratified k-fold CV (tunes on AUCPR).

Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
Last update time: 2026-04-29

Usage:
  XGB_fit_train_incremental_test.R --tuned <RDS> [-o OUTCOME] [--seed SEED] [--nfold FOLD] [--nrounds NROUND] [--early-stop EARLY] [--nthread NTHREAD] --save-model <SAVE>
  XGB_fit_train_incremental_test.R (-h | --help)

Options:
  --tuned <RDS>          Path to tuned RDS from tuning script (required).
  -o OUTCOME             Outcome column name (optional; falls back to tuned$outcome_col).
  --seed SEED            Random seed [default: 123].
  --nfold FOLD           CV folds per subset (stratified) [default: 5].
  --nrounds NROUND       Max boosting rounds for per-subset CV [default: 2000].
  --early-stop EARLY     Early stopping rounds for per-subset CV [default: 30].
  --nthread NTHREAD      Number of parallel threads for XGBoost [default: 4].
  --save-model <SAVE>    Save all trained incremental models as a single RDS (required).
" -> doc

suppressMessages(library(docopt))
suppressMessages(library(data.table))
suppressMessages(library(xgboost))

opts <- docopt(doc)

# -----------------------------
# Robust parsing helpers
# -----------------------------
# Modified get_opt to handle non-numeric/missing tuned defaults correctly
get_opt <- function(optlist, keys, default) {
  for (k in keys) {
    if (!is.null(optlist[[k]]) && length(optlist[[k]]) == 1 && nzchar(as.character(optlist[[k]]))) {
      v <- suppressWarnings(as.numeric(optlist[[k]]))
      if (!is.na(v)) return(v)
      return(optlist[[k]])
    }
  }
  default
}

# -----------------------------
# Initial Setup (Independent of Tuned File)
# -----------------------------
SEED      <- get_opt(opts, c("seed"), 123)
NFOLD_CMD <- get_opt(opts, c("nfold"), 5) # Store N_FOLD command line value
set.seed(SEED)

if (is.null(opts$tuned) || !nzchar(opts$tuned)) stop("--tuned <RDS> is required.")
save_models_path <- opts$save_model

# -----------------------------
# Load tuned object and extract defaults
# -----------------------------
tuned <- readRDS(opts$tuned)
if (is.null(tuned$best_params))
  stop("Tuned RDS missing best_params.")

best_params           <- tuned$best_params
tuned_outcome         <- tuned$outcome_col
features_from_tuned   <- tuned$features_used
tuned_scale_wt        <- tuned$scale_pos_weight
tuned_nrounds_default <- if (!is.null(tuned$nrounds) && tuned$nrounds > 0) as.numeric(tuned$nrounds) else 2000
tuned_early_stop_default <- if (!is.null(tuned$early_stop) && tuned$early_stop > 0) as.numeric(tuned$early_stop) else 30
tuned_nthread_default <- if (!is.null(best_params$nthread)) {
  as.integer(best_params$nthread)
} else if (!is.null(tuned$nthread)) {
  as.integer(tuned$nthread)
} else {
  4L
}


# -----------------------------
# HYBRID SOURCING LOGIC (Using CMD > Tuned RDS > Hardcoded Default)
# -----------------------------
MAX_ROUNDS <- get_opt(opts, c("nrounds"), tuned_nrounds_default)
EARLY_STOP <- get_opt(opts, c("early-stop", "early_stop", "early.stop"), tuned_early_stop_default)

# Use N_FOLD_CMD as input default, but still rely on class size capping
NFOLD <- NFOLD_CMD

N_THREAD   <- as.integer(get_opt(opts, c("nthread","n-thread","n.thread"), tuned_nthread_default))
if (is.na(N_THREAD) || N_THREAD < 1) N_THREAD <- 1L

# -----------------------------
# Read TRAIN (stdin) robustly
# -----------------------------
df <- tryCatch(
  suppressWarnings(data.table::fread("stdin", data.table = FALSE)),
  error = function(e) {
    suppressWarnings(read.table(file("stdin"), header=TRUE, check.names=FALSE, sep="", quote="", comment.char=""))
  }
)
if (nrow(df) == 0) stop("No input rows read from stdin.")
colnames(df) <- trimws(colnames(df))

# Outcome col (case/space-insensitive)
outcome_col <- if (!is.null(opts$o) && nzchar(opts$o)) opts$o else tuned_outcome
if (is.null(outcome_col) || !nzchar(outcome_col))
  stop("Outcome column not provided and not present in tuned RDS.")
if (!outcome_col %in% colnames(df)) {
  idx <- which(tolower(colnames(df)) == tolower(trimws(outcome_col)))
  if (length(idx) == 1) {
    outcome_col <- colnames(df)[idx]
    cat(sprintf("NOTE: mapped outcome to '%s'\n", outcome_col), file=stderr())
  } else {
    cat("Available columns:\n", paste(colnames(df), collapse=", "), "\n", file=stderr())
    stop(sprintf("Outcome column '%s' not found in input.", outcome_col))
  }
}


# Identify PGS features (order matters)
IID_col <- colnames(df)[1]
PGS_columns <- if (!is.null(features_from_tuned) && length(features_from_tuned) > 0) {
  intersect(features_from_tuned, colnames(df))
} else {
  setdiff(colnames(df), c(IID_col, outcome_col))
}
if (length(PGS_columns) < 2) stop("Need at least 2 PGS features to make incremental models.")

# Map outcome to 0/1
y_raw <- df[[outcome_col]]
if (is.numeric(y_raw)) {
  vals <- sort(unique(y_raw))
  if (all(vals %in% c(0,1))) y <- as.integer(y_raw)
  else if (length(vals) == 2) y <- as.integer(y_raw == max(vals))
  else stop("Outcome must be binary (two unique values).")
} else {
  f <- factor(y_raw)
  if (nlevels(f) != 2) stop("Outcome must have exactly two classes.")
  y <- as.integer(f) - 1L
}
tab <- table(y); if (length(tab) != 2) stop("Outcome must contain both classes.")
pos <- as.integer(tab["1"]); neg <- as.integer(tab["0"])
scale_pos_weight <- if (!is.null(tuned_scale_wt)) tuned_scale_wt else if (pos>0) neg/pos else 1
best_params$scale_pos_weight <- if (is.null(best_params$scale_pos_weight)) scale_pos_weight else best_params$scale_pos_weight

# Ensure nthread is set from CLI/tuned
best_params$nthread <- N_THREAD

# -----------------------------
# Build stratified folds ONCE (deterministic)
# -----------------------------
build_stratified_folds <- function(y, k, seed=123) {
  set.seed(seed)
  pos_idx <- sample(which(y == 1))
  neg_idx <- sample(which(y == 0))
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
# cap nfold by smallest class size (avoid empty-class folds)
NFOLD <- min(NFOLD, pos, neg)
if (NFOLD < 2) stop("nfold too small for class counts; need at least 2.")
folds_cv <- build_stratified_folds(y, NFOLD, seed = SEED)

# Ensure params have objective/eval_metric/scale_pos_weight
best_params$objective <- "binary:logistic"
best_params$eval_metric <- c("aucpr","auc")


# -----------------------------
# Prepare data & outputs
# -----------------------------
X_all <- as.matrix(df[, PGS_columns, drop=FALSE])

out <- data.frame(
  IID = df[[IID_col]],
  outcome = df[[outcome_col]],
  stringsAsFactors = FALSE
)

models <- vector("list", length(PGS_columns))
names(models) <- paste0("model_", seq_along(PGS_columns))

cat(sprintf("Incremental training with %d PGS features; NFOLD=%d, MAX_ROUNDS=%d, ES=%d, nthread=%d\n",
            length(PGS_columns), NFOLD, MAX_ROUNDS, EARLY_STOP, N_THREAD), file=stderr())


# -----------------------------
# Incremental loop: i = 2..P
# -----------------------------
for (i in 2:length(PGS_columns)) {
  sel_idx  <- 1:i
  sel_cols <- PGS_columns[sel_idx]
  dtrain_i <- xgb.DMatrix(data = X_all[, sel_idx, drop=FALSE], label = y)

  # per-subset CV to choose best nrounds (optimize AUCPR)
  cv <- try({
    xgb.cv(params = best_params,
           data   = dtrain_i,
           nrounds = MAX_ROUNDS,
           folds   = folds_cv,
           early_stopping_rounds = EARLY_STOP,
           maximize = TRUE,
           verbose  = 0)
  }, silent = TRUE)

  # fallback if early stopping crashes
  if (inherits(cv, "try-error")) {
    cv <- xgb.cv(params = best_params,
                 data   = dtrain_i,
                 nrounds = MAX_ROUNDS,
                 folds   = folds_cv,
                 maximize = TRUE,
                 verbose  = 0)
    attr(cv, "no_es") <- TRUE
  } else {
    attr(cv, "no_es") <- FALSE
  }

  # Safely extract best round using PR-AUC
  pr <- cv$evaluation_log$test_aucpr_mean
  finite_idx <- which(is.finite(pr))
  if (length(finite_idx) == 0) {
    best_nrounds_subset <- 200L # Final extreme fallback if all PR-AUCs are NA
  } else {
    best_idx <- finite_idx[ which.max(pr[finite_idx]) ]
    best_nrounds_subset <- cv$best_iteration
    if (is.null(best_nrounds_subset) || !is.finite(best_nrounds_subset) || best_nrounds_subset <= 0)
      best_nrounds_subset <- best_idx
  }

  # Train final model on full TRAIN for this subset
  final_model <- xgb.train(
    params   = best_params,
    data     = dtrain_i,
    nrounds  = best_nrounds_subset,
    watchlist= list(train = dtrain_i),
    verbose  = 0
  )


# Just use all trees in the trained model (no iterationrange / ntree_limit needed)
  out[[paste0("xgb_score_", i)]] <- as.numeric(predict(final_model, dtrain_i))



  # store model (optional)
  models[[i]] <- list(
    model = final_model,
    features = sel_cols,
    best_nrounds = best_nrounds_subset,
    params = best_params
  )

  cat(sprintf("Incremental %d/%d | features=%d | best_rounds=%d%s\n",
              i, length(PGS_columns), i, best_nrounds_subset,
              if (isTRUE(attr(cv, "no_es"))) " (no ES fallback)" else ""),
      file = stderr())
}

# -----------------------------
# Write TRAIN scores to stdout
# -----------------------------
write.table(out, "", row.names=FALSE, col.names=TRUE, quote=FALSE)

# -----------------------------
# save all models
# -----------------------------
# ---- SAVE DEBUG / PROBE ----
cat("DEBUG: save_models_path =", save_models_path, "\n", file=stderr())

if (is.null(save_models_path) || !nzchar(save_models_path)) {
  stop("❌ --save-model is required and must be a non-empty path.")
}

save_models_path <- normalizePath(save_models_path, mustWork = FALSE)
out_dirname <- dirname(save_models_path)
dir.create(out_dirname, recursive = TRUE, showWarnings = FALSE)

cat("DEBUG: getwd()       =", getwd(), "\n", file=stderr())
cat("DEBUG: dirname(path) =", out_dirname, "\n", file=stderr())
cat("DEBUG: write access? =", file.access(out_dirname, 2), "(0 means OK)\n", file=stderr())

# Probe write to catch path issues early
probe_path <- paste0(save_models_path, ".probe.rds")
ok <- tryCatch({ saveRDS(list(ok=TRUE), file=probe_path); TRUE },
               error=function(e) { cat("PROBE ERROR:", conditionMessage(e), "\n", file=stderr()); FALSE })
if (!ok) stop("Probe write failed: ", probe_path)


if (!is.null(save_models_path) && nzchar(save_models_path)) {
  save_models_path <- normalizePath(save_models_path, mustWork = FALSE)
  dir.create(dirname(save_models_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(
    models      = models,
    feature_order = PGS_columns,
    outcome_col = outcome_col,
    seed        = SEED,
    tuned_source= opts$tuned,
    nfold       = NFOLD,
    max_rounds  = MAX_ROUNDS,
    early_stop  = EARLY_STOP,
    nthread     = N_THREAD
  ), file = save_models_path)
  cat(sprintf("Saved %d incremental models to: %s\n",
              length(PGS_columns)-1, save_models_path), file=stderr())
}


