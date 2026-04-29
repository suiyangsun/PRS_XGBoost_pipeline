#!/usr/bin/env bash
set -euo pipefail
############################################
# XGBoost Integration Pipeline
# Step 1: Hyperparameter tuning
# Step 2: Incremental training on train set
# Step 3: Score test set
# Step 4: (Optional) Score full dataset
# Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
# Last update time: 2026-04-29
############################################

# -----------------------------
# Paths
# -----------------------------
script="Scripts/xgb"
out_dir="results/CAD_XGB"
train="data/train_set_seed123_0.7.txt"
test="data/test_set_seed123_0.3.txt"
full="data/full_set.txt"   # optional: for scoring full dataset

mkdir -p "$out_dir"

# -----------------------------
# Activate conda environment
# -----------------------------
source activate /path/to/conda/xgboost

# -----------------------------
# Step 1: Hyperparameter tuning
# -----------------------------
echo "========== Step 1: Hyperparameter Tuning =========="
cat "$train" | Rscript "$script/WeightedScore_XGB_train_Tuned_IB_nthread.R" \
  --save-model "$out_dir/tuned.rds" \
  -o Has_cad \
  --seed 123 \
  --n-iter 100 \
  --nfold 5 \
  --nrounds 1000 \
  --early-stop 30 \
  --nthread 20
echo "Tuning done. Saved to $out_dir/tuned.rds"

# -----------------------------
# Step 2: Incremental training on train set
# -----------------------------
echo "========== Step 2: Incremental Training =========="
cat "$train" | Rscript "$script/WeightedScore_XGB_train_incremental_test.R" \
  --tuned "$out_dir/tuned.rds" \
  -o Has_cad \
  --seed 123 \
  --nfold 5 \
  --nrounds 2000 \
  --early-stop 30 \
  --nthread 20 \
  --save-model "$out_dir/models.rds" \
  > "$out_dir/train_scores.txt" \
  2> "$out_dir/train.log"
echo "Training done. Scores saved to $out_dir/train_scores.txt"

# -----------------------------
# Step 3: Score test set
# -----------------------------
echo "========== Step 3: Test Set Scoring =========="
cat "$test" | Rscript "$script/WeightedScore_XGB_test_score.R" \
  --models "$out_dir/models.rds" \
  -o Has_cad \
  > "$out_dir/test_scores.txt" \
  2> "$out_dir/test.log"
echo "Test scoring done. Scores saved to $out_dir/test_scores.txt"

# -----------------------------
# Step 4: (Optional) Score full dataset
# -----------------------------
echo "========== Step 4: Full Dataset Scoring (Optional) =========="
if [ -f "$full" ]; then
  cat "$full" | Rscript "$script/WeightedScore_XGB_train_incremental_test.R" \
    --tuned "$out_dir/tuned.rds" \
    -o Has_cad \
    --seed 123 \
    --nfold 5 \
    --nrounds 2000 \
    --early-stop 30 \
    --nthread 20 \
    --save-model "$out_dir/full.models.rds" \
    > "$out_dir/full_scores.txt" \
    2> "$out_dir/full.log"
  echo "Full dataset scoring done. Scores saved to $out_dir/full_scores.txt"
else
  echo "Full dataset file not found ($full), skipping Step 4."
fi

# -----------------------------
# Done
# -----------------------------
conda deactivate
echo "========== All steps completed. Results in $out_dir/ =========="