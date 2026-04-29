
# Step 3: XGBoost Integration

PRS features are ranked **using the training set only** based on their individual C-statistic to avoid data leakage. This ranking is then fixed and applied to both training and test sets.
Early stopping is performed using validation folds within cross-validation and does not use the held-out test set.

## Data Leakage Prevention

- PRS ranking is performed using the training set only
- Hyperparameter tuning uses cross-validation within the training set
- The test set is strictly held out and used only for final evaluation


**Input**: Training set and test set with the following format (whitespace-separated, header required):

```
IID    Has_cad    PRS_1    PRS_2    PRS_3   ...
SAMPLE1    0    0.507311823784729    -1.3667209634        ...
SAMPLE2    1    1.515846731698480    0.841451279700       ...
```

> **Note**: the first column must be `IID`, the second column must be the outcome (e.g. `Has_cad`), and PRS features start from the third column onward. The column order of PRS features determines the order in which features are added during incremental training — order your columns from lowest to highest C-statistic (as ranked in Step 2).

**Output**: Final model (`models.rds`) and PRS scores for train and test sets

---

## 3.1. Add PRSs sequentially (low → high)

PRS features are added one at a time in order of increasing C-statistic (as ranked in Step 2). Model i uses the first i PRS features (i = 2 ... P), producing one trained model per feature subset.

---

## 3.2. Hyperparameter Tuning

[`Scripts/xgb/WeightedScore_XGB_train_Tuned_IB_nthread.R`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/xgb/WeightedScore_XGB_train_Tuned_IB_nthread.R)

Random search over XGBoost hyperparameters using stratified k-fold CV, optimizing PR-AUC. Saves the best hyperparameter set to `tuned_params.rds`.


```bash
cat $train | Rscript Scripts/xgb/WeightedScore_XGB_train_Tuned_IB_nthread.R \
  --save-model "$out_dir/tuned_params.rds" \
  -o Has_cad \
  --seed 123 \
  --n-iter 100 \
  --nfold 5 \
  --nrounds 1000 \
  --early-stop 30 \
  --nthread 20
```

| Option | Description | Default |
|---|---|---|
| `--save-model` | Output RDS path (required) | — |
| `-o` | Outcome column name (required) | — |
| `--n-iter` | Number of random hyperparameter sets | 100 |
| `--nfold` | CV folds (auto-capped by minority class size) | 5 |
| `--nrounds` | Max boosting rounds per CV run | 1000 |
| `--early-stop` | Early stopping rounds | 30 |
| `--nthread` | Parallel threads | 4 |
| `--seed` | Random seed | 123 |

**Hyperparameter search space:**

| Parameter | Values searched |
|---|---|
| `eta` (learning rate) | 0.01, 0.03, 0.05, 0.1, 0.2 |
| `max_depth` | 3, 5, 7, 9 |
| `subsample` | 0.6, 0.8, 1.0 |
| `colsample_bytree` | 0.6, 0.8, 1.0 |
| `min_child_weight` | 1, 3, 5 |
| `gamma` | 0, 1, 5 |
| `lambda` (L2 regularization) | 1, 3, 5 |
| `alpha` (L1 regularization) | 0, 1, 5 |

**Output:** `tuned_params.rds` + `tuned_params.summary.txt`

---

## 3.3. Train XGBoost model

[`Scripts/xgb/WeightedScore_XGB_train_incremental_test.R`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/xgb/WeightedScore_XGB_train_incremental_test.R)

Using the tuned parameters, trains P-1 models where model i uses the first i PRS features (i = 2 ... P). Per-subset N-fold cross-validation with early stopping determines the optimal number of boosting rounds for each model. Handles class imbalance automatically via `scale_pos_weight`.


```bash
cat $train | Rscript Scripts/xgb/WeightedScore_XGB_train_incremental_test.R \
  --tuned "$out_dir/tuned_params.rds" \
  -o Has_cad \
  --seed 123 \
  --nfold 5 \
  --nrounds 2000 \
  --early-stop 30 \
  --nthread 20 \
  --save-model "$out_dir/models.rds" \
  > "$out_dir/train_scores.txt" \
  2> "$out_dir/train.log"
```

| Option | Description | Default |
|---|---|---|
| `--tuned` | Path to tuned RDS from step 2 (required) | — |
| `--save-model` | Output RDS path for all models (required) | — |
| `-o` | Outcome column name | — |
| `--nrounds` | Max rounds for per-subset CV | 2000 |
| `--early-stop` | Early stopping per subset | 30 |
| `--nfold` | CV folds per subset | 5 |
| `--nthread` | Parallel threads | 4 |

**Output:** `models.rds` + `train_scores.txt` (`IID | outcome | xgb_score_2 | ... | xgb_score_P`)

---

## 3.4. Test Set Scoring

[`Scripts/xgb/WeightedScore_XGB_test_score.R`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/xgb/WeightedScore_XGB_test_score.R)

Applies each saved incremental model to the held-out test set.

```bash
cat $test | Rscript Scripts/xgb/WeightedScore_XGB_test_score.R \
  --models "$out_dir/models.rds" \
  -o Has_cad \
  > "$out_dir/test_scores.txt" \
  2> "$out_dir/test.log"
```

| Option | Description |
|---|---|
| `--models` | Path to models RDS from step 3 (required) |
| `-o` | Outcome column name (optional) |

**Output:** `test_scores.txt` (`IID | outcome | xgb_score_2 | ... | xgb_score_P`)

---

## Tips

**Running on the full dataset:**

If you want to apply the trained model to the full dataset (not just the test set), pass the full dataset to the incremental training script directly:

```bash
cat $full | Rscript Scripts/xgb/WeightedScore_XGB_train_incremental_test.R \
  --tuned "$out_dir/tuned_params.rds" \
  -o Has_cad \
  --seed 123 \
  --nfold 5 \
  --nrounds 2000 \
  --early-stop 30 \
  --nthread 10 \
  --save-model "$out_dir/full.models.rds" \
  > "$out_dir/xgb_scores.full.txt" \
  2> "$out_dir/xgb_full.log"
```

**Using all PRS features in a single model (no incremental):**

The incremental script always trains models from 2 features up to P features. If you only want the final model using all PRS features, simply extract the last score column from the output:

```bash
cat $train | Rscript Scripts/xgb/WeightedScore_XGB_train_incremental_test.R \
  --tuned "$out_dir/tuned_params.rds" \
  -o Has_cad ... | \
  python Scripts/Utils/wcut.py -t 'IID,outcome,xgb_score_P'
```

Where `xgb_score_P` corresponds to the model trained on all P PRS features. Replace `P` with the actual number of PRS features in your dataset (e.g. `xgb_score_10` if you have 10 PRS features).

---

## Input / Output Formats

| Stage | Input | Output |
|---|---|---|
| XGBoost tuning | IID + ranked PRS + outcome | `tuned_params.rds` |
| XGBoost training | IID + outcome + ranked PRS | `models.rds` + train scores |
| XGBoost scoring | IID + outcome + ranked PRS | Test scores per feature subset |

---

## Full Step 3 pipeline example:


```bash
#!/usr/bin/env bash
# workflow/run_xgb.sh
set -euo pipefail

script="Scripts/xgb"
out_dir="results/CAD_XGB"
train="data/train_set_seed123_0.7.txt"
test="data/test_set_seed123_0.3.txt"

mkdir -p "$out_dir"
source activate /path/to/conda/xgboost

# Step 1: tune hyperparameters
cat "$train" | Rscript "$script/WeightedScore_XGB_train_Tuned_IB_nthread.R" \
  --save-model "$out_dir/tuned.rds" \
  -o Has_cad --seed 123 --n-iter 100 --nfold 5 \
  --nrounds 1000 --early-stop 30 --nthread 20

# Step 2: incremental training
cat "$train" | Rscript "$script/WeightedScore_XGB_train_incremental_test.R" \
  --tuned "$out_dir/tuned.rds" \
  -o Has_cad --seed 123 --nfold 5 \
  --nrounds 2000 --early-stop 30 --nthread 20 \
  --save-model "$out_dir/models.rds" \
  > "$out_dir/train_scores.txt" 2> "$out_dir/train.log"

# Step 3: score test set
cat "$test" | Rscript "$script/WeightedScore_XGB_test_score.R" \
  --models "$out_dir/models.rds" -o Has_cad \
  > "$out_dir/test_scores.txt" 2> "$out_dir/test.log"

conda deactivate
echo "Done. Results in $out_dir/"
```

Full example scripts provided in [`workflow/run_xgb.sh`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/workflow/run_xgb.sh).

---

## Notes on Reproducibility

- Use the same `--seed` across all three XGBoost steps for deterministic fold assignments and training.
- `scale_pos_weight` is auto-set from class counts in the training data and stored in the tuned RDS; it is automatically reused in subsequent steps.
- `--nfold` is capped at the minority class size to prevent empty-class folds; a warning is printed when capping occurs.
- All scripts log diagnostic messages to stderr; redirect with `2> logfile.log` for HPC job tracking.
- The incremental models RDS stores the complete feature order used during training; the test scoring script validates that all required features are present before scoring.
