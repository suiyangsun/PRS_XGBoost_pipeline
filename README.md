# PRS Calculation & XGBoost Integration Pipeline

A modular pipeline for computing Polygenic Risk Scores (PRS) from imputed genotype data and combining them into an optimized prediction model using XGBoost. The pipeline covers the complete workflow: raw genotypes → PRS calculation → PC-adjustment → performance evaluation → XGBoost ensemble modeling.

---

## Table of Contents

1. [Overview](#overview)
2. [Repository Structure](#repository-structure)
3. [Dependencies](#dependencies)
4. [Full Pipeline Workflow](#full-pipeline-workflow)
5. [Part 1: PRS Calculation](#part-1-prs-calculation)
   - [Step 1: Extract SNPs from genotype file](#step-1-extract-snps-from-genotype-file)
   - [Step 2: Set variant IDs](#step-2-set-variant-ids)
   - [Step 3: Update variant IDs to match weight file](#step-3-update-variant-ids-to-match-weight-file)
   - [Step 4: Calculate PRS by chromosome](#step-4-calculate-prs-by-chromosome)
   - [Step 5: Combine chromosomes](#step-5-combine-chromosomes)
6. [Part 2: PRS Post-processing and Performance](#part-2-prs-post-processing-and-performance)
7. [Part 3: Locally Optimized PRS (Elastic Net)](#part-3-locally-optimized-prs-elastic-net)
8. [Part 4: XGBoost Ensemble Modeling](#part-4-xgboost-ensemble-modeling)
   - [Step A: Hyperparameter Tuning](#step-a-hyperparameter-tuning)
   - [Step B: Incremental Training](#step-b-incremental-training)
   - [Step C: Test Set Scoring](#step-c-test-set-scoring)
9. [Input / Output Formats](#input--output-formats)
10. [Example End-to-End Runs](#example-end-to-end-runs)
11. [Notes on Reproducibility](#notes-on-reproducibility)

---

## Overview

This pipeline is designed for large-scale GWAS / PRS analyses (e.g., UK Biobank, MGB/GSA-53K cohorts). Multiple PRS features computed from different GWAS summary statistics are first calculated per-individual, then combined using either Elastic Net regression or XGBoost to maximize prediction performance for a binary disease outcome (e.g., CAD).

**Key design choices:**

- **SUM not AVG scoring**: per-chromosome scores are summed (not averaged) to avoid biasing cross-chromosome aggregation
- **Imbalance-aware XGBoost**: `scale_pos_weight` auto-set from class counts; PR-AUC (not ROC AUC) used as the primary tuning criterion
- **Incremental feature inclusion**: XGBoost models are trained with subsets of PRS features (top 2 ... top P), so performance as a function of feature count is directly observable
- **Stratified cross-validation**: folds preserve class proportions throughout, safe for highly imbalanced disease datasets
- **stdin/stdout interface**: all R scripts read from stdin and write to stdout, making them composable with pipes and easy to integrate into HPC schedulers

---

## Repository Structure

```
PRS_cal_pipeline/
├── README.md
├── scripts/
│   ├── prs/
│   │   ├── 01.extract.bgen.sh               # Extract SNPs (supports UKB and MGB)
│   │   ├── 02.setID.sh                      # Set variant ID to chr@:#:ref:alt
│   │   ├── 03.updateID.sh                   # Reformat ID to CHR:POS:A1:A2 (alphabetical)
│   │   ├── 04.calculate.score.sh            # Per-chromosome PRS scoring with plink2
│   │   └── 05.combinechr.R                  # Combine per-chromosome scores
│   ├── utils/
│   │   ├── Residuals_YS.R                   # Regress PRS on PCs, extract residuals
│   │   ├── Scale.R                          # Normalize a column (z-score)
│   │   ├── GlmRegression_YS.R               # AUC, OR, incremental R2 for PRS
│   │   ├── WeightedScore.R                  # Elastic Net / Lasso / Ridge combination
│   │   ├── KeyMapReplacer.py                # Join two files on a key column
│   │   └── wcut.py                          # Select columns by name or index
│   └── xgb/
│       ├── WeightedScore_XGB_train_Tuned_IB_nthread.R   # XGBoost hyperparameter tuning
│       ├── WeightedScore_XGB_train_incremental_test.R   # Incremental XGBoost training
│       └── WeightedScore_XGB_test_score.R               # Apply models to test set
├── workflow/
│   ├── UKB.score.sh                         # Full example run for UK Biobank
│   ├── GSA_53K.score.sh                     # Full example run for MGB GSA-53K
│   └── run_xgb.sh                           # XGBoost stage example
├── example/
│   └── example.weight.txt                   # Example weight file format
└── envs/
    └── xgboost_env.yml                      # Conda environment specification
```

---

## Dependencies

### Tools

| Tool | Purpose |
|---|---|
| `plink2` | PRS calculation, variant filtering |
| `bgenix` / `qctool` | Extracting variants from BGEN files |
| `Python >= 3.7` | Utility scripts |
| `R >= 4.0` | Scoring, regression, XGBoost |

### R packages

```r
install.packages(c("docopt", "data.table", "xgboost", "glmnet"))
```

### Conda environment (for XGBoost steps)

```bash
conda activate /path/to/your/conda/xgboost
# or create from env file:
conda env create -f envs/xgboost_env.yml
```

---

## Full Pipeline Workflow

```
Imputed genotype files (BGEN / PLINK)
             |
             v
[Step 1]  Extract SNPs in weight file  (01.extract.bgen.sh)
             |
             v
[Step 2]  Set variant IDs              (02.setID.sh)
             |
             v
[Step 3]  Reformat IDs to match        (03.updateID.sh)
          weight file format
             |
             v
[Step 4]  Calculate PRS per chr        (04.calculate.score.sh)
             |
             v
[Step 5]  Combine chromosomes          (05.combinechr.R)
             |
             v
[Post-processing]
  Residuals_YS.R     -> regress PRS on PCs
  Scale.R            -> normalize
  GlmRegression_YS.R -> AUC / OR / incremental R2
             |
             |---- Option A: Elastic Net combination
             |     WeightedScore.R
             |
             +---- Option B: XGBoost combination
                        |
                        |-- TRAIN set (e.g. 70%)
                        |       |
                        |       v
                        |   [Step A] Hyperparameter tuning
                        |       WeightedScore_XGB_train_Tuned_IB_nthread.R
                        |       output: tuned_params.rds
                        |       |
                        |       v
                        |   [Step B] Incremental training
                        |       WeightedScore_XGB_train_incremental_test.R
                        |       output: models.rds + train_scores.txt
                        |
                        +-- TEST set (e.g. 30%)
                                |
                                v
                            [Step C] Test scoring
                                WeightedScore_XGB_test_score.R
                                output: test_scores.txt
```

---

## Part 1: PRS Calculation

### Weight file format

All weight files must have the following columns (see `example/example.weight.txt`):

```
SNP             Chr   Pos       effect_allele   other_allele   effect_weight
1:10000:A:C     1     10000     A               C              0.0031
1:20000:G:T     1     20000     G               T             -0.0012
```

> **Important**: the SNP ID format is `CHR:POS:A1:A2` where A1 and A2 are **alphabetically ordered**, regardless of which allele is the effect allele.

---

### Step 1: Extract SNPs from genotype file

A single script handles both UKB and MGB cohorts via the optional `$6` argument:

```bash
bash scripts/prs/01.extract.bgen.sh $plink $bgen $sample $extract $out [$ref]
```

| Argument | Description |
|---|---|
| `$1` | Path to plink2 executable |
| `$2` | Path to BGEN file |
| `$3` | Path to sample file |
| `$4` | Path to extract (bed0 region) file |
| `$5` | Output prefix |
| `$6` | BGEN ref allele convention: `ref-last` (MGB, **default**) or `ref-first` (UKB) |

```bash
# For MGB (GSA-53K): ref-last is the default, $6 can be omitted
bash scripts/prs/01.extract.bgen.sh $plink $bgen $sample $extract $out

# For UKB: pass ref-first explicitly
bash scripts/prs/01.extract.bgen.sh $plink $bgen $sample $extract $out 'ref-first'
```

The script runs:

```bash
REF=${6:-'ref-last'}

$1 \
  --bgen $2 $REF \
  --sample $3 \
  --extract bed0 $4 \
  --make-pgen \
  --out $5
```

Extracting only the SNPs present in the weight file dramatically reduces compute time, especially when calculating multiple scores from the same genotype data.

---

### Step 2: Set variant IDs

```bash
bash scripts/prs/02.setID.sh
```

Sets variant IDs to the intermediate format `chr@:#:ref:alt` required before the final reformat step.

---

### Step 3: Update variant IDs to match weight file

```bash
bash scripts/prs/03.updateID.sh
```

Reformats IDs to `CHR:POS:A1:A2` with A1/A2 alphabetically ordered, matching the weight file SNP column.

---

### Step 4: Calculate PRS by chromosome

```bash
bash scripts/prs/04.calculate.score.sh
```

Uses plink2 `--score` with the following key flags:

```
--score $weight_file 1 4 header-read list-variants ignore-dup-ids \
  cols='sid,nallele,dosagesum,scoresums'
--score-col-nums $col_num
```

> **Note on SUM vs AVG**: this pipeline uses `scoresums` (sum) rather than the plink2 default (average). When combining per-chromosome scores, summing is more appropriate than averaging because the denominator differs across chromosomes. With proper QC, SUM and AVG scores are effectively equivalent, but SUM avoids an arbitrary per-chromosome normalization artifact. See [this discussion](https://groups.google.com/g/prsice/c/hy-C66uo8ok?pli=1) for background.

---

### Step 5: Combine chromosomes

```bash
Rscript scripts/prs/05.combinechr.R
```

Merges per-chromosome `.sscore.gz` files into a single genome-wide PRS file per individual.

---

## Part 2: PRS Post-processing and Performance

After combining chromosomes, the raw PRS is adjusted for population stratification and evaluated. The standard post-processing pipeline:

```bash
name="BMI"
link="gaussian"   # gaussian for continuous phenotype, binomial for binary phenotype

full_model="$name~adjNormPRS+age+inferred_gender+genotyping_array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
null_model="$name~age+inferred_gender+genotyping_array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"

zcat score/sscore.gz | \
  sed 's/#IID/IID/g' | \
  python scripts/utils/KeyMapReplacer.py -k1 -a NA -p<(cat $pheno) -x | \
  sed 's/effect_weight_SUM/PRS/g' | \
  python scripts/utils/wcut.py -t 'IID,PRS,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10' | \
  Rscript scripts/utils/Residuals_YS.R \
    -f 'PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10' -t adjPRS | \
  Rscript scripts/utils/Scale.R -c adjPRS -t adjNormPRS | \
  python scripts/utils/wcut.py -t 'IID,PRS,adjPRS,adjNormPRS' | \
  python scripts/utils/KeyMapReplacer.py -k1 -a NA -p<(cat $pheno) -x | \
  python scripts/utils/wcut.py \
    -t "$name,adjNormPRS,age,inferred_gender,genotyping_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" | \
  Rscript scripts/utils/GlmRegression_YS.R \
    -f $full_model -m $link -n $null_model -r y -p y -a AUC -t $name | \
  tee result/$name.log
```

**Utility script reference:**

| Script | Purpose | Key options |
|---|---|---|
| `Residuals_YS.R` | Regress PRS on covariates (PCs), output residuals | `-f 'PRS~PC1+...'` formula; `-t` output column name |
| `Scale.R` | Z-score normalize a column | `-c` input column; `-t` output column name |
| `GlmRegression_YS.R` | Compute AUC, OR, incremental R², Pearson r | `-m gaussian` or `binomial`; `-f` full model; `-n` null model |
| `KeyMapReplacer.py` | Left-join two files on a key column | `-k1` key col; `-p` second file; `-a NA` fill missing |
| `wcut.py` | Select columns by name or index | `-t 'col1,col2'` by name; `-f4,3` by position |

---

## Part 3: Locally Optimized PRS (Elastic Net)

`WeightedScore.R` combines multiple PRS features into a single optimized score using penalized regression (`cv.glmnet`).

```bash
name="Lasso"
anum=1   # alpha: 1=Lasso, 0=Ridge, 0.5=Elastic Net

zcat $PC4.adjNorm.PRS.join.txt.gz > $unzipped_file

cat $pheno_file | \
  python scripts/utils/wcut.py -t 'IID,Has_cad' | \
  python scripts/utils/KeyMapReplacer.py -k1 -a NA -p $unzipped_file -x | \
  Rscript scripts/utils/WeightedScore.R -m binomial -a $anum \
    2>/dev/stderr > $out_dir/$name.weighted.score.txt
```

| `alpha` value | Model |
|---|---|
| `1` | Lasso (sparse, feature selection) |
| `0` | Ridge (all features retained, shrunk) |
| `0.5` | Elastic Net (compromise) |

The output — a single optimized combined score per individual — can be used directly for evaluation, or included as one of the input features in the XGBoost stage.

---

## Part 4: XGBoost Ensemble Modeling

This stage takes a table of multiple PRS features per individual (e.g., scores from different GWAS, PRS methods, or PC-adjustment strategies) and learns the optimal nonlinear combination.

**Input table format** (whitespace-separated, header required):

```
IID       PRS_A    PRS_B    PRS_C    Has_cad
SAMPLE1   0.412   -0.031    1.203    1
SAMPLE2  -0.887    0.221    0.004    0
```

> The **column order of PRS features** determines the order in which features are added during incremental training. Order your columns from most predictive to least predictive if you want the incremental performance curve to be most informative.

---

### Step A: Hyperparameter Tuning

`WeightedScore_XGB_train_Tuned_IB_nthread.R`

Random search over XGBoost hyperparameters using stratified k-fold CV, optimizing PR-AUC. Saves the best hyperparameter set to an RDS file.

```bash
cat $train | Rscript scripts/xgb/WeightedScore_XGB_train_Tuned_IB_nthread.R \
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

**Outputs:** `tuned_params.rds` (passed to Step B) + `tuned_params.summary.txt` (human-readable)

---

### Step B: Incremental Training

`WeightedScore_XGB_train_incremental_test.R`

Using the tuned parameters, trains P-1 models where model i uses the first i PRS features (i = 2 ... P). Per-subset CV determines the optimal number of boosting rounds for each model. All models are saved to a single RDS.

```bash
cat $train | Rscript scripts/xgb/WeightedScore_XGB_train_incremental_test.R \
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
| `--tuned` | Path to tuned RDS from Step A (required) | — |
| `--save-model` | Output RDS path for all models (required) | — |
| `-o` | Outcome column name (falls back to value saved in tuned RDS) | — |
| `--nrounds` | Max rounds for per-subset CV | 2000 |
| `--early-stop` | Early stopping per subset | 30 |
| `--nfold` | CV folds per subset | 5 |
| `--nthread` | Parallel threads | 4 |

**Stdout:** `IID | outcome | xgb_score_2 | xgb_score_3 | ... | xgb_score_P`

---

### Step C: Test Set Scoring

`WeightedScore_XGB_test_score.R`

Applies each saved incremental model to the held-out test set.

```bash
cat $test | Rscript scripts/xgb/WeightedScore_XGB_test_score.R \
  --models "$out_dir/models.rds" \
  -o Has_cad \
  > "$out_dir/test_scores.txt" \
  2> "$out_dir/test.log"
```

| Option | Description |
|---|---|
| `--models` | Path to incremental models RDS from Step B (required) |
| `-o` | Outcome column name (optional; included in output if found in test data) |

**Stdout:** `IID | outcome | xgb_score_2 | xgb_score_3 | ... | xgb_score_P`

---

## Input / Output Formats

| Stage | Input | Output |
|---|---|---|
| `01.extract.bgen.sh` | BGEN + sample + extract file | PLINK2 pgen/pvar/psam |
| `04.calculate.score.sh` | PLINK2 pgen + weight file | Per-chr `.sscore.gz` |
| `05.combinechr.R` | Per-chr `.sscore.gz` files | Single combined score table |
| `Residuals_YS.R` | IID + PRS + PC columns | Same table + residual column |
| `Scale.R` | Table with target column | Same table + normalized column |
| `GlmRegression_YS.R` | Phenotype + PRS + covariates | AUC / OR / R² summary |
| `WeightedScore.R` | IID + multiple PRS + outcome | Single combined score per IID |
| XGBoost Step A | IID + multiple PRS + outcome | `tuned_params.rds` |
| XGBoost Step B | IID + multiple PRS + outcome | `models.rds` + per-subset train scores |
| XGBoost Step C | IID + multiple PRS (+ outcome) | Per-subset test scores |

---

## Example End-to-End Runs

Full example scripts for UK Biobank and MGB GSA-53K are provided in `workflow/`. A minimal XGBoost-only run:

```bash
#!/usr/bin/env bash
# workflow/run_xgb.sh
set -euo pipefail

script="scripts/xgb"
out_dir="results/CAD_XGB"
train="data/train_set_seed123_0.7.txt"
test="data/test_set_seed123_0.3.txt"

mkdir -p "$out_dir"
source activate /path/to/conda/xgboost

# Step A: tune hyperparameters
cat "$train" | Rscript "$script/WeightedScore_XGB_train_Tuned_IB_nthread.R" \
  --save-model "$out_dir/tuned.rds" \
  -o Has_cad --seed 123 --n-iter 100 --nfold 5 \
  --nrounds 1000 --early-stop 30 --nthread 20

# Step B: incremental training
cat "$train" | Rscript "$script/WeightedScore_XGB_train_incremental_test.R" \
  --tuned "$out_dir/tuned.rds" \
  -o Has_cad --seed 123 --nfold 5 \
  --nrounds 2000 --early-stop 30 --nthread 20 \
  --save-model "$out_dir/models.rds" \
  > "$out_dir/train_scores.txt" 2> "$out_dir/train.log"

# Step C: score test set
cat "$test" | Rscript "$script/WeightedScore_XGB_test_score.R" \
  --models "$out_dir/models.rds" -o Has_cad \
  > "$out_dir/test_scores.txt" 2> "$out_dir/test.log"

conda deactivate
echo "Done. Results in $out_dir/"
```

---

## Notes on Reproducibility

- Use the same `--seed` across all three XGBoost steps for deterministic fold assignments and training.
- `scale_pos_weight` is auto-set from class counts in the training data and stored in the tuned RDS; it is automatically reused in Steps B and C.
- `--nfold` is capped at the minority class size to prevent empty-class folds; a warning is printed when capping occurs.
- Per-chromosome PRS uses SUM (not AVG) to avoid artifactual cross-chromosome normalization differences.
- All scripts log diagnostic messages to stderr; redirect with `2> logfile.log` for HPC job tracking.
- The incremental models RDS stores the complete feature order used during training; the test scoring script validates that all required features are present before scoring.
