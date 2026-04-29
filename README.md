# PRS Analysis Pipeline

This pipeline implements an incremental multi-PRS integration framework using XGBoost to evaluate how combining multiple polygenic risk scores improves disease risk prediction beyond individual PRSs.

From genotype data to locally optimized PRS and performance evaluation.

## Key Idea: Incremental PRS Integration

Unlike standard stacking approaches, this pipeline adds PRS features sequentially based on their individual predictive performance (C-statistic), allowing direct visualization of the incremental contribution of each PRS.

This design enables:
- transparent interpretation of multi-PRS integration
- identification of diminishing returns
- comparison with best single PRS


![Pipeline Overview](docs/pipeline_overview.png)

---


## Overview

This pipeline takes genotype data and Polygenic Risk Scores (PRS) effect sizes as input, computes per-individual PRS, evaluates their performance, and combines multiple PRS into an locally optimized prediction model using XGBoost.

The pipeline is organized into three steps:

- **Step 1**: Calculate raw PRS from genotype data (BGEN/PLINK) using effect size weight files
- **Step 2**: Merge PRS with phenotype and covariates, residualize by PCs, normalize, and evaluate performance (AUC, OR, R²)
- **Step 3**: Combine multiple PRS sequentially using XGBoost, with hyperparameter tuning and incremental feature addition

**Key design choices:**

- **Imbalance-aware XGBoost**: `scale_pos_weight` auto-set from class counts; PR-AUC used as the primary tuning criterion
- **Incremental feature inclusion**: XGBoost models trained with PRS added sequentially from lowest to highest C-statistic (ranked in Step 2), making the performance gain per feature directly observable
- **Stratified cross-validation**: folds preserve class proportions, safe for highly imbalanced disease datasets
- **stdin/stdout interface**: all R and Python scripts read from stdin and write to stdout, composable with pipes and HPC schedulers

---

## Repository Structure

```
PRS_XGBoost_pipeline/
├── README.md
├── docs/
│   └── pipeline_overview.png                # Pipeline diagram
├── Scripts/
│   ├── prs/
│   │   ├── 01.extract.bgen.sh               # Extract SNPs (supports UKB and MGB)
│   │   ├── 02.setID.sh                      # Set variant ID to chr@:#:ref:alt
│   │   ├── 03.updateID.sh                   # Reformat ID to CHR:POS:A1:A2 (alphabetical)
│   │   ├── 04.calculate.score.sh            # Per-chromosome PRS scoring with plink2
│   │   └── 05.combinechr.R                  # Combine per-chromosome scores
│   ├── Utils/
│   │   ├── KeyMapReplacer.py                # Join two files on a key column
│   │   ├── wcut.py                          # Select columns by name or index
│   │   ├── Residuals.R                      # Regress PRS on PCs, extract residuals
│   │   ├── Scale.R                          # Normalize a column (z-score)
│   │   └── Cstatic_R2_GlmRegression.R       # AUC, OR, incremental R², Nagelkerke R², Liability R²
│   └── xgb/
│       ├── WeightedScore_XGB_train_Tuned_IB_nthread.R   # XGBoost hyperparameter tuning
│       ├── WeightedScore_XGB_train_incremental_test.R   # Incremental XGBoost training
│       └── WeightedScore_XGB_test_score.R               # Apply models to test set
├── workflow/
│   ├── step1.UKB.score.sh                        # Full step1 PRS Calculation (PLINK) example run for UK Biobank
│   ├── step1.MGB.score.sh                        # Full step1 PRS Calculation (PLINK) example run for MGB
│   └── run_xgb.sh                           # XGBoost stage example
├── example/
    └── example.weight.txt                   # Example weight file format

```

---

## Dependencies

### Tools

| Tool | Purpose |
|---|---|
| `plink2` | PRS calculation, variant filtering |
| `Python >= 3.7` | Utility scripts |
| `R >= 4.0` | Scoring, regression, XGBoost |

### R packages

```r
install.packages(c("docopt", "data.table", "xgboost", "glmnet"))
```

## Documentation

Detailed usage and examples for each step:

- 📄 [Step 1: PRS Calculation](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/docs/01_PRS_Calculation.md)
- 📄 [Step 2: PRS Processing and Evaluation](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/docs/02_PRS_Processing_Evaluating.md)
- 📄 [Step 3: XGBoost Integration](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/docs/03_XGBoost.md)

---

## Citation

If you use this pipeline, please cite:

> Yang Sui et al. (2026). *[Paper title]*. *Journal*.
