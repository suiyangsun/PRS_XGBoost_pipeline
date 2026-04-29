
# Step 2: PRS Processing and Evaluation

**Input:**
- `sscore.txt`: combined PRS file from Step 1.5
- `$pheno`: phenotype and covariate file (tab-separated, first column must be named `IID`)

Example `$pheno` format:
```
IID        CAD   age   inferred_gender   genotyping_array   PC1      PC2    ...
SAMPLE1    1         55    M                 GSA                0.012   -0.003  ...
SAMPLE2    0         62    F                 GSA               -0.008    0.011  ...
```

**Output:** Performance metrics for each PRS (C-statistic, OR, R², etc.) + PRS ranking (low → high) saved as `prs_ranking.txt` for use in Step 3

---

## 2.1. Merge PRS with phenotype

[`Scripts/Utils/KeyMapReplacer.py`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/Utils/KeyMapReplacer.py)
[`Scripts/Utils/wcut.py`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/Utils/wcut.py)

```bash
cat score/sscore.txt | \
  python Scripts/Utils/KeyMapReplacer.py -k1 -a NA -p<(cat $pheno) -x | \
  python Scripts/Utils/wcut.py -t 'IID,PRS,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10'
```

---

## 2.2. Residualize PRS by PCs

Regresses PRS on top PCs and extracts residuals to remove population stratification:

[`Scripts/Utils/Residuals.R`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/Utils/Residuals.R)

```bash
Rscript Scripts/Utils/Residuals.R \
  -f 'PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10' \
  -t adjPRS
```

| Option | Description |
|---|---|
| `-f` | Regression formula (required) |
| `-t` | Output column name for residuals (required) |

Residualization is performed using linear regression with an intercept term.
---

## 2.3. Normalize PRS

[`Scripts/Utils/Scale.R`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/Utils/Scale.R)

Z-score normalizes the residualized PRS:

```bash
Rscript Scripts/Utils/Scale.R -c adjPRS -t adjNormPRS
```

| Option | Description |
|---|---|
| `-c` | Input column to normalize |
| `-t` | Output column name |

---

## 2.4. Evaluate PRS performance

[`Scripts/Utils/Cstatic_R2_GlmRegression.R`](https://github.com/suiyangsun/PRS_XGBoost_pipeline/blob/main/Scripts/Utils/Cstatic_R2_GlmRegression.R)

Fits full and null models and computes performance metrics including AUC, Delta C-statistic (with DeLong test), Nagelkerke R², and Liability R²:

```bash
name="CAD"
link="binomial"   # gaussian for continuous phenotype, binomial for binary phenotype

full_model="$name~adjNormPRS+age+inferred_gender+genotyping_array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
null_model="$name~age+inferred_gender+genotyping_array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"

Rscript Scripts/Utils/Cstatic_R2_GlmRegression.R \
  -f $full_model \
  -n $null_model \
  -m $link \
  -a y \
  -i y \
  -r y \
  -p y \
  -k 0.03 \
  -t $name
```

| Option | Description |
|---|---|
| `-f` | Full model formula (required) |
| `-n` | Null model formula |
| `-m` | Regression link: `gaussian` (continuous) or `binomial` (binary) |
| `-a` | Calculate AUC and Delta C-statistic with DeLong test (provide any value for yes) |
| `-i` | Calculate 95% CI for coefficients (provide any value for yes) |
| `-r` | Output R² metrics: variance R², Nagelkerke R², Liability R² (provide any value for yes) |
| `-p` | Output Pearson correlation between observed and predicted (provide any value for yes) |
| `-k` | Population prevalence for liability R² conversion (e.g. `0.03`); if omitted, sample prevalence is used |

**Output metrics:**

| Metric | Description |
|---|---|
| `AUC_Full` | AUC of full model with 95% CI |
| `AUC_Null` | AUC of null model with 95% CI |
| `Delta_AUC` | AUC gain over null model |
| `Delta_AUC_DeLong_p` | DeLong test p-value for Delta AUC |
| `Full_Model_Rsq` | Variance-based R² of full model |
| `Incremental_Model_Rsq` | R² gain over null model |
| `Nagelkerke_R2_Full` | Nagelkerke R² of full model |
| `Nagelkerke_R2_Incremental` | Nagelkerke R² gain over null model |
| `Liability_R2_Full` | Liability R² of full model (binomial only) |
| `Liability_R2_Incremental` | Liability R² gain over null model (binomial only) |
| `Pearson_Correlation` | Pearson r between observed and predicted with 95% CI |

The C-statistic from this step is used to **rank PRS from low to high** (`prs_ranking.txt`), which determines the order in which PRS features are added in Step 3.

---

## Full Step 2 pipeline example

```bash
name="CAD"
link="binomial"
full_model="$name~adjNormPRS+age+inferred_gender+genotyping_array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
null_model="$name~age+inferred_gender+genotyping_array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"

# Step 1: Merge PRS with phenotype and covariates
cat score/sscore.txt | \
  python Scripts/Utils/KeyMapReplacer.py -k1 -a NA -p<(cat $pheno) -x | \

# Step 2: Select IID, PRS, and top 10 PCs
  python Scripts/Utils/wcut.py -t 'IID,PRS,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10' | \

# Step 3: Residualize PRS by top 10 PCs to remove population stratification
  Rscript Scripts/Utils/Residuals.R \
    -f 'PRS~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10' -t adjPRS | \

# Step 4: Z-score normalize the residualized PRS
  Rscript Scripts/Utils/Scale.R -c adjPRS -t adjNormPRS | \

# Step 5: Keep IID, raw PRS, residualized PRS, and normalized PRS
  python Scripts/Utils/wcut.py -t 'IID,PRS,adjPRS,adjNormPRS' | \

# Step 6: Merge back with full phenotype and covariate file
  python Scripts/Utils/KeyMapReplacer.py -k1 -a NA -p<(cat $pheno) -x | \

# Step 7: Select columns needed for regression
  python Scripts/Utils/wcut.py \
    -t "$name,adjNormPRS,age,inferred_gender,genotyping_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" | \

# Step 8: Fit full and null models, compute AUC, OR, R²
  Rscript Scripts/Utils/Cstatic_R2_GlmRegression.R \
    -f $full_model -n $null_model -m $link -a y -i y -r y -p y -k 0.03 -t $name | \

# Step 9: Save output to log file
  tee result/$name.log
```

---

### Utility script reference

| Script | Purpose | Key options |
|---|---|---|
| `KeyMapReplacer.py` | Left-join two files on a key column | `-k1` key col; `-p` second file; `-a NA` fill missing |
| `wcut.py` | Select columns by name or index | `-t 'col1,col2'` by name; `-f4,3` by position |
| `Residuals.R` | Regress PRS on covariates, output residuals | `-f` formula; `-t` output column name |
| `Scale.R` | Z-score normalize a column | `-c` input column; `-t` output column name |
| `Cstatic_R2_GlmRegression.R` | Compute AUC, Delta C-statistic, R², Nagelkerke R², Liability R² | `-f` full model; `-n` null model; `-m` link function; `-k` population prevalence |
---
