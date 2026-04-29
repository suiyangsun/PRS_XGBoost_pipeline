'
Regression analysis using glm.

Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
Last update time: 2026-04-29

Notes:
    1. Read data from stdin and output results to stdout.
    2. Rows with NA (missing value) will be ignored.
    3. A glm model was fit, supported "family": https://www.statmethods.net/advstats/glm.html

Usage:
    Cstatic_R2_GlmReression.R -f regression_formula -m family_link_function [-n null_model_formula] [-a AUC] [-i beta_95%CI] [-r ModelFitRsquare] [-p PearsonCorrelation] [-k population_prevalence]
 
Options:
    -f string       Regression formular of Full model.eg:CAD~PRS+age+sex+PC.
    -m string       Family regression type that glm supports: binomial,gaussian,poisson ect.
    -n string       Regression formular of Null model.eg:CAD~age+sex+PC
    -a string       Calculate model AUC and incremental AUC (Delta C-statistic) for binary outcomes. Provide anything if "yes".
    -i string       Calculate coefficients 95% CI.
    -r string       Output R2 metrics: variance-based R2 (original), Nagelkerke R2, and logit liability R2 (binomial only). Provide anything if "yes".
    -p string       Output the pearson correlation between raw phenotype and predicted phenotype. Provide anything if "yes".
    -k string       Population prevalence (K_pop) for liability R2 conversion (e.g. 0.03). Optional; if not provided, sample prevalence (P_sample) is used instead.
    -h --help       Show this screen.
    --version       Show version.
' -> doc

## Auto-detect and install needed packages.
options(warn = -1)
suppressMessages(library(pacman))
pacman::p_load(docopt, data.table, pROC, dplyr, rsq)

opts <- docopt(doc, version="V3_2024.revised")

outAUC <- opts$a
outCI  <- opts$i
outR   <- opts$r
outP   <- opts$p
N      <- opts$n
K_pop  <- if (!is.null(opts$k)) as.numeric(opts$k) else NULL


df <- na.omit(read.table(file("stdin"), header=T, check.names=F))

outcome_col <- colnames(df)[1]
reg.formula <- opts$f

print("Regression formula:")
print(reg.formula)

glm.fit      <- glm(as.formula(reg.formula), family=opts$m, data=df)
null.glm.fit <- glm(as.formula(N),           family=opts$m, data=df)

# Output regression results
summary(glm.fit)

# ── 95% CI ─────────────────────────────────────────────────────────────────────
if (!is.null(outCI)) {
    x <- confint(glm.fit)
    y <- as.data.frame(x)
    y$CI <- '95CI'
    print(y)
}

# ── AUC and incremental AUC (Delta C-statistic) ────────────────────────────────
if (!is.null(outAUC)) {
    obs <- df[[outcome_col]]

    # Full model AUC
    pred_full <- predict(glm.fit, type='response')
    roc_full  <- roc(obs, pred_full, quiet=TRUE)
    myauc     <- auc(roc_full)
    ci_full   <- ci.auc(roc_full)
    cat('AUC_Full:\t',        round(myauc, 6),        '\n', sep='')
    cat('AUC_95%_CI_Full:\t', round(ci_full[1], 6), '-', round(ci_full[3], 6), '\n', sep='')

    # Null model AUC + Delta AUC
    if (!is.null(N)) {
        pred_null <- predict(null.glm.fit, type='response')
        roc_null  <- roc(obs, pred_null, quiet=TRUE)
        auc_null  <- auc(roc_null)
        ci_null   <- ci.auc(roc_null)
        cat('AUC_Null:\t',        round(auc_null, 6),       '\n', sep='')
        cat('AUC_95%_CI_Null:\t', round(ci_null[1], 6), '-', round(ci_null[3], 6), '\n', sep='')

        # Delta C-statistic
        delta_auc <- as.numeric(myauc) - as.numeric(auc_null)
        cat('Delta_AUC:\t', round(delta_auc, 6), '\n', sep='')

        # DeLong test
        roc_test <- roc.test(roc_full, roc_null, method='delong')
        cat('Delta_AUC_DeLong_p:\t', format(roc_test$p.value, digits=4), '\n', sep='')
    }
}

# ── R² metrics ─────────────────────────────────────────────────────────────────

# Nagelkerke R² function (built-in, consistent with rcompanion package)
# Reference: Nagelkerke NJD. Biometrika. 1991;78(3):691-692.
nagelkerke_r2 <- function(fit) {
    n           <- fit$df.null + 1
    logLik_full <- as.numeric(logLik(fit))
    fit0        <- update(fit, . ~ 1, family=family(fit))
    logLik_null <- as.numeric(logLik(fit0))
    r2_cs       <- 1 - exp(-2 * (logLik_full - logLik_null) / n)  # Cox-Snell
    r2_max      <- 1 - exp(2 * logLik_null / n)
    r2_cs / r2_max  # Nagelkerke
}

# Logit liability R² - simplified form
# Reference: Lee SH et al. Genet Epidemiol. 2012;36(3):214-224. (Eq.10)
# If K_pop is provided (population prevalence), use Lee et al. full formula with ascertainment correction
# If K_pop is not provided, use P_sample as K_pop (simplified, no ascertainment correction)
liability_r2 <- function(r2_obs, P_sample, K_pop=NULL) {
    if (is.null(K_pop)) {
        # K_pop = P_sample: simplified version (no ascertainment correction)
        t <- qnorm(1 - P_sample)
        z <- dnorm(t)
        return(r2_obs * P_sample * (1 - P_sample) / z^2)
    } else {
        # K_pop provided: use separate K_pop and P_sample (Lee et al. 2012 Eq.10)
        t <- qnorm(1 - K_pop)
        z <- dnorm(t)
        return(r2_obs * (K_pop^2 * (1 - K_pop)^2) / (z^2 * P_sample * (1 - P_sample)))
    }
}


if (!is.null(outR)) {
    is_binomial <- family(glm.fit)$family == "binomial"
 
    # ── 1. Original variance-based R² (Zhang 2017, rsq type='v') ──────────────
    # Kept for backward compatibility with downstream scripts
    model_r2 <- rsq(glm.fit, type='v', adj=F)
    model_r2 <- formatC(model_r2, digits=4, format="g")
    cat('Full_Model_Rsq:\t', model_r2, '\n', sep='')
 
    if (!is.null(N)) {
        model_r2_null <- rsq(null.glm.fit, type='v', adj=F)
        model_r2_null <- formatC(model_r2_null, digits=4, format="g")
        cat('Null_Model_Rsq:\t', model_r2_null, '\n', sep='')
        diff <- as.numeric(model_r2) - as.numeric(model_r2_null)
        diff <- formatC(diff, digits=6, format="g")
        cat('Incremental_Model_Rsq:\t', diff, '\n', sep='')
    }

    # ── 2. Nagelkerke R² ──────────────────────────────────────────────────────
    r2_full <- nagelkerke_r2(glm.fit)
    cat('Nagelkerke_R2_Full:\t', formatC(r2_full, digits=6, format="g"), '\n', sep='')
 
    if (!is.null(N)) {
        r2_null        <- nagelkerke_r2(null.glm.fit)
        r2_incremental <- r2_full - r2_null
        cat('Nagelkerke_R2_Null:\t',        formatC(r2_null,        digits=6, format="g"), '\n', sep='')
        cat('Nagelkerke_R2_Incremental:\t', formatC(r2_incremental, digits=6, format="g"), '\n', sep='')
    }
 
    # ── 3. Logit liability R² (binomial only) ─────────────────────────────────
    if (is_binomial) {
        P_sample <- mean(df[[outcome_col]])  # sample prevalence calculated from data
        cat('P_sample:\t', round(P_sample, 6), '\n', sep='')
 
        # Use K_pop if provided via -k, otherwise fall back to P_sample
        if (!is.null(K_pop)) {
            cat('K_pop:\t', K_pop, '\n', sep='')
        } else {
            cat('K_pop:\t', round(P_sample, 6), ' (using P_sample, no -k provided)\n', sep='')
        }
 
        r2_liab_full <- liability_r2(r2_full, P_sample, K_pop)
        cat('Liability_R2_Full:\t', formatC(r2_liab_full, digits=6, format="g"), '\n', sep='')
 
        if (!is.null(N)) {
            r2_liab_null        <- liability_r2(r2_null,      P_sample, K_pop)
            r2_inc_obs          <- r2_full - r2_null
            r2_liab_incremental <- liability_r2(r2_inc_obs,   P_sample, K_pop)
            cat('Liability_R2_Null:\t',        formatC(r2_liab_null,        digits=6, format="g"), '\n', sep='')
            cat('Liability_R2_Incremental:\t', formatC(r2_liab_incremental, digits=6, format="g"), '\n', sep='')
        }
    }
}

# ── Pearson correlation (predicted vs observed) ─────────────────────────────────
obs_fit  <- glm.fit$y
pred_fit <- glm.fit$fitted.values
 
if (!is.null(outP)) {
    pr <- cor.test(obs_fit, pred_fit, method="pearson")$estimate
    pr <- formatC(pr, digits=4, format="g")
    cat('Pearson_Correlation:\t', pr, '\n', sep='')
    CI <- cor.test(obs_fit, pred_fit, method="pearson")$conf.int[1:2]
    cat('95%_CI_PearsonR:\t', CI[1], '-', CI[2], '\n', sep='')
}
