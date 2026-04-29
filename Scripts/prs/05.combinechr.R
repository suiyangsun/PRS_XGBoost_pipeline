#!/usr/bin/env Rscript

############################################
# Step 5: Sum PRS across chromosomes
# Input: concatenated per-chromosome .sscore files (cat chr1.sscore chr2.sscore ... > all.sscore)
# Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
# Updated: 2026-04-29
############################################

suppressWarnings({
  library(data.table)
  library(optparse)
})

# -----------------------------
# Argument parsing
# -----------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input PRS file (combined chr or single file)"),
  make_option(c("-o", "--output"), type="character", help="Output file"),
  make_option(c("--col"), type="character", default="effect_weight_SUM",
              help="Column name to aggregate (default: effect_weight_SUM)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# -----------------------------
# Check arguments
# -----------------------------
if (is.null(opt$input) || is.null(opt$output)) {
  cat("Usage:\n")
  cat("  Rscript combinechr.R -i input.txt -o output.txt [--col effect_weight_SUM]\n")
  quit(status=1)
}

# -----------------------------
# Load data
# -----------------------------
cat("Loading data...\n")

dt <- fread(opt$input)

# -----------------------------
# Validate columns
# -----------------------------
if (!("IID" %in% colnames(dt))) {
  stop("ERROR: Input must contain column 'IID'")
}

if (!(opt$col %in% colnames(dt))) {
  stop(paste("ERROR: Column", opt$col, "not found in input"))
}

# -----------------------------
# Aggregate PRS
# -----------------------------
cat("Aggregating PRS across chromosomes...\n")

score <- dt[, .(PRS = sum(get(opt$col), na.rm = TRUE)), by = IID]

# -----------------------------
# Output
# -----------------------------
fwrite(score, file=opt$output, sep="\t")

cat("Done! Output written to:", opt$output, "\n")