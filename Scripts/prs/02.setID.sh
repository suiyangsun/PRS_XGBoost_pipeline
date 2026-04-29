#!/bin/bash

set -euo pipefail

############################################
# Step 2: Set variant ID to chr@:#:ref:alt format
# Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
# Last update time: 2026-04-29
############################################

# -----------------------------
# Usage check
# -----------------------------
if [ "$#" -lt 3 ]; then
  echo "Usage:"
  echo "  $0 <plink2_path> <pfile_prefix> <output_prefix>"
  exit 1
fi

# -----------------------------
# Input arguments
# -----------------------------
PLINK=$1
PFILE=$2
OUT=$3

# -----------------------------
# Print parameters
# -----------------------------
echo "========== Set Variant IDs =========="
echo "PLINK:   $PLINK"
echo "PFILE:   $PFILE"
echo "OUTPUT:  $OUT"
echo "====================================="

# -----------------------------
# Check input files
# -----------------------------
if [ ! -f "${PFILE}.pgen" ]; then
  echo "ERROR: ${PFILE}.pgen not found!"
  exit 1
fi

if [ ! -f "${PFILE}.pvar" ]; then
  echo "ERROR: ${PFILE}.pvar not found!"
  exit 1
fi

if [ ! -f "${PFILE}.psam" ]; then
  echo "ERROR: ${PFILE}.psam not found!"
  exit 1
fi

# -----------------------------
# Run PLINK2
# -----------------------------
echo "Setting variant IDs..."

$PLINK \
  --pfile "$PFILE" \
  --set-all-var-ids 'chr@:#:$r:$a' \
  --new-id-max-allele-len 10000 truncate \
  --make-pgen \
  --out "$OUT"

echo "Done!"