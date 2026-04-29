#!/bin/bash

set -euo pipefail

############################################
# Step 3: Reformat ID to CHR:POS:A1:A2 (alphabetical) format to match PRS weight file
# Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
# Last update time: 2026-04-29
############################################

# -----------------------------
# Usage check
# -----------------------------
if [ "$#" -lt 4 ]; then
  echo "Usage:"
  echo "  $0 <plink2_path> <pfile_prefix> <update_file> <output_prefix>"
  echo ""
  echo "IMPORTANT:"
  echo "  update_file must contain mapping to PRS weight ID format CHR:POS:A1:A2 (alphabetical)"
  exit 1
fi

# -----------------------------
# Input arguments
# -----------------------------
PLINK=$1
PFILE=$2
UPDATE_FILE=$3
OUT=$4

# -----------------------------
# Print parameters
# -----------------------------
echo "========== Update Variant IDs =========="
echo "PLINK:        $PLINK"
echo "PFILE:        $PFILE"
echo "UPDATE FILE:  $UPDATE_FILE"
echo "OUTPUT:       $OUT"
echo "========================================"

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

if [ ! -f "$UPDATE_FILE" ]; then
  echo "ERROR: update file not found!"
  exit 1
fi

# -----------------------------
# Sanity check: update file format
# -----------------------------
head -n 5 "$UPDATE_FILE" | awk '{if(NF<2){exit 1}}' || {
  echo "ERROR: update file must have at least 2 columns: oldID newID"
  exit 1
}

# -----------------------------
# Run PLINK2
# -----------------------------
echo "Updating variant IDs..."

$PLINK \
  --pfile "$PFILE" \
  --update-name "$UPDATE_FILE" \
  --make-pgen \
  --out "$OUT"

echo "Done!"