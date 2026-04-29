#!/bin/bash

set -euo pipefail

############################################
# Step 4: Calculate PRS using PLINK2 --score
# Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
# Updated: 2026-04-29
############################################

# -----------------------------
# Default parameters
# -----------------------------
SNP_COL=1
WEIGHT_COL=4
SCORE_COLS=6     
EXTRACT=""
REMOVE=""

# -----------------------------
# Usage
# -----------------------------
if [ "$#" -lt 4 ]; then
  echo "Usage:"
  echo "  $0 <plink2_path> <pfile_prefix> <score_file> <output_prefix> [options]"
  echo ""
  echo "Options:"
  echo "  --snp-col INT           (default: 1)"
  echo "  --weight-col INT        (default: 4)"
  echo "  --score-cols STR        (default: 6, e.g. 6 or 6-10)"
  echo "  --extract FILE          (optional)"
  echo "  --remove FILE           (optional)"
  echo ""
  exit 1
fi

# -----------------------------
# Required arguments
# -----------------------------
PLINK=$1
PFILE=$2
SCORE=$3
OUT=$4
shift 4

# -----------------------------
# Parse optional arguments
# -----------------------------
while [[ $# -gt 0 ]]; do
  case $1 in
    --snp-col)
      SNP_COL=$2
      shift 2
      ;;
    --weight-col)
      WEIGHT_COL=$2
      shift 2
      ;;
    --score-cols)
      SCORE_COLS=$2
      shift 2
      ;;
    --extract)
      EXTRACT=$2
      shift 2
      ;;
    --remove)
      REMOVE=$2
      shift 2
      ;;
    *)
      echo "ERROR: Unknown option: $1"
      exit 1
      ;;
  esac
done

# -----------------------------
# Print parameters
# -----------------------------
echo "========== PRS Calculation =========="
echo "PLINK:         $PLINK"
echo "PFILE:         $PFILE"
echo "SCORE file:    $SCORE"
echo "OUTPUT:        $OUT"
echo "SNP column:    $SNP_COL"
echo "Weight column: $WEIGHT_COL"
echo "Score column:  $SCORE_COLS"

if [ -n "$EXTRACT" ]; then
  echo "EXTRACT file:  $EXTRACT"
fi

if [ -n "$REMOVE" ]; then
  echo "REMOVE file:   $REMOVE"
fi

echo "====================================="

# -----------------------------
# Check input files
# -----------------------------
if [ ! -f "${PFILE}.pgen" ]; then
  echo "ERROR: ${PFILE}.pgen not found!"
  exit 1
fi

if [ ! -f "$SCORE" ]; then
  echo "ERROR: score file not found!"
  exit 1
fi

if [ -n "$EXTRACT" ] && [ ! -f "$EXTRACT" ]; then
  echo "ERROR: extract file not found!"
  exit 1
fi

if [ -n "$REMOVE" ] && [ ! -f "$REMOVE" ]; then
  echo "ERROR: remove file not found!"
  exit 1
fi

# -----------------------------
# Sanity check: score file format
# -----------------------------
head -n 5 "$SCORE" | awk '{if(NF<4){exit 1}}' || {
  echo "ERROR: score file must have at least 4 columns"
  exit 1
}

# -----------------------------
# Build PLINK command
# -----------------------------
CMD="$PLINK --pfile $PFILE"

# Optional arguments
if [ -n "$EXTRACT" ]; then
  CMD="$CMD --extract $EXTRACT"
fi

if [ -n "$REMOVE" ]; then
  CMD="$CMD --remove $REMOVE"
fi

# Core scoring
CMD="$CMD \
  --score $SCORE $SNP_COL $WEIGHT_COL header-read list-variants ignore-dup-ids cols='sid,nallele,dosagesum,scoresums' \
  --score-col-nums $SCORE_COLS \
  --out $OUT"

# -----------------------------
# Run
# -----------------------------
echo "Running command:"
echo "$CMD"

eval $CMD

echo "Done!"