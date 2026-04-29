#! /bin/bash

set -euo pipefail

############################################
# Step 1: Extract the region with SNPs present in the weight file 
# It dramatically reduces compute time, especially when calculating multiple scores from the same genotype data
# Author: Yang Sui (ysui@broadinstitute.org/ysui1@mgh.harvard.edu)
# Last update time: 2026-04-29
############################################
if [ "$#" -lt 5 ]; then
  echo "Usage:"
  echo "  $0 <plink2> <bgen> <sample> <bed_file> <output_prefix> [ref-first|ref-last]"
  echo ""
  echo "IMPORTANT:"
  echo "  bed_file MUST be in BED format (chr start end)"
  exit 1
fi

PLINK=$1
BGEN=$2
SAMPLE=$3
BEDFILE=$4
OUT=$5
REF=${6:-ref-last}

echo "========== PRS SNP Extraction (bed0) =========="
echo "BED file (chr start end): $BEDFILE"
echo "==============================================="

# -----------------------------
# sanity check: BED format
# -----------------------------
head -n 5 "$BEDFILE" | awk '{if(NF<3){exit 1}}' || {
  echo "ERROR: BED file must have at least 3 columns: chr start end"
  exit 1
}

# -----------------------------
# run plink
# -----------------------------
$PLINK \
  --bgen "$BGEN" "$REF" \
  --sample "$SAMPLE" \
  --extract bed0 "$BEDFILE" \
  --make-pgen \
  --out "$OUT"

echo "Done!"
