#!/usr/bin/env bash

set -euo pipefail

CACHE="/media/sofia/TOSHIBA EXT/raw_seqs_tumors/sra_cache"
OUTDIR="/media/sofia/TOSHIBA EXT/raw_seqs_tumors"
TMP="/media/sofia/TOSHIBA EXT/raw_seqs_tumors/tmp"

SRR_LIST=(
  #SRR27831379
  #SRR27831380
  SRR27831381
  #SRR27831382
  #SRR27831383
  #SRR27831384
  #SRR27831385
  #SRR27831386
  #SRR27831387
)

mkdir -p "$CACHE" "$TMP" "$OUTDIR"

for SRR in "${SRR_LIST[@]}"; do
  echo "========================================"
  echo "Processing $SRR"
  echo "========================================"

  prefetch \
    --output-directory "$CACHE" \
    --max-size 100G \
    "$SRR"

  fasterq-dump \
    "$CACHE/$SRR/$SRR.sra" \
    --split-files \
    --threads 8 \
    --temp "$TMP" \
    --outdir "$OUTDIR"

  echo "Finished $SRR"
done

echo "All samples completed successfully."
