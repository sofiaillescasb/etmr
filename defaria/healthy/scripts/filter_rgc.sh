#!/usr/bin/env bash
set -euo pipefail

INPUT="/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev.h5ad"
OUTPUT="/home/sofia/Projects/etmr/defaria/healthy/data/raw/radial_glia.h5ad"
CELL_TYPE="Radial glia"

python - "$INPUT" "$OUTPUT" "$CELL_TYPE" <<'PY'
import sys
import scanpy as sc

input_file, output_file, cell_type = sys.argv[1:]

adata = sc.read_h5ad(input_file, backed="r")

mask = adata.obs["CellClass"].astype(str) == cell_type
n = int(mask.sum())

print(f"Input cells: {adata.n_obs:,}")
print(f"Selected '{cell_type}': {n:,}")

if n == 0:
    raise RuntimeError(
        f"No cells found with CellClass == {cell_type!r}"
    )

adata[mask].write(
    output_file,
    compression="gzip",
    compression_opts=4,
)

print(f"Saved: {output_file}")
PY
