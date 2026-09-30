#!/usr/bin/env bash
set -euo pipefail

INPUT="/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev.h5ad"
OUTPUT="/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled.h5ad"

python - "$INPUT" "$OUTPUT" <<'PY'
import sys
import numpy as np
import scanpy as sc

input_file, output_file = sys.argv[1:]

adata = sc.read_h5ad(input_file, backed="r")

rng = np.random.default_rng(42)

donor = adata.obs["donor_id"].astype(str).to_numpy()
celltype = adata.obs["CellClass"].astype(str).to_numpy()

n_per_donor = 5000

keep = []

for d in np.unique(donor):
    donor_idx = np.flatnonzero(donor == d)

    # Cell types within this donor
    types = celltype[donor_idx]

    # Allocate approximately proportionally to cell-type abundance
    unique_types, counts = np.unique(types, return_counts=True)

    quotas = np.floor(
        counts / counts.sum() * n_per_donor
    ).astype(int)

    # Make sure rounding gives approximately n_per_donor
    remainder = n_per_donor - quotas.sum()

    if remainder > 0:
        order = np.argsort(-counts)
        quotas[order[:remainder]] += 1

    for ct, quota in zip(unique_types, quotas):
        ct_idx = donor_idx[types == ct]

        quota = min(quota, len(ct_idx))

        if quota > 0:
            selected = rng.choice(
                ct_idx,
                size=quota,
                replace=False,
            )
            keep.extend(selected)

keep = np.asarray(keep, dtype=np.int64)

print(f"Input cells: {adata.n_obs:,}")
print(f"Selected cells: {len(keep):,}")

# This is the important part:
# subset the backed AnnData and write directly.
adata[keep].write(
    output_file,
    compression="gzip",
    compression_opts=4,
)

print(f"Saved: {output_file}")
PY
