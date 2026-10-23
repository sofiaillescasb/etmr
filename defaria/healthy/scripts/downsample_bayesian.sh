#!/usr/bin/env bash
set -euo pipefail

# Require at least input and output arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_h5ad> <output_h5ad> [used_cells_file]"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"
USED_CELLS_FILE=USED_CELLS_FILE="${3:-/home/sofia/Projects/etmr/defaria/healthy/data/processed/used_cells_downsampled_bayesian.txt}" # Defaults to used_cells-qs2 if 3rd argument is omitted

python - "$INPUT" "$OUTPUT" "$USED_CELLS_FILE" <<'PY'
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc

input_file, output_file, used_cells_file = sys.argv[1:]

adata = sc.read_h5ad(input_file, backed="r")

# --- BAYESIAN PRIOR PARAMETERS ---
TARGET_TOTAL_CELLS = 100000  # Total desired reference size
MIN_CELLS_FLOOR = 300        # Floor so rare classes remain learnable
# ---------------------------------

rng = np.random.default_rng(42) # random number generator

# Decode bytes helper
def decode_bytes(arr):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr] 

donors = decode_bytes(adata.obs["donor_id"].values)
cellclasses = decode_bytes(adata.obs["CellClass"].values)

# Index corresponds to integer positions (0..n_obs-1) for backed slicing
df = pd.DataFrame({
    "cell_id": decode_bytes(adata.obs_names.values),
    "donor_id": donors,
    "CellClass": cellclasses,
}, index=np.arange(adata.n_obs))

# --- REMOVE USED CELLS ---
if os.path.exists(used_cells_file):
    with open(used_cells_file, "r") as f:
        used_cells = set(line.strip() for line in f if line.strip())
    
    initial_count = len(df)
    df = df[~df["cell_id"].isin(used_cells)]
    print(f"Removed {initial_count - len(df):,} previously used cells listed in '{used_cells_file}'.")
else:
    print(f"Warning: File '{used_cells_file}' not found. No cells removed.")
# -------------------------

# 1. Compute empirical prior distribution P(CellClass)
class_counts = df["CellClass"].value_counts()
total_input_cells = len(df)
priors = class_counts / total_input_cells

print("=== Empirical Priors & Quotas ===")
quotas = {}
for ct, count in class_counts.items():
    prior = priors[ct]
    # Proportional calculation with minimum floor
    prop_quota = int(round(prior * TARGET_TOTAL_CELLS))
    final_quota = max(prop_quota, MIN_CELLS_FLOOR)
    # Cannot exceed available total
    final_quota = min(final_quota, count)
    quotas[ct] = final_quota
    
    print(f"  {ct:<15} | Prior: {prior*100:5.2f}% | Raw: {count:10,} -> Reference Target: {final_quota:8,}")

# 2. Perform stratified sampling across donors for each class target
selected_indices = []

for celltype, ct_df in df.groupby("CellClass"):
    target_k = quotas[celltype]
    
    if len(ct_df) <= target_k:
        selected_indices.extend(ct_df.index.values)
    else:
        # Sample proportionally across donors to retain donor-level variability
        donor_counts = ct_df["donor_id"].value_counts()
        donor_priors = donor_counts / len(ct_df)
        
        ct_selected = []
        for donor_id, d_df in ct_df.groupby("donor_id"):
            d_target = int(round(donor_priors[donor_id] * target_k))
            d_target = min(d_target, len(d_df))
            if d_target > 0:
                sampled = rng.choice(d_df.index.values, size=d_target, replace=False)
                ct_selected.extend(sampled)
                
        # Fix potential rounding discrepancies
        if len(ct_selected) > target_k:
            ct_selected = rng.choice(ct_selected, size=target_k, replace=False)
        elif len(ct_selected) < target_k:
            remaining_idx = np.setdiff1d(ct_df.index.values, ct_selected)
            extra_needed = target_k - len(ct_selected)
            extra = rng.choice(remaining_idx, size=extra_needed, replace=False)
            ct_selected.extend(extra)
            
        selected_indices.extend(ct_selected)

# Sort indices sequentially for fast backed disk read/write
keep = np.sort(np.array(selected_indices, dtype=np.int64))

print("\n--- Summary ---")
print(f"Input Dataset (After Filtering): {total_input_cells:,} cells")
print(f"Reference Dataset: {len(keep):,} cells")

adata[keep].write(
    output_file,
    compression="gzip",
    compression_opts=4,
)

print(f"Saved Bayesian Prior reference to: {output_file}")
PY