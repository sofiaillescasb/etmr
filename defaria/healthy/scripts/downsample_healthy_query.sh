# INPUT="/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev.h5ad"
OUTPUT="/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_bayesian.h5ad"

python - "$INPUT" "$OUTPUT" <<'PY'
import sys
import numpy as np
import pandas as pd
import scanpy as sc

input_file, output_file = sys.argv[1:]

adata = sc.read_h5ad(input_file, backed="r")

# --- BAYESIAN PRIOR PARAMETERS ---
TARGET_TOTAL_CELLS = 50000  # Total desired reference size
MIN_CELLS_FLOOR = 300       # Floor so rare classes remain learnable
# ---------------------------------

rng = np.random.default_rng(42) #random number generator

#this is because the values are encoded as bytes, but i don't know if it's really necessary
def decode_bytes(arr):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr] 

donors = decode_bytes(adata.obs["donor_id"].values)
cellclasses = decode_bytes(adata.obs["CellClass"].values)

df = pd.DataFrame({
    "donor_id": donors,
    "CellClass": cellclasses,
}, index=np.arange(adata.n_obs)) #arange creates a sequence of integers from 0 to the argument

# 1. Compute empirical prior distribution P(CellClass) (the proportion of each class in the data)
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
print(f"Input Dataset: {total_input_cells:,} cells")
print(f"Reference Dataset: {len(keep):,} cells")

adata[keep].write(
    output_file,
    compression="gzip",
    compression_opts=4,
)

print(f"Saved Bayesian Prior reference to: {output_file}")
PY