import os
import scanpy as sc
import pandas as pd
from scipy import sparse



adata = sc.read_10x_mtx("/home/sofia/Projects/etmr/de Faria/scRNA/input/paper_data")

# --- Calculate QC metrics ---
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# --- Subset cells ---
adata = adata[
    (adata.obs['n_genes_by_counts'] > 200) &
    (adata.obs['total_counts'] < 40000) &
    (adata.obs['pct_counts_mt'] < 5),
    :
].copy()

# --- Normalize, find variable genes, scale ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat')
adata = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata, max_value=10, zero_center=False)

# --- Optional: store sparse counts matrix like in Seurat ---
adata.layers["counts"] = sparse.csr_matrix(adata.raw.X if adata.raw is not None else adata.X)

# --- Save to h5ad (Scanpy format) ---
adata.write("/home/sofia/Projects/etmr/de Faria/scRNA/input/sctumor.h5ad")
