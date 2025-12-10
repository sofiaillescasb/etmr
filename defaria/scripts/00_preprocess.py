import os
import scanpy as sc
import pandas as pd
from scipy import sparse
import anndata as ad
import sys
import hdf5plugin


print("Reading" + str(sys.argv[1]) + "...")

input_path = str(sys.argv[1])

# Load data depending on file type
if input_path.endswith(".h5ad"):
    adata = ad.read_h5ad(input_path)
else:
    # if it's a directory, assume 10x mtx format
    filenames = os.listdir(input_path)
    adatas = [sc.read_10x_mtx(os.path.join(input_path, f)) for f in filenames]
    adata = adatas[0].concatenate(adatas[1:], batch_key='orig.ident', batch_categories=filenames)

adata.write_h5ad(os.path.join(os.path.dirname(os.path.dirname(input_path)),"processed", "snrna_data.h5ad"))

print(adata)
print("adata.X type:", type(adata.X))
print("adata shape:", adata.shape)


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
sc.pp.normalize_total(adata, target_sum=1e4) #library size normalization
sc.pp.log1p(adata) #log-transformation
sc.pp.highly_variable_genes(adata, flavor='seurat') #identify highly variable genes
#adata = adata[:, adata.var['highly_variable']].copy()

adata.layers["scale.data"] = adata.X.copy()
sc.pp.scale(adata, layer="scale.data")


# --- Optional: store sparse counts matrix like in Seurat ---
adata.layers["counts"] = sparse.csr_matrix(adata.raw.X if adata.raw is not None else adata.X)

# --- Save to h5ad (Scanpy format) ---

base = os.path.splitext(os.path.basename(str(input_path)))[0]

adata.write_h5ad(os.path.join(os.path.dirname(os.path.dirname(input_path)), "processed", "snrna_data_filt.h5ad"))

print(os.path.join(os.path.dirname(os.path.dirname(input_path)), "processed","snrna_data_filt.h5ad" ))
