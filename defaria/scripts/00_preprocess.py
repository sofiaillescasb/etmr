import os
import scanpy as sc
import pandas as pd
from scipy import sparse
import anndata as ad
import sys
import hdf5plugin
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import parameters
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings("ignore")

print("Reading" + str(sys.argv[1]) + "...")

input_path = str(sys.argv[1])
base = os.path.splitext(os.path.basename(str(input_path)))[0]
snrna_path = Path(input_path).parents[3]
data_path = Path(input_path).parents[2]


adata = sc.read_10x_mtx(input_path)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)
sns.displot(adata.obs["total_counts"], bins=100, kde=False).savefig(os.path.join(snrna_path, "plots", "raw", base + "_total_counts_distribution.png"))
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color=["pct_counts_mt", "pct_counts_ribo"], show=False)
              
plt.savefig(os.path.join(snrna_path, "plots", "raw", base + "_qc_plots.png"), dpi=300, bbox_inches="tight")
plt.close()


filters = parameters.filters
adata.uns['filters'] = filters

# --- Subset cells ---
adata = adata[
    (adata.obs['n_genes_by_counts'].between(filters['min_counts'], filters['max_counts'])) &
    (adata.obs['total_counts'] < filters['max_total_counts']) &
    (adata.obs['pct_counts_mt'] < filters['max_pct_counts_mt']),
    :
].copy()

adata.layers["counts"] = adata.X.copy() #store the unnormalized data in a new layer called "counts" to use later

print("Transforming data...")

sc.pp.normalize_total(adata, target_sum=1e4) #library size normalization
sc.pp.log1p(adata)

adata.raw = adata.copy() #store the unscaled data in the .raw attribute

sc.pp.highly_variable_genes(adata, flavor='seurat_v3', layer="counts", n_top_genes=2000) #Calculate highly variable genes, using the unnormalized counts layer

adata.layers["scale.data"] = adata.X.copy() #store the unnormalized data in a new layer called "scale.data" to use for scaling
sc.pp.scale(adata, layer="scale.data")

sns.displot(adata.obs["total_counts"], bins=100, kde=False).savefig(os.path.join(snrna_path, "plots", "processed",base + "_total_counts_distribution.png"))
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color=["pct_counts_mt", "pct_counts_ribo"], show=False)
              
plt.savefig(os.path.join(snrna_path, "plots", "processed", base + "_qc_plots.png"), dpi=300, bbox_inches="tight")
plt.close()


adata.X = csr_matrix(adata.X)

print("Saving data...")
adata.write_h5ad(os.path.join(data_path, "processed", "preprocessed",base + ".h5ad"))

print("Done.")