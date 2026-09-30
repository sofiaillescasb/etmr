import os
import scanpy as sc
import pandas as pd
from scipy import sparse
import anndata as ad
import hdf5plugin
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import parameters as params



def filter_anndata(adata, input_path):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo"], percent_top=None, log1p=False, inplace=True)

    base = os.path.splitext(os.path.basename(str(input_path)))[0]

    sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
    
    plt.savefig(os.path.join(Path(input_path).parents[3], "plots", "raw", base + "_total_counts_distribution.png"), dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color=["pct_counts_mt", "pct_counts_ribo"], show=False)
              
    plt.savefig(os.path.join(Path(input_path).parents[3], "plots", "raw", base + "_qc_plots.png"), dpi=300, bbox_inches="tight")
    plt.close()



    filters = params.filters

    adata.uns['filters'] = filters

# --- Subset cells ---
    adata = adata[
        (adata.obs['n_genes_by_counts'].between(filters['min_counts'], filters['max_counts'])) &
        (adata.obs['total_counts'] < filters['max_total_counts']) &
        (adata.obs['pct_counts_mt'] < filters['max_pct_counts_mt']), 
        :
    ].copy()
    
    return adata

def preprocess_adata(adata):
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', layer="counts") #identify highly variable genes, IMPORTANT TO DO BEFORE NORMALIZATION AND LOG-TRANSFORMATION
    sc.pp.normalize_total(adata, target_sum=1e4) #library size normalization
    sc.pp.log1p(adata) #log-transformation
    adata.raw = adata.copy() #store the unscaled data in the .raw attribute
    adata.layers["scale.data"] = adata.X.copy() #store the unnormalized data in a new layer called "scale.data"
    sc.pp.scale(adata, layer="scale.data")

    return adata