from __future__ import annotations
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
import hdf5plugin #in case needed for h5ad reading (if previously saved with compression)
import sys
import os

# input_path1 = str(sys.argv[1])
# input_path2 = str(sys.argv[2])

input_path1 = "/home/sofia/Projects/etmr/defaria/scRNA/jessa/data/processed/so_etmr1_anotado_filt.h5ad"
input_path2 = "/home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/paper_data_filt.h5ad"


print("Reading" + input_path1 + "...")
etmr1 = ad.io.read_h5ad(input_path1)
print(etmr1)

print("Reading" + input_path2 + "...")
etmr2 = ad.io.read_h5ad(input_path2)
print(etmr2)

print("Finding common genes...")
common_genes = etmr1.var_names.intersection(etmr2.var_names)
etmr1 = etmr1[:, common_genes].copy()
etmr2 = etmr2[:, common_genes].copy()

print("Computing" + input_path1 + "PCA...")
sc.pp.pca(etmr1)
print("Computing" + input_path1 + "neighbors and UMAP...")
sc.pp.neighbors(etmr1, n_neighbors=10, n_pcs=50)
sc.tl.umap(etmr1)

print("Computing" + input_path1 + "Leiden clustering...")
sc.tl.leiden(etmr1, resolution=1)


print("Ingesting" + input_path2 + "...")
etmr2_int = etmr2.copy()

sc.tl.ingest(etmr2_int, etmr1, obs=['leiden', 'cell_types_etmr1', 'cell_types_2_etmr1'], embedding_method = "umap")

print("Merging datasets...")
all_data = ad.concat([etmr1, etmr2_int], label="condition", keys=["jessa", "defaria"])
all_data.obs["leiden"] = (
    all_data.obs["leiden"].astype("category").cat.reorder_categories(etmr1.obs["leiden"].cat.categories)
)



# fix category colors
for key in etmr1.uns:
    if key.endswith("_colors"):
        all_data.uns[key] = etmr1.uns[key]


print("Saving merged dataset...")

etmr1.write_h5ad(os.path.split(str(input_path1))[0] + "/merged_etmr1" + ".h5ad")
etmr2_int.write_h5ad(str(input_path2)[:-5] + "_etmrint.h5ad")
all_data.write_h5ad(os.path.split(str(input_path1))[0] + "/all_etmr_int.h5ad")



print("Computing" + input_path2 + " individual PCA...")
sc.pp.pca(etmr2)
print("Computing" + input_path2 + " individual neighbors and UMAP...")
sc.pp.neighbors(etmr2, n_neighbors=10, n_pcs=50)
sc.tl.umap(etmr2)

print("Computing" + input_path2 + " individual Leiden clustering...")
sc.tl.leiden(etmr2, resolution=1)

etmr2.write_h5ad(str(input_path2)[:-5] + "_alone.h5ad")

