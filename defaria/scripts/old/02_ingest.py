from __future__ import annotations
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
import hdf5plugin
import sys
import os

input_path1 = str(sys.argv[1])
input_path2 = str(sys.argv[2])


print("Reading" + input_path1 + "...")
ref_data = ad.io.read_h5ad(input_path1)
print(ref_data)

print("Reading" + input_path2 + "...")
experiment_data = ad.io.read_h5ad(input_path2)
print(experiment_data)

print("Computing" + input_path1 + "PCA...")
sc.pp.pca(ref_data)

print("Computing" + input_path1 + "neighbors and UMAP...")
sc.pp.neighbors(ref_data, n_neighbors=10, n_pcs=50)
sc.tl.umap(ref_data)

print("Computing" + input_path1 + "Leiden clustering...")
sc.tl.leiden(ref_data, resolution=1)

print("Ingesting" + input_path2 + "...")
sc.tl.ingest(experiment_data, ref_data, obs=["leiden", "Age", "CellClass", "Region", "Subregion", "dissection"], embedding_method = "umap")

print("X_pca" in ref_data.obsm)
print("PCs" in ref_data.varm)
print("X_pca" in ref_data.obsm)

print("Merging datasets...")
all_data = ad.concat([ref_data, experiment_data], label="condition", keys=["ref", "etmr"])
all_data.obs["leiden"] = (
    all_data.obs["leiden"].astype("category").cat.reorder_categories(ref_data.obs["leiden"].cat.categories)
)


#For some reason the Age column is read as string, convert to numeric
experiment_data.obs["Age"] = pd.to_numeric(experiment_data.obs["Age"])


# fix category colors
for key in ref_data.uns:
    if key.endswith("_colors"):
        all_data.uns[key] = ref_data.uns[key]


print("Saving merged dataset...")

ref_data.write_h5ad(os.path.split(str(input_path1))[0] + "/merged_ref_" + os.path.splitext(os.path.basename(str(input_path2)))[0] + ".h5ad")
experiment_data.write_h5ad(str(input_path2)[:-5] + "_int.h5ad")
all_data.write_h5ad(os.path.split(str(input_path1))[0] + "/all_data_int.h5ad")

print(os.path.split(str(input_path1))[0] + "/merged_ref_" + os.path.splitext(os.path.basename(str(input_path2)))[0] + ".h5ad")
print(str(input_path2)[:-5] + "_int.h5ad")
print(os.path.split(str(input_path1))[0] + "/all_data_int.h5ad")



