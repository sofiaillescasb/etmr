import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import hdf5plugin
import sys

print("Reading reference data...")
# ref_data = sc.read_h5ad("/home/sofia/Projects/etmr/defaria/healthy/data/processed/human_dev.h5ad")
ref_data = sc.read_h5ad("/home/sofia/Projects/etmr/defaria/healthy/data/processed/so_sub_1_fixed_filt.h5ad")
ref_data.obs["Age"] = pd.to_numeric(ref_data.obs["Age"])  

ref_data.X = csr_matrix(ref_data.X)

print("Reading query data...")
experiment_data = sc.read_h5ad(sys.argv[1])

print("Finding common genes...")
common_genes = ref_data.var_names.intersection(experiment_data.var_names)

common_genes.to_series().to_csv(str(sys.argv[1])[:-5] + "_comvar_" + "common_genes.txt", index=False, header=False)


print("Subsetting...")
ref_data = ref_data[:, common_genes].copy() #scanpy deepcopies when you access the data
experiment_data = experiment_data[:, common_genes].copy()
print("Done subsetting")

print("Writing h5ad files...")
ref_data.write_h5ad(str(sys.argv[1])[:-5] + "_so_sub_1_fixed_filt.h5ad", compression=hdf5plugin.FILTERS["zstd"])
experiment_data.write_h5ad(str(sys.argv[1])[:-5] + "_comvar.h5ad", compression=hdf5plugin.FILTERS["zstd"])
