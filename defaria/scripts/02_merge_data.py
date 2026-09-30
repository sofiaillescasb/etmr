import os
import scanpy as sc
import pandas as pd
import anndata as ad
import sys
import hdf5plugin
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from parameters import pc_parameters
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings("ignore")

print("Reading mouse data...")
humanized_data = sc.read_h5ad("/media/sofia/TOSHIBA EXT/data/mouse/so_12_SP_KO_WT_just_joined_humanized.h5ad")
humanized_data.X = humanized_data.layers["counts"].copy()



input_path = str(sys.argv[1])

base = os.path.splitext(os.path.basename(str(input_path)))[0]
snrna_path = Path(input_path).parents[3]
data_path = Path(input_path).parents[2]

adata = sc.read_h5ad(input_path)
adata.obs["sample_origin"] = base
adata.X = adata.layers["counts"].copy()


mask = humanized_data.var_names.intersection(adata.var_names)
mask
humanized_data_filtered = humanized_data[:, mask].copy()
tumor_filtered = adata[:, mask].copy()

print("Merging data...")
combined_adata = ad.concat(
    [humanized_data_filtered, tumor_filtered],
    label="species",
    keys=["humanized", "tumor"]
)
    
combined_adata.X = csr_matrix(combined_adata.X)
print("Saving data...")

combined_adata.write_h5ad(os.path.join("/home/sofia/Projects/etmr/defaria/snRNA/data/processed/merged", "combined_" + base + ".h5ad"), compression=hdf5plugin.FILTERS["zstd"])
# combined_adata.write_h5ad(os.path.join("/media/sofia/TOSHIBA EXT/data/merged/new", "combined_" + base + ".h5ad"))
