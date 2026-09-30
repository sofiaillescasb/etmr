import os
import scanpy as sc
import pandas as pd
import sys
import hdf5plugin
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from parameters import pc_parameters
import warnings
warnings.filterwarnings("ignore")

print("Reading" + str(sys.argv[1]) + "...")

input_path = str(sys.argv[1])

base = os.path.splitext(os.path.basename(str(input_path)))[0]
snrna_path = Path(input_path).parents[3]
data_path = Path(input_path).parents[2]

adata = sc.read_h5ad(input_path)

print("Calculating PCA...")

sc.pp.pca(adata, use_highly_variable=True,  n_comps=50, layer="scale.data") #compute PCA using the scaled data in the "scale.data" layer, and only using the highly variable genes identified earlier

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=False, show=False) #plot the variance explained by each of the first 50 PCs, to help determine how many PCs to use for downstream analyses

# Path.mkdir(os.path.join(snrna_path, "plots", "pca"), exist_ok=True)

plt.savefig(os.path.join(snrna_path, "plots", "pca", base + "_pca_variance_ratio.png"), dpi=300, bbox_inches="tight")
plt.close()

print("Calculating neighbors, UMAP, and Leiden clustering...")

sc.pp.neighbors(
                adata,
                n_neighbors=pc_parameters["n_neighbors"],
                n_pcs=pc_parameters["n_pcs"]
            ) #compute the neighborhood graph of cells using the specified number of neighbors and PCs, which will be used for UMAP and Leiden clustering

            # UMAP
sc.tl.umap(
                adata
            )

            # Leiden
sc.tl.leiden(
                adata,
                resolution=pc_parameters["leiden_resolution"]
            )


sc.pl.umap(adata, color=["leiden", "n_genes_by_counts", "pct_counts_mt"], wspace=0.4, legend_loc='on data', show = False)

# Path.mkdir(os.path.join(snrna_path, "plots", "umaps"), exist_ok=True)

plt.savefig(os.path.join(snrna_path, "plots", "umaps", base + "_initial.png"), dpi=300, bbox_inches="tight")
plt.close()

adata.uns['parameters'] = pc_parameters

print("Saving data...")

adata.write_h5ad(os.path.join(data_path, "processed","preprocessed", base + ".h5ad"))
