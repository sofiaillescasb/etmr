import os
import scanpy as sc
import pandas as pd
import sys
import hdf5plugin
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from parameters import mnn_parameters, pc_parameters
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings("ignore")
import hdf5plugin

print("Reading" + str(sys.argv[1]) + "...")

input_path = str(sys.argv[1])

base = os.path.splitext(os.path.basename(str(input_path)))[0]
snrna_path = Path(input_path).parents[3]
data_path = Path(input_path).parents[2]

adata = sc.read_h5ad(input_path)

#During preprocessing, X was set as normalized, log-transformed counts for both datasets.

print("Reading mouse data...")
humanized_data = sc.read_h5ad("/home/sofia/Projects/etmr/defaria/mouse/so_12_SP_WT_KO_joined_humanized.h5ad")

humanized_data.obs.rename({
    "nFeature_RNA":"n_genes_by_counts",
    "nCount_RNA":"total_counts",
    "sample_origin":"sample"
    }, inplace=True)

adata.obs["dataset"] = "human"
humanized_data.obs["dataset"] = "mouse"

vars_com = adata.var_names.intersection(humanized_data.var_names)

adata = adata[:, vars_com].copy()
humanized_data = humanized_data[:, vars_com].copy()

combined = adata.concatenate(humanized_data, batch_key="dataset")

sc.pp.highly_variable_genes(
    combined,
    flavor="seurat_v3",
    batch_key="dataset",
    n_top_genes=2000
)

hvg = combined.var.highly_variable

adata = adata[:, hvg].copy()
humanized_data = humanized_data[:, hvg].copy()

# 4. Run MNN
adata_mnn = sc.external.pp.mnn_correct(
    adata,
    humanized_data,
    **mnn_parameters
)[0]


sc.pp.scale(adata_mnn)
#The data has been filtered so that only the hvg are left
sc.pp.pca(adata_mnn, use_highly_variable=False, n_comps=50) #compute PCA using the highly variable genes identified earlier
sc.pl.pca_variance_ratio(adata_mnn, n_pcs=50, log=False, show=False) #plot the variance explained by the first 50 PCs

plt.savefig(os.path.join(snrna_path, "plots", "pca", base + "_pca_variance_ratio_mnn.png"), dpi=300, bbox_inches="tight")
plt.close()

print("Calculating neighbors, UMAP, and Leiden clustering...")

sc.pp.neighbors(
                adata_mnn,
                n_neighbors=pc_parameters["n_neighbors"],
                n_pcs=pc_parameters["n_pcs"]
            ) #compute the neighborhood graph of cells using the specified number of neighbors and PCs
            # UMAP
sc.tl.umap(
                adata_mnn
            )

            # Leiden
sc.tl.leiden(
                adata_mnn,
                resolution=pc_parameters["leiden_resolution"]
            )


sc.pl.umap(adata_mnn, color=["n_genes_by_counts", "pct_counts_mt","dataset", "leiden", "phase"], wspace=0.4, legend_loc='on data', show = False)

plt.savefig(os.path.join(snrna_path, "plots", "umaps", base + "_corrected.png"), dpi=300, bbox_inches="tight")
plt.close()

adata_mnn.uns['parameters'] = [pc_parameters, mnn_parameters]

adata_mnn.X = csr_matrix(adata_mnn.X)


print("Saving data...")

adata_mnn.uns["integration"] = ["mnn"]


adata_mnn.write_h5ad(os.path.join(data_path, "processed","corrected", base + ".h5ad"), compression=hdf5plugin.FILTERS["zstd"], compression_opts=9)