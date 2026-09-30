import os
import scanpy as sc
import pandas as pd
import sys
import hdf5plugin
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from parameters import pc_parameters
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings("ignore")

print("Reading" + str(sys.argv[1]) + "...")

input_path = str(sys.argv[1])
base = os.path.splitext(os.path.basename(str(input_path)))[0]
# snrna_path = Path(input_path).parents[3]
# data_path = Path(input_path).parents[2]

merged = sc.read_h5ad(input_path)
# merged.X = merged.layers["counts"].copy()

print("Calculating highly variable genes, normalizing, and log-transforming data...")

sc.pp.highly_variable_genes(merged, flavor='seurat_v3', layer="counts", batch_key="species", n_top_genes=2000) #Calculate highly variable genes, using the unnormalized counts layer
sc.pp.normalize_total(merged, target_sum=1e4) #library size normalization
sc.pp.log1p(merged)
merged.raw = merged.copy() 



# sc.pp.regress_out(merged, ["total_counts", "pct_counts_mt"], n_jobs=1)





print("Scaling data...")

merged.layers["scale.data"] = merged.X.copy() #store the unnormalized data in a new layer called "scale.data" to use for scaling
sc.pp.scale(merged, layer="scale.data")

print("Integrating data with Harmony...")

sc.pp.pca(merged, use_highly_variable=True,  n_comps=50, layer="scale.data")
sc.pl.pca_variance_ratio(merged, n_pcs=50, log=False, show=False) #plot the variance explained by each of the first 50 PCs, to help determine how many PCs to use for downstream analyses

plt.savefig(os.path.join("/home/sofia/Projects/etmr/defaria/snRNA/plots/umaps/harmony/integrated/regressed", base + "_pca_variance_ratio_regressed.png"), dpi=300, bbox_inches="tight")
plt.close()

sc.external.pp.harmony_integrate(merged, key="sample_origin") #integrate the data using Harmony, specifying the batch key, the dimensionality reduction to use as input, and the number of PCs to use

# sns.displot(merged.obs["total_counts"], bins=100, kde=False).savefig(os.path.join("/home/sofia/Projects/etmr/defaria/snRNA/plots/processed/integrated",base + "_total_counts_distribution.png"))
# sc.pl.scatter(merged, "total_counts", "n_genes_by_counts", color=["pct_counts_mt", "pct_counts_ribo"], show=False)
              
# plt.savefig(os.path.join("/home/sofia/Projects/etmr/defaria/snRNA/plots/processed/integrated", base + "_qc_plots_regressed.png"), dpi=300, bbox_inches="tight")
# plt.close()



print("Calculating neighbors, UMAP, and Leiden clustering...")

sc.pp.neighbors(
                merged,
                n_neighbors=pc_parameters["n_neighbors"],
                n_pcs=pc_parameters["n_pcs"],
                use_rep="X_pca_harmony"
            ) #compute the neighborhood graph of cells using the specified number of neighbors and PCs, which will be used for UMAP and Leiden clustering

            # UMAP
sc.tl.umap(
                merged
            )

            # Leiden
sc.tl.leiden(
                merged,
                resolution=pc_parameters["leiden_resolution"]
            )


sc.pl.umap(merged, color=["n_genes_by_counts", "pct_counts_mt", "sample_origin", "species","leiden"], wspace=0.4, show = False)

plt.savefig(os.path.join("/home/sofia/Projects/etmr/defaria/snRNA/plots/umaps/harmony/integrated/regressed", base + "_umap.png"), dpi=300, bbox_inches="tight")
plt.close()

merged.uns['parameters'] = pc_parameters

merged.X = csr_matrix(merged.X)


print("Saving data...")

merged.uns["integration"] = ["harmony"]

merged.write_h5ad(os.path.join("/home/sofia/Projects/etmr/defaria/snRNA/data/processed/harmony", base + "_integrated.h5ad"), compression=hdf5plugin.FILTERS["zstd"])
#merged.write_h5ad(os.path.join("/media/sofia/TOSHIBA EXT/data/harmony", "regressed", base + "_integrated.h5ad"), compression=hdf5plugin.FILTERS["zstd"])