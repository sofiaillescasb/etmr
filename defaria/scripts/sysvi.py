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
from scvi.external import SysVI

print("Reading" + str(sys.argv[1]) + "...")

input_path = str(sys.argv[1])

base = os.path.splitext(os.path.basename(str(input_path)))[0]
snrna_path = Path(input_path).parents[3]
data_path = Path(input_path).parents[2]

adata = sc.read_h5ad(input_path)

print("Running SysVI...")

sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], n_jobs=1) #regress out total counts and percentage of mitochondrial genes

SysVI.setup_anndata(
    adata=adata,
    batch_key="dataset"
)


model = SysVI(adata)
model.train()

embed = model.get_latent_representation(adata=adata)
embed = sc.AnnData(embed, obs=adata.obs)

sc.pp.neighbors(
                embed,
                n_neighbors=pc_parameters["n_neighbors"],
                n_pcs=pc_parameters["n_pcs"],)
sc.tl.umap(embed)
sc.tl.leiden(
                embed,
                resolution=pc_parameters["leiden_resolution"])

sc.pl.umap(embed, color=["dataset", "leiden"], show=False)
plt.savefig(os.path.join(snrna_path, "plots", "umaps", base + "_sysvi.png"), dpi=300, bbox_inches="tight")
plt.close()


sc.pp.highly_variable_genes(adata, flavor='seurat_v3', layer="counts", n_top_genes=2000) #Calculate highly variable genes, using the unnormalized counts layer


adata.obs["leiden"] = embed.obs["leiden"]

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata)

adata.uns['parameters'] = [pc_parameters, mnn_parameters]

adata.X = csr_matrix(adata.X)


print("Saving data...")

adata.uns["integration"] = ["sysvi"]


adata.write_h5ad(os.path.join(data_path, "processed","sysvi", base + ".h5ad"), compression=hdf5plugin.FILTERS["zstd"], compression_opts=9)