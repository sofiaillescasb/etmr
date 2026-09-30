##CELLCYCLE
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import hdf5plugin #in case needed for h5ad reading (if previously saved with compression)
import genes
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

print("Reading" + str(sys.argv[1]) + "...")

input_path = str(sys.argv[1])

base = os.path.splitext(os.path.basename(str(input_path)))[0]
snrna_path = Path(input_path).parents[3]
data_path = Path(input_path).parents[2]

adata = sc.read_h5ad(input_path)


s_genes = [x for x in genes.s_genes if x in adata.var_names]

g2m_genes = [x for x in genes.g2m_genes if x in adata.var_names]

cell_cycle_genes = [*s_genes, *g2m_genes]

print("Calculating cell cycle score...")

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

datacc_genes = adata[:, cell_cycle_genes].copy()
sc.tl.pca(datacc_genes, layer="scale.data", key_added="pca_cc")
sc.pl.embedding(datacc_genes, color="phase", basis="pca_cc", show=False)
plt.suptitle(base + " cell cycle")
plt.savefig(os.path.join(snrna_path, "plots", "pca", base + "_cell_cycle.png"), dpi=300, bbox_inches="tight")
plt.close()

print("Calculating gene ranking...")

sc.tl.rank_genes_groups(adata, groupby="leiden",  method="wilcoxon", use_raw=False)

print("Plotting gene expression per cluster...")

Path.mkdir(os.path.join(snrna_path, "plots", base), exist_ok=True)

sc.tl.dendrogram(adata, groupby="leiden")

for name, gene in genes.prov_ann.items(): 
    sc.pl.dotplot(adata, var_names=gene, groupby="leiden", standard_scale=None, show=False, dendrogram=True)
    plt.suptitle(base + " " + f"{name} markers")
    plt.savefig(os.path.join(snrna_path, "plots", base, name + "_expression_per_cluster.png"), dpi=300, bbox_inches="tight")
    plt.close()

for name, gene in genes.prov_ann.items(): 
    sc.pl.matrixplot(adata, gene, groupby="leiden", show=False, dendrogram=True)
    plt.suptitle(base + " " + f"{name} markers")
    plt.savefig(os.path.join(snrna_path, "plots", base, name + "_expression_per_cluster.png"), dpi=300, bbox_inches="tight")
    plt.close()


print("Saving data...")

adata.write_h5ad(os.path.join(data_path, "processed","preprocessed", base + ".h5ad"))
