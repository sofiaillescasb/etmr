import os
import scanpy as sc
import pandas as pd
import sys
import hdf5plugin
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings("ignore")

mouse = sc.read_h5ad("/media/sofia/TOSHIBA EXT/data//mouse/so_12_SP_KO_prefiltered.h5ad")
mouse.obs.drop(columns=['nCount_RNA', 'nFeature_RNA', 'nf'], inplace=True)

# mitochondrial genes
mouse.var["mt"] = mouse.var_names.str.startswith("mt")
# ribosomal genes
mouse.var["ribo"] = mouse.var_names.str.startswith("Rps", "Rpl")



mouse_wt = sc.read_h5ad("/media/sofia/TOSHIBA EXT/data/mouse/so_12_SP_WT_prefiltered.h5ad")
mouse_wt.obs.drop(columns=['nCount_RNA', 'nFeature_RNA', 'nf'], inplace=True)

# wt
mouse_wt.var["mt"] = mouse_wt.var_names.str.startswith("mt")

mouse_wt.var["ribo"] = mouse_wt.var_names.str.startswith("Rps", "Rpl")

mouse.raw = mouse.copy() 
mouse.layers["counts"] = mouse.X.copy()

adata = mouse.concatenate(mouse_wt, batch_key="condition", batch_categories=["KO", "WT"])
adata.layers["counts"] = adata.X.copy()

sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

mask = (
        adata.obs['n_genes_by_counts'].between(3000, 10000) &
        (adata.obs['total_counts'] < 100000) &
        (adata.obs['pct_counts_mt'] < 10) &
        (adata.obs['pct_counts_ribo'] < 30)
    )

adata = adata[mask, :].copy()
adata = adata[:, ~adata.var_names.duplicated()].copy()

#adata = sc.read_h5ad("/media/sofia/TOSHIBA EXT/data/mouse/so_12_SP_KO_WT_just_joined.h5ad")
# adata = adata.raw.to_adata()
adata.obs.rename(columns={"leiden":"leiden_orig", "celltype": "celltype_orig", "phase": "phase_orig", "S_score": "S_score_orig", "G2M_score": "G2M_score_orig"}, inplace=True)

adata.write_h5ad("/media/sofia/TOSHIBA EXT/data/mouse/so_12_SP_KO_WT_just_joined.h5ad")

orthologs_file = pd.read_csv("/home/sofia/Projects/etmr/defaria/mouse/orthologs_final_human_mouse.tsv", sep='\t', header=None, usecols=[0,1])
orthologs_file = pd.read_csv("/home/sofia/Projects/etmr/defaria/mouse/HOM_MouseHumanSequence.rpt", delimiter='\t').loc[:, ["DB Class Key", "NCBI Taxon ID", "Symbol"]]
orthologs_file["Symbol"].unique()

orthologs_file_wide = (
    orthologs_file
    .pivot_table(index="DB Class Key", columns="NCBI Taxon ID", values="Symbol",
        aggfunc="first" )
    .reset_index()
)


orthologs_file_wide.columns = ["class_key", "human", "mouse"]
orthologs_file_wide.drop_duplicates(subset=["human"], inplace=True)

mouse_to_human = orthologs_file_wide.set_index("mouse")["human"]

humanized_data = adata[:, adata.var_names.isin(mouse_to_human.index)].copy()
humanized_data.var_names = mouse_to_human.loc[humanized_data.var_names].values

humanized_data.write_h5ad("/media/sofia/TOSHIBA EXT/data/mouse/so_12_SP_KO_WT_just_joined_humanized.h5ad")