from __future__ import annotations
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
import hdf5plugin #in case needed for h5ad reading (if previously saved with compression)


etmr1 = ad.io.read_h5ad("/home/sofia/Projects/etmr/defaria/snRNA/data/processed/snrna_data_umap_no_integration.h5ad")


sc.tl.rank_genes_groups(etmr1, groupby="leiden_unintegrated", method="wilcoxon")
