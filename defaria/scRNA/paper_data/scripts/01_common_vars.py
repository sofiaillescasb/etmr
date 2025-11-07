import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import hdf5plugin


#For the healthy data I'll use the already processed data from the snRNA code project
ref_data = sc.read_h5ad("/home/sofia/Projects/etmr/de Faria/healthy/human_dev_cont_age.h5ad")


experiment_data = sc.read_h5ad("/home/sofia/Projects/etmr/de Faria/scRNA/input/sctumor.h5ad")

common_genes = ref_data.var_names.intersection(experiment_data.var_names)

common_genes.to_series().to_csv("/home/sofia/Projects/etmr/de Faria/scRNA/input/common_genes.txt", index=False, header=False)

