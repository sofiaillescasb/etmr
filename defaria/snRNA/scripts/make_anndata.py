from __future__ import annotations
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
import hdf5plugin #in case needed for h5ad reading (if previously saved with compression)
import os



tumor_data = sc.read_h5ad("/home/sofia/Projects/etmr/defaria/snRNA/data/processed/snrna_data.h5ad")

tumor_data.layers["counts"] = tumor_data.X
tumor_data.raw = tumor_data.copy()

tumor_data.obs["mutation"] = tumor_data.obs["orig.ident"].map({
  "GSM8058157_KK22-H-225": "C19MC",
  "GSM8058158_KK22-H-226": "C19MC",
  "GSM8058159_KK22-H-227": "Dicer",
  "GSM8058160_KK23-H-507": "C19MC",
  "GSM8058161_KK23-H-508": "C19MC",
  "GSM8058162_KK23-H-509": "C19MC",
  "GSM8058163_KK23-H-510": "C19MC",
  "GSM8058164_KK23-H-511": "C19MC",
  "GSM8058165_KK23-H-512": "C19MC"
})


tumor_data.obs["id"] = tumor_data.obs["orig.ident"].map({
  "GSM8058157_KK22-H-225": "KK22-H-225",
  "GSM8058158_KK22-H-226": "KK22-H-226",
  "GSM8058159_KK22-H-227": "KK22-H-227",
  "GSM8058160_KK23-H-507": "KK23-H-507",
  "GSM8058161_KK23-H-508": "KK23-H-508",
  "GSM8058162_KK23-H-509": "KK23-H-509",
  "GSM8058163_KK23-H-510": "KK23-H-510",
  "GSM8058164_KK23-H-511": "KK23-H-511",
  "GSM8058165_KK23-H-512": "KK23-H-512"
})

# I'm working with each tumor separately so I can see which ones integrate better with the humanized data

tumor1, tumor2, tumor3, tumor4, tumor5, tumor6, tumor7, tumor8, tumor9 = [x for _, x in tumor_data.obs.groupby(tumor_data.obs['orig.ident'])]

adata = tumor_data[adata.index].copy()