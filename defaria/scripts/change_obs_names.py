import scanpy as sc
import hdf5plugin

combined_adata_mouse = sc.read_h5ad("/home/sofia/Projects/etmr/defaria/scRNA/data/processed/merged_bothhuman_mouse.h5ad")

combined_adata_mouse.obs_names = combined_adata_mouse.obs["sample_origin"].astype(str) + "_" + combined_adata_mouse.obs_names

combined_adata_mouse.write_h5ad("/home/sofia/Projects/etmr/defaria/scRNA/data/processed/merged_bothhuman_mouse.h5ad")