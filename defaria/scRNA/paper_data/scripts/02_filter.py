import pandas as pd
import anndata as ad
import sys
import hdf5plugin #do not remove

common_genes = pd.read_csv("/home/sofia/Projects/etmr/de Faria/snRNA/data/common_genes.txt", header=None)[0].tolist()
print("Reading" + str(sys.argv[1]) + "...")
adata = ad.io.read_h5ad(str(sys.argv[1]), backed="r")

print("Subsetting...")
adata = adata[:, common_genes] #who knows if deepcopy is even possible with this much ram
print("Done subsetting")

#print("Turning into sparse...")
#adata.X = csr_matrix(adata.X)
#print("It's sparse now")

print("Writing h5ad file")
adata.write_h5ad(str(sys.argv[1])[:-5] + "_varfilt.h5ad", compression=hdf5plugin.FILTERS["zstd"])
print("Done")
