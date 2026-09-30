python - <<'PY'
import scanpy as sc

src = "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled.h5ad"
dst = "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled_seurat_genes.h5ad"

# Read only metadata initially
adata = sc.read_h5ad(src, backed="r")

# Genes from the Seurat object
with open("/home/sofia/Projects/etmr/defaria/snRNA/data/seurat_genes.txt") as f:
    genes = [x.strip() for x in f if x.strip()]

# Match against the H5AD's Gene column
gene_set = set(genes)
keep = adata.var["Gene"].astype(str).isin(gene_set).to_numpy()

print(f"Cells: {adata.n_obs:,}")
print(f"Original genes: {adata.n_vars:,}")
print(f"Genes retained: {keep.sum():,}")

if keep.sum() == 0:
    raise RuntimeError("No genes matched.")

# Preserve Seurat gene order
gene_order = [g for g in genes if g in set(adata.var.loc[keep, "Gene"].astype(str))]

# Map Gene -> H5AD column index
gene_to_idx = {
    str(g): i
    for i, g in enumerate(adata.var["Gene"])
}

indices = [gene_to_idx[g] for g in gene_order]

# Filter and write compressed
subset = adata[:, indices]

subset.write(
    dst,
    compression="gzip",
    compression_opts=4,
)

print(f"Saved: {dst}")
PY
