#!/bin/bash

python - <<'PY'
import scanpy as sc

a = sc.read_h5ad("radial_glia_seurat_genes.h5ad", backed="r")

with open("genes.tsv", "w") as f:
    f.write("Ensembl_ID\tGene\n")
    for ensembl, gene in a.var["Gene"].items():
        f.write(f"{ensembl}\t{gene}\n")

print(f"Extracted {a.n_vars} genes to genes.tsv")
PY
