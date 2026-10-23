#!/bin/bash

input_file="$1"
output_file="${2:-genes.tsv}"

python - "$input_file" "$output_file" <<'PY'
import sys
import scanpy as sc

input_file = sys.argv[1]
output_file = sys.argv[2]

a = sc.read_h5ad(input_file, backed="r")

with open(output_file, "w") as f:
    f.write("Ensembl_ID\tGene\n")
    for ensembl, gene in a.var["Gene"].items():
        f.write(f"{ensembl}\t{gene}\n")

print(f"Extracted {a.n_vars} genes to {output_file}")
a.file.close()
PY