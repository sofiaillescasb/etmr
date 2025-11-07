library(scCustomize)
library(here)
library(reticulate)

py_require("anndata")
anndata <- import("anndata")

args <- commandArgs(trailingOnly = TRUE)
print(args)

print("Reading R data...")
rdt <- readRDS(args[1])

print(rdt)

print("Converting to Anndata...")
as.anndata(x = rdt, file_path = gsub("/[^/]*\\.rds$", "", args[1]), file_name = paste0(gsub(".*/(.*)\\.rds$", "\\1", args[1]), ".h5ad"), main_layer = "counts", assay = "RNA", other_layers = NULL)
print("Conversion complete")

