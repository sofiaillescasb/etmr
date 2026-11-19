Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(R.utils)
library(data.table)
library(qs2)
library(Azimuth)
library(reticulate)
library(dplyr)
library(data.table)
library(Matrix)
library(qs2)

#Read Braun
braun_path <- "/media/sofia/TOSHIBA EXT/developing_human_brain/braun/human_dev.h5ad"

braun_mat <- open_matrix_anndata_hdf5(braun_path)

mat_counts <- convert_matrix_type(braun_mat, type = "uint32")

# Write to a compressed binary folder on disk
write_matrix_dir(mat = mat_counts, dir = "/media/sofia/TOSHIBA EXT/developing_human_brain/braun/braun_human_counts")

bp_counts <- open_matrix_dir(dir = "/media/sofia/TOSHIBA EXT/developing_human_brain/braun/braun_human_counts")

#Need to read metadata and gene symbols
# This is the only way I think this can be done on disk, 
# because reticulate converts pandas dataframes into number codes with levels as categories

ad <- import("anndata")

#note: human_dev was filtered to singlets-only
adata <- ad$read_h5ad("/media/sofia/TOSHIBA EXT/developing_human_brain/braun/human_dev.h5ad", backed = "r")

cell_meta <- py_to_r(adata$obs)
gene_meta <- py_to_r(adata$var)


braun <- CreateSeuratObject(
  counts = bp_counts,
  meta.data = cell_meta,
  project = "Braun")

braun$sex <- braun$sex_ontology_term_id
braun$postbirth_age <- NA_real_
braun$dev_stage <- "Fetal"
braun$dataset <- "Braun"

#convert names
og_rownames <- rownames(braun)
new_rownames <- setNames(gene_meta$Gene, gene_meta$Accession)
# unmatched IDs keep their original name
new_names <- new_rownames[og_rownames]
new_names[is.na(new_names)] <- og_rownames[is.na(new_names)] # Keep old ID if no symbol found


rownames(braun) <- new_names

qs_save(braun, "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/braun.qs2", nthreads = 6)