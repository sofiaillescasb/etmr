library(BPCells)
library(R.utils)
library(Azimuth)
library(reticulate)
library(dplyr)
library(data.table)
library(Matrix)
library(qs2)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)


herr_path <- "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/"

herr_files <-c("RL1612_filtered_feature_bc_matrix.h5",
               "RL1613_filtered_feature_bc_matrix.h5",
               "RL1777_filtered_feature_bc_matrix.h5",
               "RL2100_filtered_feature_bc_matrix.h5",
               "RL2103_filtered_feature_bc_matrix.h5",
               "RL2104_filtered_feature_bc_matrix.h5",
               "RL2105_filtered_feature_bc_matrix.h5",
               "RL2107_filtered_feature_bc_matrix.h5",
               "RL2108_filtered_feature_bc_matrix.h5",
               "RL2121_filtered_feature_bc_matrix.h5",
               "RL2122_filtered_feature_bc_matrix.h5",
               "RL2125_filtered_feature_bc_matrix.h5",
               "RL2129_filtered_feature_bc_matrix.h5")

for (i in seq_along(herr_files)) {
  path <- file.path(herr_path, herr_files[i])
  data <- open_matrix_10x_hdf5(path)
  matrix.dir <- file.path(herr_path, paste0(tools::file_path_sans_ext(herr_files[i]), "_BP"))
  write_matrix_dir(mat = data, dir = matrix.dir)
}

data.list <- c()

for (i in 1:length(herr_files)) {
  path <- paste0(herr_path, herr_files[i])
  matrix.dir <- file.path(herr_path, paste0(tools::file_path_sans_ext(herr_files[i]), "_BP"))
  # Load in BP matrices
  mat <- open_matrix_dir(dir = matrix.dir)
  data.list[[i]] <- mat
}

# Name layers
names(data.list) <- c("RL1612",
                      "RL1613",
                      "RL1777",
                      "RL2100",
                      "RL2103",
                      "RL2104",
                      "RL2105",
                      "RL2107",
                      "RL2108",
                      "RL2121",
                      "RL2122",
                      "RL2125",
                      "RL2129")

herr_seurat <- lapply(names(data.list), function(sample_name) {
  
  mat <- data.list[[sample_name]]
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )
  
  obj$sample <- sample_name
  
  obj <- RenameCells(
    obj,
    add.cell.id = sample_name
  )
  
  obj
})

names(herr_seurat) <- names(data.list)

herr <- merge(
  x = herr_seurat[[1]],
  y = herr_seurat[-1],
  project = "herr",
  merge.data = FALSE
)

herr <- JoinLayers(herr)

ad <- import("anndata")

#note: human_dev was filtered to singlets-only
adata <- ad$read_h5ad(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad", 
  backed = "r"
  )

cell_meta <- py_to_r(adata$obs)
colnames(cell_meta) <- make.names(colnames(cell_meta))
gene_meta <- py_to_r(adata$var)

sel_cols <- c("batch",                     
              "RL.", "age",                      
              "chem", "concat_id",                
              "numerical_age", "stage_id",                  
              "Sex", "Race",                      
              "PMI", "Brain.Regions.",            
              "Cause.of.Death", "ICD.10.Code",               
              "ICD.10.category", "Oxygen.No.Oxygen",          
              "Date.of.Collection", "Collection_year",           
              "Library.Prep.Date", "Library.Prep.Lot",          
              "n_genes_by_counts", "log1p_n_genes_by_counts",   
              "total_counts", "log1p_total_counts",        
              "doublet_score", "log10_gene_counts",        
              "log10_UMI_counts", "percent_mito",              
              "percent_ribo", "n_counts",                  
              "leiden", "mat.dev",                   
              "cell_type", "major_clust",               
              "sub_clust", "combined.leiden")


cell_meta <- cell_meta |>
  select(all_of(sel_cols)) |>
  filter(RL. %in% unique(herr@meta.data$orig.ident))

rownames(cell_meta) <- sub("^([^-]+)-([^_]+)_.*$", "\\2_\\1-1", rownames(cell_meta))

expr <- GetAssayData(herr, assay = "RNA", layer = "counts")
n_cells <- Matrix::rowSums(expr > 0)

summary(subset(herr, cells =  cells_in_common)$doublet_score)

herr <- AddMetaData(herr, metadata = cell_meta)

cells_in_common <- intersect(colnames(herr), rownames(cell_meta))

herr <- subset(herr, cells = cells_in_common)
herr$dataset <- "Herring"

og_rownames <- rownames(herr)

idx <- match(og_rownames, gene_meta$gene_ids)

new_names <- rownames(gene_meta)[idx]

# Keep original names for unmatched IDs
new_names[is.na(idx)] <- og_rownames[is.na(idx)]


rownames(herr) <- new_names

qs_save(herr, "/media/sofia/TOSHIBA EXT/developing_human_brain/data/herr.qs2", nthreads = 6)

