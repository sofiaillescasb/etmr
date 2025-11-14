

library(Seurat)
library(patchwork)
library(anndata)
library(reticulate)
library(viridis)
library(tidyverse)


convert_to_seurat <- function(adata) {

  # --- 1. Extract obsm keys ---
  obsm_list <- names(etmr_ad$obsm)

  
  # --- 2. Extract counts matrix ---
  if (!is.null(adata$layers[["counts"]])) {
    counts <- py_to_r(adata$layers[["counts"]])
  } else {
    counts <- py_to_r(adata$X)
  }
  
  counts <- Matrix::Matrix(t(counts), sparse = TRUE)
  
  # --- 3. Metadata ---
  meta <- as.data.frame(py_to_r(adata$obs))
  cell_names <- rownames(meta)
  colnames(counts) <- cell_names
  
  # --- 4. Create Seurat object ---
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta
  )
  
  # --- 5. Add embeddings from obsm ---
  for (key in obsm_list) {
    mat <- py_to_r(adata$obsm[[key]])
    mat <- as.matrix(mat)
    
    # ensure rownames match cells
    rownames(mat) <- cell_names
    
    # Clean key name
    clean_name <- gsub("^X_", "", key)
    
    if (grepl("umap", tolower(clean_name))) {
      seurat_obj[[clean_name]] <- CreateDimReducObject(
        embeddings = mat,
        key = "UMAP_",
        assay = DefaultAssay(seurat_obj)
      )
    } else if (grepl("pca", tolower(clean_name))) {
      seurat_obj[[clean_name]] <- CreateDimReducObject(
        embeddings = mat,
        key = "PC_",
        assay = DefaultAssay(seurat_obj)
      )
    } else {
      seurat_obj@misc[[clean_name]] <- mat
    }
  }
  
  # --- 6. Variable features ---
  var_df <- as.data.frame(py_to_r(adata$var))
  if ("highly_variable" %in% colnames(var_df)) {
    hvf <- rownames(var_df)[which(var_df$highly_variable)]
    VariableFeatures(seurat_obj) <- hvf
  }
  
  return(seurat_obj)
}




plot_one_over_other <- function(umapdf, x, y, colvar, condition, group1, group2) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = umapdf %>% dplyr::filter(get(condition) == group1),
               aes(x = get(x), y = get(y), color = get(colvar)),
               alpha = 1, size = 0.5) +
    ggplot2::geom_point(data = umapdf %>% dplyr::filter(get(condition) == group2),
               aes(x = get(x), y = get(y), color = get(colvar)),
               alpha = 1, size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = x, y = y, color = colvar) 
  return(p)
} 


plot_general <- function(umapdf, x, y, colvar) {
  p <- umapdf %>% 
    ggplot2::ggplot(aes(x = get(x), y = get(y), color = get(colvar))) +
    ggplot2::geom_point(alpha = 0.5, size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = x, y = y, color = colvar)
  
  return(p)
  }

plot_hist <- function(umapdf, x, colvar) {
  h <- umapdf %>%
    ggplot2::ggplot(aes(x = get(x),  fill = get(colvar))) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(colvar, "-integrated clusters")) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::labs(x = x, fill = colvar)  +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(h)
}
