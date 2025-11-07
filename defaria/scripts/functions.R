

library(Seurat)
library(patchwork)
library(anndata)
library(reticulate)
library(viridis)
library(tidyverse)



convert_to_seurat <- function(adata) {
  # Making data compatible with seurat
  counts <- Matrix::Matrix(t(adata$X), sparse = TRUE)
  meta <- adata$obs
  
  # Make sure cell names match
  cell_names <- rownames(meta)
  colnames(counts) <- cell_names
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = meta,
  )
  
  
  umap <- reticulate::py_to_r(adata$obsm$X_umap)
  
  rownames(umap) <- rownames(adata$X)
  colnames(umap) <- c("UMAP1", "UMAP2")
  
  seurat_obj[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap,
    key = "UMAP_",
    assay = Seurat::DefaultAssay(seurat_obj)
  )
  return(seurat_obj)
}

plot_both <- function(colvar) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = all_umap %>% dplyr::filter(condition == "ref"),
               aes(x = UMAP_1, y = UMAP_2, color = get(colvar)),
               alpha = 0.5, size = 0.5) +
    ggplot2::geom_point(data = all_umap %>% dplyr::filter(condition == "etmr"),
               aes(x = UMAP_1, y = UMAP_2, color = get(colvar)),
               alpha = 1, size = 0.5) +
    ggtheme::theme_classic() +
    ggplot2::guides(colour = guide_legend(override.aes = list(size=3)))
  return(p)
  
} 



