library(here)
library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(patchwork)
library(viridis)
library(tidyverse)


# args <- commandArgs(trailingOnly = TRUE)

path_j <- "/home/sofia/Projects/etmr/de Faria/scRNA/data/jessa/so_etmr1_anotado.rds"
etmr_j <- readRDS(path_j)



print("Making UMAP for plotting...")
etmr_j_umap <- cbind(Embeddings(etmr_j, "umap.unintegrated"), etmr_j@meta.data)

print("Plotting integration")

or_clus_plot_j <- ggplot(etmr_j_umap) +
  geom_point(aes(x = umapunintegrated_1, y = umapunintegrated_2, color = seurat_clusters),
             alpha = 1, size = 0.5) +
  theme_classic() 

ggsave(paste0(dirname(dirname(path_j)), "/plots/or_clus_plot_j.jpeg"), or_clus_plot_j, width = 6, height = 3, units = "in", dpi = 600)

cell_type_j <- ggplot(etmr_j_umap) +
  geom_point(aes(x = umapunintegrated_1, y = umapunintegrated_2, color = cell_types_2_etmr1),
             alpha = 1, size = 0.5) +
  theme_classic() 

ggsave(paste0(dirname(dirname(path_j)), "/plots/cell_type_j.jpeg"), cell_type_j, width = 6, height = 3, units = "in", dpi = 600)
