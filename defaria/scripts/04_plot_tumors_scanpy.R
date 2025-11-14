

source(paste0(here::here(), "/defaria/scripts/functions.R")) 

library(reticulate)
library(dplyr)
library(anndata)
library(ggplot2)
library(patchwork)

py_require(c("anndata"))
# #convert = TRUE To turn python objects into native R objects
ad <- import("anndata", convert = TRUE)


int_etmr <- ad$read_h5ad("/home/sofia/Projects/etmr/defaria/scRNA/jessa/data/processed/all_etmr_int.h5ad")
etmr1 <- readRDS("/home/sofia/Projects/etmr/defaria/scRNA/jessa/data/raw/so_etmr1_anotado.rds")
etmr2 <- ad$read_h5ad("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/paper_data_filt_alone.h5ad")

int_etmr$obs_names <- make.unique(reticulate::py_to_r(int_etmr$obs_names))

int_etmr <- convert_to_seurat(int_etmr)
etmr2 <- convert_to_seurat(etmr2)

umap_etmr1 <- cbind(Embeddings(etmr1, "umap.unintegrated"), etmr1@meta.data)

etmr1_cell_plot <- ggplot(umap_etmr1 %>% rename(CellType = cell_types_2_etmr1), aes(x = umapunintegrated_1, y = umapunintegrated_2, color = CellType)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_classic() +
  scale_color_manual(values = c("#ff9aa2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  

ggsave(etmr1_cell_plot, filename="/home/sofia/Projects/etmr/defaria/scRNA/jessa/plots/celltype_etmr1_original.jpeg", width = 8, height = 4, units = "in", dpi = 600)

etmr1_cell_hist <- umap_etmr1 %>%
  rename(CellType = cell_types_2_etmr1) %>% 
  ggplot(aes(x = seurat_clusters,  fill = CellType)) +
  geom_bar(position = "fill") +
  theme_classic() +
  scale_fill_manual(values = c("#ff9aa2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  ggtitle("Jessa original cell type by cluster")


ggsave(etmr1_cell_hist, filename="/home/sofia/Projects/etmr/defaria/scRNA/jessa/plots/celltype_etmr1_original_hist.jpeg", width = 8, height = 4, units = "in", dpi = 600)





########################################de faria tumor #############################################################################


umap_etmr2 <- cbind(Embeddings(etmr2, "umap"), etmr2@meta.data)

etmr2_cluster_plot <- ggplot(umap_etmr2, aes(x = UMAP_1, y = UMAP_2, color = leiden)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_classic()

ggsave(etmr2_cluster_plot, filename="/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/etmr2_cluster_plot.jpeg", width = 8, height = 4, units = "in", dpi = 600)



###################################################################################################################################


all_umap <- cbind(Embeddings(int_etmr, "umap"), int_etmr@meta.data)

int_etmr_plot <- ggplot(all_umap, aes(x = UMAP_1, y = UMAP_2, color = condition)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#102542", "#f87060"))  +
  theme_classic() +
  ggtitle("ETMR integrated scRNA sample (Jessa as reference)") +
  guides(colour = guide_legend(override.aes = list(size=3)))


ggsave(int_etmr_plot, filename="/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/ingest_etmr.jpeg", width = 8, height = 4, units = "in", dpi = 600)


  
celltype_etmr_plot <- ggplot(all_umap, aes(x = UMAP_1, y = UMAP_2, color = cell_types_2_etmr1)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#ff9aa2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  theme_classic() +
  ggtitle("ETMR integrated scRNA sample (Jessa as reference)") +
  guides(colour = guide_legend(override.aes = list(size=3)))


ggsave(celltype_etmr_plot, filename="/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/celltype_etmr.jpeg", width = 8, height = 4, units = "in", dpi = 600)

celltype_etmr_hist <- all_umap %>%
  ggplot(aes(x = leiden,  fill = cell_types_2_etmr1)) +
  geom_bar(position = "fill",  color = "gray") +
  theme_classic() +
  scale_fill_manual(values = c("#ff9aa2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  ggtitle("ETMR integrated scRNA sample (Jessa as reference)")

ggsave(celltype_etmr_hist, filename="/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/celltype_etmr_hist.jpeg", width = 8, height = 4, units = "in", dpi = 600)


celltype_defaria_hist_etmrint <- all_umap %>%
  filter(condition == "defaria") %>%
  rename(CellType = cell_types_2_etmr1) %>%
  ggplot(aes(x = leiden,  fill = CellType)) +
  geom_bar(position = "fill",  color = "gray") +
  theme_classic() +
  scale_fill_manual(values = c("#ff9aa2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  ggtitle("deFaria tumor cell type using Jessa as reference")

ggsave(celltype_defaria_hist_etmrint, filename="/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/ccelltype_defaria_hist_etmrint.jpeg", width = 8, height = 4, units = "in", dpi = 600)


etmr1_cell_plot_int <- ggplot(all_umap %>% rename(CellType = cell_types_2_etmr1) %>% filter(condition == "jessa"), 
                              aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_classic() +
  scale_color_manual(values = c("#ff9aa2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  

ggsave(etmr1_cell_plot_int, filename="/home/sofia/Projects/etmr/defaria/scRNA/jessa/plots/celltype_etmr1_new_filt.jpeg", width = 8, height = 4, units = "in", dpi = 600)

etmr1_cell_hist_new_filt <- all_umap %>% 
  rename(CellType = cell_types_2_etmr1) %>% 
  filter(condition == "jessa") %>% 
  ggplot(aes(x = leiden,  fill = CellType)) +
  geom_bar(position = "fill",  color = "gray") +
  theme_classic() +
  scale_fill_manual(values = c("#ff9aa2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  ggtitle("Jessa new cell type by cluster")


ggsave(etmr1_cell_hist_new_filt, filename="/home/sofia/Projects/etmr/defaria/scRNA/jessa/plots/etmr1_cell_hist_new_filt.jpeg", width = 8, height = 4, units = "in", dpi = 600)
