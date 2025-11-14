
#!/usr/bin/env Rscript

source(paste0(here::here(), "/defaria/scripts/functions.R")) 

library(reticulate)
library(dplyr)
library(anndata)
library(ggplot2)
library(patchwork)

#internet says to use package:: instead of library

# args <- commandArgs(trailingOnly = TRUE)
# 
# 
py_require(c("anndata"))
# #convert = TRUE To turn python objects into native R objects
ad <- import("anndata", convert = TRUE)
# 
# print("Reading data...")

# etmr_ad <- ad$read_h5ad(args[1])
# ref_ad <- ad$read_h5ad(args[2])
# all_data_ad <- ad$read_h5ad(args[3])

etmr_ad <- ad$read_h5ad("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/paper_data_filt_comvar_int.h5ad") 
ref_ad <- ad$read_h5ad("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/merged_ref_paper_data_filt_comvar.h5ad")
all_data_ad <- ad$read_h5ad("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/all_data_int.h5ad")


print("Converting Anndata to Seurat...")

ref_data <- convert_to_seurat(ref_ad)
etmr <- convert_to_seurat(etmr_ad)
all_data <- convert_to_seurat(all_data_ad)

all_meta <- all_data@meta.data 


print("Making UMAP for plotting...")
all_umap <- cbind(Embeddings(all_data, "umap"), all_meta)
  

print("Plotting integration")


int_plot <- plot_both("condition") +
  scale_color_manual(values = c("grey", "#c05299")) +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/integration_plot.jpeg", int_plot, width = 6, height = 3, units = "in", dpi = 600)

# ggsave(paste0(dirname(dirname(dirname(as.character(args[1])))), "/plots/integration_plot.jpeg"), int_plot, width = 6, height = 3, units = "in", dpi = 600)

print("Plotting CellType")


all_cell_plot <- plot_both("CellClass") +
  scale_color_manual(values = c("#ff9aa2", "#ffb7b2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea")) 


ggsave("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/all_data_cell_class.jpeg", all_cell_plot, width = 6, height = 3, units = "in", dpi = 600)
# ggsave(paste0(dirname(dirname(dirname(as.character(args[1])))), "/plots/all_data_cell_class.jpeg"), all_cell_plot, width = 6, height = 3, units = "in", dpi = 600)

all_cell_hist <- ggplot(all_umap, aes(x = leiden, fill = CellClass)) +
  geom_bar(position = "fill",  color = "gray") +
  scale_fill_manual(values = c("#ff9aa2", "#ffb7b2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  theme_classic() +
  ggtitle("de Faria tumor cell type according to healthy reference")

ggsave("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/all_cell_hist_ref.jpeg", all_cell_hist, width = 6, height = 3, units = "in", dpi = 600)



etmr_cell_plot <- ggplot(all_umap %>% filter(condition == "etmr"), aes(x = UMAP_1, y = UMAP_2, color = CellClass)) +
  geom_point(size = 0.6) +
  scale_color_manual(values = c("#ff9aa2", "#ffb7b2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  theme_classic() +
  ggtitle("de Faria tumor cell type according to healthy reference")

ggsave("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/etmr_cell_plot.jpeg", etmr_cell_plot, width = 6, height = 3, units = "in", dpi = 600)

etmr_cell_hist <- ggplot(all_umap %>% filter(condition == "etmr"), aes(x = leiden, fill = CellClass)) +
  geom_bar(position = "fill",  color = "gray") +
  scale_fill_manual(values = c("#ff9aa2", "#ffb7b2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  theme_classic() +
  ggtitle("de Faria tumor cell type according to healthy reference")

ggsave("/home/sofia/Projects/etmr/defaria/scRNA/paper_data/plots/etmr_cell_hist_ref.jpeg", etmr_cell_hist, width = 6, height = 3, units = "in", dpi = 600)

#Save plots individually
cell_plot <- all_cell_plot + etmr_cell_plot + plot_layout(guides = "collect")  + plot_annotation(title = "Cell type in healthy cells and ETMR snRNA data")

ggsave(paste0(dirname(dirname(dirname(as.character(args[1])))), "/plots/cell_class.jpeg"), cell_plot, width = 12, height = 4, units = "in", dpi = 600)

print("Plotting Age")

all_age_plot <- plot_both("Age") +
  geom_bar(position = "fill", alpha= .7,color = "darkgrey")  +
  scale_fill_manual(values = c("#c9cba3","#ffe1a8", "#e26d5c", "#723d46")) +
  theme(legend.position = "none")


etmr_age_plot <- ggplot(all_umap %>% filter(condition == "etmr"), aes(x = UMAP_1, y = UMAP_2, color = Age)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_bar(position = "fill", alpha= .7,color = "darkgrey")  +
  scale_fill_manual(values = c("#c9cba3","#ffe1a8", "#e26d5c", "#723d46")) +
  theme_classic() 


age_plot <- int_plot + all_age_plot + etmr_age_plot + plot_layout(guides = "collect")  + plot_annotation(title = "Age in healthy cells and ETMR snRNA data")


ggsave(paste0(dirname(dirname(dirname(as.character(args[1])))), "/plots/age.jpeg"), age_plot, width = 12, height = 4, units = "in", dpi = 600)

etmr@meta.data <- etmr@meta.data %>%
  mutate(orig.ident = ifelse("orig.ident" %in% colnames(.) & length(unique(.$orig.ident))>1, .$orig.ident, "ETMR")) %>%
  rename(donor_id = orig.ident)


ids <- rbind(ref_data@meta.data["donor_id"], etmr@meta.data["donor_id"]) %>%
  rownames_to_column("cell_id")

all_umap <- all_umap %>%
  rownames_to_column("cell_id") %>%
  left_join(ids, by = "cell_id") %>%
  column_to_rownames("cell_id")



id_plot_all <- plot_both("donor_id")

patch_id <- int_plot + id_plot_all + plot_layout(guides = "collect")  + plot_annotation(title = "Donor ID in healthy cells and ETMR")


ggsave(paste0(dirname(dirname(dirname(as.character(args[1])))), "/plots/id_plot_all.jpeg"), patch_id, width = 8, height = 4, units = "in", dpi = 600)

#I'm gonna order the ages from here, maybe later change it so i never convert them to numeric in the first place

# scale_color_viridis_d()

all_meta <- all_meta %>%
  rownames_to_column("cell_id") %>%
  left_join(ids, by = "cell_id") %>%
  column_to_rownames("cell_id")


all_meta %>%
  arrange(donor_id) %>%
  ggplot(aes(x = leiden, fill = donor_id)) +
  geom_bar(position = "fill") 


all_meta %>%
  arrange(donor_id) %>%
  ggplot(aes(x = leiden, fill = condition)) +
  geom_bar(position = "fill") 


all_meta %>%
  filter(condition== "etmr") %>%
  mutate(Age = ordered(as.character(Age), levels = c("5", "8", "11.5", "14"))) %>%
  arrange(Age) %>%
  ggplot(aes(x = leiden, fill = Age)) +
  geom_bar(position = "fill", alpha= .7,color = "darkgrey")  +
  scale_fill_manual(values = c("#c9cba3","#ffe1a8", "#e26d5c", "#723d46")) +
  theme_classic() +
  ggtitle("ETMR data: scanpy integration")

all_meta %>%
  mutate(Age = ordered(as.character(Age), levels = c("5", "8", "11.5", "14"))) %>%
  arrange(Age) %>%
  ggplot(aes(x = leiden, fill = Age)) +
  geom_bar(position = "fill", color = "darkgrey")  +
  scale_fill_manual(values = c("#c9cba3","#ffe1a8", "#e26d5c", "#723d46")) +
  theme_classic() +
  ggtitle("Reference and ETMR: scanpy integration")


all_meta %>%
  filter(condition== "etmr") %>%
  ggplot(aes(x = leiden,  fill = CellClass)) +
  geom_bar(position = "fill",  color = "gray") +
  theme_classic() +
  scale_fill_manual(values = c("#ff9aa2", "#ffb7b2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  ggtitle("ETMR data: scanpy integration")


cluster_hist <- function(colvar) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_bar(aes())
  ggtheme::theme_classic() +
    ggplot2::guides(colour = guide_legend(override.aes = list(size=3)))
  return(p)}