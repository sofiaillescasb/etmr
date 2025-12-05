
source(paste0(here::here(), "/defaria/scripts/functions.R")) 

#DO NOT LOAD TIDYVERSE, LEAVE AS IS
library(reticulate)

use_python("/home/sofia/miniconda3/bin/python", required = TRUE) #magic spell and absolute cornerstone

library(anndata)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tibble)
library(Matrix)
library(patchwork)
library(Seurat)
library(SeuratObject)

py_require(c("anndata"))
# #convert = TRUE To turn python objects into native R objects
ad <- import("anndata", convert = TRUE)

etmr_ad <- read_h5ad("/home/sofia/Projects/etmr/defaria/snRNA/data/processed/snrna_data_umap.h5ad")

patient_info <- tribble(
  ~orig.ident, ~Age, ~Sex,
  "GSM8058157_KK22-H-225", 48, "F",
  "GSM8058158_KK22-H-226", 12, "M",
  "GSM8058159_KK22-H-227", 12, "M",
  "GSM8058160_KK23-H-507", 34, "M",
  "GSM8058161_KK23-H-508", 30, "M",
  "GSM8058162_KK23-H-509", 29, "F",
  "GSM8058163_KK23-H-510", 31.2, "M",
  "GSM8058164_KK23-H-511", 35, "F",
  "GSM8058165_KK23-H-512", 39.6, "M")


etmr <- convert_to_seurat(etmr_ad)

etmr@meta.data <- etmr@meta.data %>%
  rownames_to_column() %>%
  left_join(patient_info, by = "orig.ident") %>%
  column_to_rownames()

############################################# Unintegrated integration #####################################################################################

umap_etmr <- cbind(Embeddings(etmr, "umapunintegrated"), etmr@meta.data)

############################################# Harmony integration #####################################################################################

#Use umap instead of harmony because i ran umap after harmony when i integrated the data
umap_harm <- cbind(Embeddings(etmr, "umapharmony"), etmr@meta.data)

  
############################################# umaps ###################################################################################

age_colors <- c("#d16ba5", "#86a8e7","#5ffbf1")
sex_colors <- c("#fdff9a", "#98bfff")
tt_colors <- c("#a294ff", "#a9ff95")
id_colors <-c("#33a8c7", "#52e3e1", "#a0e426", "#fdf148", "#ffab00", "#f77976", "#f050ae", "#d883ff", "#9336fd")

clus_unint <- DimPlot(etmr, reduction="umapunintegrated", group.by="leiden_unintegrated", pt.size = .6, 
                      shuffle = TRUE, label = TRUE, label.size = 6)  + 
  theme(legend.position = "none")

  
clus_harm <- DimPlot(etmr, reduction="umapharmony", group.by="leiden_harmony", pt.size = .6, 
                     shuffle = TRUE, label = TRUE, label.size = 6)  + 
  theme(legend.position = "none")

clusters_plot <- clus_unint + clus_harm

age_plots <- plot_with_patchwork("Age", "Sex","M", "F") &
  scale_color_gradientn(colors = age_colors) 
  
sex_plots <- plot_with_patchwork("Sex", "Sex","M", "F") &
  guides(colour = guide_legend(override.aes = list(size=2))) &
  scale_color_manual(values = sex_colors) 


tumor_type_plot <- plot_with_patchwork("Tumor_Type", "Tumor_Type","C19MC", "Dicer") &
  guides(colour = guide_legend(override.aes = list(size=2))) &
  scale_color_manual(values = tt_colors) 


id_plot <- plot_with_patchwork("orig.ident",  "Tumor_Type", "Dicer","C19MC") &
  guides(colour = guide_legend(override.aes = list(size=2))) &
  scale_color_manual(values =id_colors) 

all_plots <- clusters_plot / id_plot / tumor_type_plot / age_plots 

ggsave(all_plots, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/all_plots.jpeg", width = 12, height = 12, units = "in", dpi = 600)

######################################################### histograms #############################################################################

id_hist <- hist_with_patchwork("orig.ident") &
  scale_fill_manual(values = id_colors) 

sex_hist <- hist_with_patchwork("Sex") &
  scale_fill_manual(values = sex_colors)

tumor_type_hist <- hist_with_patchwork("Tumor_Type") &
  scale_fill_manual(values = tt_colors) 

umap_etmr <- umap_etmr %>% mutate(Age = ordered(factor(ifelse(Age < 20 | Age > 40 | is.factor(Age), Age, 
                                                              ifelse(Age <= 36, "< 36", "< 48"))), c("12", "< 36", "< 48", "48")))
umap_harm <- umap_harm %>% mutate(Age = ordered(factor(ifelse(Age < 20 | Age > 40 | is.factor(Age), Age, 
                                                              ifelse(Age <= 36, "< 36", "< 48"))), c("12", "< 36", "< 48", "48")))

age_colors <- c("#d16ba5", "#86a8e7", "#73d2ec","#5ffbf1")

age_hist <- hist_with_patchwork("Age") &
  scale_fill_manual(values = age_colors)

all_hists <- id_hist  / tumor_type_hist / sex_hist / age_hist 

ggsave(all_hists, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/all_hists.jpeg", width = 12, height = 12, units = "in", dpi = 600)

saveRDS(etmr, "/home/sofia/Projects/etmr/defaria/snRNA/data/processed/snrna_data.rds")


