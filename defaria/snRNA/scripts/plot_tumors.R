
source(paste0(here::here(), "/defaria/scripts/functions.R")) 

library(reticulate)
library(dplyr)
library(anndata)
library(ggplot2)
library(patchwork)

py_require(c("anndata"))
# #convert = TRUE To turn python objects into native R objects
ad <- import("anndata", convert = TRUE)

etmr_ad  <- ad$read_h5ad("/home/sofia/Projects/etmr/defaria/snRNA/data/processed/snrna_data_umap_no_integration.h5ad")

patient_info <- tribble(
  ~orig.ident, ~Age, ~Sex,
  "GSM8058157", 48, "F",
  "GSM8058158", 12, "M",
  "GSM8058159", 12, "M",
  "GSM8058160", 34, "M",
  "GSM8058161", 30, "M",
  "GSM8058162", 29, "F",
  "GSM8058163", 31.2, "M",
  "GSM8058164", 35, "F",
  "GSM8058165", 39.6, "M")

etmr <- convert_to_seurat(etmr_ad)

etmr@meta.data <- etmr@meta.data   %>%
  left_join(patient_info, by = "orig.ident")

############################################# Unintegrated integration #####################################################################################

umap_etmr <- cbind(Embeddings(etmr, "umap_unintegrated"), etmr@meta.data)

############################################# Harmony integration #####################################################################################

#Use umap instead of harmony because i ran umap after harmony when i integrated the data
umap_harm <- cbind(Embeddings(etmr, "umap_harmony"), etmr@meta.data)

############################################# Scanorama integration ###################################################################################

umap_scan <- cbind(Embeddings(etmr, "umap_scanorama"), etmr@meta.data)
  
############################################# umaps ###################################################################################

age_colors <- c("#d16ba5", "#86a8e7", "#5ffbf1")
sex_colors <- c("#fdff9a", "#98bfff")
tt_colors <- c("#a294ff", "#a9ff95")
id_colors <-c("#33a8c7", "#52e3e1", "#a0e426", "#fdf148", "#ffab00", "#f77976", "#f050ae", "#d883ff", "#9336fd")


age_plots <- plot_with_patchwork("Age", "Sex","M", "F") &
  scale_color_gradientn(colors = age_colors) 
  
sex_plots <- plot_with_patchwork("Sex", "Sex","M", "F") &
  guides(colour = guide_legend(override.aes = list(size=3))) &
  scale_color_manual(values = sex_colors) 


tumor_type_plot <- plot_with_patchwork("Tumor_Type", "Tumor_Type","C19MC", "Dicer") &
  guides(colour = guide_legend(override.aes = list(size=3))) &
  scale_color_manual(values = tt_colors) 


id_plot <- plot_with_patchwork("orig.ident", "Sex","M", "F") &
  guides(colour = guide_legend(override.aes = list(size=3))) &
  scale_color_manual(values =id_colors) 

all_plots <- id_plot / tumor_type_plot / sex_plots / age_plots 

ggsave(all_plots, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/all_plots.jpeg", width = 12, height = 12, units = "in", dpi = 600)

############################################# histograms ###################################################################################

id_hist <- hist_with_patchwork("orig.ident") &
  scale_fill_manual(values = id_colors) 

sex_hist <- hist_with_patchwork("Sex") &
  scale_fill_manual(values = sex_colors)

tumor_type_hist <- hist_with_patchwork("Tumor_Type") &
  scale_fill_manual(values = tt_colors) 

umap_etmr <- umap_etmr %>% mutate(Age = ordered(factor(ifelse(Age < 20, "12", ifelse(Age > 40, "48", "20 - 40"))), c("12", "20 - 40", "48")))
umap_harm <- umap_harm %>% mutate(Age = ordered(factor(ifelse(Age < 20, "12", ifelse(Age > 40, "48", "20 - 40"))), c("12", "20 - 40", "48")))
umap_scan <- umap_scan %>% mutate(Age = ordered(factor(ifelse(Age < 20, "12", ifelse(Age > 40, "48", "20 - 40"))), c("12", "20 - 40", "48")))

age_hist <- hist_with_patchwork("Age") &
  scale_fill_manual(values = age_colors)

all_hists <- id_hist  / tumor_type_hist / sex_hist / age_hist 

ggsave(all_hists, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/all_hists.jpeg", width = 12, height = 12, units = "in", dpi = 600)

########################################################## dotplot ###############################################################################


