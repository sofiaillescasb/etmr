
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

etmr_tumortype_plot <- plot_general(umap_etmr, "umapunintegrated_1", "umapunintegrated_2", "Tumor_Type") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values = c("#6c9a8b", "#e8998d")) 

ggsave(etmr_tumortype_plot, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_tumortype_plot_unint.jpeg", width = 8, height = 4, units = "in", dpi = 600)

etmr_unint_tumortype_hist <- plot_hist(umap_etmr, "leiden_unintegrated", "Tumor_Type") +
  scale_fill_manual(values = c("#6c9a8b", "#e8998d"))  +
  ggtitle("Unintegrated clusters")  

ggsave(etmr_unint_tumortype_hist, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_tumortype_hist_unint.jpeg", width = 8, height = 4, units = "in", dpi = 600)

############################################# Harmony integration #####################################################################################

#Use umap instead of harmony because i ran umap after harmony when i integrated the data
umap_harm <- cbind(Embeddings(etmr, "umap_harmony"), etmr@meta.data)

etmr_tumortype_harm_plot <- plot_general(umap_harm, "umapharmony_1", "umapharmony_2", "Tumor_Type") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values = c("#6c9a8b", "#e8998d")) 
  

ggsave(etmr_tumortype_harm_plot, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_tumortype_plot_harm.jpeg", width = 8, height = 4, units = "in", dpi = 600)

etmr_harm_tumortype_hist <- plot_hist(umap_harm, "leiden_harmony", "Tumor_Type") +
  scale_fill_manual(values = c("#6c9a8b", "#e8998d"))  +
  ggtitle("Harmony-integrated clusters")  


ggsave(etmr_harm_tumortype_hist, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_tumortype_hist_harm.jpeg", width = 8, height = 4, units = "in", dpi = 600)

############################################# Scanorama integration ###################################################################################

umap_scan <- cbind(Embeddings(etmr, "umap_scanorama"), etmr@meta.data)
  
etmr_tumortype_scanorama_plot <- plot_general(umap_scan, "umapscanorama_1", "umapscanorama_2", "Tumor_Type") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values = c("#6c9a8b", "#e8998d")) 

ggsave(etmr_tumortype_scanorama_plot, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_tumortype_scanorama_plot.jpeg", width = 8, height = 4, units = "in", dpi = 600)

etmr_scan_tumortype_hist <- plot_hist(umap_scan, "leiden_scanorama", "Tumor_Type") +
  scale_fill_manual(values = c("#6c9a8b", "#e8998d"))  +
  ggtitle("Scanorama-integrated clusters")  

ggsave(etmr_scan_tumortype_hist, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_scan_tumortype_hist.jpeg", width = 8, height = 4, units = "in", dpi = 600)

############################################# joined integration plots ###################################################################################

joined_etmr_int_plots <- (etmr_tumortype_plot + ggtitle("No integration")) + 
  (etmr_tumortype_harm_plot + ggtitle("Harmony integration")) + (etmr_tumortype_scanorama_plot + ggtitle("Scanorama integration")) + 
  plot_layout(guides = "collect") + plot_annotation(title = "Integration of tumor types in ETMR snRNA data")

ggsave(joined_etmr_int_plots, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/joined_etmr_int_plots.jpeg", width = 12, height = 4, units = "in", dpi = 600)

joined_etmr_int_hist <- etmr_unint_tumortype_hist + etmr_harm_tumortype_hist + etmr_scan_tumortype_hist + 
  plot_layout(guides = "collect") 

ggsave(joined_etmr_int_hist, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/joined_etmr_int_hist.jpeg", width = 12, height = 4, units = "in", dpi = 600)

############################################# Age ###################################################################################
age_plots <- plot_one_over_other(umap_etmr, "umapunintegrated_1", "umapunintegrated_2", "Age", "Sex","M", "F") +
  plot_one_over_other(umap_harm, "umapharmony_1", "umapharmony_2", "Age", "Sex","M", "F") +
  plot_one_over_other(umap_scan, "umapscanorama_1", "umapscanorama_2", "Age", "Sex","M", "F") +
  plot_layout(guides = "collect") &
  scale_color_gradientn(colors = c("#d16ba5", "#86a8e7", "#5ffbf1")) & 
  plot_annotation(title = "Age distribution accross integration methods in ETMR snRNA data") 
  

sex_plots <- plot_one_over_other(umap_etmr, "umapunintegrated_1", "umapunintegrated_2", "Sex", "Sex","M", "F") + 
  plot_one_over_other(umap_harm, "umapharmony_1", "umapharmony_2", "Sex", "Sex","M", "F") + 
  plot_one_over_other(umap_scan, "umapscanorama_1", "umapscanorama_2", "Sex", "Sex","M", "F") +
  plot_annotation(title = "Sex distribution accross integration methods in ETMR snRNA data") +
  plot_layout(guides = "collect") &
  guides(colour = guide_legend(override.aes = list(size=3))) &
  scale_color_manual(values = c("#fdff9a", "#98bfff")) 


tumor_type_plot <- plot_one_over_other(umap_etmr, "umapunintegrated_1", "umapunintegrated_2", "Tumor_Type", "Tumor_Type","C19MC", "Dicer") +
  plot_one_over_other(umap_harm, "umapharmony_1", "umapharmony_2", "Tumor_Type", "Tumor_Type","C19MC", "Dicer") +
  plot_one_over_other(umap_scan, "umapscanorama_1", "umapscanorama_2", "Tumor_Type", "Tumor_Type","C19MC", "Dicer") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Tumor type distribution accross integration methods in ETMR snRNA data") &
  guides(colour = guide_legend(override.aes = list(size=3))) &
  scale_color_manual(values = c("#a294ff", "#a9ff95")) 

  

tumor_type_plot / sex_plots / age_plots 


etmr_age_scanorama_plot <- plot_general(umap_scan, "umapscanorama_1", "umapscanorama_2", "Age") 

ggsave(etmr_tumortype_scanorama_plot, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_tumortype_scanorama_plot.jpeg", width = 8, height = 4, units = "in", dpi = 600)

etmr_scan_tumortype_hist <- plot_hist(umap_scan, "leiden_scanorama", "Tumor_Type") +
  scale_fill_manual(values = c("#6c9a8b", "#e8998d"))  +
  ggtitle("Scanorama-integrated clusters")

ggsave(etmr_scan_tumortype_hist, filename="/home/sofia/Projects/etmr/defaria/snRNA/plots/etmr_scan_tumortype_hist.jpeg", width = 8, height = 4, units = "in", dpi = 600)



