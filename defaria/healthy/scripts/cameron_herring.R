library(R.utils)
library(data.table)
library(Seurat)
library(qs2)

#https://www.sciencedirect.com/science/article/pii/S0006322322014044
#FC is frontal lobal
gz_path <- "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/cameron_2022_snRNAseq_FC_raw_count_gEX_matrix.txt.gz"

#unzip file
gunzip(gz_path, remove = FALSE, overwrite = TRUE)

txt_path <- sub("\\.gz$", "", gz_path)

#load the uncompressed count table efficiently
counts_dt <- fread(txt_path)


# Convert the first column (gene symbols/IDs) to row names
gene_names <- counts_dt[[1]]
count_matrix <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(count_matrix) <- gene_names

# 3 karyotypically normal human fetuses from the second trimester of gestation (2 of 14 and 1 of 15 postconception weeks, all female).

# Initialize the Seurat object with raw counts
prenatal <- CreateSeuratObject(
  counts = count_matrix, 
  project = "Prenatal", 
  min.cells = 3, 
  min.features = 200
)

# Read metadata
metadata <- fread("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/cameron_2022_snRNAseq_FC_metadata.txt.gz")


prenatal <- AddMetaData(prenatal, metadata = as.data.frame(metadata))

prenatal[["RNA"]] <- split(prenatal[["RNA"]], f = prenatal$sample) 

prenatal <- NormalizeData(prenatal)

prenatal <- FindVariableFeatures(prenatal)

prenatal <- ScaleData(prenatal)

prenatal$condition <- "healthy_cameron"

gc()


healthy_int <- qs_read(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/healthy_int.qs2", 
  nthreads = 6
)

#Area 8 – Includes frontal eye fields
#Area 9 – Dorsolateral prefrontal cortex

gc()

healthy_int <- DietSeurat(healthy_int,
                          assays = "RNA")

healthy_int$condition <- "healthy_herring"

her_cam <- merge(healthy_int, prenatal)
gc()

rm(healthy_int)
gc()

#get common variable features and do pca
her_cam <- FindVariableFeatures(her_cam)
her_cam <- ScaleData(her_cam)
her_cam <- RunPCA(her_cam)

her_cam_int <- IntegrateLayers(her_cam,
                               method = SeuratWrappers::FastMNNIntegration,
                               assay = "RNA")

her_cam_int <- RunUMAP(her_cam_int,
                   reduction = "integrated.mnn",
                   reduction.name = "mnn_umap",
                   return.model = TRUE,
                   dims = 1:30
                   )

DimPlot(her_cam_int, reduction = "mnn_umap", group.by = "condition")
DimPlot(her_cam_int, reduction = "mnn_umap", group.by = c("cellIDs", "major_clust"), label = TRUE) + NoLegend()

library(tidyverse)

converted_cell_type <- c("FC-InN-2" = "IN",   
"FC-IP" = "Non-Neu",
"FC-RG-1" = "Non-Neu",
"FC-ExN-1" = "PN",
"FC-ExN-4" = "PN",  
"FC-InN-1" = "IN",
"FC-ExN-3" = "PN",
"FC-InN-3" = "IN",
"FC-RG-2" = "Non-Neu",
"FC-N-undef" = "Neu-undefined",
"FC-ExN-5" = "PN",
"FC-ExN-2" = "PN",
"FC-CycPro" = "Non-Neu",
"FC-InN-4" = "IN",
"FC-Endo" = "Non-Neu",
"FC-MG" = "Non-Neu",
"FC-OPC" = "Non-Neu") 

#Clusters were annotated based on expression of known cell markers (Figures S1–S5 in Supplement 2) 
# and other differentially expressed genes (Tables S2–S6 in Supplement 1). 
# In addition to populations of radial glia (markers: GLI3, TNC), intermediate progenitors (marker: EOMES), 
# oligodendrocyte precursor cells (marker: OLIG1), endothelial cells (marker: ITM2A), 
# and microglia (marker: C3), we were able to define a variety of developing neuron populations
# appropriate for each region. Thus, in the FC, we found developing excitatory (glutamatergic)
# neuron populations differentially expressing markers of upper cortical layers, such as TLE3 
# and LHX2 (FC-ExN-1, FC-ExN-2, and FX-ExN-3), or markers of deep cortical layers, such as CRYM and FEZF2 
# (FC-ExN-4 and FC-ExN-5), as well as developing GABAergic (gamma-aminobutyric acidergic) 
# interneuron populations expressing CALB2 (FC-InN-1) or SST (FC-InN-4). 

# radial glia = c("GLI3", "TNC")
# IP = c("EOMES")
# OPC =  c("OLIG1")
# Endo = c("ITM2A")
# Microglia = c("C3")
# PN-dev_upper = c("TLE3", "LHX2") #(FC-ExN-1, FC-ExN-2, and FX-ExN-3)
# PN_dev_deep = c("CRYM", "FEZF2") 
# developing GABAergic (gamma-aminobutyric acidergic) interneuron populations expressing CALB2 (FC-InN-1) or SST (FC-InN-4)

converted_maj <- c("FC-InN-2" = "InN-2",   
                         "FC-IP" = "IP",
                         "FC-RG-1" = "RG-1",
                         "FC-ExN-1" = "PN_dev_upper-1",
                         "FC-ExN-4" = "PN_dev_deep-4",  
                         "FC-InN-1" = "IN_dev_CALB2",
                         "FC-ExN-3" = "PN_dev_upper-3",
                         "FC-InN-3" = "InN-3",
                         "FC-RG-2" = "RG_2",
                         "FC-N-undef" = "Neu-undefined",
                         "FC-ExN-5" = "PN_dev_deep-5",
                         "FC-ExN-2" = "PN_dev_upper-2",
                         "FC-CycPro" = "Cyc_Prog_EOMES",
                         "FC-InN-4" = "SST",
                         "FC-Endo" = "Vas",
                         "FC-MG" = "Micro",
                         "FC-OPC" = "OPC") 



her_cam_int@meta.data  <- her_cam_int@meta.data |>
  mutate(
    cell_type_joined = ifelse(
      is.na(cell_type), 
      converted_cell_type[cellIDs], 
      cell_type),
    major_clust_joined = ifelse(
      is.na(major_clust), 
      converted_maj[cellIDs], 
      major_clust),
    sub_clust_joined = ifelse(
      is.na(major_clust), 
      converted_maj[cellIDs], 
      sub_clust),
    stage_id = ifelse(
      !is.na(sample), 
      "Fetal", 
      stage_id),
    orig.ident = ifelse(
      !is.na(sample),
      sample,
      orig.ident),
    dataset = ifelse(
      sample %in% c("510_PFC_B1", "611_PFC_B1", "993_PFC_B2"),
      "Cameron",
      "Herring")
    )

#L5-6_TLE4

FeaturePlot(
  her_cam_int, 
  features = c("CRYM", "FEZF2"), 
  order = TRUE
  ) / 
  (DimPlot(
    her_cam_int, 
    reduction = "mnn_umap", 
    group.by = c("major_clust", "cellIDs"),
    label = TRUE)
     ) & NoLegend()

FeaturePlot(
  her_cam_int, 
  features = c("CRYM", "FEZF2"), 
  order = TRUE
) / 
  (DimPlot(
    her_cam_int, 
    reduction = "mnn_umap", 
    group.by = c("major_clust", "cellIDs"),
    label = TRUE)
  ) & NoLegend()

healthy_plot_p1 <- DimPlot(
  her_cam_int, 
  reduction = "mnn_umap", 
  group.by = c("dataset", "stage_id","cell_type_joined"),
  label = TRUE,
  shuffle = TRUE
  ) &
  NoLegend()

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/healthy_plot_p1.png", 
       healthy_plot_p1,
       units="in",
       width = 21,
       height = 8)

healthy_plot_p2 <- DimPlot(
  her_cam_int, 
  reduction = "mnn_umap", 
  group.by = c("major_clust_joined"), 
  label = TRUE
  ) & NoLegend()

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/healthy_plot_p2.png", 
       healthy_plot_p2,
       units="in",
       width = 8,
       height = 8)

qs_save(her_cam_int, "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/her_cam_int.qs2", nthreads = 6)

her_cam_int$condition <- "human_healthy"

her_cam_int <- JoinLayers(her_cam_int)

gc()

tumor <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/tumor_herr_lt.qs2", 
                 nthreads = 6)

colnames(tumor@meta.data)

tumor <- JoinLayers(tumor)

gc()

anchors.mnn <- FindTransferAnchors(
  reference = her_cam_int,
  query = tumor,
  normalization.method = "LogNormalize", 
  reduction = "pcaproject",
  reference.reduction = "integrated.mnn"
)

qs_save(anchors.mnn, "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/anchors.mnn_herr_cam_lt.qs2", 
        nthreads = 6)
gc()

full_metadata <- her_cam_int@meta.data
#later
# her_cam_int@meta.data |>
#   select("orig.ident", "nCount_RNA"             
#          "nFeature_RNA"            "percent.mt"             
#          "X"                       "batch"                  
#                     "age"                    
#          "chem"                              
#          "numerical_age"           "stage_id"               
#          "Sex"                     "Race"                   
#          "PMI"                     "Brain.Regions."         
#           "n_genes_by_counts"       "log1p_n_genes_by_counts"
#           "total_counts"            "log1p_total_counts"     
#           "doublet_score"           "log10_gene_counts"      
#           "log10_UMI_counts"        "percent_mito"           
#           "percent_ribo"            "n_counts"               
#           "stripped_name"           "RNA_snn_res.0.8"      "RNA_snn_res.1"          
#           "condition"               "cells"                  
#             "sample"                 
#           "cell_type_joined"        "major_clust_joined"     
#           "sub_clust_joined")

tumor <- MapQuery(
  anchorset = anchors.mnn,
  reference = her_cam_int,
  query = tumor,
  refdata = list(celltype_cam = "cell_type_joined",
                 major_clust_cam = "major_clust_joined",
                 sub_clust_cam = "sub_clust_joined",
                 stage_id_cam = "stage_id"), # Change this to your actual cell type column name in healthy_int
  reference.reduction = "integrated.mnn", 
  reduction.model = "mnn_umap",
  new.reduction.name = "her_cam_umap"
)

Idents(tumor) <- "orig.ident"

age <- c(
  "GSM8058157" = 31.2, 
  "GSM8058158" = 12, 
  "GSM8058159" = 30, 
  "GSM8058160" = 29, 
  "GSM8058161" = 12, 
  "GSM8058162" = 35, 
  "GSM8058163" = 48, 
  "GSM8058164" = 34, 
  "GSM8058165"= 39.6)

tumor@meta.data["age"]  <- age[Idents(tumor)]
tumor$predicted.stage_id <- factor(tumor$predicted.stage_id, levels = c("Fetal", "Neonatal", "Infancy", "Childhood"))

compare_major_clust_stage <- (DimPlot(tumor, 
        group.by = c("age", "by_group_celltypes", "predicted.stage_id", "predicted.major_clust_cam"), 
        shuffle = TRUE)) 

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/compare_major_clust_stage.png", 
       compare_major_clust_stage,
       units="in",
       width = 12,
       height = 10)

(DimPlot(tumor, group.by = c("predicted.major_clust_cam"), label = TRUE) +
          DimPlot(her_cam_int, group.by = c("major_clust_joined"), label = TRUE)) & NoLegend()



(DimPlot(tumor, group.by = c("predicted.celltype_cam"), label = TRUE) +
    DimPlot(her_cam_int, group.by = c("cell_type_joined"), label = TRUE)) 

celltype_mnn_tumor_score <- tumor@meta.data |>
  group_by(predicted.celltype_cam) |>
  filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(x = predicted.celltype_cam.score)) +
  geom_histogram(binwidth = .05, fill = "#ffc8c8", color = "#ff5858") +
  theme_classic() +
  facet_wrap(vars(predicted.celltype_cam), scales = "free_y")  +
  ggtitle("Cell type score tumor-healthy")

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/celltype_mnn_tumor_score.png", 
       celltype_mnn_tumor_score,
       units="in",
       width = 12,
       height = 6)


major_mnn_tumor_score <- tumor@meta.data |>
  group_by(predicted.major_clust_cam) |>
  filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(x = predicted.major_clust_cam.score)) +
  geom_histogram(binwidth = .05,  fill = "#83c5be", color = "#006d77") +
  theme_classic() +
  facet_wrap(vars(predicted.major_clust_cam), scales = "free_y")  +
  ggtitle("Major clust score tumor-healthy")

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/major_mnn_tumor_score.png", 
       major_mnn_tumor_score,
       units="in",
       width = 12,
       height = 8)

stage_id_mnn_tumor_score <- tumor@meta.data |>
  group_by(predicted.stage_id_cam) |>
  filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(x = predicted.stage_id_cam.score)) +
  geom_histogram(binwidth = .05,  fill = "lightgreen", color = "darkgreen") +
  theme_classic() +
  facet_wrap(vars(predicted.stage_id_cam), scales = "free_y")  +
  ggtitle("Developmental stage score tumor-healthy")

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/stage_id_mnn_tumor_score.png", 
       stage_id_mnn_tumor_score,
       units="in",
       width = 10,
       height = 8)

summary(tumor$predicted.celltype.score)
gc()
