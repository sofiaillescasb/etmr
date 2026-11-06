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

# Initialize the Seurat object with raw counts
prenatal <- CreateSeuratObject(
  counts = count_matrix, 
  project = "Prenatal", 
  min.cells = 3, 
  min.features = 200
)

# Read metadata
metadata <- fread("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/cameron_2022_snRNAseq_FC_metadata.txt.gz")

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

 
#L5-6_TLE4
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
