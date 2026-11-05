library(qs2)
library(Seurat)
library(tidyverse)


tumor <- qs_read(
  "/home/sofia/Projects/etmr/defaria/snRNA/data/processed/preprocessed/tumor_normalized.qs2", 
  nthreads = 6
)

DefaultAssay(tumor) <- "RNA"

tumor <- DietSeurat(tumor, 
                    assays = "RNA",
                    layers = c("counts", "data"))

healthy_int <- qs_read(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/healthy_int.qs2", 
  nthreads = 6
)

healthy_int <- DietSeurat(healthy_int, 
                    assays = "RNA",
                    layers = c("counts", "data"))

gc()

tumor_healthy <- merge(tumor, healthy_int)


metadata <- tumor_healthy@meta.data



rm(tumor, healthy_int)
gc()

tumor_healthy <- FindVariableFeatures(tumor_healthy)

tumor_healthy <- ScaleData(tumor_healthy)

tumor_healthy <- RunPCA(tumor_healthy)

gc()

tumor_healthy_int <- IntegrateLayers(tumor_healthy, 
                                     method = SeuratWrappers::FastMNNIntegration,
                                     assay = "RNA")

rm(tumor_healthy)

gc()

tumor_healthy_int <- IntegrateLayers(tumor_healthy_int, 
                                     method = HarmonyIntegration,
                                     assay = "RNA")


tumor_healthy_int <- JoinLayers(tumor_healthy_int)

gc()


tumor_healthy_int <- RunUMAP(tumor_healthy_int, 
                             reduction = "integrated.mnn", 
                             dims = 1:30,
                             reduction.name = "mnn_umap")


tumor_healthy_int <- RunUMAP(tumor_healthy_int, 
                             reduction = "harmony", 
                             dims = 1:30,
                             reduction.name = "harmony_umap")

tumor_healthy_int$major_clust <- coalesce(
  tumor_healthy_int$by_group_celltypes, 
  tumor_healthy_int$major_clust
  )

tumor_healthy_int$sub_clust <- coalesce(
  tumor_healthy_int$by_group_celltypes, 
  tumor_healthy_int$sub_clust
  )

tumor_healthy_int$cell_type <- coalesce(
  tumor_healthy_int$by_group_celltypes, 
  tumor_healthy_int$cell_type
  )


DimPlot(tumor_healthy_int,
        reduction = "harmony_umap",
        group.by = "orig.ident") +
  DimPlot(tumor_healthy_int,
          reduction = "mnn_umap",
          group.by = "orig.ident")

(DimPlot(tumor_healthy_int,
        reduction = "harmony_umap",
        group.by = "major_clust",
        label = TRUE) +
  DimPlot(tumor_healthy_int,
          reduction = "mnn_umap",
          group.by = "major_clust",
          label = TRUE)) &
  NoLegend()

source("~/Projects/etmr/defaria/scripts/functions.R")

plot_hist(tumor_healthy_int@meta.data, "major_clust", "orig.ident", "orig.ident")



tumor_healthy_int <- FindNeighbors(
  tumor_healthy_int,
  reduction = "integrated.mnn",
  dims = 1:30,
  graph.name = c("mnn_nn", "mnn_snn")
)

tumor_healthy_int <- FindNeighbors(
  tumor_healthy_int,
  reduction = "harmony",
  dims = 1:30,
  graph.name = c("harmony_nn", "harmony_snn")
)

gc()

qs_save(tumor_healthy_int, 
        "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/tumor_healthy_int.qs2",
        nthreads = 6)

###################################################################################################
library(qs2)
library(Seurat)
library(tidyverse)

tumor_healthy_int <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/tumor_healthy_int.qs2",
        nthreads = 6)


for (i in c(0.6, 0.8, 1, 1.2, 1.5, 1.8)) {
  tumor_healthy_int <- FindClusters(
    tumor_healthy_int,
    graph.name = "mnn_snn",
    resolution = i,
    cluster.name = paste0("mnn_clusters",i)
  )
  
  tumor_healthy_int <- FindClusters(
    tumor_healthy_int,
    graph.name = "harmony_snn",
    resolution = i,
    cluster.name = paste0("harmony_clusters",i)
  )
} #it's zero-indexed because it's forced into a factor for some reason

gc()

metadata <- tumor_healthy_int@meta.data 

tumor_healthy_int@meta.data <- tumor_healthy_int@meta.data[
  , colSums(
    is.na(
      tumor_healthy_int@meta.data
      )
    ) == 0, 
  drop = FALSE
  ]

library(clustree)


clustree_mnn_plot <- clustree(tumor_healthy_int, prefix = "mnn_clusters")

gc()

clustree_harmony_plot <- clustree(tumor_healthy_int, prefix = "harmony_clusters")

tumor_healthy_int@meta.data <- metadata

qs_save(tumor_healthy_int, 
        "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/tumor_healthy_int.qs2",
        nthreads = 6)

source("~/Projects/etmr/defaria/scripts/functions.R")

plot_hist(tumor_healthy_int@meta.data, "harmony_clusters1.5", "orig.ident", "orig.ident")
plot_hist(tumor_healthy_int@meta.data, "harmony_clusters1.5", "cell_type", "cell_type")
plot_hist(tumor_healthy_int@meta.data, "harmony_clusters1.5", "major_clust", "major_clust")

plot_hist(tumor_healthy_int@meta.data, "harmony_clusters1.5", "stage_id", "stage_id")

Idents(tumor_healthy_int) <- "mnn_clusters"

by_id_mnn <- DimPlot(tumor_healthy_int, 
        reduction = "mnn_umap", 
        group.by = "orig.ident", 
        raster = FALSE)

LabelClusters(by_id_mnn, id = 'mnn_clusters1.5')

DimPlot(tumor_healthy_int, 
        reduction = "mnn_umap", 
        group.by = "mnn_clusters1.2", 
        raster = FALSE,
        label = TRUE) +
  NoLegend()

DimPlot(tumor_healthy_int, 
        reduction = "harmony_umap", 
        group.by = "harmony_clusters1.2", 
        raster = FALSE,
        label = TRUE) +
  NoLegend()

DimPlot(tumor_healthy_int, 
        reduction = "harmony_umap", 
        group.by = "major_clust", 
        raster = FALSE,
        label = TRUE)  +
  NoLegend()

DimPlot(tumor_healthy_int, 
        reduction = "mnn_umap", 
        group.by = "major_clust", 
        raster = FALSE,
        label = TRUE)  +
  NoLegend()

DimPlot(tumor_healthy_int, 
        reduction = "mnn_umap", 
        group.by = "mnn_clusters1.2", 
        raster = FALSE,
        label = TRUE)  +
  NoLegend()

gc()

DimPlot(tumor_healthy_int, 
        reduction = "umap", 
        group.by = "stage_id", 
        raster = FALSE,
        label = TRUE)

DimPlot(tumor_healthy_int, 
        reduction = "harmony_umap", 
        group.by = "cell_type", 
        raster = FALSE,
        label = TRUE)

DimPlot(tumor_healthy_int, 
        reduction = "mnn_umap", 
        group.by = "cell_type", 
        raster = FALSE,
        label = TRUE)


DimPlot(tumor_healthy_int, 
        reduction = "harmony_umap", 
        group.by = "sub_clust", 
        raster = FALSE,
        label = TRUE) +
  NoLegend()



tumor_healthy_int_ds <- IntegrateLayers(tumor_healthy_int, 
                                     method = HarmonyIntegration,
                                     assay = "RNA")

