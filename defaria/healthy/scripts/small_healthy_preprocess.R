library(qs2)
library(BPCells)
library(Azimuth)
library(tidyverse)
library(clustree)
library(Seurat)

###############################important################################

healthy_mat <- open_matrix_anndata_hdf5(
  "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled_bayesian.h5ad")

write_matrix_dir(
  mat = healthy_mat,
  dir = "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled_bayesian_counts"
  )

healthy_mat <- open_matrix_dir(
  "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled_bayesian_counts"
  )

path = "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled_bayesian.h5ad"

metadata <- LoadH5ADobs(
  path
  )

healthy <- CreateSeuratObject(
  counts = healthy_mat, meta.data = metadata
)

#~/Projects/etmr/defaria/healthy/scripts/get_genes.sh
symbols <- read_tsv(
  "/home/sofia/Projects/etmr/defaria/healthy/data/raw/all_genes_symbol.tsv",
  show_col_types = FALSE
  )

# Match Seurat row names to Ensembl IDs in your dataframe
idx <- match(rownames(healthy), symbols$Ensembl_ID)

# Replace row names with gene names where a match exists
new_names <- symbols$Gene[idx]

# Keep original Ensembl IDs for genes without a match
new_names[is.na(new_names)] <- rownames(healthy)[is.na(new_names)]

rownames(healthy) <- new_names


healthy[["RNA"]] <- split(healthy[["RNA"]], f = healthy_int$sample_id)

healthy <- SCTransform(
  healthy, 
  variable.features.n = 3000, 
  ncells = 5000, 
  conserve.memory = TRUE
  )


healthy <- RunPCA(healthy)

healthy <- RunUMAP(
  healthy, 
  reduction = "pca", 
  dims = 1:30, 
  return.model = TRUE, 
  reduction.name = "unintegrated"
  )

healthy <- FindNeighbors(
  healthy, 
  reduction = "pca", 
  dims = 1:30,
  graph.name = c(
    "unintegrated_nn", "unintegrated_snn"
    ))

for (r in c(.6, .8, 1, 1.2, 1.5, 1.8, 2, 2.1)) {
  healthy <- FindClusters(
    healthy, 
    resolution = r,  
    algorithm = 4, 
    graph.name = "unintegrated_snn", 
    random.seed=2026
    )
} #i'm doing this to have several resolution values for clustree



qs_save(
  healthy, 
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/human_dev_downsampled_bayesian.qs2", 
  nthreads = 6
  )

used_cells <- Cells(healthy) #cells to exclude from other sampling


write.table(
  used_cells, 
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/used_cells_downsampled_bayesian.txt"
)

#############################integrate#########################################################

healthy_int <- qs_read(
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/human_dev_downsampled_bayesian.qs2", 
  nthreads = 6
)


DefaultAssay(healthy_int) <- "RNA"

healthy_int <- DietSeurat(healthy_int,
                          assays = "RNA",
                          layers = "counts",
                          dimreducs = NULL
                          )

#Chemistry seems to be the biggest difference, but the original h5ad had sample_id annotated as the batch
healthy_int[["RNA"]] <- split(healthy_int[["RNA"]], f = healthy_int$donor_id)

# healthy_int <- SCTransform(
#   healthy_int, 
#   variable.features.n = 2000, 
#   ncells = 5000, 
#   conserve.memory = FALSE,
#   return.only.var.genes = FALSE
# )
# 
# gc()
# 
# 
# 
# #try just ignoring warning Warning message:
# # In PrepDR(object = object, features = features, verbose = verbose) :
# #   The following 887 features requested have not been scaled...
# 
# healthy_int <- RunPCA(healthy_int, assay = "SCT", reduction.name = "sct.pca")
# 
# 
# healthy_int <- RunUMAP(
#   healthy_int, 
#   reduction = "sct.pca", 
#   dims = 1:30, 
#   return.model = TRUE, 
#   reduction.name = "unintegrated.sct"
# )
# 
# healthy_int <- IntegrateLayers(
#   healthy_int,
#   method = "HarmonyIntegration",
#   orig.reduction = "sct.pca",
#   assay = "SCT"
# )
# 
# healthy_int <- RunUMAP(
#   healthy_int, 
#   reduction = "harmony", 
#   dims = 1:30, 
#   return.model = TRUE, 
#   reduction.name = "integrated.sct.harmony"
# )
# 
# healthy_int <- IntegrateLayers(
#   healthy_int,
#   method = SeuratWrappers::FastMNNIntegration,
#   orig.reduction = "sct.pca",
#   assay = "SCT",
#   batch = healthy_int$donor_id
# )
# 
# healthy_int <- RunUMAP(
#   healthy_int, 
#   reduction = "integrated.mnn", 
#   dims = 1:30, 
#   return.model = TRUE, 
#   reduction.name = "integrated.sct.mnn"
# )
# 
# mnn_plot <- DimPlot(
#   healthy_int, 
#   reduction = "integrated.sct.mnn", 
#   group.by = c("donor_id", "Chemistry", "CellClass")
#   ) & 
#   labs(subtitle = "SCT MNN-integrated")
# 
# 
# 
# 
# 
# unint_sct_plot <- DimPlot(
#   healthy_int, 
#   reduction = "unintegrated.sct", 
#   group.by = c("Chemistry", "CellClass")
#   ) + 
#   ggtitle("SCT unintegrated")
# 
# 
# gc()

#I can't make MNN work with SCT
healthy_int <- NormalizeData(healthy_int, assay = "RNA")
healthy_int <- FindVariableFeatures(healthy_int, assay = "RNA")
healthy_int <- ScaleData(healthy_int, assay = "RNA")


healthy_int <- RunPCA(healthy_int, assay = "RNA", reduction.name = "pca")


healthy_int <- RunUMAP(
  healthy_int, 
  reduction = "pca", 
  dims = 1:30, 
  return.model = TRUE, 
  reduction.name = "unintegrated"
)

healthy_int <- IntegrateLayers(
  healthy_int,
  method = HarmonyIntegration,
  assay = "RNA",
  orig.reduction = "pca"
)


healthy_int <- RunUMAP(
    healthy_int,
    reduction = "harmony",
    dims = 1:30,
    return.model = TRUE,
    reduction.name = "integrated.harmony"
  )

healthy_int <- IntegrateLayers(
    healthy_int,
    method = SeuratWrappers::FastMNNIntegration,
    orig.reduction = "pca", 
    assay = "RNA",
    batch = healthy_int$donor_id
  )

healthy_int <- RunUMAP(
    healthy_int,
    reduction = "integrated.mnn",
    dims = 1:30,
    return.model = TRUE,
    reduction.name = "integrated.mnn.umap"
  )

unint_sct_plot <- DimPlot(
    healthy_int,
    reduction = "unintegrated",
    group.by = c("donor_id","Chemistry", "CellClass")
    ) +
    ggtitle("lognorm unintegrated")
  
 harmony_plot <- DimPlot(
    healthy_int,
    reduction = "integrated.harmony",
    group.by = c("donor_id", "Chemistry", "CellClass")
  ) &
    labs(subtitle = "harmony-integrated")

mnn_plot <- DimPlot(
    healthy_int,
    reduction = "integrated.mnn.umap",
    group.by = c("donor_id", "Chemistry", "CellClass")
    ) &
    labs(subtitle = "MNN-integrated")

gc()


unint_int_plot <- (
  DimPlot(
    healthy_int,
    group.by = "donor_id",
    reduction = "unintegrated") +
    NoLegend() |
  DimPlot(
    healthy_int, 
    group.by = "donor_id", 
    reduction = "integrated.harmony") +
    NoLegend() |
  DimPlot(
    healthy_int, 
    group.by = "donor_id", 
    reduction = "integrated.mnn.umap")
  ) / 
  (
    DimPlot(
      healthy_int, 
      group.by = "Chemistry", 
      reduction = "unintegrated") +
   NoLegend() |
   DimPlot(
     healthy_int, 
     group.by = "Chemistry", 
     reduction = "integrated.harmony") +
   NoLegend() |
  DimPlot(
    healthy_int, 
    group.by = "Chemistry", 
    reduction = "integrated.mnn.umap")
  ) / 
  (
    DimPlot(
      healthy_int, 
      group.by = "CellClass", 
      reduction = "unintegrated") +
   NoLegend() |
   DimPlot(
     healthy_int, 
     group.by = "CellClass", 
     reduction = "integrated.harmony") +
   NoLegend() |
  DimPlot(
    healthy_int, 
    group.by = "CellClass", 
    reduction = "integrated.mnn.umap")
  )

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognrom/unint_int_plot.png", 
  unint_int_plot, 
  width = 18, 
  height = 18, 
  units = "in")

healthy_int[["RNA"]] <- JoinLayers(healthy_int[["RNA"]])

healthy_int <- FindNeighbors(
  healthy_int, 
  reduction = "pca", 
  dims = 1:30,
  graph.name = c(
    "unintegrated_nn", "unintegrated_snn"
  ))

for (r in c(.6, .8, 1, 1.2, 1.5, 1.8, 2, 2.1)) {
  healthy_int <- FindClusters(
    healthy_int, 
    resolution = r,  
    algorithm = 4, 
    graph.name = "unintegrated_snn", 
    random.seed=2026
  )
} 

gc()


healthy_int <- FindNeighbors(
  healthy_int, 
  reduction = "harmony", 
  dims = 1:30,
  graph.name = c(
    "integrated_harmony_nn", "integrated_harmony_snn"
  ))

for (r in c(.6, .8, 1, 1.2, 1.5, 1.8, 2, 2.1)) {
  healthy_int <- FindClusters(
    healthy_int, 
    resolution = r,  
    algorithm = 4, 
    graph.name = "integrated_harmony_snn", 
    random.seed=2026
  )
} 

healthy_int <- FindNeighbors(
  healthy_int, 
  reduction = "integrated.mnn", 
  dims = 1:30,
  graph.name = c(
    "integrated_mnn_nn", "integrated_mnn_snn"
  ))

for (r in c(.6, .8, 1, 1.2, 1.5, 1.8, 2, 2.1)) {
  healthy_int <- FindClusters(
    healthy_int, 
    resolution = r,  
    algorithm = 4, 
    graph.name = "integrated_mnn_snn", 
    random.seed=2026
  )
} 

qs_save(
  healthy_int, 
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/human_dev_downsampled_bayesian_integrated.qs2", 
  nthreads = 6
)


