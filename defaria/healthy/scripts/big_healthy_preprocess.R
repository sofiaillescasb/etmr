#bash /home/sofia/Projects/etmr/defaria/healthy/scripts/downsample_bayesian.sh /home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev.h5ad /home/sofia/Projects/etmr/defaria/healthy/data/processed/human_dev_bayesian_ref.h5ad

library(qs2)
library(BPCells)
library(Azimuth)
library(tidyverse)
library(clustree)
library(Seurat)

#DO NOT DELETE
# healthy_mat <- open_matrix_anndata_hdf5(
#   "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_bayesian_ref.h5ad") #downsampled without query cells and 100k total
# 
# write_matrix_dir(
#   mat = healthy_mat,
#   dir = "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled_bayesian_counts_big"
# )

healthy_mat <- open_matrix_dir(
  "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_downsampled_bayesian_counts_big"
)

# #warning asks if matrix should be integers
# # Extract a subset into memory as a standard sparse matrix
# sub_mat <- healthy_mat[1:min(1000, nrow(healthy_mat)), 1:min(1000, ncol(healthy_mat))] |>
#   as("dgCMatrix")
# 
# # Get all non-zero values
# vals <- sub_mat@x
# 
# # Check if any non-zero value has a fractional/decimal part
# has_decimals <- any(vals %% 1 != 0)
# 
# print(has_decimals)

#i think yes it should be integers
healthy_mat <- convert_matrix_type(healthy_mat, type = "uint32")


path = "/home/sofia/Projects/etmr/defaria/healthy/data/raw/human_dev_bayesian_ref.h5ad"

metadata <- LoadH5ADobs(
  path
)

healthy <- CreateSeuratObject(
  counts = healthy_mat, meta.data = metadata
)

#~/Projects/etmr/defaria/healthy/scripts/get_genes.sh because the symbols were in another part of the anndata
symbols <- read_tsv(
  "/home/sofia/Projects/etmr/defaria/healthy/data/raw/all_genes_symbol.tsv",
  show_col_types = FALSE
)

# Match Seurat row names to Ensembl IDs in your dataframe
idx <- match(rownames(healthy), symbols$Ensembl_ID)

#if NAs bother me i can do this later healthy <- healthy[!grepl("-PAR-Y$", rownames(healthy)), ]

# rownames(healthy)[is.na(new_names)]
# [1] "ENSG00000228572.7-PAR-Y"  "ENSG00000182378.14-PAR-Y"
# [3] "ENSG00000178605.13-PAR-Y" "ENSG00000226179.6-PAR-Y" 
# [5] "ENSG00000167393.18-PAR-Y" "ENSG00000281849.3-PAR-Y" 
# [7] "ENSG00000275287.5-PAR-Y"  "ENSG00000280767.3-PAR-Y" 
# [9] "ENSG00000234958.6-PAR-Y"  "ENSG00000229232.6-PAR-Y" 
# [11] "ENSG00000185960.14-PAR-Y" "ENSG00000237531.6-PAR-Y" 
# [13] "ENSG00000225661.7-PAR-Y"  "ENSG00000205755.12-PAR-Y"
# [15] "ENSG00000198223.16-PAR-Y" "ENSG00000265658.6-PAR-Y" 
# [17] "ENSG00000223274.6-PAR-Y"  "ENSG00000185291.12-PAR-Y"
# [19] "ENSG00000169100.14-PAR-Y" "ENSG00000236871.7-PAR-Y" 
# [21] "ENSG00000236017.8-PAR-Y"  "ENSG00000169093.16-PAR-Y"
# [23] "ENSG00000182162.11-PAR-Y" "ENSG00000197976.12-PAR-Y"
# [25] "ENSG00000196433.13-PAR-Y" "ENSG00000223511.7-PAR-Y" 
# [27] "ENSG00000234622.6-PAR-Y"  "ENSG00000169084.15-PAR-Y"
# [29] "ENSG00000223571.6-PAR-Y"  "ENSG00000214717.12-PAR-Y"
# [31] "ENSG00000277120.5-PAR-Y"  "ENSG00000223773.7-PAR-Y" 
# [33] "ENSG00000230542.6-PAR-Y"  "ENSG00000002586.20-PAR-Y"
# [35] "ENSG00000237801.6-PAR-Y"  "ENSG00000237040.6-PAR-Y" 
# [37] "ENSG00000124333.16-PAR-Y" "ENSG00000228410.6-PAR-Y" 
# [39] "ENSG00000223484.7-PAR-Y"  "ENSG00000124334.17-PAR-Y"
# [41] "ENSG00000270726.6-PAR-Y"  "ENSG00000185203.12-PAR-Y"
# [43] "ENSG00000182484.15-PAR-Y" "ENSG00000227159.8-PAR-Y" 

# Replace row names with gene names where a match exists
new_names <- symbols$Gene[idx]

# Keep original Ensembl IDs for genes without a match
new_names[is.na(new_names)] <- rownames(healthy)[is.na(new_names)]

rownames(healthy) <- new_names

healthy[["RNA"]] <- split(healthy[["RNA"]], f = healthy$donor_id)


#I can't make it work with SCT
healthy <- NormalizeData(healthy)
healthy <- FindVariableFeatures(healthy)
healthy <- ScaleData(healthy, features = VariableFeatures(healthy))


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
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_transfer/integrated/human_dev_bayesian_ref.qs2", 
  nthreads = 6
)

used_cells <- Cells(healthy) #cells to exclude from other sampling


write.table(
  used_cells, 
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/used_cells_downsampled_bayesian_ref.txt"
)

#############################integrate#########################################################

healthy <- IntegrateLayers(
  healthy,
  method = HarmonyIntegration,
  assay = "RNA",
  orig.reduction = "pca"
)


healthy <- RunUMAP(
  healthy,
  reduction = "harmony",
  dims = 1:30,
  return.model = TRUE,
  reduction.name = "integrated.harmony"
)

healthy <- IntegrateLayers(
  healthy,
  method = SeuratWrappers::FastMNNIntegration,
  orig.reduction = "pca",
  assay = "RNA",
  batch = healthy$donor_id
)

healthy <- RunUMAP(
  healthy,
  reduction = "integrated.mnn",
  dims = 1:30,
  return.model = TRUE,
  reduction.name = "integrated.mnn.umap"
)


gc()


unint_int_plot <- (
  DimPlot(
    healthy,
    group.by = "donor_id",
    reduction = "unintegrated") +
    NoLegend() |
    DimPlot(
      healthy, 
      group.by = "donor_id", 
      reduction = "integrated.harmony") +
    NoLegend() |
    DimPlot(
      healthy, 
      group.by = "donor_id", 
      reduction = "integrated.mnn.umap")
) / 
  (
    DimPlot(
      healthy, 
      group.by = "Chemistry", 
      reduction = "unintegrated") +
      NoLegend() |
      DimPlot(
        healthy, 
        group.by = "Chemistry", 
        reduction = "integrated.harmony") +
      NoLegend() |
      DimPlot(
        healthy, 
        group.by = "Chemistry", 
        reduction = "integrated.mnn.umap")
  ) / 
  (
    DimPlot(
      healthy, 
      group.by = "CellClass", 
      reduction = "unintegrated") +
      NoLegend() |
      DimPlot(
        healthy, 
        group.by = "CellClass", 
        reduction = "integrated.harmony") +
      NoLegend() |
      DimPlot(
        healthy, 
        group.by = "CellClass", 
        reduction = "integrated.mnn.umap")
  )

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/unint_int_plot_big.png", 
  unint_int_plot, 
  width = 20, 
  height = 18, 
  units = "in")

healthy[["RNA"]] <- JoinLayers(healthy[["RNA"]])

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
} 

gc()


healthy <- FindNeighbors(
  healthy, 
  reduction = "harmony", 
  dims = 1:30,
  graph.name = c(
    "integrated_harmony_nn", "integrated_harmony_snn"
  ))

for (r in c(.6, .8, 1, 1.2, 1.5, 1.8, 2, 2.1)) {
  healthy <- FindClusters(
    healthy, 
    resolution = r,  
    algorithm = 4, 
    graph.name = "integrated_harmony_snn", 
    random.seed=2026
  )
} 

healthy <- FindNeighbors(
  healthy, 
  reduction = "integrated.mnn", 
  dims = 1:30,
  graph.name = c(
    "integrated_mnn_nn", "integrated_mnn_snn"
  ))

for (r in c(.6, .8, 1, 1.2, 1.5, 1.8, 2, 2.1)) {
  healthy <- FindClusters(
    healthy, 
    resolution = r,  
    algorithm = 4, 
    graph.name = "integrated_mnn_snn", 
    random.seed=2026
  )
} 


qs_save(
  healthy, 
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_transfer/integrated/human_dev_bayesian_ref.qs2", 
  nthreads = 6
)
