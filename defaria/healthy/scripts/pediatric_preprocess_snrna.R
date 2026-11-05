library(Seurat)
library(qs2)
library(Matrix)
library(tidyverse)

data_dir <- "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408"

# Get all .h5 files
h5_files <- list.files(
  path = data_dir,
  pattern = "\\.h5$",
  full.names = TRUE
)

healthy_list <- list()

for (file in h5_files) {
  id <- sub("_.*", "", basename(file))
  
  counts <- Read10X_h5(
    file, 
    use.names = TRUE, 
    unique.features = TRUE)

  healthy_list[[id]] <- CreateSeuratObject(
    counts = counts,
    project = id,
    min.cells = 5,      # Drop genes expressed in fewer than 5 nuclei
    min.features = 300  # Drop nuclei with fewer than 300 detected genes
)
}

healthy_list <- lapply(healthy_list, function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x
})

human_mito <- lapply(healthy_list, function(x) {
  FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "percent.mt")
}
)

human_mito <- patchwork::wrap_plots(human_mito)


healthy_list <- lapply(healthy_list, function(x) {
  Idents(x) <- "orig.ident"
  x
})

human_feat <- patchwork::wrap_plots(lapply(healthy_list, function(x) {
  FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
}
))

healthy_list <- lapply(healthy_list, function(x) {
  subset(
    x,
    subset = nFeature_RNA > 250 &
      nFeature_RNA < 9000 &
      nCount_RNA < 75000 &
      percent.mt < 15
  )
})

# 2. Merge the list into a single Seurat object (automatically creates layers)
healthy_merged <- merge(
  x = healthy_list[[1]], 
  y = healthy_list[-1], 
  add.cell.ids = names(healthy_list)
)

# 3. Run standard processing on the merged object (Seurat v5 handles layers natively)
healthy_merged <- NormalizeData(healthy_merged)
healthy_merged <- FindVariableFeatures(healthy_merged)
healthy_merged <- ScaleData(healthy_merged)
healthy_merged <- RunPCA(healthy_merged)

healthy_metadata <- read.csv(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/Processed_data_RNA-all_BCs-meta-data.csv"
  )

sel_cols <- c("X", "batch",                     
  "RL.", "age",                      
  "chem", "concat_id",                
  "numerical_age", "stage_id",                  
  "Sex", "Race",                      
  "PMI", "Brain.Regions.",            
  "Cause.of.Death", "ICD.10.Code",               
  "ICD.10.category", "Oxygen.No.Oxygen",          
  "Date.of.Collection", "Collection_year",           
  "Library.Prep.Date", "Library.Prep.Lot",          
  "n_genes_by_counts", "log1p_n_genes_by_counts",   
  "total_counts", "log1p_total_counts",        
  "doublet_score", "log10_gene_counts",        
  "log10_UMI_counts", "percent_mito",              
  "percent_ribo", "n_counts",                  
  "leiden", "mat.dev",                   
  "cell_type", "major_clust",               
  "sub_clust", "combined.leiden")

healthy_metadata <- healthy_metadata |>
  select(all_of(sel_cols))

healthy_metadata$stripped_name <- sub("^([^-]+)-([^_]+)_.*$", "\\2_\\1-1", healthy_metadata$X)



valid_cells_mask <- Cells(healthy_merged) %in% healthy_metadata$stripped_name

# 2. Subset the Seurat object to keep ONLY those cells
healthy_merged <- subset(healthy_merged, cells = Cells(healthy_merged)[valid_cells_mask])
matched_indices <- match(Cells(healthy_merged), healthy_metadata$stripped_name)


aligned_metadata <- healthy_metadata[matched_indices, ]


rownames(aligned_metadata) <- Cells(healthy_merged)

healthy_merged <- AddMetaData(healthy_merged, metadata = aligned_metadata)

gc()


ElbowPlot(healthy_merged, ndims = 50)

healthy_merged <- FindNeighbors(healthy_merged, dims = 1:40)
healthy_merged <- RunUMAP(healthy_merged, dims = 1:40)

DimPlot(healthy_merged, group.by = "cell_type", raster = FALSE)
DimPlot(healthy_merged, group.by = "major_clust", raster = FALSE)
FeaturePlot(healthy_merged, features = "percent.mt", raster = FALSE)
FeaturePlot(healthy_merged, features = "percent_mito", raster = FALSE) #?? filter didnt work
 gc()

healthy_merged <- FindClusters(healthy_merged, 
                               resolution = .8, 
                               algorithm = 4, 
                               random.seed = 2026)

#paper integrated by reducing data to 1000 UMI for the dimensionality reduction

gc()

healthy_int <- IntegrateLayers(healthy_merged, 
                               method = SeuratWrappers::FastMNNIntegration)

gc()

healthy_int <- FindNeighbors(healthy_int, 
                             reduction = "integrated.mnn", 
                             dims = 1:40)

healthy_int <- RunUMAP(healthy_int, 
                       reduction = "integrated.mnn", 
                       dims = 1:40, 
                       seed.use = 2026, 
                       reduction.name = "umap.integrated")

DimPlot(healthy_int, 
        group.by = "cell_type", 
        reduction = "umap.integrated")

DimPlot(healthy_int, 
        group.by = "major_clust", 
        reduction = "umap.integrated", 
        label = TRUE)

DimPlot(healthy_int, 
        split.by = "cell_type", 
        group.by = "major_clust", 
        label = TRUE, 
        reduction = "umap.integrated") + 
  NoLegend()

gc()

healthy_int <- FindClusters(healthy_int, 
                            resolution = 1, 
                            algorithm = 4, 
                            random.seed = 2026)

qs_save(healthy_int, "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/healthy_int.qs2")

subset(healthy_merged, subset = nCount_RNA) 

rgc_unint <- FeaturePlot(healthy_int, features = c(
  "PAX6","SOX2","VIM","NES",  
  "HES1","HES5","FABP7","SLC1A3",
  "PTN","HMGA2"
), order = TRUE)

rgc_unint <- FeaturePlot(healthy_int, features = c(
  "PAX6","SOX2","VIM","NES",
  "HES1","HES5","FABP7","SLC1A3",
  "PTN","HMGA2"
), order = TRUE, reduction = "umap.integrated")

FeaturePlot(healthy_int, 
            features = c("SOX2","VIM","NES", "PAX6",
                                      "ALDH1L1", "AQP4", "SLC1A3", "SLC1A2", "SOX9"), 
            order = TRUE, reduction = "umap.integrated")

FeaturePlot(subset(healthy_int, subset = major_clust == "Astro"), 
            features = c("SOX2","VIM","NES", "PAX6",
                         "ALDH1L1", "AQP4", "SLC1A3", "SLC1A2", "SOX9"), 
            order = TRUE, reduction = "umap.integrated")



DimPlot(subset(healthy_int, subset = major_clust == "Astro"), 
        group.by = "sub_clust", 
        reduction = "umap.integrated",
        split.by = "sub_clust") + NoLegend()

gc()

DimPlot(subset(healthy_int, subset = major_clust == "Astro"), 
        group.by = "sub_clust", 
        reduction = "umap.integrated", label = TRUE) + NoLegend()

DimPlot(subset(healthy_int, subset = major_clust == "Astro"), 
        group.by = "sub_clust", 
        reduction = "umap", label = TRUE) + NoLegend()
