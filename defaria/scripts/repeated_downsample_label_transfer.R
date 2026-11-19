library(Seurat)
library(dplyr)
library(tibble)
library(qs2)

merged_reference  <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/merged_reference.qs2", nthreads = 7)

n_iterations <- 5
cv_results <- list()

for (i in 1:n_iterations) {
  cat("\n--- CV Iteration", i, "of", n_iterations, "---\n")
  
  # 1. Split into training samples and test samples for this fold
  test_samples <- sample(unique(merged_reference$id), size = 3)
  ref_samples  <- setdiff(unique(merged_reference$id), test_samples)
  
  fold_query <- subset(merged_reference, subset = id %in% test_samples)
  fold_ref   <- subset(merged_reference, subset = id %in% ref_samples)
  
  # 2. Sketch the reference training set for speed
  fold_ref <- SketchData(
    object = fold_ref,
    ncells = 1000,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  DefaultAssay(fold_ref) <- "sketch"
  
  # 3. Standard processing on the sketch assay
  fold_ref <- NormalizeData(fold_ref, verbose = FALSE)
  fold_ref <- FindVariableFeatures(fold_ref, verbose = FALSE)
  fold_ref <- ScaleData(fold_ref, verbose = FALSE)
  fold_ref <- RunPCA(fold_ref, npcs = 30, verbose = FALSE)
  
  # 4. CRITICAL: Integrate the training layers to remove batch effects *inside the fold*
  fold_ref <- IntegrateLayers(
    object = fold_ref,
    method = SeuratWrappers::FastMNNIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.mnn",
    verbose = FALSE
  )
  
  # 5. Build UMAP model on the integrated reduction
  fold_ref <- RunUMAP(fold_ref, reduction = "integrated.harmony", dims = 1:30, return.model = TRUE, verbose = FALSE)
  
  # 6. Find anchors and transfer labels to the test query
  anchors <- FindTransferAnchors(
    reference = fold_ref,
    query = fold_query,
    normalization.method = "LogNormalize",
    reduction = "pcaprojection", # Map query into reference PCA space
    dims = 1:30,
    verbose = FALSE
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = fold_ref$cell_type,
    dims = 1:30,
    verbose = FALSE
  )
  
  # 7. Calculate accuracy
  fold_query <- AddMetaData(fold_query, predictions)
  valid_cells <- !is.na(fold_query$predicted.id)
  iter_acc <- sum(fold_query$cell_type[valid_cells] == fold_query$predicted.id[valid_cells]) / sum(valid_cells)
  
  cat("Iteration", i, "Accuracy:", round(iter_acc * 100, 2), "%\n")
  cv_results[[i]] <- list(accuracy = iter_acc, test_samples = test_samples)
}

#Remove original objects

#Read tumor

#Read Jessa

#Read scrna from de faria

#Normalize each

#I think label transferring each one separately makes more sense?

######loop#########

#Downsample

#save cells used

#Preprocess

#Integrate

#label transfer

