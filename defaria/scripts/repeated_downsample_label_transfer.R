library(BPCells)
library(Seurat)
library(dplyr)
library(tibble)
library(qs2)
library (ggplot2)

merged_reference  <- qs_read(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/merged_reference.qs2", 
  nthreads = 7)


#This si the final validation set, which i will not touch from here on out
holdout_id <- c("BRC2021", "RL2108")
holdout <- subset(merged_reference, subset = id %in% holdout_id)
holdout$id <- droplevels(holdout$id)

qs_save(holdout,
  "/media/sofia/TOSHIBA EXT/developing_human_brain/merged_reference_holdout.qs2", 
  nthreads = 7)

rm(holdout)
gc()

merged_reference <- subset(merged_reference, subset = !(id %in% holdout_id))
merged_reference$id <- droplevels(merged_reference$id)

gc()

####################################5-fold cross-validation#######################

set.seed(2026)  

n_folds <- 5

all_ids <- unique(merged_reference$id)

# Randomly shuffle IDs
shuffled_ids <- sample(all_ids)

# Divide IDs into 5 folds (we want to use all the samples as testing exactly once)
folds <- split(
  shuffled_ids,
  rep(1:n_folds, length.out = length(shuffled_ids))
)

cv_results <- vector("list", n_folds)

meta_df <- merged_reference@meta.data[, c("id", "dev_stage"), drop = FALSE]
all_cells <- rownames(meta_df)

for (i in seq_along(1:n_folds)) {
  cat("\n--- CV Fold", i, "of", n_folds, "---\n")
  
  test_samples <- droplevels(folds[[i]])
  ref_samples  <- droplevels(setdiff(all_ids, test_samples))
  
  query_cells <- all_cells[meta_df$id %in% test_samples]
  ref_cells   <- all_cells[meta_df$id %in% ref_samples]
  
  # 1. Subset reference and sketch
  fold_ref <- subset(merged_reference, cells = ref_cells)
  DefaultAssay(fold_ref) <- "RNA"
  
  
  fold_ref[["RNA"]] <- split(fold_ref[["RNA"]], f = fold_ref$id)
  fold_ref <- NormalizeData(fold_ref, verbose = FALSE)
  fold_ref <- FindVariableFeatures(fold_ref, nfeatures = 2000, verbose = FALSE) 
  
  
  fold_ref <- SketchData(
    object = fold_ref,
    ncells = 1000,
    method = "LeverageScore",
    sketched.assay = "sketch",
    features = VariableFeatures(fold_ref)
  )
  
  DefaultAssay(fold_ref) <- "sketch"
  
  # 2. Process sketched assay (keep nfeatures reasonable, e.g. 2000, to limit RAM)
  fold_ref <- FindVariableFeatures(fold_ref, nfeatures = 2000, verbose = FALSE)
  fold_ref <- ScaleData(fold_ref, verbose = FALSE)
  fold_ref <- RunPCA(fold_ref, npcs = 30, verbose = FALSE)
  
  # 3. Integrate batch layers inside the training fold
  fold_ref <- IntegrateLayers(
    object = fold_ref,
    method = SeuratWrappers::FastMNNIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.mnn",
    verbose = FALSE
  )
  
  fold_ref <- RunUMAP(
    fold_ref,
    reduction = "integrated.mnn",
    dims = 1:30,
    return.model = TRUE,
    verbose = FALSE
  )
  
  # 4. Subset query and normalize for transfer
  fold_query <- subset(merged_reference, cells = query_cells)
  DefaultAssay(fold_query) <- "RNA"
  fold_query <- NormalizeData(fold_query, verbose = FALSE)
  
  qs_save(fold_ref, 
          paste0("/media/sofia/TOSHIBA EXT/developing_human_brain/downsample_ref_", 
                 i, 
                 ".qs2"), nthreads = 7)
  qs_save(fold_query, 
          paste0("/media/sofia/TOSHIBA EXT/developing_human_brain/downsample_query_", 
                 i, 
                 ".qs2"), nthreads = 7)
  
  
  anchors <- FindTransferAnchors(
    reference = fold_ref,
    query = fold_query,
    reference.assay = "sketch",
    query.assay = "RNA",
    normalization.method = "LogNormalize",
    reduction = "pcaproject",
    dims = 1:30,
    verbose = FALSE
  )

  sketch_cells <- colnames(fold_ref[["sketch"]])
  ref_dev_stage <- fold_ref$dev_stage[sketch_cells]
  ref_harm_celltype <- fold_ref$harmonized_cell_type[sketch_cells]
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = list(DevStage = ref_dev_stage,
                   Celltype = ref_harm_celltype),
    dims = 1:30,
    verbose = FALSE
  )
  
  # 5. Score and store lightweight results
  pred_ids <- predictions$DevStage$predicted.id
  true_ids <- fold_query$dev_stage
  valid_cells <- !is.na(pred_ids)
  iter_acc <- sum(true_ids[valid_cells] == pred_ids[valid_cells]) / sum(valid_cells)
  
  cat("Iteration", i, "Accuracy developmental stage:", round(iter_acc * 100, 2), "%\n") 
 
  pred_cts <- predictions$Celltype$predicted.id
  true_cts <- fold_query$harmonized_cell_type
  valid_cells_ct <- !is.na(pred_cts)
  iter_acc_ct <- sum(
    true_cts[valid_cells_ct] == pred_cts[valid_cells_ct]
    ) / sum(
      valid_cells_ct
      )
  
  cat("Iteration", i, "Accuracy cell type:", round(iter_acc_ct * 100, 2), "%\n") 
  
  cv_results[[i]] <- list(
    fold = i,
    accuracy_devstage = iter_acc,
    test_samples = test_samples,
    cell = query_cells,
    predictions_devstage = data.frame(
      true_stage = true_ids,
      predicted_stage = pred_ids,
      max_score = predictions$DevStage$prediction.score.max,
      stringsAsFactors = FALSE),
    predictions_ct = data.frame(
      true_stage = true_cts,#change this name eventually
      predicted_stage = pred_cts,
      max_score = predictions$Celltype$prediction.score.max,
      stringsAsFactors = FALSE
    )
  )
  
  rm(fold_ref, fold_query, anchors, predictions)
  gc(verbose = FALSE)
}

# Save CV results separately to disk so you don't lose progress if a later fold hits a snag
qs_save(
  cv_results, 
  "/media/sofia/TOSHIBA EXT/developing_human_brain/cv_results.qs2", 
  nthreads = 7
  )
  
 #Confucion matrices
 cm_list <- vector("list", n_folds)

 for (i in seq_along(1:n_folds)){ 
   for (j in seq_along(1:2)) {
   df <- cv_results[[i]][[4+j]]
   all_classes <- union(
     as.character(df$true_stage),
     as.character(df$predicted_stage)
   )
   
   cm <- table(
     Predicted = factor(df$predicted_stage, levels = all_classes),
     Actual    = factor(df$true_stage, levels = all_classes)
   )
   
  
   prop_mat <- prop.table(cm, margin = 2)
   prop_mat[is.nan(prop_mat)] <- 0
   
   cm_norm_df <- merge(
     as.data.frame(prop_mat),
     as.data.frame(cm),
     by = c("Predicted", "Actual"),
     suffixes = c(".prop", ".count")
   )
   cm_norm_df$Predicted <- factor(cm_norm_df$Predicted, levels = rownames(cm))
   cm_norm_df$Actual    <- factor(cm_norm_df$Actual,    levels = colnames(cm))
   #separate loop later ig
   cm_list[[i]][[j]] <- ggplot(cm_norm_df, aes(x = Predicted, y = Actual, fill = Freq.prop)) +
     geom_tile(color = "white", linewidth = 0.5) +
     geom_text(
       aes(label = sprintf("%.2f\n(n=%d)", Freq.prop, Freq.count)),
       size = 3
     ) +
     scale_fill_gradient(low = "#ffeef2", high = "#d81159", limits = c(0, 1)) +
     labs(
       title = paste0(
         "Fold ", 
         i, 
         " Confusion Matrix " , 
         names(cv_results[[i]])[[4+j]],
         " (Normalized by Actual Class)"),
       x = "Predicted",
       y = "Actual",
       fill = "Recall"
     ) +
     theme_minimal(base_size = 14) +
     theme(
       panel.grid = element_blank(),
       axis.text.x = element_text(angle = 45, hjust = 1)
     ) + 
     scale_y_discrete(limits = rev)
   
  ggsave(
    paste0(
    "/media/sofia/TOSHIBA EXT/developing_human_brain/plots/confusion_matrix_", 
    i, "_", 
    names(cv_results[[i]])[[4+j]], 
    ".png"),
    width = 16,
  height = 10
  )
   } 
   }
 
 #clustering
 for (i in seq_along(1:n_folds)){ 
   for (j in seq_along(1:2)) {
     df <- cv_results[[i]][[4+j]]
     all_classes <- union(
       as.character(df$true_stage),
       as.character(df$predicted_stage)
     )
     
     cm <- table(
       Predicted = factor(df$predicted_stage, levels = all_classes),
       Actual    = factor(df$true_stage, levels = all_classes)
     )
     
     
     prop_mat <- prop.table(cm, margin = 2)
     prop_mat[is.nan(prop_mat)] <- 0
     
     cm_norm_df <- merge(
       as.data.frame(prop_mat),
       as.data.frame(cm),
       by = c("Predicted", "Actual"),
       suffixes = c(".prop", ".count")
     )
     cm_norm_df$Predicted <- factor(cm_norm_df$Predicted, levels = rownames(cm))
     cm_norm_df$Actual    <- factor(cm_norm_df$Actual,    levels = colnames(cm))
     #separate loop later ig
     cm_list[[i]][[j]] <- ggplot(cm_norm_df, aes(x = Predicted, y = Actual, fill = Freq.prop)) +
       geom_tile(color = "white", linewidth = 0.5) +
       geom_text(
         aes(label = sprintf("%.2f\n(n=%d)", Freq.prop, Freq.count)),
         size = 3
       ) +
       scale_fill_gradient(low = "#ffeef2", high = "#d81159", limits = c(0, 1)) +
       labs(
         title = paste0(
           "Fold ", 
           i, 
           " Confusion Matrix " , 
           names(cv_results[[i]])[[4+j]],
           " (Normalized by Actual Class)"),
         x = "Predicted",
         y = "Actual",
         fill = "Recall"
       ) +
       theme_minimal(base_size = 14) +
       theme(
         panel.grid = element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1)
       ) + 
       scale_y_discrete(limits = rev)
     
     ggsave(
       paste0(
         "/media/sofia/TOSHIBA EXT/developing_human_brain/plots/confusion_matrix_", 
         i, "_", 
         names(cv_results[[i]])[[4+j]], 
         ".png"),
       width = 16,
       height = 10
     )
   } 
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

