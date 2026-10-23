library(qs2)
library(Seurat)
library(tidyverse)


tumor <- qs_read("/home/sofia/Projects/etmr/defaria/snRNA/data/processed/int_human_tumor_sing.qs2", nthreads = 6)

tumor <- DietSeurat(tumor, 
                    layers = "counts",
                    dimreducs = NULL, 
                    assays = "RNA")
#to avoid potential batch effect i'm doing this by sample
tumor[["RNA"]] <- split(tumor[["RNA"]], f = tumor$orig.ident)


# tumor <- SCTransform(
#   tumor, 
#   variable.features.n = 2000, 
#   ncells = 5000, 
#   conserve.memory = TRUE
# )
# 
# gc()


#I can't make it work with SCT
tumor <- NormalizeData(tumor, assay = "RNA")
tumor <- FindVariableFeatures(tumor, assay = "RNA")
tumor <- ScaleData(tumor, features = VariableFeatures(tumor), assay = "RNA")

qs_save(
  tumor, 
  "/home/sofia/Projects/etmr/defaria/snRNA/data/processed/preprocessed/tumor_normalized.qs2", 
  nthreads = 6
)

###########################################################################################################
library(qs2)
library(Seurat)


tumor <- qs_read(
  "/home/sofia/Projects/etmr/defaria/snRNA/data/processed/preprocessed/tumor_normalized.qs2", 
  nthreads = 6
)

healthy_big <- qs_read(
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_tranfer/integrated/human_dev_bayesian_ref.qs2", 
    nthreads = 6
)

DefaultAssay(tumor) <- "RNA"

tumor <- JoinLayers(tumor)

gc()

anchors_int <-  FindTransferAnchors(
  reference = healthy_big, 
  query = tumor,
  normalization.method = "LogNormalize", 
  dims = 1:30, 
  reference.reduction = "pca", # Standard approach for projecting query onto reference features
  reference.assay = "RNA") #uses variable features from reference by default


qs_save(anchors_int,
        "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_tranfer/unintegrated/bayesian/anchors_b_MNN.qs2", 
        nthreads = 6)

gc()

tumor <- MapQuery(
  anchorset = anchors_int, reference = healthy_big, 
  query = tumor, refdata = list(
    ClusterMNN = "integrated_mnn_snn_res.1.8", #chosen from clustree
    CellTypeMNN = "CellClass",
    DevStageMNN = "development_stage_ontology_term_id"),
  reference.reduction = "pca", # Standard approach for projecting query onto reference features
  reduction.model = "integrated.mnn.umap")



summary(tumor$predicted.ClusterMNN.score)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1307  0.3989  0.5447  0.5815  0.7562  1.0000  

summary(tumor$predicted.CellTypeMNN.score)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2398  0.7779  0.9364  0.8652  0.9930  1.0000 

summary(tumor$predicted.DevStageMNN.score)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1530  0.4724  0.5804  0.5886  0.7061  1.0000


celltypes_or <- as.character(tumor$by_group_celltypes)

tumor <- AddMetaData(
  tumor, metadata = ifelse(
    celltypes_or %in% c(
      "Neurons_immature_GABA", "Neurons_immature_Glut", "Neurons_immature"
    ),
    "Neuroblast",
    ifelse(
      celltypes_or %in% c(
        "Astrocytes", "OPC"
      ), "Glioblast",
      ifelse(
        celltypes_or %in% c(
          "NEC-like_RGC", "RGC"
        ), "Radial glia",
        ifelse(
          celltypes_or %in% c(
            "Microglia"
          ), "Immune",
          ifelse(
            celltypes_or %in% c(
              "IPC"
            ), "Neuronal IPC",
            celltypes_or
          )
        )))), 
  col.name = "ComparableCellType")

# 1. Create a proportion table comparing old vs predicted labels
table_counts <- table(tumor$ComparableCellType, tumor$predicted.CellTypeMNN)
table_prop <- prop.table(table_counts, margin = 1) # Normalize by row (Cell type)

# 2. Convert table to dataframe for ggplot
df_heatmap <- as.data.frame(table_prop)
colnames(df_heatmap) <- c("Previous_CellClass", "predicted.CellTypeMNN", "Proportion")

# 3. Plot heat map
tumor_conf_mat <- ggplot(df_heatmap, aes(x = predicted.CellTypeMNN, y = Previous_CellClass, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkblue") +
  geom_text(aes(label = round(Proportion, 2)), size = 3, color = ifelse(df_heatmap$Proportion > 0.5, "white", "black")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Agreement Between Previous CellClass and Mapped Predictions",
    x = "Predicted CellTypeInt (Reference)",
    y = "Previous CellClass (Query)"
  )

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/tumor_conf_mat.png", 
  tumor_conf_mat, 
  units = "in", 
  width = 10, 
  height = 8)

qs_save(
  tumor,
  "/home/sofia/Projects/etmr/defaria/snRNA/data/processed/preprocessed/tumor_transfered.qs2", 
  nthreads = 6
)

################################################################################################################



healthy_small <- qs_read(
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_tranfer/integrated/human_dev_downsampled_bayesian_integrated.qs2",
  nthreads = 6
  )

DefaultAssay(healthy_small) <- "RNA"

gc()

anchors_healthy_int <-  FindTransferAnchors(
  reference = healthy_big, 
  query = healthy_small,
  normalization.method = "LogNormalize", 
  dims = 1:30, 
  reference.reduction = "pca") #uses variable features from reference by default


qs_save(anchors_healthy_int,
        "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_tranfer/unintegrated/bayesian/anchors_b_healthy_int.qs2", 
        nthreads = 6)

healthy_small <- MapQuery(
  anchorset = anchors_healthy_int, reference = healthy_big, 
  query = healthy_small, refdata = list(
    ClusterMNN = "integrated_mnn_snn_res.1.8", #chosen from clustree
    CellTypeMNN = "CellClass",
    DevStageMNN = "development_stage_ontology_term_id"),
  reference.reduction = "pca", 
  reduction.model = "integrated.mnn.umap")

###############################etmr score####################################

signt_etmr <- read.csv("/home/sofia/Projects/etmr/defaria/snRNA/tables/signature.csv")


signt_etmr <- as.list(signt_etmr)
signt_etmr$Positive <- signt_etmr$Positive[signt_etmr$Positive %in% Features(healthy_small)]
signt_etmr$Negative <- signt_etmr$Negative[signt_etmr$Negative %in% Features(healthy_small)]

#########filter healthy to only keep clusters that are similar to tumor###########

for (x in names(signt_etmr)) {
  healthy_small <- AddModuleScore(
    healthy_small,
    features = signt_etmr[x],
    name = x
  )
}

healthy_small@meta.data <- healthy_small@meta.data |> 
  mutate(etmr_score = Positive1 - Negative1) |>
  select(-c(Positive1, Negative1))

qs_save(healthy_small,
"~/Projects/etmr/defaria/healthy/data/processed/label_tranfer/integrated/healthy_small_transferred.qs2", 
 nthreads = 6)


for (x in names(signt_etmr)) {
  healthy_big <- AddModuleScore(
    healthy_big,
    features = signt_etmr[x],
    name = x
  )
}

healthy_big@meta.data <- healthy_big@meta.data |> 
  mutate(etmr_score = Positive1 - Negative1) |>
  select(-c(Positive1, Negative1))

qs_save(healthy_big,
        "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_tranfer/integrated/human_dev_bayesian_ref.qs2", 
        nthreads = 6)

############################################plots and diagnostics##############################################
summary(healthy_small$predicted.ClusterMNN.score)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1915  0.7801  0.9587  0.8715  1.0000  1.0000 


summary(healthy_small$predicted.CellTypeMNN.score)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3068  0.9004  0.9976  0.9224  1.0000  1.0000 

summary(healthy_small$predicted.DevStageMNN.score)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1696  0.5117  0.6745  0.6767  0.8505  1.0000 

summary(healthy_small$etmr_score)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.35134 -0.03342  0.02171  0.01704  0.07234  0.31193 

summary(healthy_big$etmr_score)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.39139 -0.03405  0.02076  0.01600  0.07099  0.29528 


# 1. Create a proportion table comparing old vs predicted labels
table_counts <- table(healthy_small$CellClass, healthy_small$predicted.CellTypeMNN)
table_prop <- prop.table(table_counts, margin = 1) # Normalize by row (CellClass)

# 2. Convert table to dataframe for ggplot
df_heatmap <- as.data.frame(table_prop)
colnames(df_heatmap) <- c("Previous_CellClass", "predicted.CellTypeMNN", "Proportion")

# 3. Plot heat map
small_healthy_conf_mat <- ggplot(df_heatmap, aes(x = predicted.CellTypeMNN, y = Previous_CellClass, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkblue") +
  geom_text(aes(label = round(Proportion, 2)), size = 3, color = ifelse(df_heatmap$Proportion > 0.5, "white", "black")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Agreement Between Previous CellClass and Mapped Predictions",
    x = "Predicted CellTypeInt (Reference)",
    y = "Previous CellClass (Query)"
  )

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/small_healthy_conf_mat.png", 
  small_healthy_conf_mat, 
  units = "in", 
  width = 10, 
  height = 8)


