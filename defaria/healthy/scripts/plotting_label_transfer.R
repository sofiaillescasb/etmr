library(qs2)
library(Seurat)
library(tidyverse)


qs_save(healthy_small,
        "~/Projects/etmr/defaria/healthy/data/processed/label_tranfer/integrated/healthy_small_transferred.qs2", 
        nthreads = 6)

qs_save(healthy_big,
        "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_tranfer/integrated/human_dev_bayesian_ref.qs2", 
        nthreads = 6)

qs_save(
  tumor, 
  "/home/sofia/Projects/etmr/defaria/snRNA/data/processed/preprocessed/tumor_normalized.qs2", 
  nthreads = 6
)

celltype_colors <-  c(
  "Endothelial" = "#ffe270",
      "Erythrocyte" = "#ff9105",
     "Fibroblast" = "#c8d5b9",
     "Glioblast" = "#8332ac",
     "Immune" = "#caa8f5",
     "Lymphocytes" = "#a8dcd1",
     "Neural crest" = "#708d81",
     "Neuroblast"= "#38a3a5",
     "Neuron"= "#22577a",
     "Neuronal IPC"= "#57cc99",
     "Neurons"= "#22577a",
     "Oligo"= "#fc6dab",
     "Oligodendrocytes" = "#fc6dab",
     "OPC" = "#f72585",
     "Pericytes" = "#c8d5b9",
     "Placodes"= "#4f5d75",
     "Radial glia" = "#80ed99",
     "Vascular" = "#832161"
     )   


tumor@meta.data$predicted.ClusterMNN <- as.factor(as.numeric(tumor@meta.data$predicted.ClusterMNN))

library(tidyverse)

cluster_int_type <- tumor@meta.data |>
  ggplot(aes(x = predicted.ClusterMNN, fill = Tumor_Type)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Cluster predicted from MNN-integrated healthy data") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#fbd960", "#d5b3ff"))

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/all_cluster_mnn_type.png", 
  cluster_int_type, 
  units = "in", 
  width = 8, 
  height = 6)

cluster_int_score <- tumor@meta.data |>
  group_by(predicted.ClusterMNN) |>
  filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(x = predicted.ClusterMNN.score)) +
  geom_histogram(binwidth = .05, fill = "#83c5be", color = "#006d77") +
  theme_classic() +
  facet_wrap(vars(predicted.ClusterMNN), scales = "free_y")  

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/all_cluster_mnn_score.png", 
  cluster_int_score, 
  units = "in", 
  width = 14, 
  height = 10
)

celltype_int_type <- tumor@meta.data |>
  ggplot(aes(x = predicted.CellTypeMNN, fill = Tumor_Type)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Cell type predicted from MNN-integrated healthy data") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#fbd960", "#d5b3ff"))

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/all_celltype_mnn_type.png", 
  celltype_int_type, 
  units = "in", 
  width = 8, 
  height = 6)


celltype_int_score <- tumor@meta.data |>
  group_by(predicted.CellTypeMNN) |>
  filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(x = predicted.CellTypeMNN.score)) +
  geom_histogram(binwidth = .05, fill = "#ffc8c8", color = "#ff5858") +
  theme_classic() +
  facet_wrap(vars(predicted.CellTypeMNN), scales = "free_y")  

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/celltype_mnn_score.png", 
  celltype_int_score, 
  units = "in", 
  width = 14, 
  height = 10
)


# Side-by-side UMAP plot
p1 <- DimPlot(
  healthy_small, 
  reduction = "ref.umap", 
  group.by = "CellClass", 
  label = TRUE,
  shuffle = TRUE, 
  cols = celltype_colors) + 
  ggtitle("Previous: CellClass") +
  NoLegend()

p2 <- DimPlot(
  healthy_small, 
  reduction = "ref.umap", 
  group.by = "predicted.CellTypeMNN", 
  label = TRUE,
  shuffle = TRUE, 
  cols = celltype_colors) + 
  ggtitle("Mapped: predicted.CellTypeMNN")

prev_vs_pred_healthy <- p1 | p2

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/prev_vs_pred_healthy.png", 
  prev_vs_pred_healthy, 
  units = "in", 
  width = 14, 
  height = 8)



p3 <- DimPlot(
  tumor, 
  reduction = "ref.umap", 
  group.by = "ComparableCellType", 
  label = FALSE,
  shuffle = TRUE, 
  cols = celltype_colors) + 
  ggtitle("Previous: CellClass")

p4 <- DimPlot(
  tumor, 
  reduction = "ref.umap", 
  group.by = "predicted.CellTypeMNN", 
  label = FALSE,
  shuffle = TRUE, 
  cols = celltype_colors) + 
  ggtitle("Mapped: predicted.CellTypeMNN")

prev_vs_pred_tumor <- p3 | p4

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/prev_vs_pred_tumor.png", 
  prev_vs_pred_tumor, 
  units = "in", 
  width = 14, 
  height = 8)

healthy_small@meta.data$predicted.ClusterMNN <- as.factor(as.numeric(healthy_small@meta.data$predicted.ClusterMNN))


cluster_mnn_healthy <- healthy_small@meta.data |>
  ggplot(aes(x = predicted.ClusterMNN, fill = Chemistry)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Cluster predicted from integrated healthy data") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#bfd8bd", "pink"))

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/cluster_mnn_healthy.png", 
  cluster_mnn_healthy, 
  units = "in", 
  width = 10, 
  height = 6)

cluster_mnn_healthy_score <- healthy_small@meta.data |>
  group_by(predicted.ClusterMNN) |>
  filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(x = predicted.ClusterMNN.score)) +
  geom_histogram(binwidth = .05, fill = "#83c5be", color = "#006d77") +
  theme_classic() +
  facet_wrap(vars(predicted.ClusterMNN), scales = "free_y")  

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/cluster_mnn_healthy_score.png", 
  cluster_mnn_healthy_score, 
  units = "in", 
  width = 18, 
  height = 18
)


celltype_mnn_healthy_type <- healthy_small@meta.data |>
  ggplot(aes(x = predicted.CellTypeMNN, fill = Chemistry)) +
  geom_bar() +
  theme_classic() +
  ggtitle("Cell type predicted from integrated healthy data") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#bfd8bd", "pink"))

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/celltype_mnn_healthy_type.png", 
  celltype_mnn_healthy_type, 
  units = "in", 
  width = 8, 
  height = 6)

celltype_mnn_healthy_score <- healthy_small@meta.data |>
  group_by(predicted.CellTypeMNN) |>
  filter(n() > 10) |>
  ungroup() |>
  ggplot(aes(x = predicted.CellTypeMNN.score)) +
  geom_histogram(binwidth = .05, fill = "#ffc8c8", color = "#ff5858") +
  theme_classic() +
  facet_wrap(vars(predicted.CellTypeMNN), scales = "free_y")  +
  ggtitle("Cell type score healthy-healthy")

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/celltype_mnn_healthy_score.png",
  celltype_mnn_healthy_score, 
  units = "in", 
  width = 14, 
  height = 10)


###############################combined plots################################################################

healthy_big$id <- 'healthy_reference'
healthy_small$id <- 'healthy_query'
tumor$id <- 'tumor'

reference_umap <- as.data.frame(Embeddings(healthy_big, "integrated.mnn.umap")) |>
  mutate(
    ID = healthy_big$id,
    Cluster = healthy_big$integrated_harmony_snn_res.1.8, 
    CellType = healthy_big$CellClass,
    DevStage = healthy_big$development_stage_ontology_term_id
  )

reference_umap$Condition <- "Healthy"

tumor_umap <- as.data.frame(Embeddings(tumor, "ref.umap")) |>
  mutate(
    ID = tumor$id,
    Cluster = tumor$predicted.ClusterMNN, 
    CellType = tumor$predicted.CellTypeMNN,
    DevStage = tumor$predicted.DevStageMNN,
    OriginalCellType = tumor$ComparableCellType
  )

tumor_umap$Condition <- "Tumor"

tumor_umap <- tumor_umap |>
  dplyr::rename(integratedmnnumap_1 = refUMAP_1, integratedmnnumap_2 = refUMAP_2)


healthy_q_umap <- as.data.frame(Embeddings(healthy_small, "ref.umap")) |>
  mutate(
    ID = healthy_small$id,
    Cluster = healthy_small$predicted.ClusterMNN, 
    CellType = healthy_small$predicted.CellTypeMNN,
    DevStage = healthy_small$predicted.DevStageMNN
  )

healthy_q_umap$Condition <- "Healthy"

healthy_q_umap <- healthy_q_umap |>
  dplyr::rename(integratedmnnumap_1 = refUMAP_1, integratedmnnumap_2 = refUMAP_2)


combined_umap <- rbind(
  reference_umap,
  healthy_q_umap,
  tumor_umap |> select(-OriginalCellType)
) |>
  mutate(factor(ID, levels = c("healthy_reference", "healthy_query", "tumor"), ordered = TRUE))

combined_plot_cond <- ggplot(combined_umap, aes(x = integratedmnnumap_1, y = integratedmnnumap_2, color = ID)) +
  geom_point(size = .6, alpha = 1)+
  theme_classic() +
  labs(
    x = "umap_1",
    y = "umap_2",
    color = "Condition"
  ) + 
  scale_color_manual(values = c(
    "healthy_reference" = "#aab2ff", 
    "healthy_query" = "#84ffc9", 
    "tumor" = "#ff0054"),
    limits = c("healthy_reference", "healthy_query", "tumor")) + 
  ggtitle("Combined MNN-integrated healthy reference, healthy query, and tumor data") +
  guides(color = guide_legend(override.aes = list(size = 4)))

# combined_plot_cond <- LabelClusters(
#   plot = combined_plot_cond,
#   id = "Cluster"
# )

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/all_combined_plot_cond.png", 
  combined_plot_cond, 
  units = "in", 
  width = 10, 
  height = 8
)

combined_plot_ct <- ggplot(combined_umap, aes(x = integratedmnnumap_1, y = integratedmnnumap_2, color = ID)) +
  geom_jitter(size = 0.6, alpha = 0.8) +
  theme_classic(base_size = 20) +
  labs(
    x = "umap_1",
    y = "umap_2",
    color = "Cell type"
  ) + 
  ggtitle("Combined unintegrated healthy and tumor data") + 
  facet_wrap(vars(CellType)) +
  scale_color_manual(
    values = c(
      "healthy_reference" = "#aab2ff", 
      "healthy_query" = "#84ffc9", 
      "tumor" = "#ff0054"
    ),
    limits = c("healthy_reference", "healthy_query", "tumor")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/all_combined_plot_ct.png", 
  combined_plot_ct, 
  units = "in", 
  width = 18, 
  height = 14
)


my_colors_healthy <- c(
  "HsapDv_0000022" = "#8f2d56", 
  "HsapDv_0000023" = "#dd2d4a", 
  "HsapDv_0000025" = "#f26a8d", 
  "HsapDv_0000027" = "#f49cbb", 
  "HsapDv_0000029" = "#f7cad0", 
  "HsapDv_0000046" = "#e9c46a", 
  "HsapDv_0000047" = "#f4a261",
  "HsapDv_0000048" = "#ee8959" ,
  "HsapDv_0000049" = "#e76f51", 
  "HsapDv_0000050" = "#a6c48a", 
  "HsapDv_0000051" = "#678d58", 
  "HsapDv_0000052" = "#74d3ae"
)

combined_plot_dev <- ggplot(
  combined_umap, 
  aes(x = integratedmnnumap_1, 
      y = integratedmnnumap_2, 
      color = DevStage)
  ) +
  geom_point(size = .6, alpha = 1)+
  theme_classic() +
  labs(
    x = "umap_1",
    y = "umap_2",
    color = "DevStage"
  ) + 
  scale_color_manual(values = my_colors_healthy) + 
  facet_wrap(vars(ID)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ggtitle("Developmental stage MNN-integrated")



ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/all_combined_plot_dev.png", 
  combined_plot_dev, 
  units = "in", 
  width = 10, 
  height = 8
)


combined_plot_dev_split <- ggplot(
  combined_umap, 
  aes(x = integratedmnnumap_1, 
      y = integratedmnnumap_2, 
      color = ID)
) +
  geom_point(size = .6, alpha = 1)+
  theme_classic() +
  labs(
    x = "umap_1",
    y = "umap_2",
    color = "ID"
  ) + 
  scale_color_manual(values = c(
    "healthy_reference" = "#aab2ff", 
    "healthy_query" = "#84ffc9", 
    "tumor" = "#ff0054"),
    limits = c("healthy_reference", "healthy_query", "tumor")) + 
  facet_wrap(vars(DevStage)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ggtitle("Developmental stage MNN-integrated")

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/all_combined_plot_dev.png", 
  combined_plot_dev, 
  units = "in", 
  width = 18, 
  height = 14
)

#get mapping score
#THIS SHOULD BE THE INTEGRATED FOLDER
anchors_int <- qs_read(
        "/home/sofia/Projects/etmr/defaria/healthy/data/processed/label_tranfer/unintegrated/bayesian/anchors_b_MNN.qs2", 
        nthreads = 6)

tumor[["mapping.score"]] <- MappingScore(anchors_int, ndim = 30)


# 1. Define a helper function for percentile ranking
to_percentile <- function(x) {
  rank(x, na.last = "keep") / sum(!is.na(x))
}

# 2. Apply to both metrics
tumor$mapping_score_percentile <- to_percentile(tumor$mapping.score)
tumor$etmr_score_percentile <- to_percentile(tumor$etmr_score)

# Define a min-max scaling function
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Apply to metadata columns
tumor$mapping_score_scaled <- min_max_scale(tumor$mapping.score)
tumor$etmr_score_scaled  <- min_max_scale(tumor$etmr_score)


whole_cor_plot <- tumor@meta.data |> 
  ggplot(aes(x = etmr_score_scaled, y =  mapping_score_scaled) ) +
  geom_point(size = .6, color = "#bfd8bd") +
  geom_smooth(color = "#2d6a4f", method = "lm") +
  theme_classic() +
  ggtitle("Min-max scaled mapping and ETMR score")  +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/whole_cor_plot.png", 
  whole_cor_plot, 
  units = "in", 
  width = 8, 
  height = 8
)


tumor@meta.data |> 
  filter(predicted.DevStageMNN %in% c("HsapDv_0000022", "HsapDv_0000027") & predicted.CellTypeMNN == "Neuroblast") |>
  ggplot(aes(x = etmr_score_scaled, y =  mapping_score_scaled)) +
  geom_smooth() +
  # geom_point(size = .6) +
  facet_grid(cols = vars(predicted.CellTypeMNN), rows = vars(predicted.DevStageMNN))

tumor@meta.data |> 
  # filter(predicted.DevStageMNN %in% c("HsapDv_0000022", "HsapDv_0000027")) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled)) +
  geom_smooth(method = "lm") +
  # geom_point(size = .6) +
  facet_grid(cols = vars(predicted.CellTypeMNN), rows = vars(predicted.DevStageMNN), scales = "free")

cts <- unique(tumor$predicted.CellTypeMNN)

by_ct_plot <- list()

for (c in cts) {
  by_ct_plot[[c]] <- tumor@meta.data |> 
    filter(predicted.CellTypeMNN == c ) |>
    ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.DevStageMNN)) +
    geom_point(size = .6) +
    facet_wrap(vars(predicted.ClusterMNN))

}

# Get proportions
ref_prop <- prop.table(table(healthy_big$development_stage_ontology_term_id))
query_prop <- prop.table(table(tumor$predicted.DevStageMNN))

# Combine into a unified data frame (handling missing categories like HsapDv_0000048)
all_cats <- union(names(ref_prop), names(query_prop))
prop_df <- data.frame(
  Category = all_cats,
  Reference = sapply(all_cats, function(x) ifelse(x %in% names(ref_prop), ref_prop[x], 0)),
  Query = sapply(all_cats, function(x) ifelse(x %in% names(query_prop), query_prop[x], 0))
)

# Calculate Log2 Ratio (Query proportion vs Reference proportion)
# Values > 0 indicate enrichment in the tumor; values < 0 indicate depletion
prop_df$Log2Ratio <- log2((prop_df$Query + 1e-4) / (prop_df$Reference + 1e-4))
print(prop_df[order(prop_df$Log2Ratio, decreasing = TRUE), ])

# Assuming your transfer method outputs a score column (e.g., predicted.DevStageMNN.score)
# If low confidence scores cluster heavily in the largest reference categories, you have a size-bias artifact.
VlnPlot(tumor, features = "predicted.DevStageMNN.score", group.by = "predicted.DevStageMNN", pt.size = 0.1) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggpubr)

# 1. Calculate correlation explicitly if you want to inspect it first
cor_by_group <- tumor@meta.data %>%
  group_by(predicted.CellTypeMNN, predicted.DevStageMNN) %>%
  summarise(
    pearson_r = cor(etmr_score_scaled, mapping_score_scaled, use = "complete.obs"),
    n = n()
  )

print(cor_by_group)

# 2. Plot with regression line and the correlation statistic built-in

tumor@meta.data |>
  filter(predicted.ClusterMNN %in% c(seq(1,10))) |>
  group_by(predicted.CellTypeMNN) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.CellTypeMNN)) +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal() +
    facet_wrap(vars(predicted.ClusterMNN))

tumor@meta.data |>
  filter(predicted.ClusterMNN %in% c(seq(11,20))) |>
  group_by(predicted.CellTypeMNN) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.CellTypeMNN)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_minimal() +
  facet_wrap(vars(predicted.ClusterMNN))

tumor@meta.data |>
  filter(predicted.ClusterMNN %in% c(seq(21,30))) |>
  group_by(predicted.CellTypeMNN) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.CellTypeMNN)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_minimal() +
  facet_wrap(vars(predicted.ClusterMNN))

tumor@meta.data |>
  filter(predicted.ClusterMNN %in% c(seq(31,40))) |>
  group_by(predicted.CellTypeMNN) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.CellTypeMNN)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_minimal() +
  facet_wrap(vars(predicted.ClusterMNN))

tumor@meta.data |>
  filter(predicted.ClusterMNN %in% c(seq(41,50))) |>
  group_by(predicted.CellTypeMNN) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.CellTypeMNN)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_minimal() +
  facet_wrap(vars(predicted.ClusterMNN))

tumor@meta.data |>
  filter(predicted.ClusterMNN %in% c(seq(51,60))) |>
  group_by(predicted.CellTypeMNN) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.CellTypeMNN)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_minimal() +
  facet_wrap(vars(predicted.ClusterMNN))

cor_etmr_map <- tumor@meta.data |>
  group_by(predicted.CellTypeMNN) |>
  ggplot(aes(x = etmr_score_scaled, y = mapping_score_scaled, color = predicted.CellTypeMNN)) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", color = "black") +
  theme_classic() +
  facet_wrap(vars(predicted.CellTypeMNN), scales = "free") +
  scale_color_manual(values = celltype_colors)

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/cor_etmr_map.png", 
  cor_etmr_map, 
  units = "in", 
  width = 14, 
  height = 10
)

FeaturePlot(tumor, features = "etmr_score_scaled", order = TRUE)
FeaturePlot(tumor, features = "mapping_score_scaled", order = TRUE)

tumor$inverse_mapping_score_scaled <- 1 - tumor$mapping_score_scaled
tumor$combined_score <- tumor$etmr_score_scaled * tumor$inverse_mapping_score_scaled

FeaturePlot(
  tumor, 
  features = c(
    "predicted.CellTypeMNN.score", 
    "predicted.ClusterMNN.score", 
    "predicted.DevStageMNN.score",
    "combined_score"
    ), 
  order = TRUE
  )

combined_score_plot <- FeaturePlot(
  tumor, 
  features = c(
    "combined_score"
  ), 
  order = TRUE
)


ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/combined_score_plot.png", 
  combined_score_plot, 
  units = "in", 
  width = 10, 
  height = 8
)


combined_score_dev_vlnplot <- VlnPlot(tumor, 
        features = "combined_score",
        group.by = "predicted.DevStageMNN")

combined_score_ct_vlnplot <- VlnPlot(tumor, 
        features = "combined_score",
        group.by = "predicted.CellTypeMNN")



ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/combined_score_dev_vlnplot.png", 
  combined_score_dev_vlnplot, 
  units = "in", 
  width = 12, 
  height = 8
)

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/lognorm/combined_score_ct_vlnplot.png", 
  combined_score_ct_vlnplot, 
  units = "in", 
  width = 12, 
  height = 8
)
