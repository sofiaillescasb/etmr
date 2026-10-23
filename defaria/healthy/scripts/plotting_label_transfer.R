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

celltype_colors <-  c("Neurons" = "#7400b8",
                      "Neuron" = "#7400b8",
                      "Neuroblast" = "#5e60ce",
                      "Neuronal IPC" = "#48bfe3",
                      "Neural crest" = "#3FAF6B",
                      "Radial glia" = "#72efdd",
                      "Pericytes" = "#ea7af4",
                      "Fibroblast" = "#ea7af4",
                      "Endothelial" = "#a7c957",
                      "Placodes" = "#6e9887",
                      "Vascular" = "#bc4749",
                      "Glioblast" = "#ffea00",
                      "Immune" = "#ffb8de",
                      "Erythrocyte" = "#ff5d8f",
                      "Oligodendrocytes" = "#fb8b24",
                      "Oligo" = "#fb8b24",
                      "OPC" = "#eeba0b",
                      "Lymphocytes" = "#d90368")


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

tumor@meta.data |> 
  ggplot(aes(x = etmr_score, y = mapping.score, color = predicted.CellTypeMNN)) +
  geom_point() +
  facet_wrap(vars(predicted.ClusterMNN))

