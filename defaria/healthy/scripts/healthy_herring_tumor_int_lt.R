library(anndataR)
library(Seurat)
library(qs2)

#from ~/Projects/etmr/defaria/healthy/scripts/kim_etal_celltypist_tumor.ipynb
tumor <- read_h5ad("/media/sofia/TOSHIBA EXT/developing_human_brain/tumor_annotated.h5ad",
                   as = "Seurat")

LayerData(tumor, assay = "RNA", layer = "counts") <- LayerData(tumor, assay = "RNA", layer = "X")


DefaultAssay(tumor) <- "RNA"

tumor <- DietSeurat(tumor, 
                    assays = "RNA",
                    layers = "counts")

tumor[["RNA"]] <- split(tumor[["RNA"]], f = tumor$orig.ident)

tumor <- NormalizeData(tumor)

tumor <- FindVariableFeatures(tumor)

tumor <- ScaleData(tumor)

healthy_int <- qs_read(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/healthy_int.qs2", 
  nthreads = 6
)

healthy_int$condition <- "human_healthy"

healthy_int <- JoinLayers(healthy_int)

gc()

healthy_int <- RunUMAP(healthy_int, 
                       reduction = "integrated.mnn", 
                       dims = 1:30,
                       reduction.name = "mnn_umap",
                       return.model = TRUE)

anchors.mnn <- FindTransferAnchors(
  reference = healthy_int,
  query = tumor,
  normalization.method = "LogNormalize", 
  reduction = "pcaproject",
  reference.reduction = "integrated.mnn"
)

qs_save(anchors.mnn, "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/anchors.mnn_herr_lt.qs2", 
        nthreads = 6)
gc()

tumor <- MapQuery(
  anchorset = anchors.mnn,
  reference = healthy_int,
  query = tumor,
  refdata = list(celltype = "cell_type",
                 major_clust = "major_clust",
                 sub_clust = "sub_clust",
                 stage_id = "stage_id",
                 region = "Brain.Regions."), # Change this to your actual cell type column name in healthy_int
  reference.reduction = "pca", 
  reduction.model = "mnn_umap"
)

gc()

qs_save(tumor, "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/tumor_herr_lt.qs2", 
        nthreads = 6)


tumor <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/processed/tumor_herr_lt.qs2", 
        nthreads = 6)

tumor[["umap"]] <- tumor[["ref.umap"]]

gc()

tumor_healthy <- merge(healthy_int, y = tumor)

rm(healthy_int)
gc()

Idents(tumor_healthy) <- "orig.ident"

age <- c(
  "GSM8058157" = 31.2, 
  "GSM8058158" = 12, 
  "GSM8058159" = 30, 
  "GSM8058160" = 29, 
  "GSM8058161" = 12, 
  "GSM8058162" = 35, 
  "GSM8058163" = 48, 
  "GSM8058164" = 34, 
  "GSM8058165"= 39.6) * 30

    DimPlot(tumor, group.by = "predicted.celltype", split.by = "predicted.stage_id", label = TRUE)

tumor@meta.data |>
  group_by(predicted.celltype) |>
  ggplot(aes(x = ))
summary(tumor$predicted.celltype.score)

tumor_healthy@meta.data["Age"]  <- age[Idents(tumor_rgc)]

ifelse(

  tumor_rgc@meta.data["Age"]  <- age[Idents(tumor_rgc)])

DimPlot(
  tumor_healthy, 
  reduction = "umap", 
  group.by = c("condition", "cell_type", "by_group_celltypes"), 
  shuffle = TRUE,      # Randomly mixes the points so one group doesn't completely hide the other
pt.size = 0.5
)
