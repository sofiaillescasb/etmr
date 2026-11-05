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
