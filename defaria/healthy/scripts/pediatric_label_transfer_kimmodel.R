library(anndataR)
library(Seurat)
library(qs2)

#from ~/Projects/etmr/defaria/healthy/scripts/kim_etal_celltypist_tumor.ipynb
tumor <- read_h5ad("/media/sofia/TOSHIBA EXT/developing_human_brain/tumor_annotated.h5ad",
                   as = "Seurat")
LayerData(tumor, assay = "RNA", layer = "counts") <- LayerData(tumor, assay = "RNA", layer = "X")

tumor <- RunUMAP(tumor, reduction = "X_pca", dims = 1:30)

kim_raw_predictions <- DimPlot(tumor, reduction = "umap", group.by = "predicted_labels") |
  DimPlot(subset(tumor, subset = conf_score > .5), reduction = "umap", group.by = "predicted_labels") |
    DimPlot(tumor, reduction = "umap", group.by = "by_group_celltypes")

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/kim_raw_predictions.png", 
       kim_raw_predictions,
       units="in",
       width = 21,
       height = 6)

celltype_clus <- c(
  "C0" = "Oligo OPALIN SEPHS1P6",
  "C1" = "Exc LINC00507 FAM163A",
  "C2" = "Astro FGFR3 LINC00500",
  "C3" = "Neuroblast LINC01210",
  "C4" = "RG HES1 WNT1",
  "C5" = "Exc RORB TRABD2A",
  "C6" = "Inh SST TNFSF15",
  "C7" = "Inh VIP TAC3",
  "C8" = "Inh PVLAB AKAIN1",
  "C9" = "OPC PDGFRA GPNMB",
  "C10" = "Exc RORB IZUMO3",
  "C11" = "Neuroblast NEUROD6",
  "C12" = "Exc Fetal STK32B",
  "C13" = "Micro TYROBP DOCK8",
  "C14" = "Exc FEZF2 MANCR",
  "C15" = "Exc LINC00507 C18orf61",
  "C16" = "Exc FEZF2 GHSR",
  "C17" = "Exc FEZF2 LINC01479",
  "C18" = "RG HES1 MRPS33P2",
  "C19" = "Inh ADARB2 LINC03069",
  "C20" = "Inh Fetal LIX1L-AS1",
  "C21" = "Inh LAMP5 CPLX3",
  "C22" = "Exc RORB REREP3",
  "C23" = "Inh Fetal LINC02802",
  "C24" = "Neuroblast GBX2",
  "C25" = "Exc FEZF2 CCDC80-1",
  "C26" = "Neuroblast ARHGAP19",
  "C27" = "Astro FGFR3 CD44",
  "C28" = "Inh LAMP5 NKX2-1",
  "C29" = "Exc Fetal TAS2R30",
  "C30" = "Endo NOSTRIN FOXC2",
  "C31" = "Inh PVALB NOG",
  "C32" = "OPC PDGFRA KIAA0040",
  "C33" = "Neuroblast NPS",
  "C34" = "Exc FEZF2 LINC00922",
  "C35" = "Exc SLC17A SNCB",
  "C36" = "Exc THEMIS IL7R",
  "C37" = "RG SOX2 NEBL-1",
  "C38" = "Oligo Neonatal ADAM33",
  "C39" = "Micro Fetal CCL3L3",
  "C40" = "Others RHD"
)

celltype_clus_gen <- c(
  "C0" = "Oligo",
  "C1" = "Exc",
  "C2" = "Astro",
  "C3" = "Neuroblast",
  "C4" = "RG",
  "C5" = "Exc",
  "C6" = "Inh",
  "C7" = "Inh",
  "C8" = "Inh",
  "C9" = "OPC",
  "C10" = "Exc",
  "C11" = "Neuroblast",
  "C12" = "Exc Fetal",
  "C13" = "Micro",
  "C14" = "Exc",
  "C15" = "Exc",
  "C16" = "Exc",
  "C17" = "Exc",
  "C18" = "RG",
  "C19" = "Inh",
  "C20" = "Inh",
  "C21" = "Inh",
  "C22" = "Exc",
  "C23" = "Inh",
  "C24" = "Neuroblast",
  "C25" = "Exc",
  "C26" = "Neuroblast",
  "C27" = "Astro",
  "C28" = "Inh",
  "C29" = "Exc Fetal",
  "C30" = "Endo",
  "C31" = "Inh",
  "C32" = "OPC",
  "C33" = "Neuroblast",
  "C34" = "Exc",
  "C35" = "Exc",
  "C36" = "Exc",
  "C37" = "RG",
  "C38" = "Oligo Neonatal",
  "C39" = "Micro Fetal",
  "C40" = "Others"
)

celltype_group <- c(
  "C0" = "Non-neuronal",
  "C1" = "Excitatory",
  "C2" = "Non-neuronal",
  "C3" = "Fetal",
  "C4" = "Fetal",
  "C5" = "Excitatory",
  "C6" = "Inhibitory",
  "C7" = "Inhibitory",
  "C8" = "Inhibitory",
  "C9" = "Non-neuronal",
  "C10" = "Excitatory",
  "C11" = "Fetal",
  "C12" = "Excitatory",
  "C13" = "Non-neuronal",
  "C14" = "Excitatory",
  "C15" = "Excitatory",
  "C16" = "Excitatory",
  "C17" = "Excitatory",
  "C18" = "Fetal",
  "C19" = "Inhibitory",
  "C20" = "Inhibitory",
  "C21" = "Inhibitory",
  "C22" = "Excitatory",
  "C23" = "Inhibitory",
  "C24" = "Fetal",
  "C25" = "Excitatory",
  "C26" = "Fetal",
  "C27" = "Non-neuronal",
  "C28" = "Inhibitory",
  "C29" = "Excitatory",
  "C30" = "Non-neuronal",
  "C31" = "Inhibitory",
  "C32" = "Non-neuronal",
  "C33" = "Fetal",
  "C34" = "Excitatory",
  "C35" = "Excitatory",
  "C36" = "Excitatory",
  "C37" = "Fetal",
  "C38" = "Non-neuronal",
  "C39" = "Non-neuronal",
  "C40" = "Non-neuronal"
)

celltype_dev_lov <- c(
  "C0" = "Adult_Cortex",
  "C1" = "Adult_Cortex",
  "C2" = "Various_Cortex",
  "C3" = "Fetal_1st_Various",
  "C4" = "Fetal_1st_Various",
  "C5" = "Adult_Cortex",
  "C6" = "Adult_Cortex",
  "C7" = "Adult_Cortex",
  "C8" = "Adult_Cortex",
  "C9" = "Various_Cortex",
  "C10" = "Adult_Cortex",
  "C11" = "Fetal_1st_Cortex",
  "C12" = "Fetal_2nd_Cortex",
  "C13" = "Adult_Cortex",
  "C14" = "Various_Cortex",
  "C15" = "Adult_Cortex",
  "C16" = "Adult_Cortex",
  "C17" = "Various_Cortex",
  "C18" = "Fetal_1st_Cortex",
  "C19" = "Adult_Cortex",
  "C20" = "Fetal_2nd_GE",
  "C21" = "Adult_Cortex",
  "C22" = "Adult_Cortex",
  "C23" = "Fetal_2nd_GE",
  "C24" = "Fetal_1st_Thalamus",
  "C25" = "Adult_Cortex",
  "C26" = "Fetal_2nd_GE",
  "C27" = "Adult_Cortex",
  "C28" = "Adult_Cortex",
  "C29" = "Fetal_2nd_Cortex",
  "C30" = "Various_Cortex",
  "C31" = "Adult_Cortex",
  "C32" = "Adult_Cortex",
  "C33" = "Fetal_1st_Cerebellum",
  "C34" = "Adult_Cortex",
  "C35" = "Adult_Cortex",
  "C36" = "Adult_Cortex",
  "C37" = "Fetal_2nd_GE",
  "C38" = "Neonatal_Cortex",
  "C39" = "Fetal_1st_Various",
  "C40" = "Fetal_1st_Various")

Idents(tumor) <- "predicted_labels"
current_ids <- as.character(Idents(tumor))

tumor@meta.data["predicted_celltype_clus"]  <- celltype_clus[current_ids]
tumor@meta.data["predicted_celltype_clus_gen"]  <- celltype_clus_gen[current_ids]
tumor@meta.data["predicted_celltype_group"]  <- celltype_group[current_ids]
tumor@meta.data["celltype_dev_loc"]  <- celltype_dev_lov[current_ids]

DimPlot(tumor, reduction = "umap", group.by = "celltype_dev_loc", label = TRUE) |
  DimPlot(subset(tumor, subset = conf_score > .5), reduction = "umap", group.by = "celltype_dev_loc") |
  DimPlot(tumor, reduction = "umap", group.by = "by_group_celltypes")

Idents(tumor) <- "predicted_celltype_group"


kim_model_gen <- DimPlot(tumor,
  reduction = "umap", 
  group.by = "predicted_celltype_clus_gen", label = TRUE)

library(tidyverse)

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/kim_model_gen.png", 
       kim_model_gen,
       units="in",
       width = 8,
       height = 6)


Idents(tumor) <- "predicted_celltype_clus_gen"

dev_loc_clus_gen_umap <- DimPlot(tumor,
  reduction = "umap", 
  split.by = "celltype_dev_loc",
  group.by = "predicted_celltype_clus_gen",
  label = TRUE,
  shuffle = TRUE
  ) +
  NoLegend()

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/dev_loc_clus_gen_umap.png", 
       dev_loc_clus_gen_umap,
       units="in",
       width = 30,
       height = 6)

qs_save(tumor, "/media/sofia/TOSHIBA EXT/developing_human_brain/tumor_annotated_seurat.qs2", nthreads = 6)
