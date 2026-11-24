library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(R.utils)
library(Azimuth)
library(reticulate)
library(dplyr)
library(data.table)
library(Matrix)
library(qs2)
library(org.Hs.eg.db)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

#Read data
cam <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/cam.qs2", nthreads = 6)
braun <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/braun/braun.qs2", nthreads = 6)
herr <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/herring_2022/herr.qs2", nthreads = 6)


##########################################################################################################
clean_name <- sub("\\.[0-9]+$", "", Features(braun))
og_names <- rownames(braun)
rownames(braun) <- clean_name

length(intersect(Features(herr), Features(braun)))

merged_reference <- merge(braun,
                          herr,
                          add.cell.ids = c("braun", "herr"))

merged_reference <- JoinLayers(merged_reference)

ensembl_ids <- rownames(merged_reference)


# Map to gene symbols #################################################################################
gene_symbols <- mapIds(org.Hs.eg.db,
                            keys = ensembl_ids,
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")

edb_all <- genes(
  EnsDb.Hsapiens.v86,
  columns = c("gene_id", "gene_name", "gene_biotype", "seqnames")
)
edb_df <- as.data.frame(edb_all)

# Filter afterwards from the single table
mirna_df <- subset(edb_df, gene_biotype == "miRNA")
ribo_df  <- subset(edb_df, gene_biotype == "protein_coding" & grepl("^RP[SL]", gene_name))
hemo_df <- subset(
  edb_df,
  gene_biotype == "protein_coding" & grepl("^HB[ABDEGMQZ][0-9]*$", gene_name)
)

gene_names <- gene_symbols

# Match against the 20,551 expressed NA-symbol IDs
is_mir_target <- ensembl_ids %in% mirna_df$gene_id
mir_hits_ids  <- ensembl_ids[is_mir_target]

mir_names <- mirna_df$gene_name[match(ensembl_ids[is_mir_target], mirna_df$gene_id)]

gene_names[is_mir_target] <- mir_names

#now mitochondrial genes
mt_edb <- genes(
  EnsDb.Hsapiens.v86,
  filter = SeqNameFilter("MT"),
  columns = c("gene_id", "gene_name", "gene_biotype")
)
mt_df <- as.data.frame(mt_edb)

is_mt_target <- ensembl_ids %in% mt_df$gene_id
mt_hits_ids  <- ensembl_ids[is_mt_target]

mt_names <- mt_df$gene_name[match(ensembl_ids[is_mt_target], mt_df$gene_id)]

gene_names[is_mt_target] <- mt_names

#now ribosomal genes
is_ribo_target <- ensembl_ids %in% ribo_df$gene_id
ribo_hits_ids  <- ensembl_ids[is_ribo_target]

ribo_names <- ribo_df$gene_name[match(ensembl_ids[is_ribo_target], ribo_df$gene_id)]

gene_names[is_ribo_target] <- ribo_names

#now hemoglobin genes
hemo_df <- subset(
  edb_df,
  gene_biotype == "protein_coding" & grepl("^HB[ABDEGMQZ][0-9]*$", gene_name)
)

is_hemo_target <- ensembl_ids %in% hemo_df$gene_id
hemo_hits_ids  <- ensembl_ids[is_hemo_target]

hemo_names <- hemo_df$gene_name[match(ensembl_ids[is_hemo_target], hemo_df$gene_id)]

gene_names[is_hemo_target] <- hemo_names


#Filter out genes with 0 expression

counts <- LayerData(merged_reference, assay = "RNA", layer = "counts")

keep <- Matrix::rowSums(counts) > 0

merged_reference <- subset(
  merged_reference,
  features = rownames(counts)[keep]
)

gene_names <- gene_names[keep]
ensembl_ids <- ensembl_ids[keep]

#This has to be TRUE
identical(na.omit(names(gene_names)), na.omit(Features(merged_reference)))

mapped <- gene_names[ensembl_ids]

# Fall back to original Ensembl ID if mapped name is NA or blank
final_names <- ifelse(is.na(mapped) | mapped == "", ensembl_ids, mapped)

# Ensure uniqueness (e.g., if two Ensembl IDs map to the same gene symbol)
final_names <- make.unique(final_names)

rownames(merged_reference) <- final_names

############################################## Cameron ###########################################################

cam$sex <- "F"
cam$region <- sub(".*_(.*?)_.*", "\\1", cam$sample)
cam$age <- "14-15ga"
cam$postbirth_age <- NA
cam$dev_stage <- "Fetal"
cam$dataset <- "Cameron"

merged_reference <- merge(merged_reference,
                          cam,
                          add.cell.ids = c("ref", "cam"))



################################################QC#############################################################


############################################################################################################


#Harmonize metadata colnames

colnames(braun@meta.data)
colnames(herr@meta.data)
colnames(cam@meta.data)


# I need donor_id, celltype, sex, region, UMIs, feature counts, age, dev_stage

donor <- as.factor(c(as.character(braun$donor_id),
           herr$sample,
           as.character(cam$orig.ident)))
sex <- as.factor(c(as.character(braun$sex),
                   as.character(herr$Sex),
         cam$sex))
counts <- c(braun$nCount_RNA,
            herr$nCount_RNA,
            cam$nCount_RNA)
features <- c(braun$nFeature_RNA,
              herr$nFeature_RNA,
              cam$nFeature_RNA)
region <- as.factor(c(as.character(droplevels(braun$dissection)),
                      as.character(herr$Brain.Regions.),
            cam$region))
gestational_age <- c(as.character(droplevels(braun$Age)),
                     ifelse(grepl("ga", as.character(herr$age)), as.character(herr$age), NA),
                     as.character(cam$age))
postbirth_age <- c(as.character(braun$postbirth_age),
                  ifelse(grepl("yr|d", as.character(herr$age)), as.character(herr$age), NA),
                  as.character(cam$postbirth_age))
dev_stage <- as.factor(c(as.character(braun$dev_stage),
               as.character(herr$stage_id),
               as.character(cam$dev_stage)))
celltype <- as.factor(c(as.character(braun$CellClass),
              ifelse(
                grepl(
                  "IN|PN", 
                  as.character(herr$cell_type)
                  ), 
                paste(
                  as.character(herr$cell_type), 
                  as.character(herr$major_clust), 
                  sep="_"
                  ),
                as.character(herr$major_clust)),
              gsub("FC-", "", ifelse(
                grepl(
                  "InN", 
                  as.character(cam$cellIDs)), 
                paste(
                  "IN", 
                  as.character(cam$cellIDs), 
                  sep = "-"
                  ), 
                ifelse(
                  grepl(
                    "ExN", 
                    as.character(cam$cellIDs)
                    ), 
                  paste("PN", as.character(cam$cellIDs), sep = "-"),
                  as.character(cam$cellIDs))))))
dataset <- c(braun$dataset,
             herr$dataset,
             cam$dataset) 

metadata_joined <- data.frame(
  "id" = donor, 
  "dataset" = dataset,
  "sex" = sex,
  "nCount_RNA" = counts,
  "nFeature_RNA" = features,
  "region" = region,
  "gestational_age" = gestational_age,
  "postbirth_age" = postbirth_age,
  "dev_stage" = dev_stage,
  "celltype" = celltype)

og_met  <- merged_reference@meta.data

merged_reference@meta.data <- metadata_joined

#percent_mt, percent_ribo i'll reprocess

merged_reference[["percent.mt"]] <- PercentageFeatureSet(merged_reference, pattern = "^MT-")

merged_reference[["percent.hb"]] <- PercentageFeatureSet(
  merged_reference,
  features = hemo_names
)

merged_reference[["percent.ribo"]] <- PercentageFeatureSet(
  merged_reference,
  features = ribo_names
)

summary(merged_reference@meta.data[, c("percent.mt", "percent.ribo", "percent.hb")])

# percent.mt       percent.ribo   
# Min.   : 0.0000   Min.   : 0.000  
# 1st Qu.: 0.8994   1st Qu.: 7.292  
# Median : 1.6128   Median :11.030  
# Mean   : 2.2659   Mean   :11.523  
# 3rd Qu.: 3.1412   3rd Qu.:15.400  
# Max.   :25.0000   Max.   :48.051  
# percent.hb      
# Min.   : 0.00000  
# 1st Qu.: 0.00000  
# Median : 0.00000  
# Mean   : 0.25189  
# 3rd Qu.: 0.01661  
# Max.   :93.95925  

# > max(braun$nCount_RNA)
# [1] 93360
# > max(herr$nCount_RNA)
# [1] 85170
# > max(cam$nCount_RNA)
# [1] 55258
# > max(cam$nFeature_RNA)
# [1] 4990
# > max(herr$nFeature_RNA)
# [1] 10821
# > max(braun$nFeature_RNA)
# [1] 11024
# > min(braun$nFeature_RNA)
# [1] 48
# > min(herr$nFeature_RNA)
# [1] 533
# > min(cam$nFeature_RNA)
# [1] 987

merged_reference <- subset(
  merged_reference,
  subset = nFeature_RNA > 250 &
    nFeature_RNA < 9000 &
    nCount_RNA < 75000 &
    percent.mt < 15
)

merged_reference <- JoinLayers(merged_reference)

# Erythrocyte    Fibroblast     Glioblast     
# Immune         Neuroblast     Neural crest  
# Neuron         Neuronal IPC   Oligo         
# Radial glia    Vascular       PN_L4_RORB    
# PN_L2-3_CUX2   IN_SST         Astro         
# PN_L5-6_TLE4   PN_L5-6_THEMIS IN_VIP        
# OPC            IN_PV          IN_PV_SCUBE3  
# IN_ID2         IN_LAMP5_NOS1  Micro         
# Vas            IN_MGE_dev     Poor-Quality  
# PN_PN_dev      IN_CGE_dev     IN-InN-2      
# IP             RG-1           PN-ExN-1      
# PN-ExN-4       IN-InN-1       PN-ExN-3      
# IN-InN-3       RG-2           N-undef       
# PN-ExN-5       PN-ExN-2       CycPro        
# IN-InN-4       Endo           MG   

char_celltype <- as.character(merged_reference$celltype)

merged_reference$harmonized_cell_type <- case_when(
  grepl("(^|[-_])(PN|ExN)([-_]|$)", char_celltype) ~ "Excitatory Neurons",
  grepl("(^|[-_])(IN|InN)([-_]|$)", char_celltype) ~ "Inhibitory Neurons",
  grepl("RG-|Radial glia", char_celltype) ~ "Radial Glia",
  grepl("Astro", char_celltype) ~ "Astrocytes",
  grepl("Micro|MG|Immune", char_celltype) ~ "Microglia",
  grepl("Vas|Endo", char_celltype) ~ "Vascular",
  grepl("Endo", char_celltype) ~ "Endothelial",
  grepl("OPC", char_celltype) ~ "OPC",
  grepl("Neuronal IPC|(^|[-_])(IP|IPC)([-_]|$)", char_celltype) ~ "IPC",
  grepl("CycPro", char_celltype) ~ "CycPro",
  grepl("Neuron|Neuroblast|Hipp-CR", char_celltype) ~ "Other Neurons / Neuroblasts",
  TRUE ~ char_celltype
)

#maybe cycpro is ipc

library(tidyverse)

filtered_qc <- og_met |>
  rownames_to_column("cell") |>
  pivot_longer(cols = c(
    nCount_RNA, 
    nFeature_RNA, 
    percent.mt, 
    percent.ribo, 
    percent.hb
    ), 
    names_to = "qc") |>
  group_by(id) |>
  ggplot(aes(x = value, fill = dataset, colour = dataset)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  facet_wrap(vars(qc), scales = "free") +
  scale_fill_discrete(palette = "Set2") +
  scale_color_discrete(palette = "Set2") +
  theme_classic()

ggsave("/media/sofia/TOSHIBA EXT/developing_human_brain/plots/qc/filtered_qc.png", 
       filtered_qc,
       width = 10,
       height = 6,
       units = "in")

qs_save(
  merged_reference, 
  "/media/sofia/TOSHIBA EXT/developing_human_brain/merged_reference.qs2", 
  nthreads = 7
  )

