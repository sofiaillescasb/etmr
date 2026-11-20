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

#Read data
cam <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/cam.qs2", nthreads = 6)
#maybe the others should be kept on disk????
braun <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/braun/data/braun.qs2", nthreads = 6)
herr <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/herr.qs2", nthreads = 6)


##########################################################################################################
#Read Herring


############################################## Cameron ###########################################################


cam$sex <- "F"
cam$region <- "frontal cortex"
cam$age <- "14-15ga"
cam$postbirth_age <- NA
cam$dev_stage <- "Fetal"
cam$dataset <- "Cameron"

###################################################################################################################



#Harmonize metadata colnames

colnames(braun@meta.data)
colnames(herr@meta.data)
colnames(cam@meta.data)


# I need donor_id, celltype, sex, region, UMIs, feature counts, age, dev_stage

donor <- as.factor(c(as.character(braun$donor_id),
           herr$sample,
           cam$sample))
sex <- as.factor(c(as.character(braun$sex),
         herr$Sex,
         cam$sex))
counts <- c(braun$nCount_RNA,
            herr$nCount_RNA,
            cam$nCount_RNA)
features <- c(braun$nFeature_RNA,
              herr$nFeature_RNA,
              cam$nFeature_RNA)
region <- as.factor(c(as.character(droplevels(braun$dissection)),
            herr$Brain.Regions.,
            cam$region))
gestational_age <- c(droplevels(braun$Age),
                     ifelse(grepl("ga", herr$age), herr$age, NA),
                     cam$age)
postbirth_age <- c(braun$postbirth_age,
                  ifelse(grepl("yr|d", herr$age), herr$age, NA),
                  cam$postbirth_age)
dev_stage <- as.factor(c(as.character(braun$dev_stage),
               herr$stage_id,
               cam$dev_stage))
celltype <- as.factor(c(as.character(braun$CellClass),
              ifelse(
                grepl(
                  "IN|PN", 
                  herr$cell_type
                  ), 
                paste(
                  herr$cell_type, 
                  herr$major_clust, 
                  sep="_"
                  ),
                herr$major_clust),
              gsub("FC-", "", ifelse(
                grepl(
                  "InN", 
                  cam$cellIDs), 
                paste(
                  "IN", 
                  cam$cellIDs, 
                  sep = "-"
                  ), 
                ifelse(
                  grepl(
                    "ExN", 
                    cam$cellIDs
                    ), 
                  paste("PN", cam$cellIDs, sep = "-"),
                  cam$cellIDs)))))
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

#percent_mt, percent_ribo i'll reprocess

#Merge 



         
library(org.Hs.eg.db)

# Get current Ensembl IDs (remove version numbers if present, e.g., .1, .2)
herr_ensembl_ids <- rownames(herr)

# Map to gene symbols
herr_gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = herr_ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

meta_lookup <- gene_meta |>
  mutate(clean_key = sub("\\.[0-9]+", "", Accession),
         clean_key = sub("_PAR_[XY]", "", clean_key)) |>
  filter(!is.na(Gene)) |>
  distinct(clean_key, .keep_all = TRUE)

# Turn it into a named vector (Name = clean Ensembl ID, Value = Gene Symbol)
symbol_dict <- setNames(meta_lookup$Gene, meta_lookup$clean_key)

# 2. Clean your matrix row names directly 
raw_names <- rownames(herr)
clean_keys <- sub("\\.[0-9]+", "", raw_names)
clean_keys <- sub("_PAR_[XY]", "", clean_keys)

# 3. Look up symbols; if none found, keep the original name
mapped_names <- symbol_dict[clean_keys]
final_names <- coalesce(mapped_names, raw_names)

rownames(herr) <- make.unique(final_names)

merged_reference <- list(braun, herr, cam)
  

merged_reference <- merge(
  x = merged_reference[[1]],
  y = merged_reference[-1],
  project = "Healthy_reference",
  merge.data = FALSE
)

merged_reference[["percent.mt"]] <- PercentageFeatureSet(merged_reference, pattern = "^MT-")

merged_reference <- subset(
  merged_reference,
  subset = nFeature_RNA > 250 &
    nFeature_RNA < 9000 &
    nCount_RNA < 75000 &
    percent.mt < 15
)


merged_reference@meta.data <- metadata_joined

merged_reference <- JoinLayers(merged_reference)

merged_reference[["RNA"]] <- split(merged_reference[["RNA"]], f = merged_reference$id)


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


merged_reference$harmonized_cell_type <- recode(
  merged_reference$cell_type,
  "Radial glia"       = "Radial Glia",
  "RG-1"                = "Radial Glia",
  "RG-2" = "Radial Glia",
  
  "Excitatory_Neurons" = "Excitatory Neurons",
  "PN-ExN-1"            = "Excitatory Neurons",
  "PN-ExN-2"            = "Excitatory Neurons",
  "PN-ExN-3"            = "Excitatory Neurons",
  "PN-ExN-4"            = "Excitatory Neurons",
  "PN-ExN-5"            = "Excitatory Neurons",
  "PN_L2-3_CUX2"            = "Excitatory Neurons",
  "PN_L4_RORB"            = "Excitatory Neurons",
  "PN_L5-6_TLE4"            = "Excitatory Neurons",
  "PN_L5-6_TLE4"            = "Excitatory Neurons",
  
  "IN_SST"       = "Inhibitory Neurons",
  "IN_ID2"       = "Inhibitory Neurons",
  "IN-InN-1"       = "Inhibitory Neurons",
  "IN-InN-2"       = "Inhibitory Neurons",
  "IN-InN-3"       = "Inhibitory Neurons",
  "IN-InN-4"       = "Inhibitory Neurons",
  "IN_PV"          = "Inhibitory Neurons",
  "IN_LAMP5_NOS1"       = "Inhibitory Neurons",
  "IN_VIP"       = "Inhibitory Neurons",
  "IN_PV_SCUBE3"          = "Inhibitory Neurons",
  
  "Astro" = "Astrocytes",
  
  "Oligo" = "Oligodendrocytes",
  
  "MG" = "Microglia",
  "Micro" = "Microglia",
  "Immune" = "Microglia",
  
  "Neuronal IPC" = "IP",
  "CycPro" = "IP",
)



qs_save(merged_reference, "/media/sofia/TOSHIBA EXT/developing_human_brain/merged_reference.qs2", nthreads = 7)

