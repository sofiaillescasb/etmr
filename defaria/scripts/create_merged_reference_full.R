library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(R.utils)
library(data.table)
library(qs2)
library(Azimuth)
library(reticulate)
library(dplyr)
library(data.table)
library(Matrix)
library(qs2)

#Read data
cam <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/cam.qs2", nthreads = 6)
#maybe the others should be kept on disk????
braun <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/braun.qs2", nthreads = 6)
herr <- qs_read("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/herr.qs2", nthreads = 6)


##########################################################################################################
#Read Herring

herr_path <- "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/"

herr_files <-c("RL1612_filtered_feature_bc_matrix.h5",
                "RL1613_filtered_feature_bc_matrix.h5",
                "RL1777_filtered_feature_bc_matrix.h5",
                "RL2100_filtered_feature_bc_matrix.h5",
                "RL2103_filtered_feature_bc_matrix.h5",
                "RL2104_filtered_feature_bc_matrix.h5",
                "RL2105_filtered_feature_bc_matrix.h5",
                "RL2107_filtered_feature_bc_matrix.h5",
                "RL2108_filtered_feature_bc_matrix.h5",
                "RL2121_filtered_feature_bc_matrix.h5",
                "RL2122_filtered_feature_bc_matrix.h5",
                "RL2125_filtered_feature_bc_matrix.h5",
                "RL2129_filtered_feature_bc_matrix.h5")

for (i in seq_along(herr_files)) {
  path <- file.path(herr_path, herr_files[i])
  data <- open_matrix_10x_hdf5(path)
  matrix.dir <- file.path(herr_path, paste0(tools::file_path_sans_ext(herr_files[i]), "_BP"))
  write_matrix_dir(mat = data, dir = matrix.dir)
}

data.list <- c()

for (i in 1:length(herr_files)) {
  path <- paste0(herr_path, herr_files[i])
  matrix.dir <- file.path(herr_path, paste0(tools::file_path_sans_ext(herr_files[i]), "_BP"))
  # Load in BP matrices
  mat <- open_matrix_dir(dir = matrix.dir)
  data.list[[i]] <- mat
}

# Name layers
names(data.list) <- c("RL1612",
                      "RL1613",
                      "RL1777",
                      "RL2100",
                      "RL2103",
                      "RL2104",
                      "RL2105",
                      "RL2107",
                      "RL2108",
                      "RL2121",
                      "RL2122",
                      "RL2125",
                      "RL2129")

herr_seurat <- lapply(names(data.list), function(sample_name) {
  
  mat <- data.list[[sample_name]]
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )
  
  obj$sample <- sample_name
  
  obj <- RenameCells(
    obj,
    add.cell.id = sample_name
  )
  
  obj
})

names(herr_seurat) <- names(data.list)

herr <- merge(
  x = herr_seurat[[1]],
  y = herr_seurat[-1],
  project = "herr",
  merge.data = FALSE
)

herr_metadata <- read.csv(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE168408/Processed_data_RNA-all_BCs-meta-data.csv"
)

sel_cols <- c("X", "batch",                     
              "RL.", "age",                      
              "chem", "concat_id",                
              "numerical_age", "stage_id",                  
              "Sex", "Race",                      
              "PMI", "Brain.Regions.",            
              "Cause.of.Death", "ICD.10.Code",               
              "ICD.10.category", "Oxygen.No.Oxygen",          
              "Date.of.Collection", "Collection_year",           
              "Library.Prep.Date", "Library.Prep.Lot",          
              "n_genes_by_counts", "log1p_n_genes_by_counts",   
              "total_counts", "log1p_total_counts",        
              "doublet_score", "log10_gene_counts",        
              "log10_UMI_counts", "percent_mito",              
              "percent_ribo", "n_counts",                  
              "leiden", "mat.dev",                   
              "cell_type", "major_clust",               
              "sub_clust", "combined.leiden")


herr_metadata <- herr_metadata |>
  tibble::column_to_rownames("X") |>
  select(all_of(sel_cols)) |>
  filter(RL. %in% unique(herr@meta.data$orig.ident))

rownames(herr_metadata) <- sub("^([^-]+)-([^_]+)_.*$", "\\2_\\1-1", rownames(herr_metadata))

herr <- AddMetaData(herr, metadata = herr_metadata)

cells_in_common <- intersect(herr$X, herr_metadata$X)

herr <- subset(herr, subset = X %in% cells_in_common)
herr$dataset <- "Herring"

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

