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


#Read Braun
braun_path <- "/media/sofia/TOSHIBA EXT/developing_human_brain/braun/human_dev.h5ad"

braun_mat <- open_matrix_anndata_hdf5(braun_path)

mat_counts <- convert_matrix_type(braun_mat, type = "uint32")

# Write to a compressed binary folder on disk
write_matrix_dir(mat = mat_counts, dir = "/media/sofia/TOSHIBA EXT/developing_human_brain/braun/braun_human_counts")


#Need to read metadata and gene symbols
# This is the only way I think this can be done on disk, 
# because reticulate converts pandas dataframes into number codes with levels as categories

ad <- import("anndata")

adata <- ad$read_h5ad("human_dev.h5ad", backed = "r")

cell_meta <- py_to_r(adata$obs)
gene_meta <- py_to_r(adata$var)

braun <- CreateSeuratObject(
  counts = mat_counts,
  meta.data = cell_meta,
  project = "Braun")

braun[["RNA"]] <- split(braun[["RNA"]], f = braun$donor_id)

braun_cortex <- subset(braun, subset = Subregion == "Cortex")
braun_cortex$sex <- braun_cortex$sex_ontology_term_id
braun_cortex$postbirth_age <- NA
braun_cortex$dev_stage <- "Fetal"
braun_cortex$dataset <- "Braun"


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
  select(all_of(sel_cols))

rownames(herr_metadata) <- sub("^([^-]+)-([^_]+)_.*$", "\\2_\\1-1", herr_metadata$X)

herr <- AddMetaData(herr, metadata = herr_metadata)

cells_in_common <- intersect(herr$X, herr_metadata$X)

herr <- subset(herr, subset = X %in% cells_in_common)
herr$dataset <- "Herring"

#Read Cameron 
#https://www.sciencedirect.com/science/article/pii/S0006322322014044
#FC is frontal lobal
gz_path <- "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/cameron_2022_snRNAseq_FC_raw_count_gEX_matrix.txt.gz"

# #unzip file
# gunzip(gz_path, remove = FALSE, overwrite = TRUE)

txt_path <- sub("\\.gz$", "", gz_path)

#load the uncompressed count table efficiently
counts_dt <- fread(txt_path)


# Convert the first column (gene symbols/IDs) to row names
gene_names <- counts_dt[[1]]
count_matrix <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(count_matrix) <- gene_names

# 3 karyotypically normal human fetuses from the second trimester of gestation (2 of 14 and 1 of 15 postconception weeks, all female).


# Initialize the Seurat object with raw counts
cam <- CreateSeuratObject(
  counts = count_matrix, 
  project = "Prenatal"
)

# Read metadata
metadata <- fread("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/cameron_2022_snRNAseq_FC_metadata.txt.gz")


cam <- AddMetaData(cam, metadata = as.data.frame(metadata))

cam$sex <- "F"
cam$region <- "frontal cortex"
cam$age <- "14-15ga"
cam$postbirth_age <- NA
cam$dev_stage <- "Fetal"
cam$dataset <- "Cameron"
#Harmonize metadata colnames

colnames(braun_cortex@meta.data)
colnames(herr@meta.data)
colnames(cam@meta.data)


# I need donor_id, celltype, sex, region, UMIs, feature counts, age, dev_stage

donor <- as.factor(c(as.character(braun_cortex$donor_id),
           herr$sample,
           cam$sample))
sex <- as.factor(c(as.character(braun_cortex$sex),
         herr$Sex,
         cam$sex))
counts <- c(braun_cortex$nCount_RNA,
            herr$nCount_RNA,
            cam$nCount_RNA)
features <- c(braun_cortex$nFeature_RNA,
              herr$nFeature_RNA,
              cam$nFeature_RNA)
region <- as.factor(c(as.character(droplevels(braun_cortex$dissection)),
            herr$Brain.Regions.,
            cam$region))
gestational_age <- c(droplevels(braun_cortex$Age),
                     ifelse(grepl("ga", herr$age), herr$age, NA),
                     cam$age)
postbirth_age <- c(braun_cortex$postbirth_age,
                  ifelse(grepl("yr|d", herr$age), herr$age, NA),
                  cam$postbirth_age)
dev_stage <- as.factor(c(as.character(braun_cortex$dev_stage),
               herr$stage_id,
               cam$dev_stage))
celltype <- as.factor(c(as.character(braun_cortex$CellClass),
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
dataset <- c(braun_cortex$dataset,
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

og_rownames <- rownames(braun_cortex)
new_rownames <- setNames(gene_meta$Gene, gene_meta$Accession)
# 3. Use the named vector to replace them instantly (unmatched IDs keep their original name)
new_names <- new_rownames[og_rownames]
new_names[is.na(new_names)] <- og_rownames[is.na(new_names)] # Keep old ID if no symbol found

# 4. Assign back
rownames(braun_cortex) <- new_names


         
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

merged_reference <- list(braun_cortex, herr, cam)
  

merged_reference <- merge(
  x = merged_reference[[1]],
  y = merged_reference[-1],
  project = "Healthy_reference",
  merge.data = FALSE
)

merged_reference@meta.data <- metadata_joined

merged_reference <- JoinLayers(merged_reference)

merged_reference[["RNA"]] <- split(merged_reference[["RNA"]], f = merged_reference$id)

qs_save(merged_reference, "/media/sofia/TOSHIBA EXT/developing_human_brain/merged_reference.qs2", nthreads = 7)

