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
library(Matrix)


#https://www.sciencedirect.com/science/article/pii/S0006322322014044

#FC is frontal cortex
gz_path <- paste0(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/counts/", 
  list.files(
    "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/counts"
  )
)



#unzip file (i think this duplicates the files)
lapply(
  gz_path,
  function(f) gunzip(
    f, 
    destname = sub("/counts/", "/counts/uncompressed/", sub("\\.gz$", "",f)),
    remove = FALSE, 
    overwrite = TRUE
  )
)

gz_metadata_path <- paste0(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/metadata/", 
  list.files(
    "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/metadata"
  )
)


#unzip file
lapply(
  gz_metadata_path,
  function(f) gunzip(
    f, 
    destname = sub("/metadata/", "/metadata/uncompressed/", sub("\\.gz$", "",f)),
    remove = FALSE, 
    overwrite = TRUE
  )
)

uncompressed_metadata_path <- paste0(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/metadata/uncompressed/", 
  list.files(
    "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/metadata/uncompressed"
  )
)

metadata_cam <- lapply(uncompressed_metadata_path, function(f) fread(f))
names(metadata_cam) <- sapply(
   uncompressed_metadata_path,
        function(x) gsub(
   "^.*/cameron_2022_snRNAseq_|_metadata\\.txt$",
   "",
  x))
metadata_cam <- lapply(metadata_cam, function(x) as.data.frame(x))
metadata_cam <- bind_rows(metadata_cam, .id = "id") 
rownames(metadata_cam) <- metadata_cam$cells

uncompressed_path <- paste0(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/counts/uncompressed/", 
  list.files(
    "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/counts/uncompressed"
  )
)

for (i in seq_along(uncompressed_path)) {
  
  # Read one file
  dt <- fread(uncompressed_path[i])
  
  # Extract gene IDs
  genes <- dt[[1]]
  
  # Convert counts to integer matrix
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "integer"
  
  rownames(mat) <- genes
  
  # Convert to sparse matrix
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Save directly to disk
  qs_save(mat, file = paste0("/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/counts/matrix/", i, ".qs2"))
  
  # Free memory
  rm(dt, mat, genes)
  gc()
}



matrix_path <- paste0(
  "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/counts/matrix/", 
  list.files(
    "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/counts/matrix"
  )
)

tissues <- c("Cer", "FC", "GE", "Hipp", "Thal")

merged_matrices <- lapply(matrix_path, function(f) {
  qs_read(f)
})

names(merged_matrices) <- tissues

all_genes <- Reduce(union, lapply(merged_matrices, rownames))

aligned_matrices <- lapply(merged_matrices, function(m) {
  m <- as(m, "dgCMatrix")
  missing_genes <- setdiff(all_genes, rownames(m))
  
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(m),
      dimnames = list(missing_genes, colnames(m)),
      sparse = TRUE
    )
    m <- rbind(m, zero_mat)
  }
  
  m[all_genes, , drop = FALSE]
})

merged_counts <- do.call(cbind, aligned_matrices)

# 3 karyotypically normal human fetuses from the second trimester of gestation (2 of 14 and 1 of 15 postconception weeks, all female).

# Initialize the Seurat object with raw counts
cam <- CreateSeuratObject(
  counts = aligned_matrices, 
  meta.data = metadata_cam,
  project = "Cameron2022"
)

cam <- JoinLayers(cam)

qs_save(cam, "/media/sofia/TOSHIBA EXT/developing_human_brain/cameron_2022/data/cam.qs2", nthreads = 6)

