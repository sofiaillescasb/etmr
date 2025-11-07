library(here)
library(dplyr)
library(tibble)
library(Seurat)
library(SeuratDisk)

ids <- c(
  "GSM8058157_KK22-H-225",
"GSM8058158_KK22-H-226", 
"GSM8058159_KK22-H-227",
"GSM8058160_KK23-H-507",
"GSM8058161_KK23-H-508",
"GSM8058162_KK23-H-509",
"GSM8058163_KK23-H-510",
"GSM8058164_KK23-H-511",
"GSM8058165_KK23-H-512"
)

path <- sapply(ids, function(x) paste0(here::here(), "/de Faria/snRNA/data/GSE254819_RAW/", x, "/"))

read_data <- Read10X(path)

snrna.data <- CreateSeuratObject(counts = read_data, project = "snfixedtumor")

sample_info <- df <- tribble(
  ~orig.ident,      ~Sample_ID,   ~Tumor_Type,
  "GSM8058157", "KK22-H-225", "C19MC",
  "GSM8058158", "KK22-H-226", "C19MC",
  "GSM8058159", "KK22-H-227", "Dicer",
  "GSM8058160", "KK23-H-507", "C19MC",
  "GSM8058161", "KK23-H-508", "C19MC",
  "GSM8058162", "KK23-H-509", "C19MC",
  "GSM8058163", "KK23-H-510", "C19MC",
  "GSM8058164", "KK23-H-511", "C19MC",
  "GSM8058165", "KK23-H-512", "C19MC")


snrna.data@meta.data <- snrna.data@meta.data %>%
  rownames_to_column(var = "rownames") %>%
  right_join(sample_info, by = "orig.ident") %>%
  column_to_rownames(var = "rownames")


saveRDS(snrna.data, paste0(here::here(), "/de Faria/snRNA/data/snrna_data.rds"))


