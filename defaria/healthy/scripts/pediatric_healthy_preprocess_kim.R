#This didn't work. I decided to try using the model they trained ~/Projects/etmr/defaria/healthy/scripts
# ~/Projects/etmr/defaria/healthy/scripts/pediatric_label_transfer_kimmodel.R


library(Seurat)
library(qs2)
library(tidyverse)

pediatric_ref <- readRDS("/media/sofia/TOSHIBA EXT/developing_human_brain/GSE280569_Steyn_etal_seurat_object.rds")
pediatric_ref_bts <- readRDS("/media/sofia/TOSHIBA EXT/developing_human_brain/BTS_atlas_seurat.rds")


gc()

table(pediatric_ref$age_group)

# Adult Paediatric 
# 49213      95225

pediatric_ref <- subset(pediatric_ref, subset = age_group == "Paediatric")
table(pediatric_ref$brain_region)

qs_save(
  pediatric_ref, 
  "/media/sofia/TOSHIBA EXT/developing_human_brain/GSE280569_Steyn_etal_pediatric.qs2", 
  nthreads = 6
  )

gc()

table(pediatric_ref_bts$`Brain Region`)

# BA8                 BA9                BA10 
# 37888              176459               11211 
# BA46      Frontal Cortex               Brain 
# 7879                9396                5888 
# Caudate+Putamen   Cerebellum   Cortex 
# 898               10288                8026 
# Cortex entorhinal Cortex temporal     Diencephalon 
# 956                2228                2114 
# Forebrain Ganglionic Eminence           Hindbrain 
# 4889                7637                1969 
# Hippocampus        Hypothalamus        M1 
# 1596                2303               74482 
# Medulla            Midbrain    Midbrain ventral 
# 5161                6399                 412 
# Pons           Subcortex            Thalamus 
# 2382                3090                9509 

unique(pediatric_ref_bts$Stage)

# [1] Fetal (1st trimester) Adult (40-60Y)        Fetal (2nd trimester)
# [4] Adult (60-80Y)        Fetal (3rd trimester) Neonatal             
# [7] Childhood (1-6Y)      Childhood (6-12Y)     Adolescence (12-20Y) 
# [10] Adult (20-40Y)        Adult (>80Y)  

pediatric_ref_bts <- subset(
  pediatric_ref_bts, 
  subset = Stage %in% c(
    "Fetal (1st trimester)",
    "Fetal (2nd trimester)",
    "Fetal (3rd trimester)",
    "Neonatal",
    "Childhood (1-6Y)", 
    "Childhood (6-12Y)"))

gc()

pediatric_ref_bts@meta.data <- pediatric_ref_bts@meta.data |>
  mutate(across(where(is.factor), ~ droplevels(.x)))

qs_save(
  pediatric_ref_bts, 
  "/media/sofia/TOSHIBA EXT/developing_human_brain/BTS_atlas_seurat_pediatric.qs2", 
  nthreads = 6
)
 



