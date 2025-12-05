

convert_to_seurat <- function(adata) {

  # 1. let's get the adata object layers  
  counts <- t(Matrix::Matrix(py_to_r(adata$layers[["counts"]]), sparse = TRUE))
  
 # 2. let's get the metadata
  meta <- as.data.frame(py_to_r(adata$obs))
  cell_names <- rownames(meta)
  colnames(counts) <- cell_names
  
  
  # 3. now we create the seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta
  )
  
  # 4. let's add embeddings from obsm 
  for (key in adata$obsm_keys()) {
    mat <- py_to_r(adata$obsm[[key]])
    mat <- as.matrix(mat)
    
    # ensure rownames match cells
    rownames(mat) <- cell_names
    
    # Clean key name
    clean_name <- gsub("^X_|_", "", key)
    
    if (grepl("umap", tolower(clean_name))) {
      seurat_obj[[clean_name]] <- CreateDimReducObject(
        embeddings = mat,
        key = paste0(clean_name,"_"),
        assay = DefaultAssay(seurat_obj)
      )
    } else if (grepl("pca", tolower(clean_name))) {
      seurat_obj[[clean_name]] <- CreateDimReducObject(
        embeddings = mat,
        key = paste0(clean_name,"_"),
        assay = DefaultAssay(seurat_obj)
      )
    } else {
      seurat_obj@misc[[clean_name]] <- mat
    }
  }
  
  # --- 6. Variable features ---
  var_df <- as.data.frame(py_to_r(adata$var))
  if ("highly_variable" %in% colnames(var_df)) {
    hvf <- rownames(var_df)[which(var_df$highly_variable)]
    VariableFeatures(seurat_obj, assay = "RNA") <- hvf
  }
  
  return(seurat_obj)
}


plot_one_over_other <- function(umapdf, x, y, colvar, condition, group1, group2, integration) {
  p <- ggplot() +
    geom_point(data = umapdf %>% dplyr::filter(get(condition) == group1),
               aes(x = get(x), y = get(y), color = get(colvar)),
               alpha = 1, size = 0.5) +
    geom_point(data = umapdf %>% dplyr::filter(get(condition) == group2),
               aes(x = get(x), y = get(y), color = get(colvar)),
               alpha = 1, size = 0.5) +
    theme_classic() +
    labs(x = x, y = y, color = colvar) +
    ggtitle(integration) 
  return(p)
} 


plot_general <- function(umapdf, x, y, colvar) {
     p <- umapdf %>% 
    ggplot(aes(x = get(x), y = get(y), color = get(colvar))) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_classic() +
    labs(x = x, y = y, color = colvar) 
  
  return(p)
}

plot_with_patchwork <- function(colvar, condition, group1, group2) {

  plot_one_over_other(umap_etmr, "umapunintegrated_1", "umapunintegrated_2", colvar, condition, group1, group2, "unintegrated") +
    plot_one_over_other(umap_harm, "umapharmony_1", "umapharmony_2", colvar, condition, group1, group2, "harmony") +
    plot_layout(guides = "collect", axes = "collect") &
    plot_annotation(title = paste0(colvar, " distribution across integration methods in ETMR snRNA data")) 
}

plot_hist <- function(umapdf, x, colvar, integration) {

  h <- umapdf %>%
    ggplot(aes(x = get(x),  fill = get(colvar))) +
    geom_bar(position = "fill") +
    theme_classic() +
    ggtitle(paste(colvar, integration)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = x, fill = colvar) 
  
  return(h)
}


hist_with_patchwork <- function(colvar) {
  plot_hist(umap_etmr, "leiden_unintegrated", colvar, "unintegrated") +
    plot_hist(umap_harm, "leiden_harmony",  colvar,  "harmony") +
    plot_layout(guides = "collect", axes = "collect") &
    plot_annotation(title = paste0(colvar, " histogram across integration methods in ETMR snRNA data")) 
}

cellmarkers <- list(
  neuroepithelial = c("NES",
                      "SOX2",
                      "NOTCH1",
                      "HES1",
                      "CDH1",
                      "OCLN",
                      "SOX10"), 
  rg = c("VIM",
         "NES",
         "PAX6",
         "HES1",
         "HES5",
         "SOX2",
         "GFAP",
         "SLC1A3",
         "FABP7",
         "TNC",
         "CDH2"), 
  progenitors = c("EOMES",
                  "ASCL1",
                  "DCX",
                  "NEUROD1",
                  "TBR1",
                  "STMN1"),
  oligo = c("PDGFRA",
            "OLIG2",
            "MOG",
            "MBP",
            "SOX10"),
  schwann = c("MPZ",
              "NCAM1",
              "GAP43",
              "S100B",
              "DHH"),
  astrocytes = c("GFAP",
                 "SLC1A3",
                 "SLC1A2",
                 "S100B",
                 "ALDH1L1"),
  microglia = c("ITGAM",
                "PTPRC",
                "AIF1",
                "CD68",
                "CD40"),
  neuron = c("RBFOX3",
             "MAP2",
             "SYP",
             "DLG4",
             "BCL11B"),
  malignant = c("MCM7",
                "NME1",
                "NONO",
                "PFDN2",
                "NHP2",
                "VBP1",
                "C1orf43",
                "CCT3",
                "TXNRD1"),
  ipc = c("PPP1R17")
  
)

etmr_sig_pos <- c("ABCB5",
                  "ACRV1",
                  "ACVR2A",
                  "ACVR2B",
                  "ADAM7",
                  "AGO1",
                  "AGPS",
                  "APELA",
                  "BLTP3A",
                  "C16orf82",
                  "C22orf42",
                  "CALHM4",
                  "CBX2",
                  "CCDC150",
                  "CCNE1",
                  "CCNJ",
                  "CCNT2",
                  "CDC25A",
                  "CDC7",
                  "CHAC2",
                  "CHRNA5",
                  "CHRNB4",
                  "CRABP1",
                  "CRYGC",
                  "CSAG2",
                  "CSAG3",
                  "CTCFL",
                  "CTPS1",
                  "CXXC4",
                  "CYP27C1",
                  "DAAM1",
                  "DEPDC1B",
                  "DIS3L2",
                  "DNA2",
                  "DNMT3B",
                  "DTX4",
                  "EPB41",
                  "EVX2",
                  "EXOSC7",
                  "FAM161A",
                  "FANCD2OS",
                  "FIGNL2",
                  "FOXB1",
                  "GABRR1",
                  "GATA3",
                  "GBX1",
                  "GJC1",
                  "GPR20",
                  "GREB1",
                  "GREB1L",
                  "H1-2",
                  "H2BC11",
                  "H2BC26",
                  "H2BC9",
                  "H3C15",
                  "HCRTR2",
                  "HELLS",
                  "HES5",
                  "HIC2",
                  "HMX2",
                  "HOXA1",
                  "HOXA3",
                  "HOXB1",
                  "HRC",
                  "IDH1",
                  "IGF2BP1",
                  "INSM2",
                  "KCNH8",
                  "KLHL23",
                  "KPNA2",
                  "LBX1",
                  "LCOR",
                  "LIN28A",
                  "LIN28B",
                  "MAGEA10",
                  "MAGEA12",
                  "MAGEA2B",
                  "MAGEA3",
                  "MAGEA6",
                  "MBTD1",
                  "MEIS1",
                  "MSH2",
                  "MTF2",
                  "MTTP",
                  "NDST4",
                  "NEUROG3",
                  "NHEJ1",
                  "NKX1-1",
                  "NKX1-2",
                  "NKX6-3",
                  "NPPC",
                  "NR6A1",
                  "NSD2",
                  "NUP210",
                  "NUP62CL",
                  "ONECUT1",
                  "ORC1",
                  "OTP",
                  "P2RX3",
                  "PAPOLG",
                  "PBX2",
                  "PDE7A",
                  "PGAP1",
                  "PHF6",
                  "PHOX2B",
                  "PIK3C2B",
                  "PLEKHA8",
                  "PMFBP1",
                  "POLQ",
                  "POTEE",
                  "POTEF",
                  "POTEI",
                  "PRDM12",
                  "PRDM13",
                  "PRDM5",
                  "PRIM1",
                  "PRTG",
                  "PTBP2",
                  "RAD54L2",
                  "RALGPS2",
                  "RARB",
                  "RASAL2",
                  "RNF152",
                  "RPS6KA6",
                  "RRM2",
                  "RTKN2",
                  "SALL4",
                  "SANBR",
                  "SCML2",
                  "SCML4",
                  "SCN11A",
                  "SCUBE3",
                  "SEMG2",
                  "SERTM2",
                  "SESTD1",
                  "SFRP2",
                  "SHISA3",
                  "SINHCAF",
                  "SMARCC1",
                  "SOX11",
                  "SOX14",
                  "SOX3",
                  "SP8",
                  "SPSB4",
                  "SRPK1",
                  "STK26",
                  "STRBP",
                  "SUV39H2",
                  "TIA1",
                  "TMC4",
                  "TRIM71",
                  "TSGA10IP",
                  "ULBP3",
                  "USP44",
                  "VSX2",
                  "WNT7A",
                  "ZBTB39",
                  "ZFHX3",
                  "ZNF439",
                  "ZNF462",
                  "ZNF48",
                  "ZNF608",
                  "ZNF620",
                  "ZNF71")

